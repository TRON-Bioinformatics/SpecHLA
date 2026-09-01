#!/usr/bin/env python3
"""Build the SpecHLA reference database from a read-only IMGT/HLA checkout.

SpecHLA does not ship a generated reference. This script turns an IMGT/HLA
release into the directory layout the pipeline expects:

    <out>/HLA/whole/HLA_<gene>.fasta   full-length allele database per gene
    <out>/HLA/exon/HLA_<gene>.fasta    coding-sequence database per gene
    <out>/HLA/hla_exons.fasta          G-group defining exons, from the XML release
    <out>/HLA/hla.ref.extend.fa        one flank-padded record per gene
    <out>/HLA/HLA_<gene>/              per-gene alignment reference and indexes
    <out>/HLA/HLA_<gene>.config.txt    tool paths consumed by the phasing steps
    <out>/ref/                         references and indexes for read extraction
    <out>/BUILD_MANIFEST.json          provenance of the produced reference

Point SPECHLA_DB at the output directory afterwards, or build straight into the
default location reported by ``--help``.
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import spechla_paths  # noqa: E402  (needs the path above to be importable)
from reference_build import construct_extended, extract_exons, split_by_gene  # noqa: E402
from reference_build.genes import GENES  # noqa: E402

# Files a usable IMGT/HLA checkout must provide. hla_gen.fasta and hla.xml are
# handled separately because IMGT ships them zipped in some releases.
REQUIRED_IMGT_FILES = ("hla_nuc.fasta", "Allelelist.txt", "wmda/hla_nom_g.txt")

# Bundled inputs that are copied verbatim into the built reference.
ASSETS_FOR_HLA_DIR = ("extend.fa", "HLA_FREQ_HLA_I_II.txt")
ASSETS_FOR_REF_DIR = (
    "DRB1_dup_extract.fasta",
    "DRB1_dup_extract_ref.fasta",
    "hla_gen.format.filter.extend.DRB.no26789.fasta",
    "hla_gen.format.filter.extend.DRB.no26789.v2.fasta",
)

# Presence of this index is used as the marker of a previously finished build.
BUILD_COMPLETE_MARKER = os.path.join("ref", "hla.ref.extend.fa.bwt")

TEMP_DIR_NAME = ".build_tmp"


def log(message):
    """Report build progress on stderr, leaving stdout free for tool output."""
    sys.stderr.write("[build_reference] %s\n" % message)


def run(command, cwd=None):
    """Run an external command, failing the build if it does not succeed."""
    log("running: %s" % " ".join(command))
    subprocess.check_call(command, cwd=cwd)


def tool_available(name):
    """Return whether ``name`` can be found on PATH."""
    return shutil.which(name) is not None


def require_tools(*names):
    """Abort with a clear message if a mandatory build tool is missing."""
    missing = [name for name in names if not tool_available(name)]
    if missing:
        raise SystemExit(
            "Required tool(s) not on PATH: %s. Activate the SpecHLA environment "
            "before building a reference." % ", ".join(missing))


# --- IMGT/HLA input handling ------------------------------------------------


def resolve_imgt_dir(requested):
    """Return the IMGT/HLA checkout to build from.

    An explicit --imgt wins; otherwise the shared resolution in spechla_paths
    is used so the builder honours SPECHLA_IMGT and the conda layout exactly
    like the runtime scripts do.
    """
    if requested:
        return requested
    try:
        return spechla_paths.get_imgt_dir()
    except RuntimeError as error:
        raise SystemExit("%s Alternatively pass --imgt." % error)


def validate_imgt_dir(imgt_dir):
    """Fail early if the checkout cannot produce a complete reference."""
    if not os.path.isdir(imgt_dir):
        raise SystemExit("IMGT/HLA directory does not exist: %s" % imgt_dir)
    missing = [name for name in REQUIRED_IMGT_FILES
               if not os.path.isfile(os.path.join(imgt_dir, name))]
    if missing:
        raise SystemExit("IMGT/HLA checkout %s is missing: %s"
                         % (imgt_dir, ", ".join(missing)))


def unpacked_imgt_file(imgt_dir, relative_path, temp_dir):
    """Return a plain-file path for an IMGT artefact that may be zipped.

    Releases ship hla_gen.fasta and xml/hla.xml either directly or as a sibling
    .zip, so the archive is expanded into the build's temporary directory when
    only the compressed form exists.
    """
    direct = os.path.join(imgt_dir, relative_path)
    if os.path.isfile(direct):
        return direct
    archive = direct + ".zip"
    if not os.path.isfile(archive):
        raise SystemExit("IMGT/HLA checkout is missing %s (or %s.zip)"
                         % (relative_path, relative_path))
    log("unpacking %s" % archive)
    os.makedirs(temp_dir, exist_ok=True)
    with zipfile.ZipFile(archive) as handle:
        handle.extractall(temp_dir)
    unpacked = os.path.join(temp_dir, os.path.basename(relative_path))
    if not os.path.isfile(unpacked):
        raise SystemExit("%s did not contain %s"
                         % (archive, os.path.basename(relative_path)))
    return unpacked


def imgt_version(imgt_dir):
    """Return the release version recorded in the checkout's Allelelist.txt."""
    with open(os.path.join(imgt_dir, "Allelelist.txt")) as handle:
        for line in handle:
            if line.startswith("# version:"):
                return line.split(":", 1)[1].strip()
            if not line.startswith("#"):
                break
    return "unknown"


# --- Indexing helpers -------------------------------------------------------


def samtools_faidx(fasta_path):
    """Create the .fai index every SpecHLA step expects next to a FASTA."""
    run(["samtools", "faidx", fasta_path])


def makeblastdb(fasta_path, prefix):
    """Create a nucleotide BLAST database used by the annotation steps."""
    run(["makeblastdb", "-in", fasta_path, "-dbtype", "nucl",
         "-parse_seqids", "-out", prefix])


# --- Build steps ------------------------------------------------------------


def copy_allele_metadata(imgt_dir, hla_dir):
    """Copy the allele list and G-group nomenclature next to the databases.

    Both files are read at runtime to translate between allele names and G
    groups, so they must match the release the reference was built from.
    """
    shutil.copy2(os.path.join(imgt_dir, "Allelelist.txt"),
                 os.path.join(hla_dir, "Allelelist.txt"))
    shutil.copy2(os.path.join(imgt_dir, "wmda", "hla_nom_g.txt"),
                 os.path.join(hla_dir, "hla_nom_g.txt"))


def build_gene_databases(gen_fasta, nuc_fasta, hla_dir):
    """Split IMGT into per-gene databases and index them for lookup.

    DRB1 typing compares candidate alleles against coding sequences, so the
    DRB1 exon database is additionally published inside ``whole/`` under the
    name the typing step looks for.
    """
    whole_dir = os.path.join(hla_dir, "whole")
    exon_dir = os.path.join(hla_dir, "exon")
    whole_fastas = split_by_gene.split_genomic(gen_fasta, whole_dir, GENES)
    exon_fastas = split_by_gene.split_coding(nuc_fasta, exon_dir, GENES)

    drb1_exon_copy = os.path.join(whole_dir, "HLA_DRB1.exon.fasta")
    shutil.copy2(os.path.join(exon_dir, "HLA_DRB1.fasta"), drb1_exon_copy)

    for fasta_path in whole_fastas + [drb1_exon_copy]:
        samtools_faidx(fasta_path)
        # The whole-gene BLAST databases drop the .fasta suffix because the
        # pipeline refers to them as <dir>/HLA_<gene>.
        makeblastdb(fasta_path, os.path.splitext(fasta_path)[0])
    for fasta_path in exon_fastas:
        samtools_faidx(fasta_path)
        makeblastdb(fasta_path, fasta_path)


def build_exon_reference(xml_path, hla_dir):
    """Derive and index the G-group defining exon reference from the XML release."""
    exon_reference = os.path.join(hla_dir, "hla_exons.fasta")
    log("extracting G-group exons from %s" % xml_path)
    extract_exons.extract_exons(xml_path, exon_reference, GENES)
    samtools_faidx(exon_reference)
    makeblastdb(exon_reference, exon_reference)
    return exon_reference


def copy_bundled_assets(assets_dir, hla_dir, ref_dir):
    """Copy the immutable inputs that cannot be derived from IMGT/HLA."""
    for name in ASSETS_FOR_HLA_DIR:
        shutil.copy2(os.path.join(assets_dir, name), os.path.join(hla_dir, name))
    for name in ASSETS_FOR_REF_DIR:
        shutil.copy2(os.path.join(assets_dir, name), os.path.join(ref_dir, name))


def build_extended_reference(gen_fasta, assets_dir, representatives, hla_dir, ref_dir):
    """Build the flank-padded per-gene reference and its alignment indexes.

    The same FASTA is published under ``HLA/`` for the typing steps and under
    ``ref/`` for read extraction, which expects a BWA index and a sequence
    dictionary alongside it.
    """
    extended = os.path.join(hla_dir, "hla.ref.extend.fa")
    selected = construct_extended.construct_extended_reference(
        gen_fasta=gen_fasta,
        extend_fasta=os.path.join(assets_dir, "extend.fa"),
        representatives_file=representatives,
        output_path=extended,
        genes=GENES,
        log=log,
    )
    samtools_faidx(extended)

    published = os.path.join(ref_dir, "hla.ref.extend.fa")
    shutil.copy2(extended, published)
    shutil.copy2(extended + ".fai", published + ".fai")
    run(["bwa", "index", published])
    with open(os.path.join(ref_dir, "hla.ref.extend.dict"), "w") as dict_file:
        subprocess.check_call(["samtools", "dict", published], stdout=dict_file)
    return extended, selected


def build_per_gene_references(extended_reference, hla_dir):
    """Split the extended reference into the per-gene alignment references.

    Phasing and assembly realign reads gene by gene and read the tool paths
    from HLA_<gene>.config.txt, which is generated here so the config always
    matches the reference that was just built.
    """
    for gene in GENES:
        gene_dir = os.path.join(hla_dir, "HLA_%s" % gene)
        os.makedirs(gene_dir, exist_ok=True)
        gene_fasta = os.path.join(gene_dir, "HLA_%s.fa" % gene)
        with open(gene_fasta, "w") as output:
            subprocess.check_call(
                ["samtools", "faidx", extended_reference, "HLA_%s" % gene],
                stdout=output)
        samtools_faidx(gene_fasta)
        run(["bwa", "index", gene_fasta])
        makeblastdb(gene_fasta, os.path.join(gene_dir, "HLA_%s" % gene))
        if tool_available("faToTwoBit"):
            run(["faToTwoBit", gene_fasta,
                 os.path.join(gene_dir, "HLA_%s.2bit" % gene)])
        else:
            log("faToTwoBit not on PATH; skipping HLA_%s.2bit (blat realignment "
                "will be unavailable)" % gene)
        with open(os.path.join(hla_dir, "HLA_%s.config.txt" % gene), "w") as config:
            config.write("bwa=%s\nfreebayes=%s\nblat=%s/\n"
                         % (gene_fasta, gene_fasta, gene_dir))


def index_extraction_references(ref_dir, threads, skip_novoindex):
    """Index the bundled references used to pull HLA reads out of a BAM.

    Novoalign is licensed software: its index is only built when both the
    binary and a license file are present, and Bowtie2 is always indexed as the
    freely available fallback aligner.
    """
    samtools_faidx(os.path.join(ref_dir, "DRB1_dup_extract.fasta"))
    duplicate_reference = os.path.join(ref_dir, "DRB1_dup_extract_ref.fasta")
    makeblastdb(duplicate_reference, duplicate_reference)

    for name in ASSETS_FOR_REF_DIR:
        if not name.startswith("hla_gen.format.filter.extend.DRB.no26789"):
            continue
        fasta_path = os.path.join(ref_dir, name)
        samtools_faidx(fasta_path)
        if not skip_novoindex and novoalign_licensed():
            run(["novoindex", "-k", "14", "-s", "1",
                 os.path.splitext(fasta_path)[0] + ".ndx", fasta_path])
        run(["bowtie2-build", "--threads", str(threads), fasta_path, fasta_path])


def novoalign_licensed():
    """Return whether a licensed novoalign installation is available."""
    novoalign = shutil.which("novoalign")
    if not novoalign:
        return False
    return os.path.isfile(os.path.join(os.path.dirname(novoalign), "novoalign.lic"))


def build_linked_read_reference(ref_dir):
    """Create the Long Ranger reference needed for 10X linked-read input."""
    if not tool_available("longranger"):
        log("longranger not on PATH; 10X input will require a manual mkref")
        return
    run(["longranger", "mkref", "hla.ref.extend.fa"], cwd=ref_dir)


# --- Provenance -------------------------------------------------------------


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def spechla_revision():
    """Return the SpecHLA commit the reference was built with, if known."""
    source_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        revision = subprocess.check_output(
            ["git", "-C", source_root, "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL)
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return revision.decode().strip()


def write_build_manifest(out_dir, imgt_dir, assets_dir, selected_alleles):
    """Record what a reference was built from so results stay traceable.

    Generated references are not version controlled, so the manifest is the
    only way to tell which IMGT release, bundled assets and representative
    alleles produced a given directory.
    """
    manifest = {
        "imgt_version": imgt_version(imgt_dir),
        "imgt_source": os.path.abspath(imgt_dir),
        "built_at": datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "spechla_ref": spechla_revision(),
        "representative_alleles": selected_alleles,
        "bundled_assets_sha256": {
            name: sha256_of(os.path.join(assets_dir, name))
            for name in sorted(os.listdir(assets_dir))
            if os.path.isfile(os.path.join(assets_dir, name))
        },
    }
    path = os.path.join(out_dir, "BUILD_MANIFEST.json")
    with open(path, "w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return path


# --- Command line -----------------------------------------------------------


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--imgt", metavar="DIR", default=None,
        help="Read-only IMGT/HLA checkout to build from "
             "(default: $SPECHLA_IMGT or $CONDA_PREFIX/share/spechla/imgt).")
    parser.add_argument(
        "--out", metavar="DIR", default=spechla_paths.get_default_db_dir() or None,
        help="Directory to write the reference to; point SPECHLA_DB here "
             "afterwards (default: $SPECHLA_DB or "
             "$CONDA_PREFIX/share/spechla/db).")
    parser.add_argument(
        "--assets", metavar="DIR", default=None,
        help="Bundled reference construction assets, i.e. the flanks and "
             "prebuilt extraction references that cannot be derived from IMGT "
             "(default: $SPECHLA_ASSETS or share/reference_assets).")
    parser.add_argument(
        "--representative-alleles", metavar="FILE", default=None,
        help="TSV pinning one representative allele per gene, used as the "
             "backbone of the extended reference "
             "(default: representative_alleles.tsv from the assets directory).")
    parser.add_argument(
        "--threads", type=int, default=4, metavar="N",
        help="Threads to use for bowtie2-build (default: %(default)s).")
    parser.add_argument(
        "--skip-novoindex", action="store_true",
        help="Do not build the Novoalign index even if a licensed novoalign "
             "is installed; Bowtie2 is used for read extraction instead.")
    parser.add_argument(
        "--force", action="store_true",
        help="Rebuild into an output directory that already holds a reference.")
    args = parser.parse_args(argv)
    if not args.out:
        parser.error("--out is required when neither SPECHLA_DB nor "
                     "CONDA_PREFIX is set")
    return args


def prepare_output_dir(out_dir, force):
    """Create the output layout, refusing to silently overwrite a build."""
    if os.path.exists(os.path.join(out_dir, BUILD_COMPLETE_MARKER)) and not force:
        raise SystemExit("Reference already exists at %s; use --force to rebuild"
                         % out_dir)
    hla_dir = os.path.join(out_dir, "HLA")
    ref_dir = os.path.join(out_dir, "ref")
    temp_dir = os.path.join(out_dir, TEMP_DIR_NAME)
    for path in (hla_dir, ref_dir, temp_dir):
        os.makedirs(path, exist_ok=True)
    if force:
        # Drop the previous reference so a failed rebuild cannot be mistaken
        # for a complete one.
        for stale in (os.path.join(hla_dir, "hla.ref.extend.fa"),
                      os.path.join(ref_dir, "hla.ref.extend.fa")):
            if os.path.exists(stale):
                os.remove(stale)
    return hla_dir, ref_dir, temp_dir


def main(argv=None):
    args = parse_arguments(argv)
    require_tools("samtools", "makeblastdb", "bwa", "bowtie2-build")

    imgt_dir = resolve_imgt_dir(args.imgt)
    validate_imgt_dir(imgt_dir)
    assets_dir = args.assets or spechla_paths.get_assets_dir()
    representatives = (args.representative_alleles
                       or os.path.join(assets_dir, "representative_alleles.tsv"))

    hla_dir, ref_dir, temp_dir = prepare_output_dir(args.out, args.force)
    gen_fasta = unpacked_imgt_file(imgt_dir, "hla_gen.fasta", temp_dir)
    xml_path = unpacked_imgt_file(imgt_dir, os.path.join("xml", "hla.xml"), temp_dir)

    log("building reference from %s into %s" % (imgt_dir, args.out))
    copy_allele_metadata(imgt_dir, hla_dir)
    build_gene_databases(gen_fasta, os.path.join(imgt_dir, "hla_nuc.fasta"), hla_dir)
    copy_bundled_assets(assets_dir, hla_dir, ref_dir)
    build_exon_reference(xml_path, hla_dir)
    extended, selected = build_extended_reference(
        gen_fasta, assets_dir, representatives, hla_dir, ref_dir)
    build_per_gene_references(extended, hla_dir)
    index_extraction_references(ref_dir, args.threads, args.skip_novoindex)
    build_linked_read_reference(ref_dir)
    write_build_manifest(args.out, imgt_dir, assets_dir, selected)

    shutil.rmtree(temp_dir, ignore_errors=True)
    log("reference ready; run: export SPECHLA_DB=%s" % os.path.abspath(args.out))


if __name__ == "__main__":
    main()
