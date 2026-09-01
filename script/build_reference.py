#!/usr/bin/env python3
"""Build the SpecHLA reference database from a read-only IMGT/HLA checkout.

SpecHLA does not ship a generated reference. This script turns an IMGT/HLA
release into the directory layout the pipeline expects:

    <out>/HLA/whole/HLA_<gene>.fasta   full-length allele database per gene
    <out>/HLA/exon/HLA_<gene>.fasta    coding-sequence database per gene

Further stages of the reference (the exon reference, the extended alignment
reference and the read extraction indexes) are added on top of this layout.

Point SPECHLA_DB at the output directory afterwards, or build straight into the
default location reported by ``--help``.
"""

import argparse
import os
import shutil
import subprocess
import sys
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import spechla_paths  # noqa: E402  (needs the path above to be importable)
from reference_build import split_by_gene  # noqa: E402
from reference_build.genes import GENES  # noqa: E402

# Files a usable IMGT/HLA checkout must provide. hla_gen.fasta and hla.xml are
# handled separately because IMGT ships them zipped in some releases.
REQUIRED_IMGT_FILES = ("hla_nuc.fasta", "Allelelist.txt", "wmda/hla_nom_g.txt")


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
    args = parser.parse_args(argv)
    if not args.out:
        parser.error("--out is required when neither SPECHLA_DB nor "
                     "CONDA_PREFIX is set")
    return args


def prepare_output_dir(out_dir):
    """Create the directory layout the build stages write into."""
    hla_dir = os.path.join(out_dir, "HLA")
    ref_dir = os.path.join(out_dir, "ref")
    temp_dir = os.path.join(out_dir, TEMP_DIR_NAME)
    for path in (hla_dir, ref_dir, temp_dir):
        os.makedirs(path, exist_ok=True)
    return hla_dir, ref_dir, temp_dir


def main(argv=None):
    args = parse_arguments(argv)
    require_tools("samtools", "makeblastdb")

    imgt_dir = resolve_imgt_dir(args.imgt)
    validate_imgt_dir(imgt_dir)

    hla_dir, _ref_dir, temp_dir = prepare_output_dir(args.out)
    gen_fasta = unpacked_imgt_file(imgt_dir, "hla_gen.fasta", temp_dir)

    log("building reference from %s into %s" % (imgt_dir, args.out))
    copy_allele_metadata(imgt_dir, hla_dir)
    build_gene_databases(gen_fasta, os.path.join(imgt_dir, "hla_nuc.fasta"), hla_dir)

    shutil.rmtree(temp_dir, ignore_errors=True)
    log("reference ready; run: export SPECHLA_DB=%s" % os.path.abspath(args.out))


if __name__ == "__main__":
    main()
