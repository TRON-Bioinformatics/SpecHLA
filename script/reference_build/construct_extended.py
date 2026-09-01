"""Construct the extended, one-record-per-gene HLA reference.

Reads are first aligned against a single representative allele per gene. That
allele is padded with the bundled upstream and downstream flanks (extend.fa) so
reads overlapping the gene boundaries still align, which is what makes this
reference "extended".
"""

import re

from . import fasta

# Suffixes that mark alleles which are not normally expressed (null, questionable,
# low, secreted, cytoplasmic, aberrant). They make poor alignment references.
NON_EXPRESSED_SUFFIXES = ("N", "Q", "L", "S", "C", "A")

# Reference alleles are pinned to four fields so a build stays reproducible
# even when IMGT adds new sub-variants of the requested allele.
FOUR_FIELD_PATTERN = r"^%s\*(?:\d+):(?:\d+):(?:\d+):(?:\d+)$"


def read_representatives(path):
    """Read the ``gene<TAB>allele`` file pinning one reference allele per gene."""
    requested = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            gene, allele = line.split("\t")[:2]
            requested[gene] = allele
    return requested


def select_representative(gene, requested, available):
    """Return the allele to use for ``gene``.

    The pinned allele is preferred. IMGT occasionally renames or withdraws an
    allele, so the build falls back to the first expressed four-field allele of
    the gene instead of failing the whole release upgrade.
    """
    if requested in available:
        return requested
    pattern = re.compile(FOUR_FIELD_PATTERN % re.escape(gene))
    for allele in sorted(available):
        if pattern.match(allele) and not allele.endswith(NON_EXPRESSED_SUFFIXES):
            return allele
    raise ValueError("No expressed four-field allele found for %s" % gene)


def construct_extended_reference(gen_fasta, extend_fasta, representatives_file,
                                 output_path, genes, log=None):
    """Write the extended reference and return the ``{gene: allele}`` selection.

    The returned mapping records the alleles actually used, including fallbacks,
    so the build manifest describes the reference that was produced rather than
    the one that was requested.
    """
    alleles = fasta.read_alleles(gen_fasta)
    flanks = fasta.read_alleles(extend_fasta)
    requested = read_representatives(representatives_file)
    selected = {}
    with open(output_path, "w") as output:
        for gene in genes:
            available = {name for name in alleles if name.startswith(gene + "*")}
            allele = select_representative(gene, requested.get(gene, ""), available)
            if log is not None and allele != requested.get(gene):
                log("representative %s is unavailable; using %s instead"
                    % (requested.get(gene) or "<unset>", allele))
            upstream = flanks.get("HLA_%s_1" % gene)
            downstream = flanks.get("HLA_%s_2" % gene)
            if upstream is None or downstream is None:
                raise ValueError("Missing flanking sequences for HLA_%s" % gene)
            selected[gene] = allele
            fasta.write_record(output, "HLA_%s" % gene,
                               upstream + alleles[allele] + downstream)
    return selected
