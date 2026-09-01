"""Split the IMGT/HLA allele FASTA files into the per-gene SpecHLA layout.

SpecHLA aligns and types one gene at a time, so the monolithic IMGT/HLA
``hla_gen.fasta`` (genomic sequences) and ``hla_nuc.fasta`` (coding sequences)
are split into ``<db>/HLA/whole/HLA_<gene>.fasta`` and
``<db>/HLA/exon/HLA_<gene>.fasta`` respectively.
"""

import os

from . import fasta


def _collect_by_gene(path, genes, key):
    """Group the records of a FASTA file by HLA gene.

    ``key`` selects what is written as the record name: the plain allele for
    the genomic database, or the full IMGT header for the exon database, whose
    accession and length fields are parsed by the annotation steps.
    """
    records = {gene: [] for gene in genes}
    for header, sequence in fasta.iter_records(path):
        allele = fasta.allele_of(header)
        gene = fasta.gene_of(allele)
        if gene in records:
            records[gene].append((key(header, allele), sequence))
    return records


def _write_per_gene(records, genes, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    written = []
    for gene in genes:
        path = os.path.join(output_dir, "HLA_%s.fasta" % gene)
        with open(path, "w") as output:
            for name, sequence in sorted(records[gene]):
                fasta.write_record(output, name, sequence)
        written.append(path)
    return written


def split_genomic(gen_fasta, output_dir, genes):
    """Write the per-gene genomic (``whole``) databases and return their paths."""
    records = _collect_by_gene(gen_fasta, genes, key=lambda header, allele: allele)
    return _write_per_gene(records, genes, output_dir)


def split_coding(nuc_fasta, output_dir, genes):
    """Write the per-gene coding (``exon``) databases and return their paths."""
    records = _collect_by_gene(nuc_fasta, genes, key=lambda header, allele: header)
    return _write_per_gene(records, genes, output_dir)
