#!/usr/bin/env python
"""Split IMGT/HLA FASTA files into the per-gene SpecHLA layout."""

import argparse
import os


GENES = ("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")


def read_fasta(path):
    name = None
    sequence = []
    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(sequence)
                name = line[1:]
                sequence = []
            else:
                if name is None:
                    raise ValueError("FASTA sequence encountered before a header")
                sequence.append(line)
    if name is not None:
        yield name, "".join(sequence)


def gene_from_allele(allele):
    return allele.split("*", 1)[0] if "*" in allele else None


def split_gen(path, output_dir, genes=GENES):
    records = {gene: [] for gene in genes}
    for header, sequence in read_fasta(path):
        fields = header.split()
        if len(fields) < 2:
            continue
        allele = fields[1]
        gene = gene_from_allele(allele)
        if gene in records:
            records[gene].append((allele, sequence))
    for gene in genes:
        with open(os.path.join(output_dir, "whole", "HLA_%s.fasta" % gene), "w") as out:
            for allele, sequence in sorted(records[gene]):
                out.write(">%s\n%s\n" % (allele, sequence))


def split_nuc(path, output_dir, genes=GENES):
    records = {gene: [] for gene in genes}
    for header, sequence in read_fasta(path):
        fields = header.split()
        if len(fields) < 2:
            continue
        allele = fields[1]
        gene = gene_from_allele(allele)
        if gene in records:
            records[gene].append((header, sequence))
    for gene in genes:
        with open(os.path.join(output_dir, "exon", "HLA_%s.fasta" % gene), "w") as out:
            for header, sequence in sorted(records[gene]):
                out.write(">%s\n%s\n" % (header, sequence))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gen", required=True)
    parser.add_argument("--nuc", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    os.makedirs(os.path.join(args.out, "whole"), exist_ok=True)
    os.makedirs(os.path.join(args.out, "exon"), exist_ok=True)
    split_gen(args.gen, args.out)
    split_nuc(args.nuc, args.out)


if __name__ == "__main__":
    main()
