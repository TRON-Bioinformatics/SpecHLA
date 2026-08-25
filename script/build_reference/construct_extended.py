#!/usr/bin/env python
"""Construct the one-record-per-gene extended reference."""

import argparse
import re


GENES = ("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")
NON_NULL_SUFFIXES = ("N", "Q", "L", "S", "C", "A")


def read_fasta(path):
    records = {}
    name = None
    sequence = []
    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(sequence)
                fields = line[1:].split()
                name = fields[1] if len(fields) > 1 else fields[0]
                sequence = []
            else:
                if name is None:
                    raise ValueError("FASTA sequence encountered before a header")
                sequence.append(line)
    if name is not None:
        records[name] = "".join(sequence)
    return records


def load_representatives(path):
    selected = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            gene, allele = line.split("\t")[:2]
            selected[gene] = allele
    return selected


def select_allele(gene, requested, alleles):
    if requested in alleles:
        return requested
    pattern = re.compile(r"^%s\*(?:\d+):(?:\d+):(?:\d+):(?:\d+)$" % re.escape(gene))
    for allele in sorted(alleles):
        if pattern.match(allele) and not allele.endswith(NON_NULL_SUFFIXES):
            return allele
    raise ValueError("No non-null four-field allele found for %s" % gene)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gen", required=True)
    parser.add_argument("--extend", required=True)
    parser.add_argument("--representatives", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    alleles = read_fasta(args.gen)
    flanks = read_fasta(args.extend)
    requested = load_representatives(args.representatives)
    chosen = {}
    with open(args.out, "w") as output:
        for gene in GENES:
            available = [name for name in alleles if name.startswith(gene + "*")]
            allele = select_allele(gene, requested.get(gene, ""), set(available))
            if allele != requested.get(gene):
                print("warning: %s missing; using %s" % (requested.get(gene), allele))
            for suffix in ("1", "2"):
                flank_name = "HLA_%s_%s" % (gene, suffix)
                if flank_name not in flanks:
                    raise ValueError("Missing flank %s" % flank_name)
            chosen[gene] = allele
            sequence = flanks["HLA_%s_1" % gene] + alleles[allele] + flanks["HLA_%s_2" % gene]
            output.write(">HLA_%s\n%s\n" % (gene, sequence))
    return chosen


if __name__ == "__main__":
    main()
