#!/usr/bin/env python
"""Construct the one-record-per-gene extended reference.

Besides the FASTA, this step writes ``reference.json``: the description of the
gene content and the coordinates that were actually produced. Every downstream
step reads that file instead of hardcoding gene names or intervals.
"""

import argparse
import json
import os
import re

from _spec import add_spec_arguments, spec_from_args

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
    parser.add_argument("--out", required=True)
    parser.add_argument("--reference-json", required=True)
    parser.add_argument("--selected-out")
    add_spec_arguments(parser)
    args = parser.parse_args()

    spec = spec_from_args(args)
    alleles = read_fasta(args.gen)
    flanks = read_fasta(args.extend)
    chosen = {}
    described = []
    with open(args.out, "w") as output:
        for gene, hla_class, requested in spec:
            available = [name for name in alleles if name.startswith(gene + "*")]
            allele = select_allele(gene, requested, set(available))
            if allele != requested:
                print("warning: %s missing; using %s" % (requested, allele))
            flank5_name = "HLA_%s_1" % gene
            flank3_name = "HLA_%s_2" % gene
            for flank_name in (flank5_name, flank3_name):
                if flank_name not in flanks:
                    raise ValueError("Missing flank %s" % flank_name)
            chosen[gene] = allele
            flank5 = flanks[flank5_name]
            flank3 = flanks[flank3_name]
            sequence = flank5 + alleles[allele] + flank3
            output.write(">HLA_%s\n%s\n" % (gene, sequence))
            # The gene body starts where the 5' flank ends and spans the
            # representative allele; downstream steps restrict calling and
            # typing to this interval.
            described.append({
                "gene": gene,
                "name": "HLA_%s" % gene,
                "hla_class": hla_class,
                "representative_allele": allele,
                "requested_allele": requested,
                "flank_5": len(flank5),
                "flank_3": len(flank3),
                "allele_length": len(alleles[allele]),
                "start": len(flank5),
                "end": len(flank5) + len(alleles[allele]),
                "length": len(sequence),
            })

    directory = os.path.dirname(os.path.abspath(args.reference_json))
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(args.reference_json, "w") as handle:
        json.dump({"version": 1, "genes": described}, handle, indent=2)
        handle.write("\n")

    if args.selected_out:
        with open(args.selected_out, "w") as selected_file:
            json.dump(chosen, selected_file, indent=2, sort_keys=True)
            selected_file.write("\n")
    return chosen


if __name__ == "__main__":
    main()
