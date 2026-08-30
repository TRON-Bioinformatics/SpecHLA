#!/usr/bin/env python
"""Materialise the gene-derived files of a reference.

Run in two passes by ``build_reference.sh``:

* from the gene spec, before anything is built, to publish the gene list and
  the curated region files restricted to the selected genes;
* from ``reference.json``, once the extended reference exists, to publish the
  per-gene intervals that were actually produced.

Everything written here lives inside the built reference, so downstream shell
steps never need a hardcoded gene name or coordinate.
"""

import argparse
import json
import os

from _spec import (
    add_spec_arguments,
    load_assembly_regions,
    load_exon_extent,
    spec_from_args,
)


def _write_lines(path, lines):
    with open(path, "w") as handle:
        for line in lines:
            handle.write("%s\n" % line)


def emit_from_spec(args):
    spec = spec_from_args(args)
    genes = [gene for gene, _, _ in spec]
    os.makedirs(args.out, exist_ok=True)
    _write_lines(os.path.join(args.out, "gene_list.txt"), genes)

    assets = args.assets
    whole, exon = load_assembly_regions(
        os.path.join(assets, "assembly_regions.tsv"), genes
    )
    _write_lines(os.path.join(args.out, "select.region.txt"), whole)
    _write_lines(os.path.join(args.out, "select.region.exon.txt"), exon)
    _write_lines(
        os.path.join(args.out, "exon_extent.bed"),
        load_exon_extent(os.path.join(assets, "exon_extent.bed"), genes),
    )


def emit_from_reference(args):
    with open(args.reference_json) as handle:
        reference = json.load(handle)
    os.makedirs(args.out, exist_ok=True)
    _write_lines(
        os.path.join(args.out, "gene_regions.txt"),
        [
            "%s:%s-%s" % (entry["name"], entry["start"], entry["end"])
            for entry in reference["genes"]
        ],
    )
    _write_lines(
        os.path.join(args.out, "gene_list.txt"),
        [entry["gene"] for entry in reference["genes"]],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, help="<db>/HLA directory")
    parser.add_argument("--assets", help="Bundled reference assets directory")
    parser.add_argument("--reference-json", help="Built <db>/reference.json")
    add_spec_arguments(parser)
    args = parser.parse_args()
    if args.reference_json:
        emit_from_reference(args)
    elif args.assets:
        emit_from_spec(args)
    else:
        parser.error("either --assets or --reference-json is required")


if __name__ == "__main__":
    main()
