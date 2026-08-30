#!/usr/bin/env python
"""Extract IMGT/HLA exon sequences from the XML database."""

import argparse
import re
import xml.etree.ElementTree as ET

from _spec import G_GROUP_EXONS, add_spec_arguments, spec_from_args


def gene_of(name, genes):
    for gene in genes:
        if re.match(r"^HLA-%s\*" % gene, name):
            return gene
    return None


def local_name(tag):
    return tag.rsplit("}", 1)[-1]


def allele_name(element):
    name = element.attrib.get("name", "")
    if not name:
        return ""
    if name.startswith("HLA-"):
        return name
    return "HLA-" + name


def sequence_text(element):
    for child in element.iter():
        if local_name(child.tag).lower() == "sequence":
            return "".join(child.itertext()).replace(" ", "").replace("\n", "")
    return ""


def extract(xml_path, output_path, spec):
    # G groups are defined by exon 2 for every gene plus exon 3 for class I
    # only. Emitting further exons would skew the identity sums computed by
    # g_group_annotation.py.
    wanted_by_gene = dict(
        (gene, G_GROUP_EXONS[hla_class]) for gene, hla_class, _ in spec
    )
    genes = [gene for gene, _, _ in spec]
    with open(output_path, "w") as output:
        for _, element in ET.iterparse(xml_path, events=("end",)):
            if local_name(element.tag).lower() != "allele":
                continue
            name = allele_name(element)
            gene = gene_of(name, genes)
            if gene is None:
                element.clear()
                continue
            wanted_exons = wanted_by_gene[gene]
            sequence = sequence_text(element)
            if not sequence:
                element.clear()
                continue
            for feature in element.iter():
                if feature.attrib.get("featuretype", "").lower() != "exon":
                    continue
                exon_name = feature.attrib.get("name", "")
                match = re.search(r"(\d+)", exon_name)
                coordinates = next(
                    (child for child in feature.iter()
                     if local_name(child.tag).lower() == "sequencecoordinates"),
                    None,
                )
                if not match or coordinates is None:
                    continue
                if int(match.group(1)) not in wanted_exons:
                    continue
                try:
                    start = int(coordinates.attrib["start"])
                    end = int(coordinates.attrib["end"])
                except (KeyError, ValueError):
                    continue
                if start < 1 or end < start or end > len(sequence):
                    raise ValueError(
                        "Invalid coordinates for %s exon %s: %s-%s"
                        % (name, match.group(1), start, end)
                    )
                output.write(
                    ">%s|Exon|%s\n%s\n"
                    % (name, match.group(1), sequence[start - 1:end])
                )
            element.clear()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--xml", required=True)
    parser.add_argument("--out", required=True)
    add_spec_arguments(parser)
    args = parser.parse_args()
    extract(args.xml, args.out, spec_from_args(args))


if __name__ == "__main__":
    main()
