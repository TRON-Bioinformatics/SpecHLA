"""Extract the G-group defining exon sequences from the IMGT/HLA XML database.

The XML release is the only IMGT/HLA artefact that carries exon boundaries for
every allele, so the exon reference used for G-group annotation is derived from
it instead of from a checked-in, quickly outdated FASTA file.
"""

import re
import xml.etree.ElementTree as ET

from . import fasta
from .genes import g_group_exons


def _local_name(tag):
    """Strip the XML namespace from ``tag``, which varies between releases."""
    return tag.rsplit("}", 1)[-1]


def _gene_of(allele_name, genes):
    for gene in genes:
        if re.match(r"^HLA-%s\*" % re.escape(gene), allele_name):
            return gene
    return None


def _allele_name(element):
    """Return the ``HLA-``-prefixed allele name of an ``<allele>`` element."""
    name = element.attrib.get("name", "")
    if not name:
        return ""
    return name if name.startswith("HLA-") else "HLA-" + name


def _full_sequence(element):
    """Return the genomic sequence of an allele, without whitespace."""
    for child in element.iter():
        if _local_name(child.tag).lower() == "sequence":
            return "".join(child.itertext()).replace(" ", "").replace("\n", "")
    return ""


def _iter_exons(element):
    """Yield ``(exon_number, start, end)`` for the exon features of an allele.

    Coordinates are 1-based and inclusive, as stored in the XML. Features
    without a parsable number or without coordinates are skipped, because some
    alleles are only partially characterised.
    """
    for feature in element.iter():
        if feature.attrib.get("featuretype", "").lower() != "exon":
            continue
        number = re.search(r"(\d+)", feature.attrib.get("name", ""))
        coordinates = next(
            (child for child in feature.iter()
             if _local_name(child.tag).lower() == "sequencecoordinates"),
            None,
        )
        if not number or coordinates is None:
            continue
        try:
            yield (int(number.group(1)),
                   int(coordinates.attrib["start"]),
                   int(coordinates.attrib["end"]))
        except (KeyError, ValueError):
            continue


def extract_exons(xml_path, output_path, genes):
    """Write ``><allele>|Exon|<number>`` records for the G-group exons.

    The XML is streamed with ``iterparse`` and each allele is cleared after use
    because the full IMGT/HLA XML does not fit comfortably in memory.
    """
    with open(output_path, "w") as output:
        for _, element in ET.iterparse(xml_path, events=("end",)):
            if _local_name(element.tag).lower() != "allele":
                continue
            name = _allele_name(element)
            gene = _gene_of(name, genes)
            sequence = _full_sequence(element) if gene else ""
            if not sequence:
                element.clear()
                continue
            wanted = g_group_exons(gene)
            for number, start, end in _iter_exons(element):
                if number not in wanted:
                    continue
                if start < 1 or end < start or end > len(sequence):
                    raise ValueError(
                        "Invalid coordinates for %s exon %s: %s-%s"
                        % (name, number, start, end))
                fasta.write_record(
                    output, "%s|Exon|%s" % (name, number), sequence[start - 1:end])
            element.clear()
    return output_path
