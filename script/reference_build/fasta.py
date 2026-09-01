"""Minimal FASTA helpers shared by the reference build steps.

The reference build runs before any SpecHLA environment is guaranteed to be
present, so these helpers deliberately avoid third-party dependencies such as
Biopython and only implement what the build steps need.
"""


def iter_records(path):
    """Yield ``(header, sequence)`` pairs from a FASTA file.

    The header is returned without the leading ``>`` and with the original
    whitespace-separated fields intact, because IMGT/HLA encodes the allele
    name in the second field (e.g. ``>HLA:HLA00001 A*01:01:01:01 3503 bp``).
    """
    header = None
    sequence = []
    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:]
                sequence = []
            else:
                if header is None:
                    raise ValueError(
                        "FASTA sequence encountered before a header in %s" % path)
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def allele_of(header):
    """Return the allele name encoded in an IMGT/HLA FASTA header.

    IMGT/HLA puts its internal accession first and the allele name second, so
    the allele is read from field two and falls back to field one for FASTA
    files that only carry a plain name (such as the bundled flank sequences).
    """
    fields = header.split()
    if not fields:
        return ""
    return fields[1] if len(fields) > 1 else fields[0]


def gene_of(allele):
    """Return the HLA gene an allele belongs to, e.g. ``A*01:01`` -> ``A``."""
    return allele.split("*", 1)[0] if "*" in allele else None


def read_alleles(path):
    """Read a FASTA file into an ``{allele_name: sequence}`` mapping."""
    return {allele_of(header): sequence for header, sequence in iter_records(path)}


def write_record(handle, name, sequence):
    """Write a single unwrapped FASTA record."""
    handle.write(">%s\n%s\n" % (name, sequence))
