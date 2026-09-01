"""The HLA genes SpecHLA types, and the exons that define their G groups."""

# Class I and class II genes supported by the SpecHLA pipeline. Every build
# step derives its per-gene outputs from this list, so adding a gene here is
# the single change needed to extend the reference.
CLASS_I_GENES = ("A", "B", "C")
CLASS_II_GENES = ("DPA1", "DPB1", "DQA1", "DQB1", "DRB1")
GENES = CLASS_I_GENES + CLASS_II_GENES


def g_group_exons(gene):
    """Return the exons that define the G group of ``gene``.

    G groups are defined by exon 2 for every gene, plus exon 3 for class I.
    Emitting further exons would skew the identity sums computed by
    g_group_annotation.py, so the exon reference is restricted to these.
    """
    return (2, 3) if gene in CLASS_I_GENES else (2,)
