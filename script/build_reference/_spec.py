"""Import helper so build steps can use script/reference_config.py."""

import os
import sys

sys.path.insert(
    0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)

from reference_config import (  # noqa: E402,F401
    G_GROUP_EXONS,
    EXON_EXTENT_EXONS,
    load_assembly_regions,
    load_exon_extent,
    load_gene_spec,
)


def add_spec_arguments(parser):
    """Add the arguments every build step uses to locate the gene spec."""
    parser.add_argument(
        "--genes-config",
        help="Path to genes.tsv (defaults to the bundled reference assets)",
    )
    parser.add_argument(
        "--genes",
        help="Comma separated subset of genes to build (default: all)",
    )
    return parser


def spec_from_args(args):
    genes = None
    if getattr(args, "genes", None):
        genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    return load_gene_spec(getattr(args, "genes_config", None), genes)
