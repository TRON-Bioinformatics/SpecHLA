"""
Single source of truth for the gene content and coordinates of a SpecHLA
reference.

Two kinds of information are handled here:

1. The *build input specification*, kept next to the other bundled assets in
   ``share/reference_assets``:

     genes.tsv             which genes exist, their HLA class and the
                           representative allele used for the extended
                           reference
     assembly_regions.tsv  curated local-reassembly regions
     exon_extent.bed       curated exon extents on the extended reference

2. The *built reference description*, ``<db>/reference.json``, written by
   ``script/build_reference.sh``. It records the gene list together with the
   coordinates that were actually produced, so no downstream step has to
   hardcode either gene names or intervals.

Runtime code should only ever use :func:`load_reference`.
"""

import json
import os

CLASS_I = "I"
CLASS_II = "II"

# G groups are defined by exon 2 for every gene plus exon 3 for class I only.
G_GROUP_EXONS = {CLASS_I: (2, 3), CLASS_II: (2,)}

# Number of leading exons recorded in exon_extent.bed, i.e. the exons that
# exon-only typing considers.
EXON_EXTENT_EXONS = {CLASS_I: 5, CLASS_II: 4}

REFERENCE_JSON = "reference.json"
GENE_LIST_FILE = os.path.join("HLA", "gene_list.txt")
GENE_REGIONS_FILE = os.path.join("HLA", "gene_regions.txt")


class Gene(object):
    """One gene of a built reference."""

    def __init__(self, gene, hla_class, representative_allele, start, end,
                 length):
        self.gene = gene
        self.name = "HLA_%s" % gene
        self.hla_class = hla_class
        self.representative_allele = representative_allele
        self.start = start
        self.end = end
        self.length = length

    @property
    def region(self):
        """Region string on the extended reference, e.g. ``HLA_A:1000-4503``."""
        return "%s:%s-%s" % (self.name, self.start, self.end)

    @property
    def g_group_exons(self):
        return G_GROUP_EXONS[self.hla_class]

    def as_dict(self):
        return {
            "gene": self.gene,
            "name": self.name,
            "hla_class": self.hla_class,
            "representative_allele": self.representative_allele,
            "start": self.start,
            "end": self.end,
            "length": self.length,
        }


class Reference(object):
    """The gene content and coordinates of one built reference."""

    def __init__(self, genes, db=None):
        self._genes = list(genes)
        self.db = db
        self._by_gene = dict((g.gene, g) for g in self._genes)
        self._by_name = dict((g.name, g) for g in self._genes)

    def __iter__(self):
        return iter(self._genes)

    def __len__(self):
        return len(self._genes)

    def __getitem__(self, key):
        """Look up by bare gene symbol (``A``) or reference name (``HLA_A``)."""
        if key in self._by_gene:
            return self._by_gene[key]
        return self._by_name[key]

    def __contains__(self, key):
        return key in self._by_gene or key in self._by_name

    @property
    def gene_list(self):
        """Bare gene symbols, e.g. ``['A', 'B', ...]``."""
        return [g.gene for g in self._genes]

    @property
    def gene_names(self):
        """Reference record names, e.g. ``['HLA_A', 'HLA_B', ...]``."""
        return [g.name for g in self._genes]

    @property
    def regions(self):
        """Region strings, e.g. ``['HLA_A:1000-4503', ...]``."""
        return [g.region for g in self._genes]

    def interval_dict(self):
        """``{'A': 'HLA_A:1000-4503', ...}`` keyed by bare gene symbol."""
        return dict((g.gene, g.region) for g in self._genes)

    def focus_region(self):
        """``{'HLA_A': [1000, 4503], ...}`` keyed by reference record name."""
        return dict((g.name, [g.start, g.end]) for g in self._genes)

    def result_header(self, prefix="Sample"):
        """Tab separated typing-result header for the genes in this reference."""
        columns = [prefix]
        for gene in self._genes:
            columns.append("%s_1" % gene.name)
            columns.append("%s_2" % gene.name)
        return "\t".join(columns)


def _rows(path):
    """Yield tab separated, comment stripped rows of a spec file."""
    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            yield line.split("\t")


def default_assets_dir():
    env = os.environ.get("SPECHLA_ASSETS", "")
    if env:
        return env
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(os.path.dirname(here), "share", "reference_assets")


def load_gene_spec(path=None, genes=None):
    """Read ``genes.tsv``.

    Returns a list of ``(gene, hla_class, representative_allele)`` tuples in
    file order, optionally restricted to (and ordered by the file order of)
    ``genes``.
    """
    if path is None:
        path = os.path.join(default_assets_dir(), "genes.tsv")
    wanted = set(genes) if genes else None
    spec = []
    seen = set()
    for fields in _rows(path):
        if len(fields) < 3:
            raise ValueError("%s: expected 3 columns, got %r" % (path, fields))
        gene, hla_class, allele = (f.strip() for f in fields[:3])
        if hla_class not in (CLASS_I, CLASS_II):
            raise ValueError(
                "%s: gene %s has unknown HLA class %r (expected I or II)"
                % (path, gene, hla_class)
            )
        if gene in seen:
            raise ValueError("%s: gene %s listed more than once" % (path, gene))
        seen.add(gene)
        if wanted is not None and gene not in wanted:
            continue
        spec.append((gene, hla_class, allele))
    if wanted is not None:
        missing = sorted(wanted - seen)
        if missing:
            raise ValueError(
                "%s: requested gene(s) not defined: %s" % (path, ", ".join(missing))
            )
    if not spec:
        raise ValueError("%s: no genes defined" % path)
    return spec


def load_assembly_regions(path=None, genes=None):
    """Read ``assembly_regions.tsv``.

    Returns ``(whole_regions, exon_regions)``, each a list of region strings
    for the genes in ``genes`` (all genes when ``genes`` is None).
    """
    if path is None:
        path = os.path.join(default_assets_dir(), "assembly_regions.tsv")
    wanted = set(genes) if genes else None
    whole, exon = [], []
    for fields in _rows(path):
        if len(fields) < 4:
            raise ValueError("%s: expected 4 columns, got %r" % (path, fields))
        gene, start, end, mode = (f.strip() for f in fields[:4])
        if mode not in ("whole", "exon", "both"):
            raise ValueError("%s: unknown mode %r" % (path, mode))
        if wanted is not None and gene not in wanted:
            continue
        region = "HLA_%s:%s-%s" % (gene, int(start), int(end))
        if mode in ("whole", "both"):
            whole.append(region)
        if mode in ("exon", "both"):
            exon.append(region)
    return whole, exon


def load_exon_extent(path=None, genes=None):
    """Read ``exon_extent.bed``, returning BED lines for the wanted genes."""
    if path is None:
        path = os.path.join(default_assets_dir(), "exon_extent.bed")
    wanted = set("HLA_%s" % g for g in genes) if genes else None
    lines = []
    for fields in _rows(path):
        if len(fields) < 3:
            raise ValueError("%s: expected 3 columns, got %r" % (path, fields))
        name = fields[0].strip()
        if wanted is not None and name not in wanted:
            continue
        lines.append("\t".join(f.strip() for f in fields[:3]))
    return lines


def load_reference(db=None):
    """Load ``<db>/reference.json``.

    ``db`` defaults to the database resolved by :mod:`spechla_paths`.
    """
    if db is None:
        from spechla_paths import get_db_dir
        db = get_db_dir()
    path = os.path.join(db, REFERENCE_JSON)
    if not os.path.isfile(path):
        raise RuntimeError(
            "%s not found. The reference at %s predates the centralized gene "
            "configuration; rebuild it with script/build_reference.sh."
            % (REFERENCE_JSON, db)
        )
    return _parse_reference(path, db)


_CACHE = {}


def reference(db=None):
    """Cached :func:`load_reference`, for use at the point of need."""
    key = db or ""
    if key not in _CACHE:
        _CACHE[key] = load_reference(db)
    return _CACHE[key]


def _parse_reference(path, db):
    with open(path) as handle:
        payload = json.load(handle)
    genes = [
        Gene(
            gene=entry["gene"],
            hla_class=entry["hla_class"],
            representative_allele=entry.get("representative_allele", ""),
            start=int(entry["start"]),
            end=int(entry["end"]),
            length=int(entry["length"]),
        )
        for entry in payload["genes"]
    ]
    return Reference(genes, db=db)
