# SpecHLA reference assets

These files are the bundled inputs for `script/build_reference.sh`. The
version-specific HLA FASTA files are generated from a local IMGT/HLA clone;
these assets are either curated data or compatibility inputs retained from
the historical SpecHLA reference.

## Gene configuration

- `genes.tsv` is the single source of truth for the reference: which genes are
  covered, their HLA class, and the representative allele used for the
  extended per-gene reference. Adding, removing or reordering a gene here is
  the only change needed to change the gene content of a reference; no script
  contains a gene list.
- `assembly_regions.tsv` lists the curated regions that are re-assembled
  locally before variant calling, tagged for full-length typing, exon typing
  or both. The build writes the matching, gene-filtered subsets to
  `<db>/HLA/select.region.txt` and `<db>/HLA/select.region.exon.txt`.
- `exon_extent.bed` holds the curated exon extents on the extended reference.
  The build writes the gene-filtered subset to `<db>/HLA/exon_extent.bed`.

Coordinates in `assembly_regions.tsv` and `exon_extent.bed` refer to the
extended reference and therefore depend on the representative alleles in
`genes.tsv`; changing a representative allele means these regions should be
re-curated.

## Sequence assets

- `extend.fa` contains the 1 kb flanks used around representative alleles. It
  must provide `HLA_<gene>_1` and `HLA_<gene>_2` for every gene in
  `genes.tsv`.
- `HLA_FREQ_HLA_I_II.txt` contains population frequencies from the historical
  SpecHLA database.
- `DRB1_dup_extract*.fasta` contains the curated DRB1 duplicated-region
  references used by phasing.
- `hla_gen.format.filter.extend.DRB.no26789*.fasta` are bundled read-binning
  references built from IPD-IMGT/HLA 3.51.0. With a newer IMGT release,
  novel alleles are covered only to gene-level accuracy by this step.

The bundled read-binning references and curated duplication references are
not regenerated from IMGT by the initial build workflow.

## What ends up in a built reference

`script/build_reference.sh` makes the build self describing:

- `<db>/reference.json` records every gene together with its HLA class,
  the representative allele that was actually used, the flank lengths and the
  interval of the gene body on `HLA_<gene>`. All runtime code reads gene names
  and coordinates from here via `script/reference_config.py`.
- `<db>/HLA/gene_list.txt` and `<db>/HLA/gene_regions.txt` expose the same
  information to the shell steps.
- `<db>/spec/` keeps a copy of `genes.tsv`, `assembly_regions.tsv` and
  `exon_extent.bed`, so a built database can be inspected, adjusted and
  rebuilt without the source tree.
- `<db>/BUILD_MANIFEST.json` records the IMGT version, the resolved gene list
  and regions, and SHA-256 hashes of every bundled asset used.
