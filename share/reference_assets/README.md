# SpecHLA reference assets

These files are bundled inputs for `script/build_reference.py`. The
version-specific HLA FASTA files are generated from a local IMGT/HLA clone;
these assets are either curated data or compatibility inputs retained from
the historical SpecHLA reference.

- `extend.fa` contains the 1 kb flanks used around representative alleles.
- `HLA_FREQ_HLA_I_II.txt` contains population frequencies from the historical
  SpecHLA database.
- `DRB1_dup_extract*.fasta` contains the curated DRB1 duplicated-region
  references used by phasing.
- `hla_gen.format.filter.extend.DRB.no26789*.fasta` are bundled read-binning
  references built from IPD-IMGT/HLA 3.51.0. With a newer IMGT release,
  novel alleles are covered only to gene-level accuracy by this step.
- `representative_alleles.tsv` selects the allele used for each gene-level
  extended reference and may be edited for a particular IMGT release.

The bundled read-binning references and curated duplication references are
not regenerated from IMGT by the initial build workflow.
