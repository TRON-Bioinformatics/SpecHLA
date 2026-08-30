#!/bin/bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash script/build_reference.sh --imgt <IMGTHLA_dir> --out <build_dir>

Required:
  --imgt PATH                 Read-only IMGT/HLA clone
  --out PATH                  Output directory

Optional:
  --genes-config FILE         Gene definition table (default: bundled
                              share/reference_assets/genes.tsv). This file is
                              the single source of truth for which genes the
                              reference covers, their HLA class and their
                              representative allele.
  --genes LIST                Comma separated subset of the genes in
                              --genes-config to build (default: all)
  --threads N
  --skip-novoindex
  --force
EOF
}

IMGT="${SPECHLA_IMGT:-}"
OUT=""
GENES_CONFIG=""
GENES=""
THREADS=4
SKIP_NOVOINDEX=0
FORCE=0
while [ "$#" -gt 0 ]; do
    case "$1" in
        --imgt) IMGT="${2:-}"; shift 2 ;;
        --out) OUT="${2:-}"; shift 2 ;;
        --genes-config|--representative-alleles) GENES_CONFIG="${2:-}"; shift 2 ;;
        --genes) GENES="${2:-}"; shift 2 ;;
        --threads) THREADS="${2:-}"; shift 2 ;;
        --skip-novoindex) SKIP_NOVOINDEX=1; shift ;;
        --force) FORCE=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [ -z "$IMGT" ] || [ -z "$OUT" ]; then
    echo "--imgt and --out are required" >&2
    usage >&2
    exit 2
fi
if [ -e "$OUT/ref/hla.ref.extend.fa.bwt" ] && [ "$FORCE" = 0 ]; then
    echo "Reference already exists at $OUT; use --force to rebuild" >&2
    exit 2
fi
if [ ! -d "$IMGT" ]; then
    echo "IMGT/HLA directory does not exist: $IMGT" >&2
    exit 2
fi
if [ ! -f "$IMGT/hla_nuc.fasta" ] || [ ! -f "$IMGT/Allelelist.txt" ] \
    || [ ! -f "$IMGT/wmda/hla_nom_g.txt" ] \
    || { [ ! -f "$IMGT/xml/hla.xml" ] && [ ! -f "$IMGT/xml/hla.xml.zip" ]; }; then
    echo "IMGT/HLA clone is missing required files" >&2
    exit 2
fi

ROOT=$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/.." && pwd)
PYTHON="${PYTHON:-python}"
ASSETS="${SPECHLA_ASSETS:-$ROOT/share/reference_assets}"
GENES_CONFIG="${GENES_CONFIG:-$ASSETS/genes.tsv}"
SPEC_ARGS=(--genes-config "$GENES_CONFIG")
if [ -n "$GENES" ]; then
    SPEC_ARGS+=(--genes "$GENES")
fi
GEN="$IMGT/hla_gen.fasta"
if [ ! -f "$GEN" ] && [ -f "$IMGT/hla_gen.fasta.zip" ]; then
    mkdir -p "$OUT/.build_tmp"
    unzip -oq "$IMGT/hla_gen.fasta.zip" -d "$OUT/.build_tmp"
    GEN="$OUT/.build_tmp/hla_gen.fasta"
fi
if [ ! -f "$GEN" ]; then
    echo "IMGT/HLA clone is missing hla_gen.fasta or hla_gen.fasta.zip" >&2
    exit 2
fi

mkdir -p "$OUT/HLA/whole" "$OUT/HLA/exon"
mkdir -p "$OUT/.build_tmp"
if [ "$FORCE" = 1 ]; then
    rm -f "$OUT/HLA/hla.ref.extend.fa" "$OUT/ref/hla.ref.extend.fa"
fi
cp "$IMGT/Allelelist.txt" "$OUT/HLA/Allelelist.txt"
cp "$IMGT/wmda/hla_nom_g.txt" "$OUT/HLA/hla_nom_g.txt"

# Resolve the gene spec once and materialise everything that is derived from
# it inside the reference, so that no downstream step has to know the gene
# list, the curated regions or their coordinates.
"$PYTHON" "$ROOT/script/build_reference/emit_gene_files.py" \
    "${SPEC_ARGS[@]}" --assets "$ASSETS" --out "$OUT/HLA"
mapfile -t HLA_GENES < "$OUT/HLA/gene_list.txt"
if [ "${#HLA_GENES[@]}" -eq 0 ]; then
    echo "No genes resolved from $GENES_CONFIG" >&2
    exit 2
fi
echo "[build_reference] building genes: ${HLA_GENES[*]}"

"$PYTHON" "$ROOT/script/build_reference/split_by_gene.py" \
    "${SPEC_ARGS[@]}" \
    --gen "$GEN" --nuc "$IMGT/hla_nuc.fasta" --out "$OUT/HLA"

if [ -f "$OUT/HLA/exon/HLA_DRB1.fasta" ]; then
    cp "$OUT/HLA/exon/HLA_DRB1.fasta" "$OUT/HLA/whole/HLA_DRB1.exon.fasta"
fi

if command -v samtools >/dev/null 2>&1; then
    for fasta in "$OUT"/HLA/whole/*.fasta "$OUT"/HLA/exon/*.fasta; do
        samtools faidx "$fasta"
    done
else
    echo "samtools is required to index the generated FASTA files" >&2
    exit 127
fi
if command -v makeblastdb >/dev/null 2>&1; then
    for fasta in "$OUT"/HLA/whole/*.fasta "$OUT"/HLA/exon/*.fasta; do
        if [[ "$fasta" == */exon/* ]]; then
            blast_prefix="$fasta"
        else
            blast_prefix="${fasta%.fasta}"
        fi
        makeblastdb -in "$fasta" -dbtype nucl -parse_seqids -out "$blast_prefix"
    done
else
    echo "makeblastdb is required to index the generated FASTA files" >&2
    exit 127
fi

mkdir -p "$OUT"/{ref,ref/refdata-hla.ref.extend}
cp "$ASSETS/extend.fa" "$OUT/HLA/extend.fa"
cp "$ASSETS/HLA_FREQ_HLA_I_II.txt" "$OUT/HLA/"
cp "$ASSETS/DRB1_dup_extract.fasta" "$OUT/ref/"
cp "$ASSETS/DRB1_dup_extract_ref.fasta" "$OUT/ref/"
cp "$ASSETS/hla_gen.format.filter.extend.DRB.no26789.fasta" "$OUT/ref/"
cp "$ASSETS/hla_gen.format.filter.extend.DRB.no26789.v2.fasta" "$OUT/ref/"

# Keep the inputs that produced this reference inside the reference itself, so
# a built database is self describing and can be adjusted and rebuilt without
# the source tree.
mkdir -p "$OUT/spec"
cp "$GENES_CONFIG" "$OUT/spec/genes.tsv"
cp "$ASSETS/assembly_regions.tsv" "$OUT/spec/assembly_regions.tsv"
cp "$ASSETS/exon_extent.bed" "$OUT/spec/exon_extent.bed"

XML="$IMGT/xml/hla.xml"
if [ ! -f "$XML" ]; then
    mkdir -p "$OUT/.build_tmp"
    unzip -oq "$IMGT/xml/hla.xml.zip" -d "$OUT/.build_tmp"
    XML="$OUT/.build_tmp/hla.xml"
fi
"$PYTHON" "$ROOT/script/build_reference/extract_exons_from_xml.py" \
    "${SPEC_ARGS[@]}" \
    --xml "$XML" --out "$OUT/HLA/hla_exons.fasta"
samtools faidx "$OUT/HLA/hla_exons.fasta"
makeblastdb -in "$OUT/HLA/hla_exons.fasta" -dbtype nucl -parse_seqids \
    -out "$OUT/HLA/hla_exons.fasta"

"$PYTHON" "$ROOT/script/build_reference/construct_extended.py" \
    "${SPEC_ARGS[@]}" \
    --gen "$GEN" --extend "$ASSETS/extend.fa" \
    --out "$OUT/HLA/hla.ref.extend.fa" \
    --reference-json "$OUT/reference.json" \
    --selected-out "$OUT/.build_tmp/selected_representatives.json"
samtools faidx "$OUT/HLA/hla.ref.extend.fa"
# reference.json now carries the real per-gene coordinates; publish them for
# the shell steps as well.
"$PYTHON" "$ROOT/script/build_reference/emit_gene_files.py" \
    --reference-json "$OUT/reference.json" --out "$OUT/HLA"
cp "$OUT/HLA/hla.ref.extend.fa" "$OUT/ref/hla.ref.extend.fa"
cp "$OUT/HLA/hla.ref.extend.fa.fai" "$OUT/ref/hla.ref.extend.fa.fai"
bwa index "$OUT/ref/hla.ref.extend.fa"
samtools dict "$OUT/ref/hla.ref.extend.fa" > "$OUT/ref/hla.ref.extend.dict"

for hla in "${HLA_GENES[@]}"; do
    gene_dir="$OUT/HLA/HLA_${hla}"
    mkdir -p "$gene_dir"
    samtools faidx "$OUT/HLA/hla.ref.extend.fa" "HLA_${hla}" > "$gene_dir/HLA_${hla}.fa"
    samtools faidx "$gene_dir/HLA_${hla}.fa"
    bwa index "$gene_dir/HLA_${hla}.fa"
    makeblastdb -in "$gene_dir/HLA_${hla}.fa" -dbtype nucl -parse_seqids \
        -out "$gene_dir/HLA_${hla}"
    if command -v faToTwoBit >/dev/null 2>&1; then
        faToTwoBit "$gene_dir/HLA_${hla}.fa" "$gene_dir/HLA_${hla}.2bit"
    else
        echo "[build_reference] faToTwoBit not on PATH; skipping ${hla}.2bit" >&2
    fi
    cfg="$OUT/HLA/HLA_${hla}.config.txt"
    printf 'bwa=%s\nfreebayes=%s\nblat=%s/\n' \
        "$gene_dir/HLA_${hla}.fa" "$gene_dir/HLA_${hla}.fa" "$gene_dir" > "$cfg"
done

# The per-gene reference is also the one phasing loads by gene name.
for hla in "${HLA_GENES[@]}"; do
    cp "$OUT/HLA/HLA_${hla}/HLA_${hla}.fa" "$OUT/ref/HLA_${hla}.fa"
    cp "$OUT/HLA/HLA_${hla}/HLA_${hla}.fa.fai" "$OUT/ref/HLA_${hla}.fa.fai"
done

samtools faidx "$OUT/ref/DRB1_dup_extract.fasta"
makeblastdb -in "$OUT/ref/DRB1_dup_extract_ref.fasta" -dbtype nucl -parse_seqids \
    -out "$OUT/ref/DRB1_dup_extract_ref.fasta"
for fasta in "$OUT"/ref/hla_gen.format.filter.extend.DRB.no26789*.fasta; do
    samtools faidx "$fasta"
    if [ "$SKIP_NOVOINDEX" = 0 ] && command -v novoalign >/dev/null 2>&1 \
        && [ -f "$(dirname "$(command -v novoalign)")/novoalign.lic" ]; then
        novoindex -k 14 -s 1 "${fasta%.fasta}.ndx" "$fasta"
    fi
    bowtie2-build --threads "$THREADS" "$fasta" "$fasta"
done

if command -v longranger >/dev/null 2>&1; then
    (cd "$OUT/ref" && longranger mkref hla.ref.extend.fa)
else
    echo "[build_reference] longranger not on PATH; 10X data will require manual mkref" >&2
fi

IMGT_VERSION=$(awk -F': ' '/^# version:/{print $2; exit}' "$IMGT/Allelelist.txt")
SPECHLA_REF=$(git -C "$ROOT" rev-parse HEAD 2>/dev/null || printf 'unknown')
"$PYTHON" - "$OUT/BUILD_MANIFEST.json" "$IMGT_VERSION" "$IMGT" "$SPECHLA_REF" \
    "$OUT/.build_tmp/selected_representatives.json" "$ASSETS" "$OUT/reference.json" <<'PY'
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

manifest, version, source, ref, selected_path, assets, reference_path = sys.argv[1:]
with open(selected_path) as handle:
    selected = json.load(handle)
with open(reference_path) as handle:
    reference = json.load(handle)
asset_hashes = {}
for name in os.listdir(assets):
    path = os.path.join(assets, name)
    if os.path.isfile(path):
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        asset_hashes[name] = digest.hexdigest()
payload = {
    "imgt_version": version,
    "imgt_source": os.path.abspath(source),
    "built_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
    "spechla_ref": ref,
    "representative_alleles": selected,
    "genes": [entry["gene"] for entry in reference["genes"]],
    "gene_regions": {
        entry["name"]: [entry["start"], entry["end"]] for entry in reference["genes"]
    },
    "bundled_assets_sha256": asset_hashes,
}
with open(manifest, "w") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY

rm -rf "$OUT/.build_tmp"
