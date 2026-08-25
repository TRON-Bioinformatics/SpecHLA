#!/bin/bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash script/build_reference.sh --imgt <IMGTHLA_dir> --out <build_dir>

Required:
  --imgt PATH                 Read-only IMGT/HLA clone
  --out PATH                  Output directory
EOF
}

IMGT="${SPECHLA_IMGT:-}"
OUT=""
while [ "$#" -gt 0 ]; do
    case "$1" in
        --imgt) IMGT="${2:-}"; shift 2 ;;
        --out) OUT="${2:-}"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [ -z "$IMGT" ] || [ -z "$OUT" ]; then
    echo "--imgt and --out are required" >&2
    usage >&2
    exit 2
fi
if [ ! -d "$IMGT" ]; then
    echo "IMGT/HLA directory does not exist: $IMGT" >&2
    exit 2
fi
if [ ! -f "$IMGT/hla_nuc.fasta" ] || [ ! -f "$IMGT/Allelelist.txt" ] \
    || [ ! -f "$IMGT/wmda/hla_nom_g.txt" ]; then
    echo "IMGT/HLA clone is missing required files" >&2
    exit 2
fi

ROOT=$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/.." && pwd)
PYTHON="${PYTHON:-python}"
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
cp "$IMGT/Allelelist.txt" "$OUT/HLA/Allelelist.txt"
cp "$IMGT/wmda/hla_nom_g.txt" "$OUT/HLA/hla_nom_g.txt"
"$PYTHON" "$ROOT/script/build_reference/split_by_gene.py" \
    --gen "$GEN" --nuc "$IMGT/hla_nuc.fasta" --out "$OUT/HLA"

cp "$OUT/HLA/exon/HLA_DRB1.fasta" "$OUT/HLA/whole/HLA_DRB1.exon.fasta"

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
        makeblastdb -in "$fasta" -dbtype nucl -parse_seqids -out "${fasta%.fasta}"
    done
else
    echo "makeblastdb is required to index the generated FASTA files" >&2
    exit 127
fi
