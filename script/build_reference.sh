#!/bin/bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash script/build_reference.sh --imgt <IMGTHLA_dir> --out <build_dir>

Required:
  --imgt PATH                 Read-only IMGT/HLA clone
  --out PATH                  Output directory
  --representative-alleles FILE
  --threads N
  --skip-novoindex
  --force
EOF
}

IMGT="${SPECHLA_IMGT:-}"
OUT=""
REPRESENTATIVES=""
THREADS=4
SKIP_NOVOINDEX=0
FORCE=0
while [ "$#" -gt 0 ]; do
    case "$1" in
        --imgt) IMGT="${2:-}"; shift 2 ;;
        --out) OUT="${2:-}"; shift 2 ;;
        --representative-alleles) REPRESENTATIVES="${2:-}"; shift 2 ;;
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
REPRESENTATIVES="${REPRESENTATIVES:-$ASSETS/representative_alleles.tsv}"
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
if [ "$FORCE" = 1 ]; then
    rm -f "$OUT/HLA/hla.ref.extend.fa" "$OUT/ref/hla.ref.extend.fa"
fi
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

mkdir -p "$OUT"/{ref,ref/refdata-hla.ref.extend}
cp "$ASSETS/extend.fa" "$OUT/HLA/extend.fa"
cp "$ASSETS/HLA_FREQ_HLA_I_II.txt" "$OUT/HLA/"
cp "$ASSETS/DRB1_dup_extract.fasta" "$OUT/ref/"
cp "$ASSETS/DRB1_dup_extract_ref.fasta" "$OUT/ref/"
cp "$ASSETS/hla_gen.format.filter.extend.DRB.no26789.fasta" "$OUT/ref/"
cp "$ASSETS/hla_gen.format.filter.extend.DRB.no26789.v2.fasta" "$OUT/ref/"

XML="$IMGT/xml/hla.xml"
if [ ! -f "$XML" ]; then
    mkdir -p "$OUT/.build_tmp"
    unzip -oq "$IMGT/xml/hla.xml.zip" -d "$OUT/.build_tmp"
    XML="$OUT/.build_tmp/hla.xml"
fi
"$PYTHON" "$ROOT/script/build_reference/extract_exons_from_xml.py" \
    --xml "$XML" --out "$OUT/HLA/hla_exons.fasta"
samtools faidx "$OUT/HLA/hla_exons.fasta"
makeblastdb -in "$OUT/HLA/hla_exons.fasta" -dbtype nucl -parse_seqids \
    -out "$OUT/HLA/hla_exons.fasta"

"$PYTHON" "$ROOT/script/build_reference/construct_extended.py" \
    --gen "$GEN" --extend "$ASSETS/extend.fa" \
    --representatives "$REPRESENTATIVES" \
    --out "$OUT/HLA/hla.ref.extend.fa"
samtools faidx "$OUT/HLA/hla.ref.extend.fa"
cp "$OUT/HLA/hla.ref.extend.fa" "$OUT/ref/hla.ref.extend.fa"
cp "$OUT/HLA/hla.ref.extend.fa.fai" "$OUT/ref/hla.ref.extend.fa.fai"
bwa index "$OUT/ref/hla.ref.extend.fa"
samtools dict "$OUT/ref/hla.ref.extend.fa" > "$OUT/ref/hla.ref.extend.dict"

for hla in A B C DPA1 DPB1 DQA1 DQB1 DRB1; do
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
