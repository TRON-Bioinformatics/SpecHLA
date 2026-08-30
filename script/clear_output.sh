#!/bin/bash

# remove temporary files

outdir=$1

echo "Clean output dir."

rm -f  $outdir/HLA_*spechap.vcf.gz*
rm -f  $outdir/HLA_*.fasta
rm -f  $outdir/assembly*
rm -f  $outdir/*fa
rm -f  $outdir/sample*
rm -f  $outdir/*hap*.fasta.out
rm -f  $outdir/rematch.*
rm -f  $outdir/*uniq.name.R*.gz
rm -f  $outdir/header.sam
rm -f $outdir/newref_insertion*
rm -f $outdir/*map_database.bam*
rm -f $outdir/*merge.bam*
rm -f $outdir/*realign.bam*
rm -f $outdir/*realign.vcf.gz*
rm -f $outdir/fragment*file
rm -f  $outdir/*.R*.fq.gz
# Gene-specific leftovers; the gene list comes from the reference in use.
gene_list_file="${SPECHLA_DB:-}/HLA/gene_list.txt"
if [ -s "$gene_list_file" ]; then
    while read -r hla; do
        [ -n "$hla" ] || continue
        rm -f "$outdir/$hla.bam"*
        rm -f "$outdir/$hla.hla.count"
    done < "$gene_list_file"
fi
rm -f  $outdir/extract.read.blast
rm -rf  $outdir/tmp/







