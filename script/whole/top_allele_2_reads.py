import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from spechla_paths import get_db_dir
from reference_config import reference
from Bio import SeqIO

outdir = sys.argv[1]
detail_anno = outdir + "/hla.result.details.txt"
top_allele_fasta = outdir + "/top_allele.fasta"
top_allele_fastq = outdir + "/top_allele.fastq"
samtools = "samtools"
db_dir = get_db_dir() + "/HLA/whole/"
f=open(top_allele_fasta, 'w')
f.close()

def read():
    dict = {gene: '' for gene in reference().gene_list}
    top_allele = []
    f = open(detail_anno, 'r')
    for line in f:
        array = line.strip().split()
        if array[2] == "allele":
            continue
        allele_list = array[2].split(";")
        for allele in allele_list:
            if allele == '':
                continue
            allele_array = allele.split("*")
            gene = allele_array[0]
            dict[gene] += " " + allele
            top_allele.append(allele)
    # print (dict)
    for gene in dict.keys():
        fasta = db_dir + "HLA_" + gene + ".fasta"
        command = samtools + " faidx " + fasta + dict[gene] + ">>" + top_allele_fasta
        os.system(command)
        # print (command)

def fasta2fastq():
    f = open(top_allele_fastq, 'w')
    for r in SeqIO.parse(top_allele_fasta, "fasta"):
        r.letter_annotations["solexa_quality"] = [40] * len(r)
        print(r.format("fastq"), end='', file = f)
    f.close()

read()
fasta2fastq()
