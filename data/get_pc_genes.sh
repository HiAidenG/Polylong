#!/bin/bash
set -euxo pipefail
#This script aims to get the DNA sequences for all protein coding genes 200 bp upstream and 500 bp downstream of TSS.
mkdir ../data/
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz | gunzip -d > ../data/gencode.v40.annotation.gtf

#take all entries that are protein coding and do not have a transcript ID
#these should be the coordinates for all protein coding genes
grep 'gene_type "protein_coding"' ../data/gencode.v40.annotation.gtf | grep -v 'transcript_id' > ../data/gencode.v40.annotation.pc.gtf

#remove any potential duplicate gene names
cat ../data/gencode.v40.annotation.pc.gtf | awk 'BEGIN { FS = ";" } ; !seen[$3]++' > ../data/gencode.v40.annotation.pc.unique.gtf

wc -l ../data/gencode.v40.annotation.pc.unique.gtf

#get sequences from primary assembly fasta for protein coding genes; flat format
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz | gunzip -d > ../data/GRCh38.primary_assembly.genome.fa
bedtools getfasta -fi ../data/GRCh38.primary_assembly.genome.fa -bed ../data/gencode.v40.annotation.pc.unique.gtf -s -name \
-fo ../data/gencode.v40.annotation.pc.unique.fa

#should match above wc call for the .gtf file
grep -c ">" ../data/gencode.v40.annotation.pc.unique.fa

#trim each sequence to first 500 bases or full sequence if smaller than 500 bases and get 200 bp upstream.
cat ../data/gencode.v40.annotation.pc.unique.fa | seqkit subseq -u 200 -r 1:500 > ../data/gencode.v40.annotation.pc.unique.trim.fa

#double check above command worked
seqkit stats ../data/gencode.v40.annotation.pc.unique.trim.fa

#determine line number of 80th percent sequence
#split fasta at line 159249 for training (80%) and predicting (20%)
#rename files
# SL=$(grep -n ">" ../data/gencode.v38.annotation.pc.unique.fa | sed '15944q;d' | cut -f 1 -d:)
# csplit -f positive. ../data/gencode.v38.annotation.pc.unique.fa "$SL"
# mv ../data/positive.00 ../data/positive_pc_train_80percent.fa
# mv ../data/positive.01 ../data/positive_pc_test_20percent.fa

