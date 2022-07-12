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
#WARNING: This is inefficent but it works and I don't know a better way ATM.
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz | gunzip -d > ../data/GRCh38.primary_assembly.genome.fa
seqkit subseq --gtf ../data/gencode.v40.annotation.pc.unique.gtf ../data/GRCh38.primary_assembly.genome.fa -u 200 -d 500 > ../data/gencode.v40.annotation.pc.unique.fa

#should match above wc call for the .gtf file
grep -c ">" ../data/gencode.v40.annotation.pc.unique.fa

#Get 500bp downstream and 200 bp upstream of TSS
cat ../data/gencode.v40.annotation.pc.unique.fa | seqkit subseq -r 1:700 > ../data/raw_data.fa

#double check above command worked
seqkit stats ../data/raw_data.fa
