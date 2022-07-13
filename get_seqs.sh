#!/bin/bash
set -euxo pipefail
#This script collects all protein coding sequences and unstable sequences -200bp and +312bp from TSS.

#Download reference genome
mkdir data
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz | gunzip -d > data/gencode.v40.annotation.gtf

#take all entries that are protein coding and do not have a transcript ID
#these should be the coordinates for all protein coding genes
grep 'gene_type "protein_coding"' data/gencode.v40.annotation.gtf | grep -v 'transcript_id' > data/gencode.v40.annotation.pc.gtf

#remove any potential duplicate gene names
cat data/gencode.v40.annotation.pc.gtf | awk 'BEGIN { FS = ";" } ; !seen[$3]++' > data/gencode.v40.annotation.pc.unique.gtf

wc -l data/gencode.v40.annotation.pc.unique.gtf

#get sequences from primary assembly fasta for protein coding genes; flat format
#WARNING: This is inefficent but it works and I don't know a better way ATM.
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz | gunzip -d > data/GRCh38.primary_assembly.genome.fa
seqkit subseq --gtf data/gencode.v40.annotation.pc.unique.gtf data/GRCh38.primary_assembly.genome.fa -u 200 -d 312 > data/gencode.v40.annotation.pc.unique.fa

#should match above wc call for the .gtf file
grep -c ">" data/gencode.v40.annotation.pc.unique.fa

#Get 312bp downstream and 200 bp upstream of TSS
cat data/gencode.v40.annotation.pc.unique.fa | seqkit subseq -r 1:512 > PC_all_merged_trim.fa

#double check above command worked
seqkit stats PC_all_merged_trim.fa

#Now for the unstable sequences; get the bed files from the Core et al. paper, these are the coordinates for unstable transcripts (UU files)
#link to Core et al. paper https://www.nature.com/articles/ng.3142#MOESM77

#get Core et al. paper data
cd data
wget https://static-content.springer.com/esm/art%3A10.1038%2Fng.3142/MediaObjects/41588_2014_BFng3142_MOESM78_ESM.zip
unzip 41588_2014_BFng3142_MOESM78_ESM.zip

#get hg19 chrom size data for slopBed to extend sequences without going over chromosome ends
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

#get chain file for liftover to hg38
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O hg19ToHg38.over.chain.gz

#extend sequences to -200bp/+312bp; do for both strands
slopBed -i tss_UU_gm12878_plus.bed -g hg19.chrom.sizes -l 0 -r 512 -s > UU_gm_plus_ext.bed
slopBed -i tss_UU_gm12878_minus.bed -g hg19.chrom.sizes -l 0 -r 512 -s > UU_gm_minus_ext.bed

#repeat above for k562 cell type
slopBed -i tss_UU_k562_minus.bed -g hg19.chrom.sizes -l 0 -r 512 -s > UU_k_minus_ext.bed
slopBed -i tss_UU_k562_plus.bed -g hg19.chrom.sizes -l 0 -r 512 -s > UU_k_plus_ext.bed

liftOver UU_gm_plus_ext.bed hg19ToHg38.over.chain.gz UU_gm_plus_ext_hg38.bed UU_gm_plus_ext_um.bed
liftOver UU_gm_minus_ext.bed hg19ToHg38.over.chain.gz UU_gm_minus_ext_hg38.bed UU_gm_minus_ext_um.bed
liftOver UU_k_minus_ext.bed hg19ToHg38.over.chain.gz UU_k_minus_ext_hg38.bed UU_k_minus_ext_um.bed
liftOver UU_k_plus_ext.bed hg19ToHg38.over.chain.gz UU_k_plus_ext_hg38.bed UU_k_plus_ext_um.bed

#merge all bed files
cat ./*ext_hg38.bed > UU_all_merged.bed

#get fasta files from bed files
bedtools getfasta -fi GRCh38.primary_assembly.genome.fa -bed UU_all_merged.bed -s -name -fo UU_all_merged.fa

#trim all fasta sequences to exactly 512 bases
cat UU_all_merged.fa | seqkit subseq -r 1:512 > ../UU_all_merged_trim.fa

cd ..
seqkit stats UU_all_merged_trim.fa

echo "Done!"
echo "Now run preprocess.py to get final training/test data."

