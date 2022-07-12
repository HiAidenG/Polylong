#!/bin/bash

set -euxo
#This script gets the bed files from the Core et al. paper, gets the coordinates for unstable transcripts (UU files)\
#extends them to 500 bp, then lifts them over to hg38
#link to Core et al. paper https://www.nature.com/articles/ng.3142#MOESM77

#get Core et al. paper data
wget https://static-content.springer.com/esm/art%3A10.1038%2Fng.3142/MediaObjects/41588_2014_BFng3142_MOESM78_ESM.zip
unzip 41588_2014_BFng3142_MOESM78_ESM.zip

#get hg19 chrom size data for slopBed to extend sequences without going over chromosome ends
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

#get chain file for liftover to hg38
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O hg19ToHg38.over.chain.gz

#extend sequences 700 bases; do for both strands
slopBed -i tss_UU_gm12878_plus.bed -g hg19.chrom.sizes -l 0 -r 700 -s > UU_gm_plus_ext.bed
slopBed -i tss_UU_gm12878_minus.bed -g hg19.chrom.sizes -l 0 -r 700 -s > UU_gm_minus_ext.bed

#repeat above for k562 cell type
slopBed -i tss_UU_k562_minus.bed -g hg19.chrom.sizes -l 0 -r 700 -s > UU_k_minus_ext.bed
slopBed -i tss_UU_k562_plus.bed -g hg19.chrom.sizes -l 0 -r 700 -s > UU_k_plus_ext.bed

liftOver UU_gm_plus_ext.bed hg19ToHg38.over.chain.gz UU_gm_plus_ext_hg38.bed UU_gm_plus_ext_um.bed
liftOver UU_gm_minus_ext.bed hg19ToHg38.over.chain.gz UU_gm_minus_ext_hg38.bed UU_gm_minus_ext_um.bed
liftOver UU_k_minus_ext.bed hg19ToHg38.over.chain.gz UU_k_minus_ext_hg38.bed UU_k_minus_ext_um.bed
liftOver UU_k_plus_ext.bed hg19ToHg38.over.chain.gz UU_k_plus_ext_hg38.bed UU_k_plus_ext_um.bed

#merge all bed files
cat ./*ext_hg38.bed > UU_all_merged.bed

#get fasta files from bed files
bedtools getfasta -fi GRCh38.primary_assembly.genome.fa -bed UU_all_merged.bed -s -name -fo UU_all_merged.fa

#trim all fasta sequences to exactly 500 bases
cat UU_all_merged.fa | seqkit subseq -r 1:700 > UU_all_merged_trim.fa

