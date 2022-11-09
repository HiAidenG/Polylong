#!/bin/bash
set -euxo pipefail

#Download reference genome
mkdir data && cd data
cd data
#curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz | gunzip -d > gencode.v40.annotation.gtf

#take all entries that are protein coding and do not have a transcript ID
#these should be the coordinates for all protein coding genes
grep 'gene_type "protein_coding"' gencode.v40.annotation.gtf | grep -v 'transcript_id' > gencode.v40.annotation.pc.gtf

#remove any potential duplicate gene names
cat gencode.v40.annotation.pc.gtf | awk 'BEGIN { FS = ";" } ; !seen[$3]++' > gencode.v40.annotation.pc.unique.gtf

wc -l gencode.v40.annotation.pc.unique.gtf

#get sequences from primary assembly fasta for protein coding genes; flat format
#WARNING: This is inefficent but it works and I don't know a better way ATM.
#curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz | gunzip -d > GRCh38.primary_assembly.genome.fa
seqkit subseq --gtf gencode.v40.annotation.pc.unique.gtf GRCh38.primary_assembly.genome.fa -u 1024 -d 1024 > gencode.v40.annotation.pc.unique.fa

#should match above wc call for the .gtf file
grep -c ">" gencode.v40.annotation.pc.unique.fa

#temp
cat gencode.v40.annotation.pc.unique.fa | seqkit subseq -r 1025:1536 > posd512.fa

#double check above command worked
seqkit stats posd512.fa

#Now for the unstable sequences; get the bed files from the Core et al. paper, these are the coordinates for unstable transcripts (UU files)
#link to Core et al. paper https://www.nature.com/articles/ng.3142#MOESM77

#get Core et al. paper data
wget https://static-content.springer.com/esm/art%3A10.1038%2Fng.3142/MediaObjects/41588_2014_BFng3142_MOESM78_ESM.zip
unzip 41588_2014_BFng3142_MOESM78_ESM.zip

#get hg19 chrom size data for slopBed to extend sequences without going over chromosome ends
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

#get chain file for liftover to hg38
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O hg19ToHg38.over.chain.gz

slopBed -i tss_UU_gm12878_plus.bed -g hg19.chrom.sizes -b 1024 > UU_gm_plus_ext.bed
slopBed -i tss_UU_gm12878_minus.bed -g hg19.chrom.sizes -b 1024 > UU_gm_minus_ext.bed

#repeat above for k562 cell type
slopBed -i tss_UU_k562_minus.bed -g hg19.chrom.sizes -b 1024 > UU_k_minus_ext.bed
slopBed -i tss_UU_k562_plus.bed -g hg19.chrom.sizes -b 1024 > UU_k_plus_ext.bed

liftOver UU_gm_plus_ext.bed hg19ToHg38.over.chain.gz UU_gm_plus_ext_hg38.bed UU_gm_plus_ext_um.bed
liftOver UU_gm_minus_ext.bed hg19ToHg38.over.chain.gz UU_gm_minus_ext_hg38.bed UU_gm_minus_ext_um.bed
liftOver UU_k_minus_ext.bed hg19ToHg38.over.chain.gz UU_k_minus_ext_hg38.bed UU_k_minus_ext_um.bed
liftOver UU_k_plus_ext.bed hg19ToHg38.over.chain.gz UU_k_plus_ext_hg38.bed UU_k_plus_ext_um.bed

#merge all bed files
cat ./*ext_hg38.bed > UU_all_merged.bed

#get fasta files from bed files
bedtools getfasta -fi GRCh38.primary_assembly.genome.fa -bed UU_all_merged.bed -s -name -fo UU_all_merged.fa

#trim all fasta sequences to exactly 512 bases
cat UU_all_merged.fa | seqkit subseq -r 1025:1536 > negd512.fa

seqkit stats negd512.fa

echo "Running python script to preprocess data..."

python ../utils/preprocess.py

echo "Done!"
