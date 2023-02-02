#!/bin/bash
set -euxo pipefail

#Download reference genome
mkdir data && cd ./data
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz | gunzip -d > gencode.v40.annotation.gtf

#Download primary assembly fasta
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz | gunzip -d > GRCh38.primary_assembly.genome.fa

#take all entries that are protein coding and do not have a transcript ID
#these should be the coordinates for all protein coding genes
grep 'gene_type "protein_coding"' gencode.v40.annotation.gtf | grep -v 'transcript_id' | awk 'BEGIN { FS = ";" } ; !seen[$3]++' > gencode.v40.annotation.pc.gtf

#check that the above command worked
wc -l gencode.v40.annotation.pc.gtf

#TODO: This is inefficent but it works and I don't know a better way ATM.
seqkit subseq --gtf gencode.v40.annotation.pc.gtf GRCh38.primary_assembly.genome.fa -u 512 -d 512 > pos_TSSflanking_1024.fa

#should match above wc call for the .gtf file
grep -c ">" pos_TSSflanking_1024.fa

#double check above command worked
seqkit stats pos_TSSflanking_1024.fa

# clean up:
rm -f gencode.v40.annotation.pc.gtf
rm -f gencode.v40.annotation.gtf

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

#Again, change this to get desired subsequence for negatives. 0 is 1024 bp upstream of TSS, 2048 is 1024 bp downstream of TSS.
cat UU_all_merged.fa | seqkit subseq -r 512:1024 > neg_TSSflanking_1024.fa

seqkit stats neg_TSSflanking_1024.fa

#clean up, TOOD: should probably just pipe everything instead of having to do this
rm -f UU_all_merged.bed
rm -f UU_all_merged.fa
rm -f UU_gm_minus_ext.bed
rm -f UU_gm_minus_ext_hg38.bed
rm -f UU_gm_minus_ext_um.bed
rm -f UU_gm_plus_ext.bed
rm -f UU_gm_plus_ext_hg38.bed
rm -f UU_gm_plus_ext_um.bed
rm -f UU_k_minus_ext.bed
rm -f UU_k_minus_ext_hg38.bed
rm -f UU_k_minus_ext_um.bed
rm -f UU_k_plus_ext.bed
rm -f UU_k_plus_ext_hg38.bed
rm -f UU_k_plus_ext_um.bed
rm -f tss_UU_gm12878_minus.bed
rm -f tss_UU_gm12878_plus.bed
rm -f tss_UU_k562_minus.bed
rm -f tss_UU_k562_plus.bed
rm -f hg19.chrom.sizes
rm -f hg19ToHg38.over.chain.gz


echo "Running python script to preprocess data..."

# python ../../utils/preprocessdata.py

echo "Done!"
