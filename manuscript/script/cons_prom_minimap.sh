#!/bin/bash

# Jan 2022
# Script to défine conserved promoters bw dog and human
# Map onto human genome GRCh38
#########################################################

# activate minimap bedtools samtools
source /local/env/envconda.sh
source /local/env/envminimap2-2.15.sh
source /local/env/envsamtools-1.6.sh
source /local/env/envbedtools-2.25.0.sh

# Dog sequence promoters
query_fa="./human_dog_8x1024_dog.fa"
# Genome human
target_fa=~tderrien/DATA/hg38_GRCh38/sequence/softmasked/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa


# Map with minmap2 
# following heng li and CK post
# ./minimap2 -ax map-ont ref.fa ont.fq.gz > aln.sam         # Oxford Nanopore genomic reads
# https://github.com/lh3/minimap2/issues/238
############################################

output_pref=$(basename $query_fa .fasta)
minimap2 -ax splice -t16 $target_fa $query_fa >  ${output_pref}.sam
#real	1m8.147s
#user	2m51.029s
#sys	0m58.472s


# Get uniq match
# with blast defined seq identity 
# -F 260 = read unmapped & not primary alignment criteria 3 & 9 are selected for exclusion
#                       bit 3 + bit 9 = 4 + 256 = 260
# http://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
samtools view -F 260 ${output_pref}.sam  | awk '{seen[$1]++; if (seen[$1] == 1){print $0}}' | perl -ane 'if(/NM:i:(\d+)/){$n=$1;$l=0;$l+=$1 while/(\d+)[MID]/g;print(join( "\t", @F ),"\t",($l-$n)/$l,"\n")}'  > ${output_pref}.primAlign.blastId.sam

# COnversion in .bed12
paftools.js splice2bed ${output_pref}.primAlign.blastId.sam > ${output_pref}.primAlign.blastId.bed12


# Overlap with human gene promoter 
# see README
target_bed=./human_dog_8x1024_human_1191_Biomart_prom.bed 


# Check if same genes in tageted regions
# if $NF==0 	=> not in target (no overlap)
# $16==ar[1]  	=> in target (same genes from the dog sequence to the human promoter)
more human_dog_8x1024_dog.fa.primAlign.blastId.bed12 | awk '{print "chr"$0}' | \
bedtools intersect -a - -b human_dog_8x1024_human_1191_Biomart_prom.bed -wao  | \
awk '{split($4,ar,"_"); if ($NF==0){print $0"\t""notInTarget"};if ($16==ar[1]){print $0"\t""Target"}}' > human_dog_8x1024_dog.fa.primAlign.blastId.Target.bed12

# Add %id from previsou analysis
awk 'BEGIN{while(getline<"human_dog_8x1024_dog.fa.primAlign.blastId.sam">0){id[$1]=$NF}}{if($4 in id){print $0"\t"id[$4]}}' human_dog_8x1024_dog.fa.primAlign.blastId.Target.bed12 > human_dog_8x1024_dog.fa.primAlign.blastId.Target.blastId.bed12

## Add alignement size x2
# genomic
# by blocks
more human_dog_8x1024_dog.fa.primAlign.blastId.Target.blastId.bed12 | awk '{size=0; n=split($11,ar,","); for (i in ar){size+=ar[i]}; print $0, size,$3-$2}' | colt
######################################@@@@@@
#CCL
# From the initial 1,267 dog promoter sequences
#	- 917 (72%) have an unique alignement of GRCH38
#	- 448 (35%) have an unique alignement of GRCH38 that overlap the corresponding human gene promoter
