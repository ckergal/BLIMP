# Basenji Usage

Hear, we briefly describe how we managed the Basenji tool suite to train and use a prediction model of canine gene expression.



## Data integration

The first step to establish prediction model is the data integration. The Basenji framework consists in joining transcriptomics CAGE profiles to a reference genome.
As CAGE profiles are in a `.bam` format, the use of Basenji requires a conversion into BigWig format. This can be done with the next command line, applied on all CAGE profiles you want to integrate into your prediction model.

`bam_cov.py -o CAGE/BW/ CAGE/BAM/cage_profile.bam CAGE/BW/cage_profile.bw`


In our case, we applied this command line :

`basenji_data.py -d 1 -l 131072 --local -t .1 -v .1 -w 128 -c 0 -o data/cf4_basenji_data data/canFam4.fa data/targets.txt`
