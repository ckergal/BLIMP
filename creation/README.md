# Basenji Usage

Hear, we briefly describe how we managed the [Basenji](https://github.com/calico/basenji) tool suite to train and use a prediction model of canine gene expression.



## Data integration

The first step to establish prediction model is the data integration. The Basenji framework consists in joining transcriptomics CAGE profiles to a reference genome.
As CAGE profiles are in a `.bam` format, the use of Basenji requires a conversion into BigWig format. This can be done with the next command line, applied on all CAGE profiles you want to integrate into your prediction model.

```
bam_cov.py \
-o CAGE/BW/ \
CAGE/BAM/cage_profile.bam \
CAGE/BW/cage_profile.bw
```


The -o option allows to create the folder where BigWig files are generated. Then, `.bw` files can be used in the next step of data integration.



In our case, we applied this command line :

```shell
basenji_data.py \
-d 1 -l 131072 \ 
--local \
-t .1 -v .1 -w 128 \
-c 0 \
-o data/cf4_basenji_data \
data/canFam4.fa \
data/targets.txt
```

`-d`is the genome part you want to use (here 100%) and `-l` is the sequence length. `-t` and `-v` are the part of test and validation set respectively and `-w` represents the window size of base pairs aggregation. `-c` is a cropping factor and `-o` is the directory where output files are generated. The first argument is the genome Fasta file and the second argument is a text file specifying BigWig files to integrate, as presented [here](https://github.com/ckergal/BLIMP/blob/main/manuscript/input_data/models/cf4_targets.txt).


### Output

Then, in the directory **cf4_basenji_data** you have different files and sub-directories. The file `contigs.bed` contains genomic regions (contigs) in which sequences of length `-l` are picked and `sequences.bed` contains genomic coordinates of 131kb sequences, with their distribution between train, validation and test sets. The file `targets.txt` is diplicated here and the file `statistics.json` contains informations about this data integration step.


The folder **seqs_cov** contains `.err` and `.h5` files numbered from 0 to *number_of_CAGE_profiles* according to their position in `targets.txt`. `.err` files should be empty and `.h5` files contain expression profiles of CAGE sequenced tissues, as used in the algorithm. So he continuous signal of transcription level is aggregate in 128 bp windows.


The folder **tfrecords** contains `.err` and `.tfr` files. Hopefully `.err` files are empty too and each `.tfr` file contains 256 sequences of 131kb and the associated expression level for all CAGE samples.




## Training algorithm
When data integration is achieved, the second step consists in training the prediction model of gene expression. 
