![Logo BLIMP](./logo-BLIMP.png)

# BLIMP : Basenji Like IMpacting variants Predictor

Based on the [Basenji](https://github.com/calico/basenji) tool suite, we propose :

- two prediciton models of canine gene expression (related on both canFam3 and canFam4 geneome assemblies)
- a user-friendly way to perform *in silico* saturated mutagenesis.

Specifying genomic regions of interest (max size = 1024 bp) and tissues of interest (among those composing prediciton models), you can obtain predicitons impact of all possible mutations on gene expression.

In order to evaluate our prediction model, we redesigned some Basenji scripts to realise and assess cross-species and within-species predictions. This process is presented [here](https://github.com/ckergal/BLIMP/blob/main/bin/README.md).

## Installation

See `basenji` tutorial for the general installation, dataset creation and model training:  <https://github.com/calico/basenji#installation>

For specific scripts used to predict impacting mutations :

```bash

git clone https://github.com/ckergal/BLIMP
```

## Data requirement

If you're interested in canine genomics, we trained prediction models based on canFam3 and canFam4 verison assemblies you can download just below.

```sh
curl -o cf4_pred_model.h5 \
  http://tools.genouest.org/data/tderrien/cf4_pred_model.h5
  
curl -o cf3_pred_model.h5 \
  http://tools.genouest.org/data/tderrien/cf3_pred_model.h5
```

Otherwise, you can train a prediciton model of gene expression with your own CAGE data following different steps :

#### 1. Data Integration

<details>
<summary>
<i>Click to see instructions.</i>
</summary>

The first step to establish prediction model is the data integration. The Basenji framework consists in joining sequence biochemical activities (*e.g.* transcriptomics CAGE profiles) to a reference genome. As these data (e.g. CAGE profiles) are in a `.bam` format, the use of Basenji requires a conversion into BigWig format. This can be done with the next command line, applied on all NGS profiles you want to integrate into your prediction model.
  
 **Scripts**

```
bam_cov.py \
  -o CAGE/BW/ \
  CAGE/BAM/cage_profile.bam \
  CAGE/BW/cage_profile.bw
```

The `-o` option allows to create the folder where BigWig files are generated. Then, `.bw` files can be used in the next step of data integration.

In our case, we applied this command line :

```sh
basenji_data.py \
  -d 1 -l 131072 \ 
  --local \
  -t .1 -v .1 -w 128 \
  -c 0 \
  -o data/cf4_basenji_data \
  data/canFam4.fa \
  data/targets.txt
```

`-d` is the genome part you want to use (here 100%) and `-l` is the sequence length. `-t` and `-v` are the part of test and validation set respectively and `-w` represents the window size of base pairs aggregation. `-c` is a cropping factor and `-o` is the directory where output files are generated. The first argument is the genome Fasta file and the second argument is a text file specifying BigWig files to integrate, as presented [here](https://github.com/ckergal/BLIMP/blob/main/manuscript/input_data/models/cf4_targets.txt).

**Output**

Then, in the directory **cf4_basenji_data** you have different files and sub-directories. The file `contigs.bed` contains genomic regions (contigs) in which sequences of length `-l` are picked and `sequences.bed` contains genomic coordinates of 131kb sequences, with their distribution between **train**, **validation** and **test** sets. The file `targets.txt` is duplicated here and the file `statistics.json` contains informations about this data integration step.

The folder **seqs_cov** contains `.err` and `.h5` files numbered from 0 to *number_of_CAGE_profiles* according to their position in `targets.txt`. `.err` files should be empty and `.h5` files contain expression profiles of CAGE sequenced tissues, as used in the algorithm. So he continuous signal of transcription level is aggregate in 128 bp windows.

The folder **tfrecords** contains `.err` and `.tfr` files. Hopefully `.err` files are empty too and each `.tfr` file contains 256 sequences of 131kb and the associated expression level for all CAGE samples.

</details>

#### 2. Model Training

<details>
<summary> <i>Click to see instructions.</i> </summary>

When data integration is achieved, the second step consists in training the prediction model of gene expression. The algorithm learns with the training data set and assessment during training is established on the validation data set. A model obtained during an epoch is used to predict expression level of genomic sequences from the validation set and predictions are compared to the real expression level with a Pearson correlation coefficient. Sothe aim is to maximize this value.

**Script**

```shell
basenji_train.py \
  -o models/cf4_basenji_train \
  models/cf4_params.json \
  data/cf4_basenji_data
```

`-o` parameter allows to create the output directory. The first argument is the file of hyperparamters and the second one is the path to the data directory.

**Output**
The folder **cf4_basenji_train** contains the file `params.json`, duplicating the one in argument. `model_check.h5` saves the model corresponding to the last epoch and `model_best.h5` corresponds to the model leading to the best performance (e.g. correlation coefficient).
</details>

#### 3. Model Test

<details>
<summary>
<i>Click to see instructions.</i>
</summary>

The last step to completely train a prediction model consists in assessing it. As with the validation data set, it's about predicting the expression level of genomic sequences from the test data set and comparing predictions with real data calculating a Pearson correlation coefficient.

**Script**

```
basenji_test.py \
  -o output/cf4_basenji_test \
  --save \
  models/cf4_basenji_train/params.json \
  models/cf4_basenji_train/model_best.h5 \
  data/cf4_basenji_data
```

`-o` allows to create the output folder and `--save` means we want to store values of predicitons and real expression levels of test genomic sequences. The first argument is the path to the params used to train the prediction model, the second is the pathto the prediciton model we want to test and the third is the path to the folder of data we used during training.

**Output**
In the folder indicated as output, here **cf4_basenji_test**, you can find a file called `acc.txt`. It lists Pearson correlation coefficients between predited expression level and the real one for test sequences in each CAGE samples composing the model.
</details>

## Usage

### *In silico* saturated mutagenesis

The tool ***in silico* saturated mutagenesis** allows to predict the impact of all possible mutations occuring in a given genomic sequence on gene expression.

**Script**

```
basenji_sat.py \
  -f data/canFam4.fa \
  -o output/sat_mut \
  -l 1024 \
  -t data/targets.txt \
  models/cf4_basenji_train/params.json \
  models/cf4_basenji_train/model_best.h5 \
  data/seq_1024bp.bed
```

Element specified with `-f` option is the genome FASTA file you want to extract reference alleles and `-o` allows to create the output folder. `-l` is to precise the length of genomic sequences of interest and
`-t` is the text file specifying CAGE samples of interest. You can select fewer than those composing the model.

**Output**

The output of this step is a file named *scores.h5*. This format is totally suitable to store the result of *in silico* saturated mutagenesis. It gives predictions of impact of `-l` possible mutations from the number of regions described in `.bed` file in the different tissues presented in `-t` options.

### Impacting variants detection

In order to detect impacting mutations from the set resulting from the saturated mutagenesis, we established a method implemented in Python. The aim is to detect outliers from the dataset that we then count as the most impating variants. We also propose a more user-friendly output format.

**Script**

```
outliers.py \
  --sd 4 \
  output/sat_mut/scores.h5
```

`--sd` allows to precise the threshold of standard deviation to select the most impacting variants and the argument corresponds to the path to the `scores.h5` file from the saturated mutagenesis.

**Output**

Results are stored in `outliers.tsv` file. Each line describes a mutation predicted as impacting gene expression. We find the genomic location and gene name, reference and mutate allele, differential expression level implied (`outlier`). `Distance` column indicates how impacting is the mutation (the *sd* value). We also inform about the tissue where the mutation is tested and how far from the genomic coordinates of the bed file.
