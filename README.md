# BLIMP : Basenji Like IMpacting variants Predictor

            **(Note: This repo is in construction)**

Based on the [Basenji](https://github.com/calico/basenji) tool suite, we propose two prediciton models of canine gene expression (related to canFam3 and canFam4) and an easy way to perform *in silico* saturated mutagenesis.

Specifying genes of interest (max size = 1024 bp) and tissues of interest (among those composing prediciton models), you can obtain  predicitons impact of all possible mutations on gene expression.

## Data requirement

Follow the [import.ipynb](https://github.com/ckergal/BLIMP/blob/main/import.ipynb) to get models files.

## Installation

See `basenji` tutorial for the general installation:  <https://github.com/calico/basenji#installation>

For specific scripts used to create dog dataset, train the model and predict impacting mutations :

```bash

git clone https://github.com/ckergal/BLIMP
```
