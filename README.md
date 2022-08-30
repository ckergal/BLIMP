# BLIMP : Basenji Like IMpacting variants Predictor


Based on the [Basenji](https://github.com/calico/basenji) tool suite, we propose two prediciton models of canine gene expression (related to canFam3 and canFam4) and an easy way to perform *in silico* saturated mutagenesis.

Specifying genomic regions of interest (max size = 1024 bp) and tissues of interest (among those composing prediciton models), you can obtain predicitons impact of all possible mutations on gene expression.

## Data requirement

Python script to get model files

```python
import os, subprocess

if not os.path.isfile('BLIMP/input_data/models/cf4_pred_model.h5'):
  subprocess.call('curl -o BLIMP/manuscript/input_data/models/cf4_pred_model.h5 http://tools.genouest.org/data/tderrien/cf4_pred_model.h5', shell=True)
  
if not os.path.isfile('BLIMP/input_data/models/cf4_pred_model.h5'):
  subprocess.call('curl -o BLIMP/manuscript/input_data/models/cf3_pred_model.h5 http://tools.genouest.org/data/tderrien/cf3_pred_model.h5', shell=True)
```


## Installation

See `basenji` tutorial for the general installation, dataset creation and model training:  <https://github.com/calico/basenji#installation>


For specific scripts used to predict impacting mutations :

```bash

git clone https://github.com/ckergal/BLIMP
```

## Usage

