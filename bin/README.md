# Scripts info
In the case you realised your own prediction model of gene expression, you might be interesseted in comparing it with models already established. 
In our work, we used this approach in order to assess the interest of developping such models. In fact, the human prediciton model developped by Kelley et *.al*
can be employed in a cross-species approach so we wanted to compare both methods.


## Usage

### Data
In order to conduct cross-species study, prediction model of human gene expression is necessary. You can find it [here](https://github.com/calico/basenji/tree/master/manuscripts/cross2020). 


### Script
```
basenji_predict_par.py \
  -o output/cross_species \
  -t data/targets_human.txt \
  -f data/canFam4.fa \
  --expr data/cf4_basenji_data \
  models/human_model/params.json \
  models/human_model/model.h5 \
  data/human.bed \
  data/dog.bed
```

The option `-o` is the path where you want to find outputs. `-t` is the `targets.txt` file where you specify samples of interest among those composing the prediction model (human here). `-f` is the genome you want to predict the expresison level. `--expr` leads to the folder created by the [data integration](https://github.com/ckergal/BLIMP#1-data-integration) step. The first argument is the parameters file of the prediction model and the second one the file of the model (here human ones). Finally, the first `.bed` file represents genomic regions of interest from the species of the model (here human) and the second `.bed` file is composed of the orthologous genomic regions.


### Output
