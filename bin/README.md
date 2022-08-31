# Scripts info
In the case you realised your own prediction model of gene expression, you might be interesseted in comparing it with models already established. 
In our work, we used this approach in order to assess the interest of developping such models. In fact, the human prediciton model developped by Kelley et *.al*
can be employed in a cross-species approach so we wanted to compare both methods.


## Usage

### Script


```
basenji_predict_par.py \
  -o output/cross_species \
  -t data/targets.txt \
  -f data/canFam4.fa \
  -p 6 \
  --expr data/cf4_basenji_data \
  -l 40 \
  --length 131072 \
  models/cf4_basenji_train/
```
