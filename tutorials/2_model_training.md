# Training algorithm
When data integration is achieved, the second step consists in training the prediction model of gene expression. The algorithm learns with the training data set and assessment during training is established on the validation data set. A model obtained during an epoch is used to predict expression level of genomic sequences from the validation set and predictions are compared to the real expression level with a Pearson correlation coefficient. Sothe aim is to maximize this value.



## Script
```shell
basenji_train.py \
  -o models/cf4_basenji_train \
  models/cf4_params.json \
  data/cf4_basenji_data
```

`-o` parameter allows to create the output directory. The first argument is the file of hyperparamters and the second one is the path to the data directory.


## Output
The folder **cf4_basenji_train** contains the file `params.json`, duplicating the one in argument. `model_check.h5` saves the model corresponding to the last epoch and `model_best.h5` corresponds to the model leading to the best correlation coefficient.
