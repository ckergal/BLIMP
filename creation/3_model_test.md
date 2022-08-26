# Model test
The last step to completely train a prediction model consists in assessing it. As with the validation data set, it's about predicting the expression level of genomic sequences from the test data set and comparing predictions with real data calculating a Pearson correlation coefficient.



## Script
```
basenji_test.py \
-o output/cf4_basenji_test \
--save \
models/cf4_basenji_train/params.json \
models/cf4_basenji_train/model_best.h5 \
data/cf4_basenji_data
```

`-o` allows to create the output folder and `--save` means we want to store values of predicitons and real expression levels of test genomic sequences. The first argument is the path to the params used to train the prediction model, the second is the pathto the prediciton model we want to test and the third is the path to the folder of data we used during training.



## Output
In the folder indicated as output, here **cf4_basenji_test**, you can find a file called `acc.txt`. It lists Pearson correlation coefficients between predited expression level and the real one for test sequences in each CAGE samples composing the model.
