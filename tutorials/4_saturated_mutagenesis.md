# *In silico* saturated mutagenesis
The tool of *in silico* saturated mutagenesis allows to perform a prediction of gene expression impact induced by all possible mutations occuring in genomic sequences.


## Script
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
