# Directory storing all informations for Camille Kergal et al. manuscript

For most of the analyses, the human reference genome was `GRCh38` and the dog reference genome was `canFam4` and `canFam3`.

The structure of this repository is as following:

- `input_data` : data files used in the manuscript:
  - `coordinates/` : important files in `.bed`
    - `orthologous_cancer_genes.bed` : files of the 1,317 promoters of the orthologous cancer genes panel on canFam4 and GRCh38 genomes assembles.
    Columns of this file describe : Gene NAME, coordinates on dog genome (D prefix) and coordinates on human genome (H prefix).
  - `models/` : dog models of gene expression prediction for canFam3 and canFam4 with `.json` for NN architecture and `.txt` for correspondance with tissue files.
    - `canFam4.h5` = can be reached here <http://tools.genouest.org/data/tderrien/cf4_pred_model.h5>
    - `canFam3.h5` = can be reached here <http://tools.genouest.org/data/tderrien/cf3_pred_model.h5>
- `scripts` = scripts to generate figures/data of the paper
  - `conservation/` : Analysis of promoter evolutionary conservation
    - `cons_prom_minimap.sh` : shell script using minimap2/bedtools to anlayze conservation of dog promoters sequences on the human genome. 
     - `human_dog_8x1024_dog.fa.primAlign.blastId.Target.blastId.bed12` : output files of canine promoters on human genome.
  - `transposable_elements/` : analysis of cross-species (human) versus within-species (dog) on dog test sequences wrt to the content in trnasposable elements (TE).
    - `dog_human_cross_1119.txt` : predicted expression with human and dog models for cancer genes.
    - `cf4_1024prom_CancerGenes.vrai_OVER_TranspEl.txt`: TE content per cancer gene.
    - `eda_prom_HumModel_clean_tissu_grpTE.R` : Rscript to generate figureXX
  - `GC_content/` : 
    - [ ] **TODO** influence of GC content on gene expression prediction
  - `sat_mut/` : in silico saturated mutagenesis of the given test sequences
    - `cancer_panel_genes/dog/` and `cancer_panel_genes/human/` : as tested sequences for dog and human respectively
      - `outliers.tsv` : with standard deviation >=8
      - `outliers_sup4.tsv ` : with standard deviation >=4
      - [ ] **TODO** for orthologous human genes with same SD and matching tissues
    - `all_genes/` 
      - [ ] **TODO** tested sequences all genes cf4
