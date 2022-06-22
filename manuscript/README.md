# Directory storing all informations for Camille Kergal et al. manuscript

For most of the analyses, the dog reference genome was `canFam4` and the human reference genome was `GRCh38`.

The structure of this repository is as following:


- `input_data` : data files used in the manuscript:
  - `coordinates/` : important files in `.bed`
    - `cf3_cancer_1024prom.bed` : files of the 1,332 promoters of the cancer gene panel on canFam3 genome assembly.
    - `cf4_cancer_1024prom.bed` : files of the 1,332 promoters of the cancer gene panel on canFam4 genome assembly.
    - `human_dog_8x1024_dog.fa.primAlign.blastId.Target.blastId.bed12` : conservation analysis of canine promoters on human genome as given by the script `../../script/cons_prom_minimap.sh`.
  - `models/` : dog models of gene expression prediction
- `scripts` = scripts to generate figures/data of the paper
  - `cons_prom_minimap.sh` : shel script using minimap2/bedtools to anlayze conservation of dog promoters sequences on the human genome (output file is `../input_data/coordinates/human_dog_8x1024_dog.fa.primAlign.blastId.Target.blastId.bed12`)
