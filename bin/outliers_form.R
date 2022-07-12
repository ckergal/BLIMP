#!/usr/local/bin/Rscript
# Packages ----------------------------------------------------------------
library(tidyverse)
library(showtext)
library(scales) ## show_col
library(ggsci) ## pal_uchicago
library(stringr)



# Options -----------------------------------------------------------------
args = commandArgs(TRUE)

dog_out_file = args[1]
min_sd = args[2]
out_dir = args[3]
# min_sd = 8
# out_dir = "~/outliers_results/output"


# Style -------------------------------------------------------------------
font_add_google("Poppins", "Poppins")
showtext_auto()
theme_set(theme_light(base_size = 13, base_family = "Poppins"))



# Import ------------------------------------------------------------------
dog_outliers = read.delim("~/outliers_results/dog_outliers_sup4.tsv")
hum_outliers = read.table("~/outliers_results/hum_outliers_sup4.tsv", sep = "\t", header = T)
genes_cf4 = read.table("~/BLIMP/manuscript/input_data/coordinates/cross_1024_prom.bed", sep = "\t", header = T)
targets = read.table("~/BLIMP/manuscript/input_data/models/targets.txt", sep = "\t", header = T)



# Manip -------------------------------------------------------------------
dog_outliers$species = "Dog"
dog_outliers$chrom = str_remove_all(str_remove(dog_outliers$chrom, "b"), "'")

hum_outliers$species = "Human"
hum_outliers$chrom = str_remove_all(str_remove(hum_outliers$chrom, "b"), "'")



# Manage positions --------------------------------------------------------
position = 1023:0

### Human
hum_outliers = hum_outliers %>% 
  left_join(genes_cf4[,c("NAME", "HPROMSTART", "HSTRAND")], c("start" = "HPROMSTART")) %>% 
  left_join(targets[,c("index", "description")], c("tissue" = "index")) %>%
  select(-tissue, -gene, -strand) %>% 
  rename(gene = NAME, tissue = description, strand = HSTRAND)

hum_outliers$distTSS = ifelse(
  hum_outliers$strand == "1", 
  position[hum_outliers$position + 1], 
  hum_outliers$position)

hum_outliers$genomic_position = 
  ifelse(hum_outliers$strand == 1,
         hum_outliers$start + hum_outliers$position + 1,
         hum_outliers$start + hum_outliers$position)

write.table(hum_outliers, 
            "~/outliers_results/hum_4_cancer_cross_name.tsv", 
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F)


### Dog
dog_outliers = dog_outliers %>%
  left_join(genes_cf4[, c("NAME", "PROMSTART", "STRAND")], c("start" = "PROMSTART")) %>%
  left_join(targets[,c("index", "description")], c("tissue" = "index")) %>%
  select(-tissue, -gene, -strand) %>%
  rename(gene = NAME, tissue = description, strand = STRAND)

dog_outliers$distTSS = ifelse(dog_outliers$strand == "+",
       position[dog_outliers$position + 1],
       dog_outliers$position)

dog_outliers$genomic_position = ifelse(dog_outliers$strand == "+",
       dog_outliers$start + dog_outliers$position + 1,
       dog_outliers$start + dog_outliers$position)

write.table(dog_outliers, 
            "~/outliers_results/dog_4_cancer_cross_name.tsv", 
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F)



# Outliers filter ---------------------------------------------------------
outliers = rbind(hum_outliers, dog_outliers)
outliers = outliers %>% filter(distance >= min_sd)

com_outliers = outliers %>%
  group_by(species, distTSS, tissue, gene) %>% 
  summarise(n()) %>%
  group_by(distTSS, tissue, gene) %>% summarise(nb_species = n()) %>%
  filter(nb_species == 2)



# Descriptive stats -------------------------------------------------------

## Number of outliers per tissues and species
c1 = as.data.frame(unique(cbind(distTSS = outliers$distTSS, 
                           gene = outliers$gene, 
                           tissue = outliers$tissue, 
                           species = outliers$species))) %>%
  group_by(species, tissue) %>% 
  summarise(nb_outliers = n()) %>% 
  ggplot() + 
  aes(x = reorder(tissue, nb_outliers), y = nb_outliers, fill = species) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c( "#155F83FF", "#FFA319FF")) +
  theme(axis.text.y = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14)) +
  labs(title = "Outliers per tissue") +
  coord_flip()
ggsave(paste(out_dir, "out_tissues.png", sep = "/"), c1, dpi = 1000)



### Pourcentage outliers
basic_stats = rbind(
  cbind(
    "Impacting outliers by position and canine gene: ",
    paste(round(100 * nrow(unique(dog_outliers[,c("position", "gene")])) / (1317*1024), 2), "%")),
  
  cbind(
    "Impacting outliers by tissue, position and canine gene: ",
    paste(round(100 * nrow(unique(dog_outliers[,c("position", "gene", "tissue")])) / (1317*1024*19), 2), "%")),
  
  cbind(
    "Impacting outliers by position and human gene: ",
    paste(round(100 * nrow(unique(hum_outliers[,c("position", "gene")])) / (1317*1024), 2), "%")),
  
  cbind(
    "Impacting outliers by tissue, position and human gene: ",
    paste(round(100 * nrow(unique(hum_outliers[,c("position", "gene", "tissue")])) / (1317*1024*19), 2), "%")),
  
  cbind(
    "Common outliers by gene and TSS distance: ",
    nrow(unique(com_outliers[,c("distTSS", "gene")])))
)

write.table(basic_stats,
            paste(out_dir, "stats.tsv", sep="/"),
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F)




### Distance to TSS
c2 = outliers %>%
  ggplot() +
  aes(x = distTSS, fill = species) +
  geom_histogram(
    binwidth = 30, 
    col = "grey87", 
    alpha = .9) +
  labs(x="TSS Distance (nt)",
       y= "Impacting variants") +
  facet_grid(~species) + 
  scale_fill_manual(values = c("#155F83FF", "#FFA319FF")) +
  theme(strip.text = element_text(size = 10, color = "gray30"),
        strip.background = element_rect(fill = "gray96"),
        legend.position = "None")

ggsave(paste(out_dir, "tss_dist.png", sep = "/"), 
       c2, 
       dpi = 1000, 
       width = 9, 
       height = 3)
  
