options(encoding = "UTF-8")
library(tidyverse)
library(ggpubr)
library(scales)
library(matrixStats)
library(ggrepel)
library(RColorBrewer)


# Input data within-species cor
predexp       = read.table("dog_human_cross_1119.txt", h=T, stringsAsFactors = F)
# rename human col
colnames(predexp)[10:18] = 
  c("chr_hum","start_hum", "stop_hum", "preds_hum", "expr_dog","tissu_hum","model_hum","NAME_hum", "WINID_hum")

# Input = TE metrics
te            = read.table("cf4_1024prom_CancerGenes.vrai_OVER_TranspEl_Carnivora.txt", h=T, stringsAsFactors = F)
colnames(te)[4] = "NAME" # rename gene_name to allow merge

######################################################
# Compute correlation Pe and Spearman for all 8 windows by TE group
corall_te = left_join(predexp, te, by = c("NAME")) %>% drop_na()

corall_te_allwin = corall_te %>% 
  mutate(TE_content=cut(SINE_cov, breaks=c(-Inf, 1,  Inf), labels=c("TE-Low","TE-High"))) %>%
  group_by(tissu,TE_content) %>%
  mutate(cor_sp        = cor(preds,      expr, method = "sp"),
            cor_pe        = cor(preds,      expr, method = "pe"),
            cor_sp_cross  = cor(preds_hum,  expr, method = "sp"),
            cor_pe_cross  = cor(preds_hum,  expr, method = "pe")
  ) %>%
  select(tissu, cor_pe, cor_sp, TE_content, cor_sp_cross,cor_pe_cross ) %>% 
  distinct()



# All in one with arrow
 within_plt=ggplot(data=corall_te_allwin, aes(x=cor_pe, y=cor_pe_cross, color=TE_content) ) +  
  geom_abline(intercept = 0, slope = 1) +
  #geom_line(aes(group=tissu),  size=0.1) + 
  geom_path(aes(x = cor_pe, y = cor_pe_cross, group = tissu ), 
#            linetype ="dotted", size=0.5, color="black",
              arrow = arrow(length = unit(0.45, "cm"))) +
  geom_point(size=3) +
  geom_text_repel(
   data = subset(corall_te_allwin, TE_content=="TE-High"),
    aes(label = tissu), size = 5) + 
  scale_color_brewer(palette="Paired") +
  #stat_ellipse() +
  # xlim(0,1) +
  # ylim(0,1) +
  xlab("PearsonR - DOG trained (e.g. within-species)")+
  ylab("PearsonR - HUMAN trained (e.g. cross-species)")+
    theme_bw(base_size = 15) 

 pdf("eda_prom_HumModel_clean_tissu_grpTE.pdf")
 within_plt
 dev.off()
 
 
# Facet by TE content

  # Mean for facet
  means_within <- aggregate(cor_pe ~  TE_content      , corall_te_allwin, mean)
  means_across <- aggregate(cor_pe_cross ~  TE_content, corall_te_allwin, mean)

    ggplot(data=corall_te_allwin, aes(x=cor_pe, y=cor_pe_cross, color=TE_content) ) +  
    geom_abline(intercept = 0, slope = 1) +
    geom_line(aes(group=tissu),  size=0.1) + 
    geom_point(size=3) +
    facet_grid(. ~ TE_content )+
    geom_text(data=means_within, aes( label=
                                        paste("mean", round(cor_pe, digits = 2),sep=" = " )),size=5, color="black",
              x = Inf,  y = -Inf, hjust = 1, vjust = -0.5,) +
    geom_text(data=means_across, aes(label=paste("mean",round(cor_pe_cross, digits = 2),sep=" = ")), size=5, color="black",
               x = -Inf, y = Inf, hjust = 0, vjust = 1, inherit.aes = FALSE) +
    geom_text_repel(
      data = subset(corall_te_allwin, cor_pe_cross<cor_pe),
      aes(label = tissu), size = 5) + 
    scale_color_brewer(palette="Paired") +
    #stat_ellipse() +
      # xlim(0,1) +
      # ylim(0,1) +

      xlab("PearsonR - DOG trained (e.g. within-species)")+
      ylab("PearsonR - HUMAN trained (e.g. cross-species)")+
    theme_bw(base_size = 20) +
    theme(legend.position="none") + theme(legend.title= element_blank())
    
  
#################################################################################
# Stop here
#################################################################################
    
# Try with new file 

predexp       = read.table("atester_clean.bed", h=T, stringsAsFactors = F)
# rename columns as in the first file
predexp <- predexp %>% rename( preds = pred_dog,  preds_hum = preds_cross)


    
#  data summarizing of the 8 window in sum or max and then compute correlations sp and pe
predexp_cor_wind = predexp %>% group_by(NAME,tissu) %>%
  summarize(
    preds_sum     = sum(preds),
    preds_max    = max(preds),
    expr_sum      = sum(expr),
    expr_max     = max(expr),
    preds_sum_hum = sum(preds_hum),
    preds_max_hum= max(preds_hum)
    ) 


# Then merged Correlation with TE content and drop NA
corall_te_wind = left_join(predexp_cor_wind, te, by = c("NAME")) %>% drop_na()

corall_te_wind_TEcat = corall_te_wind %>% 
  #mutate(TE_content=cut(SINE_cov, breaks=c(-Inf, 1, 10, Inf), labels=c("Low","Medium","High"))) %>%
  mutate(TE_content=cut(SINE_cov, breaks=c(-Inf, 0, Inf), labels=c("Low","High"))) %>% 
  filter(expr_max>0 & preds_max>0) %>%  
  group_by(tissu,TE_content) %>%  #   summarise(n = n())
  summarise(cor_sum_sp  = cor(preds_sum,  expr_sum, method = "sp"),
         cor_max_sp = cor(preds_max, expr_max, method = "sp"),
         cor_sum_pe  = cor(preds_sum,  expr_sum, method = "pe"),
         cor_max_pe = cor(preds_max, expr_max, method = "pe"),
         cor_sum_sp_cross  = cor(preds_sum_hum,  expr_sum, method = "sp"),
         cor_max_sp_cross = cor(preds_max_hum, expr_max, method = "sp"),
         cor_sum_pe_cross  = cor(preds_sum_hum,  expr_sum, method = "pe"),
         cor_max_pe_cross = cor(preds_max_hum, expr_max, method = "pe")
  )  %>% group_by(tissu,TE_content) # %>%     summarise(n = n())  #%>%
#  filter(tissu=="thalamus" & TE_content =="High") 
  
View(corall_te_wind_TEcat)
# plot   
ggplot(corall_te_wind_TEcat,  aes(x=cor_max_sp, y=cor_max_sp_cross, color=TE_content) ) +
   geom_abline(intercept = 0, slope = 1) +
   geom_line(aes(group=tissu),  size=0.1) + 
   geom_point(size=6) +
   geom_text_repel(
     data = subset(  corall_te_wind_TEcat, TE_content=="High"  ),
     aes(label = tissu), size = 5) + 
     facet_wrap( tissu ~. )+
  #   facet_wrap( TE_content ~. , nrow = 1)+
   # xlim(0.15,1) +
   # ylim(0.40,1) +
   scale_color_brewer(palette="Paired") +
   theme_bw(base_size = 16) 

