#=============Stacked barplot=================
library(ggplot2)
library(cowplot)
library(ggpubr)
library(remotes)
library(immunedeconv)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(RColorBrewer)
#This script details the proteomic deconvolution steps
#============Deconvolution method============================
df = read.table("proteinMatrix.txt",header = T,sep = '\t',row.names = 1)
df

#cibersort deconvolution
set_cibersort_binary("C:/Users/Javan/Desktop/sars-cov-2Project/data/human_sarscov_db/cybersortFiles/CIBERSORT.R")
set_cibersort_mat("C:/Users/Javan/Desktop/sars-cov-2Project/data/human_sarscov_db/cybersortFiles/LM22.txt")

#check the available deconvolution methods
deconvolution_methods

# 1. quantiseq method
res_quantiseq = deconvolute(df, "quantiseq", tumor = F)

res_quantiseq %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_brewer(palette="Paired") +
  scale_x_discrete(limits = rev(levels(res_quantiseq)))


#test

df = read.csv("testData.csv",header = T,sep = ',')

colourCount = length(unique(df$cell_type))

df %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  scale_x_discrete(limits = rev(levels(df))) +
  scale_fill_manual(values = get_palette(palette = "Paired", 22)) +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),plot.background = element_rect(fill = "white")
        ,strip.text = element_text(colour = "white",face = "bold")) +
  labs(x = NULL, y = "Immune cell fraction in %")









