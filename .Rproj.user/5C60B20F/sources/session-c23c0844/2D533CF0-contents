##########################################
### Train unified CRC ml model
# Load packages

library(dplyr)
library(tidyverse)
library(stringr)
library(ranger)
library(mlr3tuning)
library(mlr3extralearners)
library(ggembl)
library(mlr3)
library(SIAMCAT)
library(here)

# Load functions
source(here('src','utils.R'))

# Generate a list to save ml model
models.all.rf<- list()

# Use balanced data and metadata 
all.meta<- read_tsv(here('data',"Metadata.all.samples.balanced.tsv")) %>% filter(Condition=='CRC'| Condition=='CTR') %>% 
  as.data.frame() %>% column_to_rownames('Sample_ID')

all.data <- read.table(here('data','Relab.all.samples.balanced.tsv'),sep='\t', check.names = F) %>%
  rownames_to_column('genus') %>% filter(genus!='unassigned') %>% column_to_rownames('genus')


models.all.rf<-train_model_rf(meta_df = all.meta,mat = all.data,label_column = 'Condition',case_label = 'CRC',control_label = 'CTR',
                              num_trees = 200, seed = 2000,prev.cutoff = 0.1)


save(models.all.rf, 
     file='/g/scb/zeller/pekel/meta_analysis/src/analysis/Rdata_updated/Training.all.data.rf.Rdata')


# Plot unified CRC model

models <- list(models.all.rf)
labels <- c("Classifier cross validated on CRC")
colours <- c( "black" )

crc_model_auc_plot<- plot_roc_siamcat_models(models, labels, colours, alpha=0.8)


