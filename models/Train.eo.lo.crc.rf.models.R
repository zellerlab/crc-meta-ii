##########################################
### Train EO-CRC and LO-CRC ml models
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

# Generate lists to save models
models.eo.rf <- models.lo.rf<- list()

# Balance EO and LO-CRC cases to train models

all.meta<- read_tsv(here('data',"Metadata.all.samples.tsv")) %>% filter(Condition=='CRC'| Condition=='CTR') %>% 
  as.data.frame() %>% column_to_rownames('Sample_ID')

all.data <- read.table(here('data','Relab.all.samples.tsv'),sep='\t', check.names = F) %>%
  rownames_to_column('genus') %>% filter(genus!='unassigned') %>% column_to_rownames('genus')


balanced_data_metadata_eocrc<-balance_data_by_groups(meta_df = all.meta %>%
                                                       rownames_to_column('Sample_ID') %>%
                                                       filter(Group %in% c('EO-CRC', 'EO-CTR')) %>%
                                                       group_by(Cohort) %>%
                                                       filter(all(c('EO-CRC', 'EO-CTR') %in% Group)) %>%
                                                       ungroup(), data_df=all.data, group_column = 'Group', groups  = c('EO-CRC','EO-CTR'), seed = 2002)

balanced_data_metadata_locrc<-balance_data_by_groups(meta_df = all.meta %>%
                                                       rownames_to_column('Sample_ID') %>%
                                                       filter(Group %in% c('LO-CRC', 'LO-CTR')) %>%
                                                       group_by(Cohort) %>%
                                                       filter(all(c('LO-CRC', 'LO-CTR') %in% Group)) %>%
                                                       ungroup(), data_df=all.data, group_column = 'Group', groups  = c('LO-CRC','LO-CTR'), seed = 2002)


models.eo.rf<-train_model_rf(meta_df = balanced_data_metadata_eocrc$meta,mat = balanced_data_metadata_eocrc$data,label_column = 'Group',case_label = 'EO-CRC',
                             control_label = 'EO-CTR', block_label = 'Cohort', num_trees = 200, seed = 2000,prev.cutoff = 0.1)

models.lo.rf<-train_model_rf(meta_df = balanced_data_metadata_locrc$meta,mat = balanced_data_metadata_locrc$data,label_column = 'Group',case_label = 'LO-CRC',
                             control_label = 'LO-CTR', block_label = 'Cohort', num_trees = 200, seed = 2000,prev.cutoff = 0.1)


# LO-CRC hold out testing using EO-CRC trained models

label<-create.label(meta = balanced_data_metadata_locrc$meta,label='Group',case='LO-CRC',control = 'LO-CTR')

siamcat.test <- siamcat(feat= balanced_data_metadata_locrc$data, meta= balanced_data_metadata_locrc$meta,
                       label=label, case='LO-CRC')

siamcat.test.predicted<- make.predictions(models.eo.rf , siamcat.test)

siamcat.test.evaluated.lo.holdout.rf<-evaluate.predictions(siamcat.test.predicted)
print(siamcat.test.evaluated.lo.holdout.rf)

# EO-CRC hold out testing using LO-CRC trained models

label<-create.label(meta = balanced_data_metadata_eocrc$meta, label='Group',case='EO-CRC',control = 'EO-CTR')

siamcat.test <- siamcat(feat=balanced_data_metadata_eocrc$data, meta=balanced_data_metadata_eocrc$meta,
                       label=label, case='EO-CRC')

siamcat.test.predicted<- make.predictions(models.lo.rf , siamcat.test)
siamcat.test.evaluated.eo.holdout.rf <- evaluate.predictions(siamcat.test.predicted)
print(siamcat.test.evaluated.eo.holdout.rf)

save(models.eo.rf,models.lo.rf,siamcat.test.evaluated.eo.holdout.rf,siamcat.test.evaluated.lo.holdout.rf,
     file='/g/scb/zeller/pekel/meta_analysis/src/analysis/Rdata_updated/Training.eolo.all.data.Rdata')


# Plot ROC for the EO-CRC model

models <- list(models.eo.rf, siamcat.test.evaluated.eo.holdout.rf)
labels <- c("Classifier cross validated on EO-CRC","Classifier trained on LO-CRC and tested on EO-CRC")
trained_on <- list(NULL, models.lo.rf)
colours <- c( "darkorange", "black" )
eo_models_auc_plot <- plot_roc_siamcat_models(models, labels, colours, trained_on, alpha=0.7)


# Plot ROC for the LO-CRC model

models <- list(models.lo.rf, siamcat.test.evaluated.lo.holdout.rf)
labels <- c("Classifier cross validated on LO-CRC", 'Classifier trained on EO-CRC and tested on LO-CRC"')
trained_on <- list(NULL, models.eo.rf)
colours <- c( "darkred" , 'black')

lo_models_auc_plot<- plot_roc_siamcat_models(models, labels, colours, trained_on,alpha=0.8)

