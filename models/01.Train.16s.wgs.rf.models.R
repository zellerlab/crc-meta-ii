##########################################
### Train 16S and WGS ml models

# Load functions
source(here('src','utils.R'))

# Load balanced metadata and data
all.meta<- read_tsv(here('data',"Metadata.all.samples.balanced.tsv")) %>% filter(Condition=='CRC'| Condition=='CTR') %>% 
  as.data.frame() %>% column_to_rownames('Sample_ID')

all.data <- read.table(here('data','Relab.all.samples.balanced.tsv'),sep='\t', check.names = F) %>%
  rownames_to_column('genus') %>% filter(genus!='unassigned') %>% column_to_rownames('genus')


# Extract 16S and WGS data
meta_df_16s  = all.meta %>%  filter(Assay=='16S')

models.rf.16S<-train_model_rf(meta_df = meta_df_16s, mat = all.data, label_column = 'Condition', case_label = 'CRC',control_label = 'CTR',
                              seed = 2025,block_label ='Cohort',num_trees = 200, prev.cutoff = 0.1)

meta_df_wgs  = all.meta %>%  filter(Assay=='WGS')

models.rf.wgs<-train_model_rf(meta_df = meta_df_wgs, mat = all.data,label_column = 'Condition',case_label = 'CRC',control_label = 'CTR',
                        seed = 2000,block_label ='Cohort',num_trees = 200, prev.cutoff = 0.1)


# 16S hold out testing using WGS trained model

label<-create.label(meta = meta_df_16s,label='Condition',case='CRC',control = 'CTR')

siamcat.test <- siamcat(feat=all.data, meta=meta_df_16s,
                        label=label, case='CRC')


siamcat.test <- normalize.features(siamcat.test,
                                 norm.param=norm_params(models.rf.wgs),
                                 feature.type='original',
                                 verbose = 2)

siamcat.test.predicted<- make.predictions(models.rf.wgs , siamcat.test)

siamcat.test.evaluated.16s.holdout.rf<-evaluate.predictions(siamcat.test.predicted)

print(siamcat.test.evaluated.16s.holdout.rf)

# WGS hold out testing using 16s trained model

label<-create.label(meta = meta_df_wgs, label='Condition',case='CRC',control = 'CTR')

siamcat.test <- siamcat(feat=all.data, meta=meta_df_wgs,
                        label=label, case='CRC')


siamcat.test <- normalize.features(siamcat.test,
                                   norm.param=norm_params(models.rf.16S),
                                   feature.type='original',
                                   verbose = 2)

siamcat.test.predicted<- make.predictions(models.rf.16S , siamcat.test)
siamcat.test.evaluated.wgs.holdout.rf <- evaluate.predictions(siamcat.test.predicted)
print(siamcat.test.evaluated.wgs.holdout.rf)

# Save models
save(models.rf.16S, models.rf.wgs, siamcat.test.evaluated.16s.holdout.rf,siamcat.test.evaluated.wgs.holdout.rf,
     file= here('data','results','Training.16s.wgs.rf.models.Rdata'))
    

##########################################
### Train genus and species level models on only WGS dataset

#Load motus level data
all.meta.wgs <- read_tsv(here('data',"Metadata.all.samples.tsv")) %>% filter(Condition=='CRC'| Condition=='CTR') %>% 
  as.data.frame() %>%
  filter(Assay=='WGS')

all.data.motus <- read.table(here('data','Raw.counts.wgs.motus.all.samples.tsv'), sep='\t', header = T, check.names = F)

all.data.motus<- all.data.motus  %>% mutate(across(everything(), ~ . / sum(.)))  %>%
  rownames_to_column('motus') %>%
  filter(motus!='unassigned') %>%
  column_to_rownames('motus')


all.data.motus<- all.data.motus %>% mutate(across(everything(), ~ . / sum(.)))


# Balance it

balanced_motus<-balance_data_by_groups(meta_df = all.meta.wgs, data_df = all.data.motus, group_column = 'Condition',groups = c('CRC','CTR'),seed = 300)


model.rf.wgs.motus<- train_model_rf(meta_df = balanced_motus$meta , mat = balanced_motus$data,label_column = 'Condition',
                                              case_label = 'CRC',control_label = 'CTR',
                                              num_trees = 200,seed = 2025 ,block_label = 'Cohort',prev.cutoff = 0.1)





#  Run a ml model with same samples but on genus level 



model.rf.wgs.genus<- train_model_rf(meta_df = balanced_motus$meta, mat = all.data,label_column = 'Condition',
                                                  case_label = 'CRC',control_label = 'CTR',
                                                  num_trees = 200,seed = 2025 ,block_label = 'Cohort',prev.cutoff = 0.1)

save(model.rf.wgs.genus,model.rf.wgs.motus, file=here('data','results','Training.wgs.rf.models.Rdata'))




