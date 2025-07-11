##############################
# LOSO Training & Evaluation
##############################

# Load libraries

# The script run on bash with running together with requirem

# Script runs loso model for each datasets

commmandArgs <- commandArgs(trailingOnly = TRUE)

source(here('src','utils.R'))

# Load balanced data/metadata and keep only CRC and CTR samples

meta_balanced<- read_tsv(here('data','Metadata.all.samples.balanced.tsv')) %>%
  filter(Condition=='CRC'| Condition=='CTR') %>% as.data.frame() %>% column_to_rownames('Sample_ID')

data_balanced <- read.table(here('data', 'Relab.all.samples.balanced.tsv'),sep='\t', check.names = F) %>%
  rownames_to_column('genus') %>% filter(genus!='unassigned') %>% column_to_rownames('genus')


# Train a LOSO (Leave-One-Study-Out) model for each cohort  

study<- commmandArgs[1]

print(study)

models.loso<-list()

meta.train <- meta_balanced %>%
  filter(Cohort != study)

feat.train <- data_balanced[,rownames(meta.train)]

siamcat<-train_model_rf(meta_df = meta.train, mat = feat.train,label_column = 'Condition',case_label = 'CRC',control_label = 'CTR',
                         seed = 2000,block_label ='Cohort',num_trees = 200, prev.cutoff = 0.1)


models.loso[[paste0(study, '_LOSO')]] <- siamcat

cat('Successfully trained a LOSO model for study', study, '\n')

# Store trained model
save(models.loso, file=paste0(here('data','results','scv.loso','crc.loso.train/'),'Model.',study,'.CRC.LOSO.training.Rdata'))


# Evaluate LOSO Model on CRC Samples

loso.eval.crc<-list()

# load model
siamcat <- models.loso[[paste0(study, '_LOSO')]]

meta.test <- meta_balanced %>%
  filter(Cohort == study)

feat.test <- data_balanced[,(rownames(meta.test))]

label<-create.label(meta = meta.test,label='Condition',case='CRC',control = 'CTR')

siamcat.test <- siamcat(feat=feat.test,meta=meta.test,label=label, case='CRC')

siamcat.test <- make.predictions(siamcat, siamcat.holdout = siamcat.test)

siamcat.test<- evaluate.predictions(siamcat.test)

loso.eval.crc[[study]]<-siamcat.test

cat('LOSO model successfully tested on CRC for the study', study, '\n')

# Store and save results

save(loso.eval.crc, file=paste0(here('data','results','scv.loso','crc.loso.test/'),'Model.',study,'.CRĆ.LOSO.test.Rdata'))


# Evaluate LOSO Model on AD Samples

# Load complete metadata including all AD classes
meta.ad <- read_tsv(here('data','Metadata.all.samples.tsv')) %>%
  filter(Condition %in% c('AD', 'AdvAD', 'smallAD','CTR')) %>%
  as.data.frame() %>%
  column_to_rownames('Sample_ID') %>%
  mutate(Condition_general = case_when(
    Condition == 'AdvAD' ~ 'AD',
    Condition == 'smallAD' ~ 'AD',
    TRUE ~ Condition
  ))


ad.studies<- meta.ad %>% filter(Condition_general=='AD') %>% pull(Cohort) %>% unique()

all.data <- read.table(here('data','Relab.all.samples.tsv'),sep='\t',check.names = F)  %>%
  rownames_to_column('genus') %>%
  filter(genus!='unassigned') %>%
  column_to_rownames('genus')

all.data <- all.data[which(rowSums(all.data > 0) / ncol(all.data) > 0.1),]

loso.eval.ad<- loso.eval.small.ad <- loso.eval.adv.ad<- list()

if (study %in% ad.studies) {
  
  # Load model
  siamcat <- models.loso[[paste0(study, '_LOSO')]]
  
  # Filter test dataset
  meta.test <- meta.ad %>% filter(Cohort == study)
  
  # Process feature matrix
  feat.test <- all.data[, rownames(meta.test), drop = FALSE]
  
  # Ensure SIAMCAT does not fail due to missing training features
  tmp <- lapply(colnames(feat.test), function(x) {
    missing_feats <- norm_params(siamcat)$retained.feat[!norm_params(siamcat)$retained.feat %in% rownames(feat.test)]
    tmp <- data.frame(dummyColumnName = rep(0, length(missing_feats)))
    rownames(tmp) <- missing_feats
    colnames(tmp)[1] <- x
    return(tmp)
  })
  
  feat.test <- rbind(feat.test, as.data.frame(tmp, check.names = FALSE))
  
  # Evaluate on the AD dataset
  label <- create.label(meta = meta.test, label = 'Condition_general', case = 'AD', control = 'CTR')
  siamcat.test <- siamcat(feat = feat.test, meta = meta.test, label = label, case = 'AD')
  siamcat.test <- make.predictions(siamcat, siamcat.holdout = siamcat.test)
  siamcat.test <- evaluate.predictions(siamcat.test)
  
  # Store results
  loso.eval.ad[[study]] <- siamcat.test
  cat('LOSO model successfully tested on AD for the study', study, '\n')
  
  # Save results
  save(loso.eval.ad,  file = paste0(here('data','results','scv.loso','ad.loso.test/'),'Model.', study, '.AD.LOSO.test.Rdata'))
}



# Evaluate on smallAD Subtype

if ('smallAD' %in% meta.test$Condition) {
  
  # Filter and check dataset
  meta.test.small <- meta.test %>% filter(Condition %in% c('smallAD', 'CTR'))
  if (length(unique(meta.test.small$Condition)) == 2) {
    
    # Create label
    label <- create.label(meta = meta.test.small, label = 'Condition', case = 'smallAD', control = 'CTR')
    
    # Create SIAMCAT test object
    siamcat.test <- siamcat(feat = feat.test, meta = meta.test.small, label = label, case = 'smallAD')
    siamcat.test <- make.predictions(siamcat, siamcat.holdout = siamcat.test)
    siamcat.test <- evaluate.predictions(siamcat.test)
    
    # Store results
    loso.eval.small.ad[[study]] <- siamcat.test
    cat('LOSO model successfully tested on small AD for the study', study, '\n')
    
    # Save results
    save(loso.eval.small.ad,  file = paste0(here('data','results','scv.loso','ad.loso.test/'),'Model.', study, '.smallAD.LOSO.test.Rdata'))
    
  } else {
    cat('Skipping study', study, 'due to missing class in smallAD evaluation\n')
  }
}

# Evaluate on advanced AD Subtype

if ('AdvAD' %in% meta.test$Condition) {
  
  # Filter and check dataset
  meta.test.adv <- meta.test %>% filter(Condition %in% c('AdvAD', 'CTR'))
  if (length(unique(meta.test.adv$Condition)) == 2) {
    
    # Create label
    label <- create.label(meta = meta.test.adv, label = 'Condition', case = 'AdvAD', control = 'CTR')
    
    # Create SIAMCAT test object
    siamcat.test <- siamcat(feat = feat.test, meta = meta.test.adv, label = label, case = 'AdvAD')
    siamcat.test <- make.predictions(siamcat, siamcat.holdout = siamcat.test)
    siamcat.test <- evaluate.predictions(siamcat.test)
    
    # Store results
    loso.eval.adv.ad[[study]] <- siamcat.test
    cat('LOSO model successfully tested on AdvAD for the study', study, '\n')
    
    # Save results
    save(loso.eval.adv.ad,  file = paste0(here('data','results','scv.loso','ad.loso.test/'),'Model.', study, '.advAD.LOSO.test.Rdata'))
    
  } else {
    cat('Skipping study', study, 'due to missing class in advAD evaluation\n')
  }
}