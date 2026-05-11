##########################################
# Train single study models for CRC vs CTR and test them across studies (SCV-SST)

# Load functions
source(here('src','requirements.R'))
source(here('src','utils.R'))

# Load metadata and data
all.meta<- read_tsv(here('data','Metadata.all.samples.balanced.tsv')) %>% filter(Condition=='CRC'| Condition=='CTR') %>% as.data.frame() %>% column_to_rownames('Sample_ID')

all.data <- read.table(here('data', 'Relab.all.samples.balanced.tsv'),sep='\t', check.names = F) %>%
  rownames_to_column('genus') %>% filter(genus!='unassigned') %>% column_to_rownames('genus')


########################################################################
# Train a single study models for each cohort  

studies<- unique(all.meta$Cohort)

models.scv <- list()

for (study in (studies)){
  # single study model
  meta.train <- all.meta %>%
    filter(Cohort == study)
  
  feat.train <- all.data[,(rownames(meta.train))]
  
  siamcat<-train_model_rf(meta_df = meta.train, mat = feat.train,label_column = 'Condition',case_label = 'CRC',control_label = 'CTR',
                          seed = 2000,block_label = NULL,num_trees = 200, prev.cutoff = 0.1)
  
  models.scv[[study]] <- siamcat
  
  cat('Successfully trained a single study model for study', study, '\n')
  

}

# Make study to study transfer predictions 

holdout.evaluation.crc<- evaluation<- list()

for (study in (studies)){
  
  # load model 
  siamcat <- models.scv[[study]]
  
  for (i in (studies)){
    
    if (i!= study) {
      
      meta.test <- all.meta %>%
        filter(Cohort == i)
      feat.test <- all.data[,(rownames(meta.test))]
      
      
      label<-create.label(meta = meta.test,label='Condition',case='CRC',control = 'CTR')
      
      siamcat.test <- siamcat(feat=feat.test, meta=meta.test,
                              label=label, case='CRC')
      
      siamcat.test.predicted<- make.predictions(siamcat , siamcat.test)
      
      siamcat.test.evaluated<-evaluate.predictions(siamcat.test.predicted)
      
      evaluation[i] <- siamcat.test.evaluated
      
      
      cat('Successfully' , paste0(study), 'tested using', paste0(i) ,'model', '\n')
    }
    else
    {
      evaluation[i]<-NA
    }
    
    holdout.evaluation.crc[[study]]<-evaluation
    
  }
}

holdout.evaluation.crc  <- holdout.evaluation.crc %>%
  map(~ .x[!is.na(.x)]) %>%
  compact()


save(holdout.evaluation.crc, models.scv, file= here('data', 'results','scv.loso','crc.scv.sst', 'Models.SCV.SST.Rdata'))

