######################
# Train Adenoma ml model
######################

# Load functions
source(here('src','utils.R'))


# Load metadata and data

meta.ad <- read_tsv(here('data','Metadata.all.samples.tsv')) %>% 
  filter(Condition %in% c('AD', 'AdvAD', 'smallAD','CTR')) %>% 
  as.data.frame() %>% 
  column_to_rownames('Sample_ID') %>% 
  mutate(Condition_general = case_when(
    Condition == 'AdvAD' ~ 'AD', 
    Condition == 'smallAD' ~ 'AD',
    TRUE ~ Condition
  )) 

selected_cohort<- meta.ad  %>% filter(Condition_general=='AD') %>% pull(Cohort) %>% unique()

meta.ad <- meta.ad %>% filter(Cohort %in% selected_cohort) %>% rownames_to_column('Sample_ID')

all.data <- read.table(here('data','Relab.all.samples.tsv'),sep='\t', check.names = F) %>%
  rownames_to_column('genus') %>% 
  filter(genus!='unassigned') %>%
  column_to_rownames('genus')


models.ad <-train_model_rf(meta_df = meta.ad %>% column_to_rownames('Sample_ID')  ,mat = all.data ,label_column = 'Condition_general',
                           case_label = 'AD',control_label = 'CTR')

save(models.ad, file= here('data','results','Training.ad.ctr.rf.model.Rdata'))


