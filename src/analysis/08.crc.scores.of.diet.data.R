######################
# CRC scores of diet data
######################
# This script calculates the CRC scores for the datasets used in figure 5 and supplementary figure 5

siamcat_lib <- '/g/scb2/zeller/karcher/mambaforge/envs/siamcat/lib/R/library/'

load_pkg <- function(pkg, fallback_lib = NULL, strict = TRUE) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    library(pkg, character.only = TRUE)
    return(invisible(TRUE))
  }
  if (!is.null(fallback_lib) && dir.exists(fallback_lib) &&
      requireNamespace(pkg, quietly = TRUE, lib.loc = fallback_lib)) {
    library(pkg, character.only = TRUE, lib.loc = fallback_lib)
    return(invisible(TRUE))
  }
  if (strict) {
    stop(sprintf("Package '%s' could not be loaded from active libraries or fallback path.", pkg))
  }
  return(invisible(FALSE))
}

load_pkg('here')
load_pkg('ggplot2')
load_pkg('dplyr')
load_pkg('stringr')
load_pkg('pROC')
load_pkg('ggembl', fallback_lib = '/g/scb2/zeller/pekel/R/CRC.meta.libPath/')
load_pkg('tidyverse')
load_pkg('ranger')
load_pkg('ggforce')
load_pkg('ggbeeswarm')
load_pkg('mlr3', fallback_lib = siamcat_lib)
load_pkg('mlr3learners', fallback_lib = siamcat_lib)
load_pkg('mlr3extralearners', fallback_lib = siamcat_lib)
load_pkg('SIAMCAT', fallback_lib = siamcat_lib)
source(here('src','cancerness.utils.R'))
load_pkg('readxl')
load_pkg('rstatix') 
load_pkg('ggpubr')

params <- yaml::read_yaml(here("src", "parameters.yml"))
plotting <- params$plotting

mOTUS.3.1.metadata<- read_tsv(here('data','mOTUs3.1.0.genome_metadata_edited.tsv'))
# Load the universal CRC model

load(here('data','results','Training.unified.crc.model.Rdata'))

# Load SowahSA_2022 data
Sowah_relab <-readRDS(here('data','fiber_data','SowahSA_2022_res_IDTaxa.rds')) %>% 
  as.data.frame() %>%
  rownames_to_column('taxon') %>%
  pivot_longer(-taxon) %>%
  mutate(genus = str_split_fixed(taxon, "\\|", n = 7)[, 6],
         genus = str_replace(genus, "g__", "")) %>%
  dplyr::rename(sampleID = name) %>%
  mutate(genus = ifelse(genus == "", "unassigned", genus)) %>% 
  group_by(sampleID,genus) %>%
  summarise(value=sum(value)) %>% 
  pivot_wider(names_from = sampleID,values_from = value) %>%
  column_to_rownames('genus') %>%
  mutate(across(everything(), ~ . / sum(.)))


# Load metadata
Sowah_meta <- read.table(here('data','fiber_data','SraRunTable_Sowah.csv'), sep=',', header = T) %>%
  select(Run,Sample.Name, Sample_name) %>%
  mutate(Time_points=  str_split_fixed(Sample_name,'_',4)[,2],
         Batch=  str_split_fixed(Sample_name,'_',4)[,3],
         Patient_ID = str_split_fixed(Sample_name,'_',4)[,1],
         Identifier = str_split_fixed(Sample_name,'_',4)[,4]) 

# assign CRC scores to the Sowah samples
dummyControls <- data.frame(dummySample = rpois(n = dim(Sowah_relab)[1], lambda = 10))


rownames(dummyControls) <- rownames(Sowah_relab)
testProfile <- cbind(Sowah_relab, dummyControls)

# Pad missing features in testProfile to match model's retained features
tmp <- lapply(colnames(testProfile), function(x) {
  miss <- norm_params(models.all.rf)$retained.feat[!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile)]
  tmp <- data.frame(dummyColumnName = rep(0, length(miss)))
  rownames(tmp) <- miss
  colnames(tmp)[1] <- x
  return(tmp)
})

testProfile <- rbind(testProfile, as.data.frame(tmp, check.names = F))
dummyMeta <- data.frame(label = c(rep("CRC", dim(testProfile)[2] - 1 ), "Control"))
rownames(dummyMeta) <- colnames(testProfile)

Sowah_score_plot <- get_cancerness(siam = models.all.rf, pro = testProfile, meta = dummyMeta, evaluated.label='SowahSA_2022')

Sowah_scores <- Sowah_score_plot[[2]] %>%  
  filter(type=='evaluated') 

# DeFilippis_2016 data

DeFilippis_relab<- readRDS(here('data','fiber_data','DeFilippis_2016-red_IDTaxa.rds')) %>% 
  as.data.frame() %>%
  rownames_to_column('taxon') %>%
  pivot_longer(-taxon) %>%
  mutate(genus = str_split_fixed(taxon, "\\|", n = 7)[, 6],
         genus = str_replace(genus, "g__", "")) %>%
  dplyr::rename(sampleID = name) %>%
  mutate(sampleID= str_remove(sampleID, '.singles')) %>%
  mutate(genus = ifelse(genus == "", "unassigned", genus)) %>% 
  group_by(sampleID,genus) %>%
  summarise(value=sum(value)) %>% 
  pivot_wider(names_from = sampleID,values_from = value) %>%
  column_to_rownames('genus') %>%
  mutate(across(everything(), ~ . / sum(.)))


# assign CRC scores to the Baxter samples
dummyControls <- data.frame(dummySample = rpois(n = dim(DeFilippis_relab)[1], lambda = 10))


rownames(dummyControls) <- rownames(DeFilippis_relab)
testProfile <- cbind(DeFilippis_relab, dummyControls)

tmp <- lapply(colnames(testProfile), function(x) {
  tmp <- data.frame(dummyColumnName = rep(0, sum(!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile))))
  rownames(tmp) <- norm_params(models.all.rf)$retained.feat[!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile)]
  colnames(tmp)[1] <- x
  return(tmp)
})

testProfile <- rbind(testProfile, as.data.frame(tmp, check.names = F))
dummyMeta <- data.frame(label = c(rep("CRC", dim(testProfile)[2] - 1 ), "Control"))
rownames(dummyMeta) <- colnames(testProfile)

DeFilippis_score_plot <- get_cancerness(siam = models.all.rf, pro = testProfile, meta = dummyMeta, evaluated.label='DeFilippis_2016')

DeFilippis_meta<- read.table(here('data','fiber_data','SraRunTable-DeFilippis_2016.csv'),sep=',', header = T) %>% 
  filter(env_biome=='feces') %>% 
  select(Run, Label, Sample.Name) %>% 
  mutate(host_diet = str_split_fixed(Label, " subject", 2)[,1],
         host_diet = str_split_fixed(host_diet, "from ",2) [,2]) %>% 
  mutate(host_diet_general = ifelse(host_diet=='omnivorous','omnivorous','vegetarian/vegan'))


DeFilippis_scores <- DeFilippis_score_plot[[2]] %>%  
  filter(type == "evaluated") 

# BarberC_2021 data

Barber_relab <- readRDS(here('data','fiber_data','Barber_2021_res_mOTUs.rds')) %>%
  as.data.frame() %>%
  rownames_to_column('taxon') %>%
  pivot_longer(-taxon) %>%
  group_by(name) %>%
  mutate(mOTUs = str_split_fixed(taxon, '\\|', n = 8)[, 8]) %>%
  mutate(mOTUs= ifelse(taxon=='unassigned','unassigned',mOTUs)) %>% 
  dplyr::rename(sampleID = name) %>%
  group_by(sampleID, mOTUs) %>%
  filter(mOTUs != '') %>%
  pivot_wider(id_cols = mOTUs, names_from = sampleID, values_from = value) %>%
  as.data.frame %>%
  pivot_longer(-mOTUs) %>% 
  dplyr::rename('sampleID'= 'name', 'count'='value') %>% 
  left_join(mOTUS.3.1.metadata, by= c('mOTUs' = 'GENOME_SOURCE')) %>% 
  filter(DOMAIN!='Unclassified') %>% 
  dplyr::rename('genus'= 'GTDB_genus') %>% 
  mutate(genus=ifelse(mOTUs=='unassigned', 'unassigned',genus)) %>% 
  filter(mOTUs!='meta_mOTU_v31_12968') %>% # these specific motus don't have any entry in mOTUs3.1.0.genome_metadata.tsv 
  select(-DOMAIN) %>% 
  group_by(sampleID, genus) %>% 
  dplyr::summarise(count=sum(count)) %>% 
  pivot_wider(id_cols = genus , names_from = sampleID, values_from = count, values_fill = 0)  %>% 
  filter(genus!='unassigned') %>% 
  column_to_rownames('genus')  %>%
  mutate(across(everything(), ~ . / sum(.))) 

# Calculate CRC scores

dummyControls <- data.frame(dummySample = rpois(n = dim(Barber_relab)[1], lambda = 10))

rownames(dummyControls) <- rownames(Barber_relab)
testProfile <- cbind(Barber_relab, dummyControls)

tmp <- lapply(colnames(testProfile), function(x) {
  tmp <- data.frame(dummyColumnName = rep(0, sum(!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile))))
  rownames(tmp) <- norm_params(models.all.rf)$retained.feat[!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile)]
  colnames(tmp)[1] <- x
  return(tmp)
})

testProfile <- rbind(testProfile, as.data.frame(tmp, check.names = F))
dummyMeta <- data.frame(label = c(rep("CRC", dim(testProfile)[2] - 1 ), "Control"))
rownames(dummyMeta) <- colnames(testProfile)

Barber_score_plot <- get_cancerness(siam = models.all.rf, pro = testProfile, meta = dummyMeta, evaluated.label='BarberC_2021')

Barber_scores<- Barber_score_plot[[2]] %>% filter(type=='evaluated')

# Delannoy_Bruno_2021 data

Delannoy_relab <- readRDS(here('data','fiber_data','Delannoy-Bruno_2021_res_mOTUs.rds')) %>%
  as.data.frame() %>%
  rownames_to_column('taxon') %>%
  pivot_longer(-taxon) %>%
  group_by(name) %>%
  mutate(mOTUs = str_split_fixed(taxon, '\\|', n = 8)[, 8]) %>%
  mutate(mOTUs= ifelse(taxon=='unassigned','unassigned',mOTUs)) %>% 
  dplyr::rename(sampleID = name) %>%
  group_by(sampleID, mOTUs) %>%
  filter(mOTUs != '') %>%
  pivot_wider(id_cols = mOTUs, names_from = sampleID, values_from = value) %>%
  as.data.frame %>%
  pivot_longer(-mOTUs) %>% 
  dplyr::rename('sampleID'= 'name', 'count'='value') %>% 
  left_join(mOTUS.3.1.metadata, by= c('mOTUs' = 'GENOME_SOURCE')) %>% 
  filter(DOMAIN!='Unclassified') %>% 
  dplyr::rename('genus'= 'GTDB_genus') %>% 
  mutate(genus=ifelse(mOTUs=='unassigned', 'unassigned',genus)) %>% 
  filter(mOTUs!='meta_mOTU_v31_12968') %>% # these specific motus don't have any entry in mOTUs3.1.0.genome_metadata.tsv 
  select(-DOMAIN) %>% 
  group_by(sampleID, genus) %>% 
  dplyr::summarise(count=sum(count)) %>% 
  pivot_wider(id_cols = genus , names_from = sampleID, values_from = count, values_fill = 0)  %>% 
  filter(genus!='unassigned') %>% 
  column_to_rownames('genus')  %>%
  mutate(across(everything(), ~ . / sum(.))) 


# Calculate CRC scores
dummyControls <- data.frame(dummySample = rpois(n = dim(Delannoy_relab)[1], lambda = 10))
rownames(dummyControls) <- rownames(Delannoy_relab)
testProfile <- cbind(Delannoy_relab, dummyControls)

tmp <- lapply(colnames(testProfile), function(x) {
  tmp <- data.frame(dummyColumnName = rep(0, sum(!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile))))
  rownames(tmp) <- norm_params(models.all.rf)$retained.feat[!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile)]
  colnames(tmp)[1] <- x
  return(tmp)
})


testProfile <- rbind(testProfile, as.data.frame(tmp, check.names = F))
dummyMeta <- data.frame(label = c(rep("CRC", dim(testProfile)[2] - 1 ), "Control"))
rownames(dummyMeta) <- colnames(testProfile)

Delannoy_score_plot <- get_cancerness(siam = models.all.rf, pro = testProfile, meta = dummyMeta, evaluated.label='Delannoy-Bruno_2021')

Delannoy_scores<- Delannoy_score_plot[[2]] %>% 
  filter(type=='evaluated')

# HealeyG_2018 data

Healey_meta <- read.table(here("data", "fiber_data", "HealeyG_2018_metadata.txt"), sep='\t', header = T)  %>% 
  mutate(Run=str_split_fixed(sampleid, '_',2)[,2])

Healey_sra <- read.table(here("data", "fiber_data", "SraRunTable-Healey.csv"), sep=',', header = T)  %>% select(Run,  Host_Diet)

Healey_meta<- Healey_meta %>% left_join(Healey_sra, by='Run')

Healey_relab <-readRDS(here('data','fiber_data','HealeyG_2018_res_IDTaxa.rds')) %>% 
  as.data.frame() %>%
  rownames_to_column('taxon') %>%
  pivot_longer(-taxon) %>%
  mutate(genus = str_split_fixed(taxon, "\\|", n = 7)[, 6],
         genus = str_replace(genus, "g__", "")) %>%
  dplyr::rename(sampleID = name) %>%
  mutate(genus = ifelse(genus == "", "unassigned", genus)) %>% 
  group_by(sampleID,genus) %>%
  summarise(value=sum(value)) %>% 
  pivot_wider(names_from = sampleID,values_from = value) %>%
  column_to_rownames('genus') %>%
  mutate(across(everything(), ~ . / sum(.)))


dummyControls <- data.frame(dummySample = rpois(n = dim(Healey_relab)[1], lambda = 10))

rownames(dummyControls) <- rownames(Healey_relab)
testProfile <- cbind(Healey_relab, dummyControls)

tmp <- lapply(colnames(testProfile), function(x) {
  tmp <- data.frame(dummyColumnName = rep(0, sum(!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile))))
  rownames(tmp) <- norm_params(models.all.rf)$retained.feat[!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile)]
  colnames(tmp)[1] <- x
  return(tmp)
})

testProfile <- rbind(testProfile, as.data.frame(tmp, check.names = F))
dummyMeta <- data.frame(label = c(rep("CRC", dim(testProfile)[2] - 1 ), "Control"))
rownames(dummyMeta) <- colnames(testProfile)

Healey_score_plot <- get_cancerness(siam = models.all.rf, pro = testProfile, meta = dummyMeta, evaluated.label='HealeyG_2018')

Healey_score <- Healey_score_plot[[2]] %>%  
  filter(type == "evaluated") 

# Delannoy_Bruno_2022 data
Delannoy_2022 <- readRDS(here('data','fiber_data','Delannoy-Bruno_2022_res_mOTUs.rds')) %>%
  as.data.frame() %>%
  rownames_to_column('taxon') %>%
  pivot_longer(-taxon) %>%
  group_by(name) %>%
  mutate(motus = str_split_fixed(taxon, '\\|', n = 8)[, 8]) %>%
  mutate(motus= ifelse(taxon=='unassigned','unassigned',motus)) %>% 
  dplyr::rename(sampleID = name) %>%
  group_by(sampleID, motus) %>%
  filter(motus != '') %>%
  pivot_wider(id_cols = motus, names_from = sampleID, values_from = value) %>%
  pivot_longer(-motus) %>% 
  dplyr::rename('sampleID'= 'name', 'count'='value') %>% 
  left_join(mOTUS.3.1.metadata, by= c('motus' = 'GENOME_SOURCE')) %>% 
  filter(DOMAIN!='Unclassified') %>% 
  dplyr::rename('genus'= 'GTDB_genus') %>% 
  mutate(genus=ifelse(motus=='unassigned', 'unassigned',genus)) %>% 
  filter(motus!='meta_mOTU_v31_12968') %>% # these specific motus don't have any entry in mOTUs3.1.0.genome_metadata.tsv 
  select(-DOMAIN) %>% 
  group_by(sampleID, genus) %>% 
  dplyr::summarise(count=sum(count)) %>% 
  pivot_wider(id_cols = genus , names_from = sampleID, values_from = count, values_fill = 0)  %>% 
  filter(genus!='unassigned') %>% 
  column_to_rownames('genus')  %>%
  mutate(across(everything(), ~ . / sum(.)))  %>% 
     as.data.frame() %>%
    select(where(~!all(is.na(.))))
 
# Calculate CRC scores
dummyControls <- data.frame(dummySample = rpois(n = dim(Delannoy_2022)[1], lambda = 10))
rownames(dummyControls) <- rownames(Delannoy_2022)
testProfile <- cbind(Delannoy_2022, dummyControls)

tmp <- lapply(colnames(testProfile), function(x) {
  tmp <- data.frame(dummyColumnName = rep(0, sum(!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile))))
  rownames(tmp) <- norm_params(models.all.rf)$retained.feat[!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile)]
  colnames(tmp)[1] <- x
  return(tmp)
})

testProfile <- rbind(testProfile, as.data.frame(tmp, check.names = F))
dummyMeta <- data.frame(label = c(rep("CRC", dim(testProfile)[2] - 1 ), "Control"))
rownames(dummyMeta) <- colnames(testProfile)

Delannoy_score_plot2 <- get_cancerness(siam = models.all.rf, pro = testProfile, meta = dummyMeta, evaluated.label='Delannoy-Bruno_2022')

Delannoy_scores2<- Delannoy_score_plot2[[2]] %>% 
  filter(type=='evaluated')

# Now combine all scores into one table
all_scores <- bind_rows(Sowah_scores, DeFilippis_scores, Barber_scores, Delannoy_scores, Healey_score, Delannoy_scores2)

# Save the scores
write_tsv(all_scores, here('data','results','CRC_scores_fiber_datasets.tsv'))
