##########################################
### Train EO-CRC and LO-CRC ml models for different cutoffs

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
library(ggembl)

# Load functions
source(here("src", "utils.R"))
load(here('data','results','Training.eo.lo.crc.rf.models.updated.Rdata'))

# Generate lists to save models and evaluations
models.eo.rf.list <- list()

# Parameters
age_cutoffs <- c(40, 45, 55)
age_col <- "Age"   

# Load metadata and data
all.meta <- read_tsv(here("data", "Metadata.all.samples.tsv")) %>%
  filter(Condition %in% c("CRC", "CTR")) %>%
  as.data.frame() %>%
  column_to_rownames("Sample_ID")

all.data <- read.table(here("data", "Relab.all.samples.tsv"),
                       sep = "\t", check.names = FALSE) %>%
  rownames_to_column("genus") %>%
  filter(genus != "unassigned") %>%
  column_to_rownames("genus")

for (cutoff in age_cutoffs) {
  message("=== Processing age cutoff: ", cutoff, " ===")
  
  ## 1) Define EO/LO groups based on age cutoff
  if(cutoff == 55){
    meta_age = all.meta %>%
      rownames_to_column("Sample_ID") %>%
      mutate(
        !!age_col := as.numeric(.data[[age_col]]),
        Onset = case_when( .data[[age_col]] < cutoff ~ "EO",
                           .data[[age_col]] >= cutoff ~ "LO",
                           Age_status == 'Early onset' ~ 'EO'),
        Group_age = paste0(Onset, "-", Condition)
      )
  }else{
    
    meta_age <- all.meta %>%
      rownames_to_column("Sample_ID") %>%
      mutate(
        !!age_col := as.numeric(.data[[age_col]]),
        Onset = if_else(.data[[age_col]] < cutoff, "EO", "LO"),
        Group_age = paste0(Onset, "-", Condition) # EO-CRC, EO-CTR, LO-CRC, LO-CTR
      )
    
  }
  
  # Balance EO and LO datasets per cutoff
  # EO: EO-CRC vs EO-CTR
  balanced_data_metadata_eocrc <- balance_data_by_groups(
    meta_df = meta_age %>%
      filter(Group_age %in% c("EO-CRC", "EO-CTR")) %>%
      group_by(Cohort) %>%
      filter(all(c("EO-CRC", "EO-CTR") %in% Group_age)) %>%
      ungroup(),
    data_df = all.data,
    group_column = "Group_age",
    groups = c("EO-CRC", "EO-CTR"),
    seed = 2002
  )
  
  # LO: LO-CRC vs LO-CTR
  balanced_data_metadata_locrc <- balance_data_by_groups(
    meta_df = meta_age %>%
      filter(Group_age %in% c("LO-CRC", "LO-CTR")) %>%
      group_by(Cohort) %>%
      filter(all(c("LO-CRC", "LO-CTR") %in% Group_age)) %>%
      ungroup() ,
    data_df = all.data,
    group_column = "Group_age",
    groups = c("LO-CRC", "LO-CTR"),
    seed = 2002
  )
  
  print(str_c(cutoff,':',balanced_data_metadata_eocrc$meta %>% count(Group_age)))
  ## 3) Train EO and LO models for this cutoff
  
  # EO model: EO-CRC vs EO-CTR
  models.eo.rf.list[[as.character(cutoff)]] <- train_model_rf(
    meta_df       = balanced_data_metadata_eocrc$meta,
    mat           = balanced_data_metadata_eocrc$data,
    label_column  = "Group",
    case_label    = "EO-CRC",
    control_label = "EO-CTR",
    block_label   = "Cohort",
    num_trees     = 200,
    seed          = 2000,
    prev.cutoff   = 0.1
  )
}

models.eo.rf.list[[as.character(50)]]<- models.eo.rf

##########################
# Test EO patients separated with different cutoffs using Lo data

siamcat_eo_tested <- list()
eval_eo_tested    <- list()

for (cutoff in age_cutoffs) {
  message("=== EO holdout using LO model, cutoff: ", cutoff, " ===")
  
  ## 1) Recreate EO/LO grouping for this cutoff
  if (cutoff == 55) {
    meta_age <- all.meta %>%
      rownames_to_column("Sample_ID") %>%
      mutate(
        !!age_col := as.numeric(.data[[age_col]]),
        Onset = case_when(
          Age_status == "Early onset"        ~ "EO",
          .data[[age_col]] <  cutoff         ~ "EO",
          .data[[age_col]] >= cutoff         ~ "LO",
          TRUE                               ~ NA_character_
        ),
        Group_age = paste0(Onset, "-", Condition)
      )
  } else {
    meta_age <- all.meta %>%
      rownames_to_column("Sample_ID") %>%
      mutate(
        !!age_col := as.numeric(.data[[age_col]]),
        Onset    = if_else(.data[[age_col]] < cutoff, "EO", "LO"),
        Group_age = paste0(Onset, "-", Condition)
      )
  }
  
  # Get EO-CRC vs EO-CTR for this cutoff
  balanced_data_metadata_eocrc <- balance_data_by_groups(
    meta_df = meta_age %>%
      filter(Group_age %in% c("EO-CRC", "EO-CTR")) %>%
      group_by(Cohort) %>%
      filter(all(c("EO-CRC", "EO-CTR") %in% Group_age)) %>%
      ungroup(),
    data_df     = all.data,
    group_column = "Group_age",
    groups       = c("EO-CRC", "EO-CTR"),
    seed         = 2002
  )
  
  # Sanity check
  print(balanced_data_metadata_eocrc$meta %>% count(Group_age))
  
  # Build siamcat object with EO samples
  label_eo <- create.label(
    meta    = balanced_data_metadata_eocrc$meta,
    label   = "Group_age",      # EO-CRC / EO-CTR are in Group_age
    case    = "EO-CRC",
    control = "EO-CTR"
  )
  
  siamcat_eo <- siamcat(
    feat  = balanced_data_metadata_eocrc$data,
    meta  = balanced_data_metadata_eocrc$meta,
    label = label_eo,
    case  = "EO-CRC"
  )
  
  # Apply the LO model to these EO samples
  siamcat_eo_tested[[as.character(cutoff)]] <- make.predictions(
    models.lo.rf,          # single LO model (50y cutoff)
    siamcat_eo             # EO test data for this cutoff
  )
  
  eval_eo_tested[[as.character(cutoff)]] <- evaluate.predictions(
    siamcat_eo_tested[[as.character(cutoff)]]
  )
  
  print(eval_eo_tested[[as.character(cutoff)]])
}

# Same evaluation without balancing 
for (cutoff in age_cutoffs) {
  message("=== EO holdout using LO model, cutoff: ", cutoff, " ===")
  
  ## 1) Recreate EO/LO grouping for this cutoff
  if (cutoff == 55) {
    meta_age <- all.meta %>%
      rownames_to_column("Sample_ID") %>%
      mutate(
        !!age_col := as.numeric(.data[[age_col]]),
        Onset = case_when(
          Age_status == "Early onset"        ~ "EO",
          .data[[age_col]] <  cutoff         ~ "EO",
          .data[[age_col]] >= cutoff         ~ "LO",
          TRUE                               ~ NA_character_
        ),
        Group_age = paste0(Onset, "-", Condition)
      )
  } else {
    meta_age <- all.meta %>%
      rownames_to_column("Sample_ID") %>%
      mutate(
        !!age_col := as.numeric(.data[[age_col]]),
        Onset    = if_else(.data[[age_col]] < cutoff, "EO", "LO"),
        Group_age = paste0(Onset, "-", Condition)
      )
  }
  
  # Get EO-CRC vs EO-CTR for this cutoff
  balanced_data_metadata_eocrc <- balance_data_by_groups(
    meta_df = meta_age %>%
      filter(Group_age %in% c("EO-CRC", "EO-CTR")) %>%
      group_by(Cohort) %>%
      filter(all(c("EO-CRC", "EO-CTR") %in% Group_age)) %>%
      ungroup(),
    data_df     = all.data,
    group_column = "Group_age",
    groups       = c("EO-CRC", "EO-CTR"),
    seed         = 2002
  )
  
  # Build siamcat object with EO samples
  label_eo <- create.label(
    meta    = balanced_data_metadata_eocrc$meta,
    label   = "Group_age",      # EO-CRC / EO-CTR are in Group_age
    case    = "EO-CRC",
    control = "EO-CTR"
  )
  
  siamcat_eo <- siamcat(
    feat  = balanced_data_metadata_eocrc$data,
    meta  = balanced_data_metadata_eocrc$meta,
    label = label_eo,
    case  = "EO-CRC"
  )
  
  # Apply the LO model to these EO samples
  siamcat_eo_tested[[as.character(cutoff)]] <- make.predictions(
    models.lo.rf,          # single LO model (50y cutoff)
    siamcat_eo             # EO test data for this cutoff
  )
  
  eval_eo_tested[[as.character(cutoff)]] <- evaluate.predictions(
    siamcat_eo_tested[[as.character(cutoff)]]
  )
  
  print(eval_eo_tested[[as.character(cutoff)]])
}

eval_eo_tested[[as.character(50)]]<- siamcat.test.evaluated.eo.holdout.rf

# Save models and evaluations
save( 
  models.eo.rf.list,
  eval_eo_tested,
  file = here("data", "results", "Testing.eo.crc.by.age.cutoffs.updated.Rdata"))

