##########################################
### Train 16S and WGS ml models
# Load packages

library(dplyr)
library(tidyverse)
library(stringr)
library(ranger)
library(mlr3)
library(mlr3tuning)
library(mlr3extralearners)
library(ggembl)
library(SIAMCAT)
library(here)

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
     file='/g/scb/zeller/pekel/meta_analysis/src/analysis/Rdata/Training.all.16s.wgs.data.rf.Rdata')

# Plot ROC Curve for the models 
plot_roc_siamcat_models <- function(models, labels, colours, linetypes, trained_on = NULL, alpha = NULL) {
  
  determine_tpr_fpr_auc <- function(eval_data, auroc) {
    tpr_list <- list()
    fpr_list <- list()
    
    for (i in seq_along(eval_data)) {
      tp <- eval_data[[i]]$tp
      fp <- eval_data[[i]]$fp
      tn <- eval_data[[i]]$tn
      fn <- eval_data[[i]]$fn
      
      if (is.null(tp) || is.null(fp) || is.null(tn) || is.null(fn)) {
        next
      }
      
      tpr <- tp / (tp + fn)
      fpr <- fp / (fp + tn)
      
      
      tpr_list[[i]] <- tpr
      fpr_list[[i]] <- fpr
    }
    
    if (length(tpr_list) == 0 || length(fpr_list) == 0) {
      stop("No valid evaluation data found.")
    }
    
    max_len <- max(sapply(tpr_list, length), sapply(fpr_list, length))
    tpr_matrix <- do.call(rbind, lapply(tpr_list, function(x) c(x, rep(NA, max_len - length(x)))))
    fpr_matrix <- do.call(rbind, lapply(fpr_list, function(x) c(x, rep(NA, max_len - length(x)))))
    
    mean_tpr <- apply(tpr_matrix, 2, mean, na.rm = TRUE)
    mean_fpr <- apply(fpr_matrix, 2, mean, na.rm = TRUE)
    
    roc_data <- data.frame(FPR = mean_fpr, TPR = mean_tpr)
    return(list(roc_data = roc_data, auc = auroc))
  }
  
  all_roc_data <- data.frame(FPR = numeric(), TPR = numeric(), model = character())
  auc_list <- list()
  
  for (i in seq_along(models)) {
    cat("Processing model:", labels[i], "\n")
    
    eval_data <- models[[i]]@eval_data$ev.all
    auroc <- models[[i]]@eval_data$auroc
    
    if (is.null(eval_data) || is.null(auroc)) {
      warning(paste("Model", labels[i], "has no valid evaluation data or AUC. Skipping."))
      next
    }
    
    if (!is.null(trained_on) && !is.null(trained_on[[i]])) {
      training_sizes <- unlist(sapply(trained_on[[i]]@data_split$training.folds, function(fold) sapply(fold, length)))
      mean_train_size <- ifelse(length(training_sizes) > 0, round(mean(training_sizes, na.rm = TRUE)), NA)
      mean_test_size <- round(dim(models[[i]]@pred_matrix)[1])
    } else {
      training_sizes <- unlist(sapply(models[[i]]@data_split$training.folds, function(fold) sapply(fold, length)))
      test_sizes <- unlist(sapply(models[[i]]@data_split$test.folds, function(fold) sapply(fold, length)))
      mean_train_size <- ifelse(length(training_sizes) > 0, round(mean(training_sizes, na.rm = TRUE)), NA)
      mean_test_size <- ifelse(length(test_sizes) > 0, round(mean(test_sizes, na.rm = TRUE)), NA)
    }
    
    result <- determine_tpr_fpr_auc(eval_data, auroc)
    roc_data <- result$roc_data
    auc <- result$auc
    
    roc_data$model <- labels[i]
    all_roc_data <- rbind(all_roc_data, roc_data)
    
    auc_list[[i]] <- paste0(labels[i], " (AUC = ", round(auc, 2), ", Train N = ", mean_train_size, ", Test N = ", mean_test_size, ")")
  }
  
  if (nrow(all_roc_data) == 0) {
    stop("No valid data for ROC plot. Please check your models.")
  }
  
  # Ensure `alpha` is set
  if (is.null(alpha)) {
    alpha <- 1
  }
  
  # Ensure `linetypes` has the correct length
  if (length(linetypes) != length(labels)) {
    stop("Error: `linetypes` must have the same length as `labels`.")
  }
  
  # Ensure `colours` matches `labels`
  if (length(colours) != length(labels)) {
    stop("Error: `colours` must have the same length as `labels`.")
  }
  
  # Ensure 'model' is a factor with the correct level order
  all_roc_data$model <- factor(all_roc_data$model, levels = labels)
  
  labels <- unique(all_roc_data$model) 
  
  # Format labels for better readability
  wrapped_labels <- sub(", ", ",\n", unlist(auc_list))
  
  #correct_labels <- unique(all_roc_data$model)
  
  # Assign linetypes to models
  #linetype_mapping <- setNames(linetypes, correct_labels)
  

  
  # Plot
  p<-ggplot(all_roc_data, aes(x = FPR, y = TPR, color = model, linetype = model)) +
    geom_line(size = 1.8, alpha = alpha) +
    scale_color_manual(values = setNames(colours, labels), labels = wrapped_labels, name = "Model") +  
    scale_linetype_manual(values = setNames(linetypes, labels), labels = wrapped_labels, name = "Model") +  
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey66", size = 0.8) +
    labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", color = "Model", linetype = "Model") +  
    theme_paper +
    theme(
      legend.position = c(0.25, 0.05),
      legend.justification = c(0, 0),
      legend.box = "horizontal",
      legend.box.background = element_rect(color = "black", size = 0.8),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_text(size = 11)
    ) +
    coord_fixed() +
    guides(color = guide_legend(override.aes = list(linetype = linetypes)))  # Combine legends
  
  
  return(p)
}


models <- list(models.rf.16S, siamcat.test.evaluated.16s.holdout.rf ,models.rf.wgs,  siamcat.test.evaluated.wgs.holdout.rf)
labels <- c("Classifier cross validated on 16S",
            "Classifier trained on WGS and tested on 16S",
            "Classifier cross validated on WGS",
            "Classifier trained on 16S and tested on WGS")
trained_on <- list(NULL,models.rf.wgs, NULL,  models.rf.16S)
colours <- c( 'brown','chocolate4', 'cyan3', 'cyan4' )
linetypes<- c('solid','solid','solid','solid')


cross_tech_models_auc_plot<- plot_roc_siamcat_models(models, labels, colours, trained_on, alpha=0.6, linetypes=linetypes)











