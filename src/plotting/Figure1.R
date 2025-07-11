######################
# Figure 1 
######################

source(here('src','utils.R'))
params <- yaml::read_yaml(here("src", "parameters.yml"))
plotting <- params$plotting

# Load LMM tables for 16S and WGS
load(here('data','results', 'lmm.tables.16S.WGS.Rdata'))


######################
# Figure 1b: Volcano plot for WGS

volcano_16S <- plot_volcano(
  plot_df = lmm.table.16S %>% select(Taxon, P.val, P.adj, Effect.size, pr.shift, pr.CRC, pr.CTR, n.CTR, n.CRC),
  fdr_thresh = 0.05,
  group_case = 'CRC',
  group_control = 'CTR',
  min_segment_length= 0.2,
  max.overlaps=20,
  feature_column_name ='Taxon',
  color_vector = c(group_case =  plotting$condition_colors$CRC, group_control = plotting$condition_colors$CTR, "n.s." = "white")) +
  xlab('16S enrichment effect size')


ggsave(volcano_16S, file=here('figures','figure1','Figure1b.pdf'), height = 5, width = 5)

######################
# Figure 1c: Volcano plot for WGS

volcano_wgs <- plot_volcano(
  plot_df = lmm.table.wgs %>% select(Taxon, P.val, P.adj, Effect.size, pr.shift,pr.CRC,pr.CTR,  n.CTR, n.CRC),
  fdr_thresh = 0.05, 
  group_case = 'CRC', 
  group_control = 'CTR',
  min_segment_length=0.2,
  feature_column_name = 'Taxon') +
  xlab('WGS enrichment effect size')

ggsave(volcano_wgs, file=here('figures','figure1','Figure1c.pdf'), height = 5, width = 5)


# Combine LMM tables for scatter plot 
lmm.table.all<- left_join(lmm.table.16S %>% select(Taxon, P.val, P.adj, Effect.size, pr.shift,pr.CRC,pr.CTR) ,
                           lmm.table.wgs %>% select(Taxon, P.val, P.adj, Effect.size, pr.shift,pr.CRC,pr.CTR), 
                           by='Taxon', suffix = c('.16S','.WGS')) 

corr <- cor(
  lmm.table.all$Effect.size.WGS,
  lmm.table.all$Effect.size.16S,
  method = "pearson",
  use = "complete.obs"
)

# Format the label
corr_label <- paste0("Pearson r = ", round(corr, 2))

# Compute jaccard index

jaccard_table <- lmm.table.all  %>%
  mutate(
    enriched_in = case_when(
      # Both tests show significant enrichment
      (P.adj.16S < 0.05 & Effect.size.16S > 0) & 
        (P.adj.WGS < 0.05 & Effect.size.WGS > 0) ~ 'both enriched',
      
      # Both tests show significant depletion
      (P.adj.16S < 0.05 & Effect.size.16S < 0) & 
        (P.adj.WGS < 0.05 & Effect.size.WGS < 0) ~ 'both depleted',
      
      # Both tests non-significant
      P.adj.16S >= 0.05 & P.adj.WGS >= 0.05 ~ 'both n.s.',
      
      # disagreements or partial significance
      TRUE ~ "one sig"
    ),
    label = ifelse(P.adj.16S < 0.05 | P.adj.WGS < 0.05, Taxon, ""),
    font = ifelse(P.adj.16S < 0.05 | P.adj.WGS < 0.05, "bold.italic", "italic")
  ) %>% 
  select(Taxon, P.adj.16S, P.adj.WGS, Effect.size.WGS, Effect.size.16S, enriched_in)

jaccard_index <- jaccard_table %>% 
  summarise(
    numerator = sum(enriched_in %in% c("both enriched", "both depleted")),
    denominator = n() - sum(enriched_in == "both n.s."),  
    jaccard = numerator / denominator
  ) %>%
  pull(jaccard)

######################
# Figure 1d: Scatter plot
scatter <-plot_comparison_scatter(
  data = lmm.table.all,
  x_col = "Effect.size.WGS",
  y_col = "Effect.size.16S",
  x_label = "Effect Size (16S)",
  y_label = "Effect Size (WGS)",
  feature_column_name = "Taxon") +
  scale_x_continuous(limits = c(-0.5,1),breaks = c(-0.5, -0.25, 0, 0.25, 0.5,0.75,1)) +
  scale_y_continuous(limits=c(-0.5,1),breaks = c(-0.5, -0.25, 0, 0.25, 0.5,0.75,1)) 

ggsave(scatter, file=here('figures','figure1','Figure1d.pdf'), height=5, width=5)

######################
# Figure 1e: AUROC for assay comparison
load(here('data','results','Training.16s.wgs.rf.models.Rdata'))


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
  
  if (is.null(alpha)) {
    alpha <- 1
  }
  
  if (length(linetypes) != length(labels)) {
    stop("Error: `linetypes` must have the same length as `labels`.")
  }

  if (length(colours) != length(labels)) {
    stop("Error: `colours` must have the same length as `labels`.")
  }
  
 
  all_roc_data$model <- factor(all_roc_data$model, levels = labels)
  
  labels <- unique(all_roc_data$model) 
  
  wrapped_labels <- sub(", ", ",\n", unlist(auc_list))
  
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

ggsave(cross_tech_models_auc_plot, file=here('figures','figure1','Figure1e.pdf'), height=7, width=7)


