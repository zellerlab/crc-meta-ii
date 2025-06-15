
.libPaths(c(
  "/g/scb/zeller/pekel/R/CRC.meta.libPath",
  "/g/easybuild/x86_64/Rocky/8/haswell/software/R-bundle-Bioconductor/3.19-foss-2023b-R-4.4.1",
  "/g/easybuild/x86_64/Rocky/8/haswell/software/arrow-R/16.1.0-foss-2023b-R-4.4.1",
  "/g/ea§sybuild/x86_64/Rocky/8/haswell/software/R/4.4.1-gfbf-2023b",
  "/g/scb/zeller/pekel/Conda/miniconda3/envs/crc_meta/lib/R/library",
  "/g/scb/zeller/pekel/R/4.4.1/"
))



#Load libraries 
require(tidyverse)
require(ape)
require(RColorBrewer)
require(ggrepel)
require(labdsv)
require(dplyr)
require(tidyverse)
require(vegan)
require(ggembl)
require(lmerTest)
require(lme4)
require(here)
library(ggplot2)
library(stringr)



theme_paper <- ggembl::theme_presentation() +
  theme(
    axis.title = element_text(face = "bold", size = 12), # Bold axis titles
    panel.border = element_rect(fill=NA, colour='black', size=1.5),
    #axis.line = element_line(size = 1.5), # Change line thickness of axes
    
    axis.text = element_text(face = "bold", size = 12),
  )



compute_PCoA <- function(meta_df, mat, method = "bray",  prevalence_threshold = 0.1, grouping, transformed) {
  # Takes a metadata dataframe, a genus relative abundance matrix and computes a PCoA
  # Filters rel. abundance matrix to keep only feature with at least prevalence_threshold prevalence
  # A feature is considered "prevalent" in a sample if its relative abundance is above threshold_for_prevalence
  
  stopifnot(all(c("Sample_ID", grouping) %in% colnames(meta_df)))
  mat_clean <- mat[, meta_df$Sample_ID]
  
  if (transformed == 'true') {
    transformed_relAB_sel <- mat_clean
    message('Matrix is already transformed: number of features in the matrix are ',nrow(mat))
  } else {
    
    # Filter matrix to keep only features with at least 10% prevalence
    relAB_sel <- mat_clean[rowSums(mat_clean > 0) / ncol(mat_clean) > prevalence_threshold, ]
    
    message(paste0(nrow(relAB_sel), " out of ", nrow(mat_clean), " features passed prevalence cutoff"))
    
    # # remove empty columns
    # if (method == "bray") {
    #   relAB_sel <- relAB_sel[, colSums(relAB_sel) > 0]
    #   # Compute Bray-Curtis distances
    #   transformed_relAB_sel <- microbiome::transform(relAB_sel+1e-05, "log10")
    #}
  }
  
  dist <- vegan::vegdist(x = t(relAB_sel), method = method)
  
  
  pco <- labdsv::pco(dis = dist, k = 2)
  res <- as_tibble(pco$points, rownames = "Sample_ID")
  colnames(res)[2:ncol(res)] <- paste0("PCoA", 2:ncol(res) - 1)
  
  res <- res %>% left_join(., meta_df %>% dplyr::select(grouping, Sample_ID))
  #* Calclualte proportion of variance explained
  var_expl <- data.frame(matrix(ncol = 3, nrow = 1))
  
  for (i in seq(1, 3)) {
    var_expl[1, i] <- round((pco$eig[i] / sum(pco$eig)) * 100, 1)
  }
  colnames(var_expl) <- paste0("PCo", seq(1:ncol(var_expl)), "[%]")
  
  # PERMANOVA test
  set.seed(1)
  formula_str <- as.formula(paste("dist ~", grouping))
  adonis_result <- adonis2(formula_str, data = meta_df, permutations = 999, parallel = 32, sqrt.dist = F)
  perm_p_value <- adonis_result$`Pr(>F)`[1]
  perm_r_squared <- adonis_result$R2[1]
  
  return(list(PCO_res = res, var_expl_df = var_expl,perm_p_value = perm_p_value, perm_r_squared = perm_r_squared ))
}



plot_PCoA <- function(meta_df, mat, method = "bray", prevalence_threshold = 0.1, point_alpha = 0.6, grouping, transformed) {
  # Takes a metadata dataframe, a species relative abundance matrix and computes a PCoA and returns a ggplot object.
  # Calls the f_compute PCoA function for computation
  
  stopifnot(all(c("Sample_ID", grouping) %in% colnames(meta_df)))
  
  pcoa_list <- compute_PCoA(meta_df, mat, method, prevalence_threshold, grouping, transformed)
  #load('/g/scb/zeller/pekel/meta_analysis/src/analysis/Rdata/pcoa_list.Rdata')
  
  print(pcoa_list)
  
  x_lab <- paste0("PCoA1 [", pcoa_list$var_expl_df[1, 1], "%]")
  y_lab <- paste0("PCoA2 [", pcoa_list$var_expl_df[1, 2], "%]")
  z_lab <- paste0("PCoA3 [", pcoa_list$var_expl_df[1, 3], "%]")
  asp_ratio <- pcoa_list$var_expl_df$`PCo2[%]` / pcoa_list$var_expl_df$`PCo1[%]`
  
  # Creating the text for the annotation
  annotation_text <- paste0('PERMANOVA p-value: ', round(pcoa_list$perm_p_value, 3),
                            ', R²: ', round(pcoa_list$perm_r_squared, 3))
  #save(pcoa_list, file = '/g/scb/zeller/pekel/meta_analysis/src/analysis/Rdata/pcoa_list.Rdata')
  
  pt <- pcoa_list$PCO_res %>%
    ggplot(aes(x = PCoA1, y = PCoA2, fill = !!sym(grouping))) +
    geom_point(size = 1.8, alpha = point_alpha, color = "black", pch = 21) +
    xlab(x_lab) +
    ylab(y_lab) +
    ggtitle(paste0("PCoA (",method,')'), subtitle = annotation_text) +  
    coord_fixed(ratio = asp_ratio) +
    theme_paper
  
  
  return(pt)
}


f_convert_qval_pval = function(qvalues,fdr_threshold) {
  # This function is used to indicate in the volcano plots the p-value that would correspond to a given FDR threshold (e.g. 0.2)    
  
  #Function to approximate p-value from a distribution of q-values
  #copy-pasted from here: https://stats.stackexchange.com/questions/51070/how-can-i-convert-a-q-value-distribution-to-a-p-value-distribution
  
  # you need to know the estimate of pi0 used to create the q-value
  # that's the maximum q-value (or very, very close to it)
  qvalues_extended <- c(qvalues,fdr_threshold)
  pi0 = max(qvalues_extended)
  # compute m0, the estimated number of true nulls
  m0 = length(qvalues_extended) * pi0
  # then you multiply each q-value by the proportion of true nulls
  # expected to be under it (the inverse of how you get there from
  # the p-value):
  p_values_reconstructed <- qvalues_extended * rank(qvalues_extended) / m0
  #return the last element of the reconstructed p-values (corresponds to the fdr_threshold - p-value)
  return(p_values_reconstructed[length(p_values_reconstructed)])
}



# Generate a function run ml models
train_model_rf <- function(meta_df, mat, label_column, case_label, control_label, block_label = NULL, num_trees = 200, seed = 2025, prev.cutoff=0.1) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Create label
  label <- create.label(meta = meta_df, label = label_column, case = case_label, control = control_label)
  
  # Initialize SIAMCAT object
  siamcat <- siamcat(feat = mat, meta = meta_df, label = label)
  
  # Filter features by prevalence and abundance
  siamcat <- filter.features(siamcat, filter.method = 'prevalence', cutoff = prev.cutoff)
  
  # Normalize features
  siamcat <- normalize.features(siamcat, norm.method = "log.std", 
                                norm.param = list(log.n0 = 1e-05, sd.min.q = .1), 
                                feature.type = "filtered")
  
  # Create data splits for cross-validation
  if (is.null(block_label)) {
    siamcat <- create.data.split(siamcat, num.folds = 10, num.resample = 10)
  } else {
    siamcat <- create.data.split(siamcat, num.folds = 10, num.resample = 10, inseparable = block_label)
  }
  
  # Train Random Forest model
  siamcat <- train.model(siamcat, method = "randomForest", param.set = list('num.trees' = num_trees))
  
  # Make predictions
  siamcat <- make.predictions(siamcat)
  
  # Evaluate predictions
  evaluation <- evaluate.predictions(siamcat)
  
  # Return the evaluation object
  return(evaluation)
}

normalize_data <- function(data_df, offset = 1e-5) {
  # Apply log transformation and standardize each feature
  normalized_df <- data_df %>%
    mutate(across(everything(), ~ scale(log10(. + offset), center = TRUE, scale = TRUE)[, 1]))
  return(normalized_df)
}


run_lmem <- function(data_df, meta_df, ref_group, column_name) {
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lme4)
  
  # Normalize data
  data_df <- normalize_data(data_df, offset = 1e-5)
  
  # Initialize result containers
  bacteria <- p.val <- NULL
  effect.size <- list()
  conf.int <- list()
  pr.shift <- list()
  sample.counts <- list()  # Container for sample counts
  
  # Check if `Cohort` column exists
  if (!"Cohort" %in% colnames(meta_df)) {
    stop("The `meta_df` dataset must contain a `Cohort` column. Please ensure your metadata is correctly formatted.")
  }
  
  # Check if the group column has more than two levels
  unique_groups <- unique(meta_df[[column_name]])
  if (length(unique_groups) > 2) {
    stop(paste(
      "The column", column_name, "has more than two levels:",
      paste(unique_groups, collapse = ", "),
      "\nPlease reduce the levels to two before running the function."
    ))
  }
  
  # Reshape and join metadata
  data <- data_df %>%
    rownames_to_column('Bacteria') %>%
    pivot_longer(-Bacteria, names_to = "Sample_ID", values_to = "value") %>%
    left_join(meta_df %>% select(column_name, Cohort, Sample_ID), by = "Sample_ID")
  
  Features <- rownames(data_df)
  
  # Iterate over each feature
  for (i in seq_along(Features)) {
    feature_name <- Features[i]
    tmp <- data %>% filter(Bacteria == feature_name)
    
    # Skip if subset is empty or group has fewer than 2 levels
    if (nrow(tmp) == 0 || length(unique(tmp[[column_name]])) < 2) {
      next
    }
    
    # Set reference group
    tmp[[column_name]] <- as.factor(tmp[[column_name]])
    tmp[[column_name]] <- relevel(tmp[[column_name]], ref = ref_group)
    
    # Extract non-reference group
    non_ref_group <- setdiff(levels(tmp[[column_name]]), ref_group)
    
    # Count samples in each group
    n.group1 <- sum(tmp[[column_name]] == ref_group)
    n.group2 <- sum(tmp[[column_name]] == non_ref_group)
    sample.counts[[i]] <- tibble(n.group1 = n.group1, n.group2 = n.group2)
    
    # Per-cohort estimates
    datasetWiseEstimates <- list()
    uq <- unique(tmp$Cohort)
    for (dsIndex in seq_along(uq)) {
      dataset <- tmp %>% filter(Cohort == uq[dsIndex])
      if (length(unique(dataset[[column_name]])) < 2) {
        datasetWiseEstimates[[length(datasetWiseEstimates) + 1]] <- NA
      } else {
        datasetWiseEstimates[[length(datasetWiseEstimates) + 1]] <- coefficients(summary(lm(value ~ get(column_name), data = dataset)))[2, 1]
      }
    }
    names(datasetWiseEstimates) <- uq
    
    # General linear mixed-effects model
    lmem <- lmer(value ~ get(column_name) + (1 | Cohort), data = tmp)
    bacteria[i] <- unique(tmp$Bacteria)
    p.val[i] <- coefficients(summary(lmem))[2, "Pr(>|t|)"]
    effect.size[[i]] <- c(coefficients(summary(lmem))[2, "Estimate"], unlist(datasetWiseEstimates))
    conf.int[[i]] <- tryCatch({
      confint(lmem, level = 0.95)[4, ]
    }, error = function(e) {
      c(NA, NA)  # Fallback if confidence intervals cannot be computed
    })
    
    # Prevalence shift
    x.pos <- tmp %>% filter(get(column_name) == non_ref_group)
    x.neg <- tmp %>% filter(get(column_name) == ref_group)
    pr.n <- mean(x.neg$value >= 1e-04)
    pr.p <- mean(x.pos$value >= 1e-04)
    pr.shift[[i]] <- c(pr.p - pr.n, pr.n, pr.p)
  }
  
  # Combine results into a table
  result_table <- tibble(
    Bacteria = bacteria,
    P.val = p.val,
    Effect.size = effect.size,
    conf.int = conf.int,
    pr.shift = pr.shift,
    sample.counts = sample.counts
  )
  
  # Expand nested results
  expanded_table <- result_table %>%
    unnest_wider(c(Effect.size, conf.int, pr.shift, sample.counts), names_sep = "_") %>%
    dplyr::rename(
      Effect.size = Effect.size_1,
      pr.shift = pr.shift_1,
      !!paste0("pr.", ref_group) := pr.shift_2,
      !!paste0("pr.", non_ref_group) := pr.shift_3,
      !!paste0("n.", ref_group) := sample.counts_n.group1,
      !!paste0("n.", non_ref_group) := sample.counts_n.group2
    ) %>%
    mutate(P.adj = p.adjust(P.val, method = "BH"))
  
  # Return the final table
  return(expanded_table)
}


plot_volcano <- function(plot_df, xBreaks, fdr_thresh = 0.05, man_y_breaks = NULL,nudge_y= NULL,min_segment_length=NULL,
                         add_to_y_axis = 0.25, group_control, group_case, color_vector = NULL,max.overlaps=NULL) {
  
  if (is.null(color_vector)) {
    color_vector <- c(
      "#C44600",        # Default color for group_case
      "dodgerblue3",    # Default color for group_control
      "white"           # Default color for non-significant
    )
    names(color_vector) <- c(group_case, group_control, "n.s.")
  }
  else{
    names(color_vector) <- c(group_case, group_control, "n.s.")
  }
  
  
  if (is.null(max.overlaps)){
    max.overlaps=10
  }
  
  
  if (is.null(min_segment_length)){
    min_segment_length= 0.1
  }
  
  if (is.null(nudge_y)){
    nudge_y= 3
  }
  
  plot_df <- plot_df %>%
    mutate(
      enriched_in =
        case_when(
          P.adj < 0.05 & Effect.size > 0 ~ group_case,
          P.adj < 0.05 & Effect.size < 0 ~ group_control,
          TRUE ~ "n.s."
        ),
      fdr_sig = ifelse(P.adj < 0.05, TRUE, FALSE),
      enriched_in = factor(enriched_in, levels = c(group_case, group_control, "n.s.")),
      font = ifelse(fdr_sig, "bold.italic", "italic"), # for labeling
      group_prev = case_when(
        Effect.size > 0 ~ !!sym(paste0("pr.", group_case)),
        Effect.size < 0 ~ !!sym(paste0("pr.", group_control)),
        TRUE ~ 0
      ),
      lab = case_when( P.adj < fdr_thresh ~ Bacteria, TRUE ~ ""),
      prev_size = case_when(
        group_prev >= 0.75 ~ 5,
        group_prev >= 0.50 ~ 3,
        group_prev >= 0.25 ~ 2,
        TRUE ~ 1.5 # Default size for points
      )
    )
  
  
  # Extract sample sizes
  sample_sizes <- plot_df %>%
    select(starts_with("n.")) %>%
    distinct() %>%
    pivot_longer(cols = everything(), names_to = "group", values_to = "sample_count") %>%
    mutate(group = gsub("^n\\.", "", group)) # Remove "n." prefix to get group names
  
  # Compute the q-value threshold
  q_value_threshold <- f_convert_qval_pval(qvalues = plot_df$P.adj, fdr_threshold = fdr_thresh)
  
  # Create enriched_in labels with sample counts
  enriched_in_labels <- sample_sizes %>%
    mutate(label = paste0(group, " (N:", sample_count, ")")) %>%
    select(group, label) %>%
    deframe() # Convert to a named vector for use in scale_fill_manual()
  
  # Add "n.s." label separately for non-significant
  enriched_in_labels <- c(enriched_in_labels, `n.s.` = "n.s")
  
  # Create the volcano plot
  pt <-plot_df %>%
    ggplot(aes(x = Effect.size, y = -log10(P.adj))) +
    theme_paper +
    geom_hline(yintercept = -log10(0.05), color = "gray", lty = "dashed", lwd = 0.5) +
    #geom_hline(yintercept = -log10(q_value_threshold), color = "darkgrey", lty = "dashed", lwd = 0.75) +
    geom_point(
      aes(
        size = prev_size,
        fill = enriched_in
      ),
      shape = 21, color = "black", alpha = 0.8
    ) +
    scale_size(
      range = c(1.5, 6),
      breaks = c(1.5, 2, 3, 5), # Sizes corresponding to thresholds
      labels = c("0.10", "0.25", "0.50", "0.75"), # Threshold labels
      guide = guide_legend(reverse = TRUE)
    ) +
    ggrepel::geom_text_repel(
      aes(label = lab, fontface = font),
      color = "black",
      segment.color = "black",
      max.overlaps = max.overlaps,
      size = 3,
      seed = 420, 
      force = 1,
      min.segment.length = min_segment_length, 
      nudge_y = nudge_y) +
    theme(
      legend.box = "horizontal",
      legend.spacing.y = unit(0.1, "cm"),
      legend.position = c(0.01, .99),
      legend.key.size = unit(0.75, "lines"),
      legend.justification = c(0, 1)
    ) +
    labs(
      fill = paste0("P adj. < 0.05"), # Legend title
      size = "Prevalence\n"
    ) +
    guides(
      color = "none",
      size = guide_legend(order = 1, reverse = T, title = "Prevalence"),
      fill = guide_legend(order = 2, override.aes = list(size = 3.5))
    ) +
    scale_fill_manual(
      values = color_vector,
      labels = enriched_in_labels # Use the enriched_in labels with sample sizes
    ) +
    xlab("Enrichment effect size")
  
  
      
  
  return(pt)
}


# To plot ROC curve of ml models 

plot_roc_siamcat_models <- function(models, labels, colours, trained_on = NULL,alpha=NULL) {
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
  
  # if(is.null(linetype)){
  #   linetype='solid'
  # }
  
  if(is.null(alpha)){
    alpha=1
  }
  wrapped_labels <- str_wrap(unlist(auc_list), width = 55)
  
  
  p <- ggplot(all_roc_data, aes(x = FPR, y = TPR, color = model)) +
  geom_line(size = 1.8, alpha = alpha) +
  scale_color_manual(values = colours, labels = wrapped_labels) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey66", size = 0.8) +
  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", color = "Model") +
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
  coord_fixed()

  return(p)
}


# Balanced number of samples across selected groups fore ML models 
balance_data_by_groups <- function(meta_df, data_df, group_column = NULL, groups = NULL ,seed =2000) {
  
  set.seed(seed)
  # Ask for the column containing group labels if not provided
  if (is.null(group_column)) {
    group_column <- readline(prompt = "Enter the column name that contains the groups to balance (e.g., 'Condition'): ")
  }
  
  # Ask for the specific groups to balance if not provided
  if (is.null(groups)) {
    groups <- unlist(strsplit(readline(prompt = "Enter the two groups to balance, separated by a comma (e.g., 'CRC,CTR'): "), ","))
    groups <- trimws(groups)  # Remove extra spaces
  }
  
  # Ensure correct input
  if (length(groups) != 2) {
    stop("Please provide exactly two groups to balance.")
  }
  
  group1 <- groups[1]
  group2 <- groups[2]
  
  # Group by Cohort and selected group column, count occurrences
  counts <- meta_df %>% 
    filter(.data[[group_column]] %in% groups) %>% 
    group_by(Cohort, .data[[group_column]]) %>% 
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = group_column, values_from = n, id_cols = Cohort) %>%
    mutate(downsampled_to = pmin(.data[[group1]], .data[[group2]], na.rm = TRUE))  # Ensure min is properly handled
  
  # Generate a balanced dataset
  data_balanced_list <- list()
  
  for (i in seq_along(counts$Cohort)) {
    study <- counts$Cohort[i]
    message("Processing Cohort: ", study)
    
    data <- meta_df %>% filter(Cohort == study)
    sampling_size <- counts$downsampled_to[i]
    
    # Sample the larger group down to the smaller group's size
    if (counts[[group1]][i] > counts[[group2]][i]) {
      data_group1 <- data %>% filter(.data[[group_column]] == group1) %>% sample_n(size = sampling_size, replace = FALSE)
      data_group2 <- data %>% filter(.data[[group_column]] == group2)
    } else {
      data_group1 <- data %>% filter(.data[[group_column]] == group1)
      data_group2 <- data %>% filter(.data[[group_column]] == group2) %>% sample_n(size = sampling_size, replace = FALSE)
    }
    
    # Combine balanced data for this cohort
    data_balanced_list[[i]] <- bind_rows(data_group1, data_group2)
  }
  
  # Combine all cohorts into a single balanced dataset
  meta <- bind_rows(data_balanced_list) %>%
    as.data.frame() %>%
    column_to_rownames('Sample_ID')
  
  # Subset the main dataset based on the balanced metadata
  data <- data_df[, rownames(meta), drop = FALSE]
  
  return(list(data = data, meta = meta))
}























