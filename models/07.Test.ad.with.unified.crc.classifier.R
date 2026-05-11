##########################################
# Test unified CRC model on AD vs CTR samples (Chaniging test folds to include AD samples in test set)

library(SIAMCAT)
library(here)
library(progress)
library(dplyr)
library(phyloseq)
library(ranger)
library(pROC)
library(readr)
library(tibble)
library(stringr)

load(here('data','results','Training.unified.crc.model.Rdata'))

# Load training and test folds od CRC unified model
training.folds <- models.all.rf@data_split$training.folds
test_folds     <- models.all.rf@data_split$test.folds

lab_train <- models.all.rf@label$label  # Label of the samples (1,-1)

crc_ids <- names(lab_train)[lab_train ==  1]
ctr_ids <- names(lab_train)[lab_train == -1]

# SIAMCAT-trained normalized feature matrix 
otu_train_normalized <- models.all.rf@norm_feat$norm.feat

# metadata from trained phyloseq 
metadata_train <- models.all.rf@phyloseq@sam_data %>% data.frame()

cat("Training matrix dim (features x samples):\n")
print(dim(otu_train_normalized))

# Load AD samples
meta.ad <- read_tsv(here("data", "Metadata.all.samples.tsv")) %>%
  filter(Condition %in% c("AD", "AdvAD", "smallAD", "CTR")) %>%
  as.data.frame() %>%
  column_to_rownames("Sample_ID") %>%
  mutate(
    Condition_general = case_when(
      Condition %in% c("AdvAD", "smallAD") ~ "AD",
      TRUE ~ Condition
    )
  )

all.data <- read.table(here("data", "Relab.all.samples.tsv"),
                       sep = "\t", check.names = FALSE) %>%
  rownames_to_column("genus") %>%
  filter(genus != "unassigned") %>%
  column_to_rownames("genus")

# Restrict raw table to samples present in meta.ad
data.ad <- all.data[, rownames(meta.ad), drop = FALSE]

# AD sample IDs
ad_samples <- rownames(meta.ad)[meta.ad$Condition_general == "AD"]

# Build ad siamcat object
label <- create.label(meta = meta.ad, label = 'Condition_general', case = "AD", control = "CTR")
sc.ad <- siamcat(feat = data.ad, meta = meta.ad, label = label)

# extract features 
train_feats <- rownames(models.all.rf@norm_feat$norm.feat)

project_to_training_space <- function(raw_feat_mat, train_feat_names, fill_value = 0) {
  # raw_feat_mat: features x samples
  # returns: features x samples with exact train_feat_names (same order)
  raw_feat_mat <- as.matrix(raw_feat_mat)
  
  out <- matrix(
    fill_value,
    nrow = length(train_feat_names),
    ncol = ncol(raw_feat_mat),
    dimnames = list(train_feat_names, colnames(raw_feat_mat))
  )
  
  common <- intersect(train_feat_names, rownames(raw_feat_mat))
  if (length(common) > 0) {
    out[common, ] <- raw_feat_mat[common, , drop = FALSE]
  }
  out
}

Xte_proj_raw <- project_to_training_space(
  raw_feat_mat     = data.ad,     
  train_feat_names = train_feats,
  fill_value       = 0
)

# Normalization parameters used in unified CRC classifier
norm_param_train <- models.all.rf@norm_feat$norm.param

# get normalization method
norm_method_train <- norm_param_train$norm.method

label_tmp <- create.label(meta = meta.ad, label = "Condition_general", case = "AD", control = "CTR")

sc.tmp <- siamcat(feat = Xte_proj_raw, meta = meta.ad, label = label_tmp)

sc.tmp <- normalize.features(
  sc.tmp,
  norm.method  = norm_method_train,
  norm.param   = norm_param_train,
  feature.type = "original"
)

Test_norm_frozen <- sc.tmp@norm_feat$norm.feat  
stopifnot(identical(rownames(Test_norm_frozen), rownames(otu_train_normalized)))

# Check parameters to see if they are same
# Training normalization parameters
# Test normalization parameters (from temp object)
norm_test  <- sc.tmp@norm_feat$norm.param

identical(norm_param_train, norm_test)

# Check number of AD samples 
cat("\n# AD samples in meta.ad:\n")
print(length(ad_samples))
cat("# AD samples present in data.ad columns:\n")
print(sum(ad_samples %in% colnames(data.ad)))

#  Build new test folds: (AD split per fold) + (previous CTR test samples)
make_nested_test_folds_ad_plus_prevctr <- function(ad_ids,
                                                   orig_test_folds,
                                                   ctr_ids,
                                                   n_runs = 10,
                                                   n_folds = 10,
                                                   seed = 123) {
  stopifnot(length(orig_test_folds) == n_runs)
  stopifnot(all(vapply(orig_test_folds, length, integer(1)) == n_folds))
  
  set.seed(seed)
  
  out <- vector("list", n_runs)
  
  for (r in seq_len(n_runs)) {
    shuffled_ad <- sample(ad_ids)
    ad_groups <- split(shuffled_ad, rep(seq_len(n_folds), length.out = length(shuffled_ad)))
    
    out[[r]] <- vector("list", n_folds)
    for (k in seq_len(n_folds)) {
      prev_test <- orig_test_folds[[r]][[k]]
      prev_ctr  <- intersect(prev_test, ctr_ids)   # keep only CTR from previous test fold
      out[[r]][[k]] <- unique(c(ad_groups[[k]], prev_ctr))
    }
  }
  out
}

new_test_folds <- make_nested_test_folds_ad_plus_prevctr(
  ad_ids = ad_samples,
  orig_test_folds = test_folds,
  ctr_ids = ctr_ids,
  n_runs = 10,
  n_folds = 10,
  seed = 123
)

# Sanity check
stopifnot(length(new_test_folds) == 10)
stopifnot(all(vapply(new_test_folds, length, integer(1)) == 10))

# sanity check to see if CTR samples used in test folds are same in models.all.rf test folds

orig_test_folds <- models.all.rf@data_split$test.folds

extract_ctr <- function(folds, ctr_ids) {
  lapply(folds, function(run) {
    lapply(run, function(fold) {
      sort(intersect(fold, ctr_ids))
    })
  })
}

ctr_orig <- extract_ctr(orig_test_folds, ctr_ids)
ctr_new  <- extract_ctr(new_test_folds, ctr_ids)

ctr_same_matrix <- matrix(FALSE, nrow = length(ctr_orig), ncol = length(ctr_orig[[1]]))

for (r in seq_along(ctr_orig)) {
  for (k in seq_along(ctr_orig[[r]])) {
    ctr_same_matrix[r, k] <- identical(
      ctr_orig[[r]][[k]],
      ctr_new[[r]][[k]]
    )
  }
}

# Sanity check
ctr_same_matrix

# Train per run/fold on fixed training.folds and test on new_test_folds
# Trains on otu_train_normalized object 
# Tests on data.ad 
-
train_predict_nested <- function(X_train_norm, y_train_labels,
                                 training_folds, test_folds,
                                 Test_norm, meta_test,
                                 num_trees = 200, seed = 2025,
                                 verbose_counts = TRUE) {
  
  # X_train_norm: features x samples (normalized, SIAMCAT space)
  X_train_norm <- as.matrix(X_train_norm)
  
  # y_train_labels: named vector over training samples (1/-1)
  y_train_labels <- y_train_labels[colnames(X_train_norm)]
  stopifnot(!any(is.na(y_train_labels)))
  
  train_feats <- rownames(X_train_norm)
  
  n_runs  <- length(training_folds)
  n_folds <- length(training_folds[[1]])
  n_total <- n_runs * n_folds
  
  pb <- progress_bar$new(
    total  = n_total,
    format = "  [:bar] :current/:total (:percent) | run :run fold :fold",
    clear  = FALSE,
    width  = 60
  )
  
  models <- vector("list", n_runs)
  preds  <- vector("list", n_runs)
  
  for (r in seq_len(n_runs)) {
    models[[r]] <- vector("list", n_folds)
    preds[[r]]  <- vector("list", n_folds)
    
    for (k in seq_len(n_folds)) {
      set.seed(seed + 100*r + k)
      
      tr_ids <- intersect(training_folds[[r]][[k]], colnames(X_train_norm))
      te_ids <- test_folds[[r]][[k]]
      
      # ---- Train ----
      Xtr <- t(X_train_norm[, tr_ids, drop = FALSE])  # samples x features
      ytr <- factor(ifelse(y_train_labels[tr_ids] == 1, "CRC", "CTR"))
      
      if (verbose_counts) {
        cat("\n", str_c("training folds ", r, "-", k), "\n", sep = "")
        # join with training metadata if available
        tmp <- tibble(Sample_ID = rownames(Xtr)) %>%
          left_join(metadata_train %>%
                      rownames_to_column("Sample_ID") %>%
                      select(Sample_ID, Condition),
                    by = "Sample_ID") %>%
          count(Condition)
        print(tmp)
      }
      
      fit <- ranger(
        x = as.data.frame(Xtr),
        y = ytr,
        num.trees = num_trees,
        probability = TRUE
      )
      models[[r]][[k]] <- fit
      
      # ---- Test ----
      te_ids <- intersect(te_ids, colnames(Test_norm))  # safety
      if (length(te_ids) > 0) {
        
        if (verbose_counts) {
          cat(str_c("test folds ", r, "-", k), "\n")
          tmp2 <- tibble(Sample_ID = te_ids) %>%
            left_join(meta_test %>%
                        rownames_to_column("Sample_ID") %>%
                        select(Sample_ID, Condition_general),
                      by = "Sample_ID") %>%
            count(Condition_general)
          print(tmp2)
        }
        
        # raw_test_table is features x samples. First subset to te_ids:
        Xte_raw <- Test_norm[, te_ids, drop = FALSE]
        
        # Project to SAME features as training:
        Xte_proj <- project_to_training_space(
          raw_feat_mat = Xte_raw,
          train_feat_names = train_feats,
          fill_value = 0
        )
      
        
        Xte <- as.data.frame(t(Xte_proj))  # samples x features
        
        pr <- predict(fit, data = Xte)$predictions
        preds[[r]][[k]] <- setNames(pr[, "CRC"], rownames(Xte))
        
      } else {
        preds[[r]][[k]] <- numeric(0)
      }
      
      pb$tick(tokens = list(run = r, fold = k))
    }
  }
  
  list(models = models, preds = preds)
}

res <- train_predict_nested(
  X_train_norm   = otu_train_normalized,
  y_train_labels = lab_train,
  training_folds = training.folds,
  test_folds     = new_test_folds,
  Test_norm      = as.data.frame(Test_norm_frozen),
  meta_test      = meta.ad,
  num_trees      = 200,
  seed           = 2025,
  verbose_counts = TRUE
)

#  label AD as positive (1) and CTR as negative (-1)
y_test <- rep(NA_integer_, length(unique(c(ad_samples, ctr_ids))))
names(y_test) <- unique(c(ad_samples, ctr_ids))
y_test[ad_samples] <-  1L
y_test[ctr_ids]    <- -1L

auc_by_run <- sapply(seq_along(res$preds), function(r) {
  pr <- unlist(res$preds[[r]], use.names = TRUE)
  pr <- pr[!is.na(pr)]
  yy <- y_test[names(pr)]
  
  keep <- !is.na(yy)
  pr <- as.numeric(pr[keep])
  yy <- yy[keep]
  
  roc_obj <- pROC::roc(
    response  = factor(yy, levels = c(-1, 1)),
    predictor = pr,
    quiet = TRUE
  )
  as.numeric(pROC::auc(roc_obj))
})

cat("\nAUC per run:\n")
print(auc_by_run)
cat("\nMean AUC ± SD:\n")
cat(mean(auc_by_run), "±", sd(auc_by_run), "\n")


pr_all <- unlist(unlist(res$preds, recursive = FALSE), use.names = TRUE)
pr_all <- pr_all[!is.na(pr_all)]
yy_all <- y_test[names(pr_all)]
keep_all <- !is.na(yy_all)

roc_all <- pROC::roc(
  response  = factor(yy_all[keep_all], levels = c(-1, 1)),
  predictor = as.numeric(pr_all[keep_all]),
  quiet = TRUE
)

cat("\nPooled AUC:\n")
print(as.numeric(pROC::auc(roc_all)))

# Sanity checks
cat("\n#AD samples in meta.ad:\n")
print(length(ad_samples))

cat("#AD samples included in new_test_folds (unique):\n")
print(length(unique(unlist(new_test_folds))))


# Test also CRC samples as a sanity check
train_predict_nested_internal <- function(X_train_norm, y_train_labels,
                                          training_folds, test_folds,
                                          num_trees = 200, seed = 2025,
                                          verbose_counts = FALSE) {
  
  X_train_norm <- as.matrix(X_train_norm)
  y_train_labels <- y_train_labels[colnames(X_train_norm)]
  stopifnot(!any(is.na(y_train_labels)))
  
  n_runs  <- length(training_folds)
  n_folds <- length(training_folds[[1]])
  n_total <- n_runs * n_folds
  
  pb <- progress_bar$new(
    total  = n_total,
    format = "  [:bar] :current/:total (:percent) | run :run fold :fold",
    clear  = FALSE,
    width  = 60
  )
  
  models <- vector("list", n_runs)
  preds  <- vector("list", n_runs)
  
  for (r in seq_len(n_runs)) {
    models[[r]] <- vector("list", n_folds)
    preds[[r]]  <- vector("list", n_folds)
    
    for (k in seq_len(n_folds)) {
      set.seed(seed + 100*r + k)
      
      tr_ids <- intersect(training_folds[[r]][[k]], colnames(X_train_norm))
      te_ids <- test_folds[[r]][[k]]
      
      # Train
      Xtr <- t(X_train_norm[, tr_ids, drop = FALSE])      # samples x features
      ytr <- factor(ifelse(y_train_labels[tr_ids] == 1, "CRC", "CTR"))
      
      if (verbose_counts) {
        cat("\n", str_c("internal training folds ", r, "-", k), "\n", sep = "")
        print(table(ytr))
      }
      
      fit <- ranger(
        x = as.data.frame(Xtr),
        y = ytr,
        num.trees = num_trees,
        probability = TRUE
      )
      models[[r]][[k]] <- fit
      
      # Test (same feature space, same normalization)
      if (length(te_ids) > 0) {
        Xte <- as.data.frame(t(X_train_norm[, te_ids, drop = FALSE]))  # samples x features
        pr  <- predict(fit, data = Xte)$predictions
        preds[[r]][[k]] <- setNames(pr[, "CRC"], rownames(Xte))
      } else {
        preds[[r]][[k]] <- numeric(0)
      }
      
      pb$tick(tokens = list(run = r, fold = k))
    }
  }
  
  list(models = models, preds = preds)
}

test_folds <- models.all.rf@data_split$test.folds

#  Run internal sanity CV using the ORIGINAL test_folds from models.all.rf
res_internal <- train_predict_nested_internal(
  X_train_norm   = otu_train_normalized,
  y_train_labels = lab_train,
  training_folds = training.folds,
  test_folds     = test_folds,     
  num_trees      = 200,
  seed           = 2025,
  verbose_counts = FALSE
)


# Compute AUCs for CRC vs CTR on the original test folds
# lab_train already encodes CRC=1, CTR=-1 for ALL samples in the model 
pr_all_crc <- unlist(unlist(res_internal$preds, recursive = FALSE), use.names = TRUE)
pr_all_crc <- pr_all_crc[!is.na(pr_all_crc)]


yy_all_crc <- models.all.rf@label$label[names(pr_all_crc)]
keep_all_crc <- !is.na(yy_all_crc)

roc_all_crc <- pROC::roc(
  response  = factor(yy_all_crc[keep_all_crc], levels = c(-1, 1)),
  predictor = as.numeric(pr_all_crc[keep_all_crc]),
  quiet = TRUE
)

cat("\nPooled AUC:\n")
print(as.numeric(pROC::auc(roc_all_crc)))

y_internal <- lab_train  # named vector: sample -> 1/-1

auc_by_run_internal <- sapply(seq_along(res_internal$preds), function(r) {
  
  pr <- unlist(res_internal$preds[[r]], use.names = TRUE)
  
  # align by sample IDs (names)
  ids <- intersect(names(pr), names(y_internal))
  pr2 <- as.numeric(pr[ids])
  yy2 <- as.integer(y_internal[ids])
  
  # sanity: need both classes
  if (length(yy2) == 0L) return(NA_real_)
  if (length(unique(yy2)) < 2L) return(NA_real_)
  
  roc_obj <- pROC::roc(
    response  = factor(yy2, levels = c(-1, 1)),
    predictor = pr2,
    quiet = TRUE
  )
  as.numeric(pROC::auc(roc_obj))
})

auc_by_run_internal


auc_by_run_internal <- sapply(seq_along(res_internal$preds), function(r) {
  pr <- unlist(res_internal$preds[[r]], use.names = TRUE)
  pr <- pr[!is.na(pr)]
  yy <- y_internal[names(pr)]
  
  keep <- !is.na(yy)
  pr <- as.numeric(pr[keep])
  yy <- yy[keep]
  
  roc_obj <- pROC::roc(
    response  = factor(yy, levels = c(-1, 1)),
    predictor = pr,
    quiet = TRUE
  )
  as.numeric(pROC::auc(roc_obj))
})

cat("\n[Sanity check] Internal CRC-vs-CTR AUC per run (using original test_folds):\n")
print(auc_by_run_internal)
cat("\n[Sanity check] Internal Mean AUC ± SD:\n")
cat(mean(auc_by_run_internal), "±", sd(auc_by_run_internal), "\n")

# Optional pooled internal AUC
pr_int_all <- unlist(unlist(res_internal$preds, recursive = FALSE), use.names = TRUE)
pr_int_all <- pr_int_all[!is.na(pr_int_all)]
yy_int_all <- y_internal[names(pr_int_all)]
keep_int   <- !is.na(yy_int_all)

roc_int_all <- pROC::roc(
  response  = factor(yy_int_all[keep_int], levels = c(-1, 1)),
  predictor = as.numeric(pr_int_all[keep_int]),
  quiet = TRUE
)

cat("\n[Sanity check] Internal pooled AUC:\n")
print(as.numeric(pROC::auc(roc_int_all)))


# Plot now 

load(here('data','results','Training.ad.ctr.rf.model.Rdata'))

models <- list(models.ad)
labels <- c("AD Classifier CV:")
trained_on <- list(NULL)
colours <- c( "#FFBA08")

ad_auc_plot<- plot_roc_siamcat_models(models, labels, colours, trained_on,alpha=0.8)

ggsave(ad_auc_plot, file=here('figures','extended.data.figure4','Extended.Data.Figure4f.pdf'),width = 6, height = 6)


# Plot pooled ROC curve (AD=case, CTR=control)
auc_val <- as.numeric(pROC::auc(roc_all))



plot_roc_siamcat_models <- function(models, labels, colours,
                                    trained_on = NULL, alpha = NULL,
                                    extra_rocs = NULL, extra_labels = NULL,
                                    extra_colours = NULL, extra_alpha = NULL) {
  
  determine_tpr_fpr_auc <- function(eval_data, auroc) {
    tpr_list <- list(); fpr_list <- list()
    for (i in seq_along(eval_data)) {
      tp <- eval_data[[i]]$tp; fp <- eval_data[[i]]$fp
      tn <- eval_data[[i]]$tn; fn <- eval_data[[i]]$fn
      if (is.null(tp) || is.null(fp) || is.null(tn) || is.null(fn)) next
      tpr_list[[i]] <- tp / (tp + fn)
      fpr_list[[i]] <- fp / (fp + tn)
    }
    if (length(tpr_list) == 0 || length(fpr_list) == 0) stop("No valid evaluation data found.")
    
    max_len <- max(sapply(tpr_list, length), sapply(fpr_list, length))
    tpr_matrix <- do.call(rbind, lapply(tpr_list, function(x) c(x, rep(NA, max_len - length(x)))))
    fpr_matrix <- do.call(rbind, lapply(fpr_list, function(x) c(x, rep(NA, max_len - length(x)))))
    
    roc_data <- data.frame(
      FPR = apply(fpr_matrix, 2, mean, na.rm = TRUE),
      TPR = apply(tpr_matrix, 2, mean, na.rm = TRUE)
    )
    list(roc_data = roc_data, auc = auroc)
  }
  
  if (is.null(alpha)) alpha <- 1
  
  all_roc_data <- data.frame(FPR = numeric(), TPR = numeric(), model = character())
  auc_list <- list()
  
  # ---- SIAMCAT models ----
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
      mean_test_size  <- round(dim(models[[i]]@pred_matrix)[1])
    } else {
      training_sizes <- unlist(sapply(models[[i]]@data_split$training.folds, function(fold) sapply(fold, length)))
      test_sizes     <- unlist(sapply(models[[i]]@data_split$test.folds, function(fold) sapply(fold, length)))
      mean_train_size <- ifelse(length(training_sizes) > 0, round(mean(training_sizes, na.rm = TRUE)), NA)
      mean_test_size  <- ifelse(length(test_sizes) > 0, round(mean(test_sizes, na.rm = TRUE)), NA)
    }
    
    result <- determine_tpr_fpr_auc(eval_data, auroc)
    roc_data <- result$roc_data
    auc      <- result$auc
    
    roc_data$model <- labels[i]
    all_roc_data <- rbind(all_roc_data, roc_data)
    
    auc_list[[length(auc_list) + 1]] <-
      paste0(labels[i], " (AUC = ", round(auc, 2),
             ", Train N = ", mean_train_size,
             ", Test N = ", mean_test_size, ")")
  }
  
  # ---- Extra ROC curves (e.g., pooled roc_all) ----
  if (!is.null(extra_rocs)) {
    if (is.null(extra_labels)) extra_labels <- paste0("Extra ", seq_along(extra_rocs))
    if (is.null(extra_colours)) extra_colours <- rep("black", length(extra_rocs))
    if (is.null(extra_alpha)) extra_alpha <- rep(1, length(extra_rocs))
    
    for (j in seq_along(extra_rocs)) {
      ro <- extra_rocs[[j]]
      
      # pROC roc object: specificities + sensitivities
      roc_df <- data.frame(
        FPR = 1 - ro$specificities,
        TPR = ro$sensitivities
      )
      
      roc_df$model <- extra_labels[j]
      all_roc_data <- rbind(all_roc_data, roc_df)
      
      auc_list[[length(auc_list) + 1]] <-
        paste0(extra_labels[j], " (AUC = ", round(as.numeric(pROC::auc(ro)), 2), ")")
    }
    
    # extend colours vector with extra
    colours <- c(colours, extra_colours)
    alpha_map <- c(rep(alpha, length(labels)), extra_alpha)
  } else {
    alpha_map <- rep(alpha, length(labels))
  }
  
  if (nrow(all_roc_data) == 0) stop("No valid data for ROC plot. Please check your models.")
  
  wrapped_labels <- str_wrap(unlist(auc_list), width = 55)
  
  # Ensure factor order so legend matches label order
  all_roc_data$model <- factor(all_roc_data$model, levels = c(labels, extra_labels))
  
  p <- ggplot(all_roc_data, aes(x = FPR, y = TPR, color = model)) +
    geom_line(aes(alpha = model), size = 1.8) +
    scale_color_manual(values = colours, labels = wrapped_labels) +
    scale_alpha_manual(values = setNames(alpha_map, c(labels, extra_labels)), guide = "none") +
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
  
  p
}


models  <- list(models.ad)
labels  <- c("AD Classifier CV")
colours <- c("#FFBA08")

ad_auc_plot <- plot_roc_siamcat_models(
  models   = models,
  labels   = labels,
  colours  = colours,
  trained_on = list(NULL),
  alpha = 0.8,
  extra_rocs    = list(roc_all),
  extra_labels  = c("Unified CRC model evaluated on AD vs CTR"),
  extra_colours = c("darkgrey"),
  extra_alpha   = c(0.9)
)

ggsave(ad_auc_plot, file=here('figures','extended.data.figure4','Extended.Data.Figure4g.pdf'),width = 6, height = 6)














