######################
# Linear Mixed Effect Model for CRC tissue and adj. normal tissue comparisons
######################
# This script performs linear mixed effects modeling (LMM) to compare microbial taxonomic differences
# between CRC tissue and adjacent normal tissue samples.

# Load libraries and setup
source(here('src','utils.R'))

meta_df <- read_tsv(here('data', 'tissue_data','Amplicon_meta-analysis_meta_clean_2501007_filtered.tsv')) %>%
  mutate(sample_type = factor(sample_type, levels = c("Primary tumor", "Adj. non-tumor")))

count_df <- read_tsv(here('data', 'tissue_data','Amplicon_tissue_meta-analysis_data_251007.tsv')) %>% rename(bac=genus)
count_mat <- count_df %>%
  column_to_rownames("bac") %>%
  as.matrix()
stopifnot(!("Bacteria" %in% rownames(count_mat))) # make sure that "Bacteria" is not part of the genus names

count_rel_mat <- prop.table(count_mat, margin = 2)

# Run LMM for CRC tissue vs adjacent normal tissue comparison
count_rel_mat2<-count_rel_mat[which(rowSums(count_rel_mat > 0) / ncol(count_rel_mat) > 0.1),] # 10 % prev threshold

meta_df2<-meta_df %>% 
  filter(sample_type=='Adj. non-tumor' | sample_type=='Primary tumor')  %>% rename(Cohort='study_name') 

count_rel_mat2 <- count_rel_mat2[,meta_df2$Sample_ID]
count_rel_mat2 <- count_rel_mat2 %>% 
  as.data.frame() %>%
  rownames_to_column('Taxon') %>% 
  filter(Taxon!='unassigned') %>% 
  column_to_rownames('Taxon')

# Define the function to run linear mixed effects model accounting for cohort and patient_id as random effects
run_lmem_tissue <- function(data_df, meta_df,ref_group,column_name, feature_column_name = "Bacteria" ) {
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lme4)
  
  # Normalize data
  data_df <- normalize_data(data_df, offset = 1e-5)
  
  # Initialize result containers
  features <- p.val <- NULL
  effect.size <- list()
  conf.int <- list()
  pr.shift <- list()
  sample.counts <- list()
  
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
    tibble::rownames_to_column(feature_column_name) %>%
    pivot_longer(-all_of(feature_column_name), names_to = "Sample_ID", values_to = "value") %>%
    left_join(meta_df %>% select(all_of(column_name), Cohort, Sample_ID, patient_id), by = "Sample_ID")
  
  Features <- rownames(data_df)
  
  # Iterate over each feature
  for (i in seq_along(Features)) {
    feature_name <- Features[i]
    tmp <- data %>% filter(.data[[feature_column_name]] == feature_name)
    
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
        # Build formula safely
        form <- as.formula(paste("value ~", column_name))
        datasetWiseEstimates[[length(datasetWiseEstimates) + 1]] <- coefficients(summary(lm(form, data = dataset)))[2, 1]
      }
    }
    names(datasetWiseEstimates) <- uq
    
    # General linear mixed-effects model
    lmem_form <- as.formula(paste("value ~", column_name, "+ (1 | Cohort)", '+ (1 | patient_id)'))
    lmem <- lmer(lmem_form, data = tmp)
    features[i] <- unique(tmp[[feature_column_name]])
    p.val[i] <- coefficients(summary(lmem))[2, "Pr(>|t|)"]
    effect.size[[i]] <- c(coefficients(summary(lmem))[2, "Estimate"], unlist(datasetWiseEstimates))
    conf.int[[i]] <- tryCatch({
      confint(lmem, level = 0.95)[4, ]
    }, error = function(e) {
      c(NA, NA)  # Fallback if confidence intervals cannot be computed
    })
    
    # Prevalence shift
    x.pos <- tmp %>% filter(.data[[column_name]] == non_ref_group)
    x.neg <- tmp %>% filter(.data[[column_name]] == ref_group)
    pr.n <- mean(x.neg$value >= 1e-04)
    pr.p <- mean(x.pos$value >= 1e-04)
    pr.shift[[i]] <- c(pr.p - pr.n, pr.n, pr.p)
  }
  
  # Combine results into a table
  result_table <- tibble(
    !!feature_column_name := features,
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


lmm.tissue <- run_lmem_tissue(data_df = count_rel_mat2, meta_df = meta_df2, ref_group = 'Adj. non-tumor',
               column_name = 'sample_type', feature_column_name = 'Taxon')

write_tsv(lmm.tissue$table, file = here('data','results','lmm.tables.tissue.tsv'))

