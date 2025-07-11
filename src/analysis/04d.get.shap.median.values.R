######################
# SHAP value extraction and annotation
######################
# This script loads trained machine learning models, SHAP value objects, and microbial 
# abundance profiles for a given dataset (e.g., Alldata(Unified CRC classifier), EO-CRC, LO-CRC). 
# It extracts metadata (fold, resampling round, model type), aligns SHAP values with profile data, 
# and computes the median SHAP values per feature/sample. 
#
# For selected model types (e.g., RF), the script:
#   - Ranks features by mean absolute SHAP value
#   - Calculates the fraction of positive SHAP values
#   - Computes Spearman correlations between SHAP values and original feature abundances
#
# The final annotated table of SHAP values is written to a TSV file for downstream plotting.


get_shap_median <- function(dataset, label_case, model_types_to_evaluate) {
  
  # load models and clean up
  modelPaths <- list.files(here('data','results','shap.analysis' ,'models'), pattern = ".rds", full.names = TRUE)
  modelPaths <- modelPaths[str_detect(modelPaths, str_c(dataset, "__"))]
  models <- map(modelPaths, \(x) {
    readRDS(x)
  })
  names(models) <- modelPaths
  models <- enframe(models, name = 'raw_path', value = "model_object")
  models <- models %>%
    mutate(path = str_split_fixed(raw_path, paste0(here()), n = 2)[, 2]) %>%
    mutate(fold = str_split_fixed(path, "__", n = 3)[, 2]) %>%
    mutate(fold = str_replace(fold, "fold_id_", "")) %>%
    mutate(fold = as.numeric(fold)) %>%
    mutate(resampling = str_split_fixed(path, "__", n = 4)[, 3]) %>%
    mutate(resampling = str_replace(resampling, "repeat_", "")) %>%
    mutate(resampling = as.numeric(resampling)) %>%
    mutate(model_type = str_split_fixed(path, "__", n = 6)[, 5]) %>%
    mutate(model_type = str_replace(model_type, ".rds", "")) %>%
    select(-raw_path) %>%
    mutate(beta_values = map(model_object, \(x) {
      lam_min_index <- x$model$index['1se', ]
      betas <- x$model$glmnet.fit$beta[, lam_min_index]
      return(as.data.frame(betas) %>% rownames_to_column('feature'))
    }))
  
  print(head(models))
  
  
  # load shap values and clean up
  shapPaths <- list.files(here("data","results","shap.analysis", 'kernelshap_objects'), pattern = ".rds", full.names = TRUE)
  shapPaths <- shapPaths[str_detect(shapPaths, str_c(dataset, ".rds"))]
  shap <- map(shapPaths, \(x) {
    readRDS(x)
  })
  names(shap) <- shapPaths
  shap <- enframe(shap, name = 'raw_path', value = "shap_values")
  shap <- shap %>%
    mutate(path = str_split_fixed(raw_path, paste0(here()),n=2)[, 2]) %>%
    mutate(fold = str_split_fixed(path, "__", n = 4)[, 2]) %>%
    mutate(fold = str_replace(fold, "fold_id_", "")) %>%
    mutate(fold = as.numeric(fold)) %>%
    mutate(resampling = str_split_fixed(path, "__", n = 4)[, 1]) %>%
    mutate(resampling = str_split_fixed(resampling, "resamp_id_", 2)[,2]) %>%
    mutate(resampling = as.numeric(resampling)) %>%
    mutate(on = str_split_fixed(path, "__", n = 5)[, 3]) %>%
    mutate(on = str_replace(on, "on_", "")) %>%
    mutate(model_type = str_split_fixed(path, "__", n = 5)[, 4]) %>%
    mutate(model_type = str_replace(model_type, ".rds", "")) %>%
    mutate(model_type = str_replace(model_type, "which_model_", "")) %>%
    select(-raw_path) %>%
    identity()
  
  print(head(shap))
 
  # Load training/testing data
  profilePaths <- list.files(here("data","results","shap.analysis", 'fold_info'), pattern = ".tsv", full.names = TRUE)
  profilePaths <- profilePaths[str_detect(profilePaths, str_c(dataset, ".tsv"))]
  profiles <- map(profilePaths, \(x) {
    prof <- read.table(x, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) %>% as_tibble()
    colnames(prof) <- map_chr(colnames(prof), \(x) str_replace(x, "-", '.'))
    colnames(prof) <- ifelse(colnames(prof) == "51.20", "X51.20", colnames(prof))        
    return(prof)
  })
  names(profiles) <- profilePaths
  profiles <- enframe(profiles, name = 'raw_path', value = "profile")
  profiles <- profiles %>%
    mutate(path = str_split_fixed(raw_path, paste0(here()),n=2)[, 2]) %>%
    mutate(fold = str_split_fixed(path, "__", n = 4)[, 1]) %>%
    mutate(fold = str_replace(fold, ".*fold", "")) %>%
    mutate(fold = as.numeric(fold)) %>%
    mutate(resampling = str_split_fixed(path, "__", n = 4)[, 2]) %>%
    mutate(resampling = str_replace(resampling, "repeat_", "")) %>%
    mutate(resampling = str_replace(resampling, ".tsv", "")) %>%
    mutate(resampling = as.numeric(resampling)) %>%
    mutate(on = str_split_fixed(path, "_fold", n = 4)[, 1]) %>%
    mutate(on = sapply(str_split(on, "/"), function(x) tail(x, 1)))%>%
    mutate(on = case_when(
      on == "test_data" ~ "testing",
      on == "training_data" ~ "training"
    )) %>%
    select(-path, -raw_path)
  
  profiles <- profiles %>%
    mutate(profile = map(profile, ~ .x %>%
                           mutate(sampleID = as.character(sampleID))))  # ensure sampleID is character
  
  print(head(profiles))
  
  shap <- shap %>%
    left_join(profiles, by = c("fold", "resampling", 'on'))
  
  shap2 <- shap %>%
    mutate(shap_values_long = pmap(
      list(i = 1:n(), x = shap_values, p = profile),
      function(i, x, p) {
        tmp <- as.data.frame(x$S)
        tmp$sampleID <- p$sampleID
        tmp %>%
          pivot_longer(-sampleID, names_to = "feature", values_to = "shap_value")
      }
    ))
  
  print(head(shap2))
  
  shap2 <- shap2 %>%
    mutate(on = ifelse(on == "training", "training\nset", "test\nset")) %>%
    mutate(shap_values_long = pmap(
      list(i = 1:n(), x = shap_values, p = profile),
      function(i, x, p) {
        tmp <- as.data.frame(x$S)
        tmp$sampleID <- as.character(p$sampleID)  # ← convert here
        tmp %>%
          pivot_longer(-sampleID, names_to = "feature", values_to = "shap_value")
      }
    ))
  
  # Boxplots of shap values
  shap_tmp <- shap2 %>%
    select(fold, resampling, on, model_type, shap_values_long) %>%
    unnest(shap_values_long) %>%
    # I evalaute shap on training and testing
    # For testing folds, I get one shap value per model and resampling
    # for training folds, I get  4 shap values (in 5x cv) per model and resampling
    # In any case, take the median shap value for each sampleID
    group_by(sampleID, on, model_type, feature)
  
  # Understand mean/variance relationship of test-sample shap values over resampling rounds
  
  lel <- shap_tmp %>%
    filter(on == "test\nset") %>%
    group_by(sampleID, feature, model_type) %>%
    summarize(`mean(shap)\n(over resampled models)` = mean(shap_value), `var(shap)\n(over resampled models)` = var(shap_value))
  
  
  # IMPORTANT...
  shap_tmp <- shap_tmp %>% summarize(shap_value = median(shap_value))

  
  for (mt in model_types_to_evaluate) {
    
    shap_tmp2 <- shap_tmp %>%
      filter(model_type == mt)
    shap_tmp2$feature <- factor(shap_tmp2$feature, levels = shap_tmp2 %>%
                                  group_by(feature) %>%
                                  summarize(n = mean(abs(shap_value))) %>%
                                  arrange(desc(n)) %>%
                                  pull(feature))
    
    shap_tmp2 <- shap_tmp2 %>%
      inner_join(shap_tmp2 %>%
                   group_by(feature) %>%
                   summarize(n = mean(abs(shap_value))) %>%
                   arrange(desc(n)) %>%
                   head(10), by = "feature")
    

   
    shap_tmp3 <- shap_tmp %>%
      group_by(feature) %>%
      summarize(frac_pos = mean(shap_value > 0)) %>%
      arrange(desc(frac_pos))

  }
  
  
  # For simplicity, let's move on with RF on testing data
  mt <- "RF"
  more_plot_data_all <- shap_tmp %>%
    filter(model_type == mt, on == "test\nset")
  
  more_plot_data_all$feature <- factor(more_plot_data_all$feature, levels = more_plot_data_all %>%
                                         group_by(feature) %>%
                                         summarize(n = mean(abs(shap_value))) %>%
                                         arrange(desc(n)) %>%
                                         pull(feature))
  l <- levels(more_plot_data_all$feature)        
  
  more_plot_data_all <- more_plot_data_all %>%
    left_join(
      profiles %>%
        select(profile) %>%
        unnest() %>%
        distinct() %>%
        pivot_longer(-c(sampleID, Condition)) %>%
        rename(feature = name, feature_value = value), by = c('sampleID', "feature"))
  
  # Get spearman cors between genus abundance and shap to pimp the mean(abs(shap)) summary metric
  more_plot_data_all <- more_plot_data_all %>%
    left_join(more_plot_data_all %>%
                group_by(feature) %>%
                summarize(
                  spearman = cor(shap_value, feature_value, method = "spearman")
                ), by = 'feature') %>%
    mutate(spearman_sign = ifelse(spearman > 0, 1, -1))
  
  more_plot_data_all_to_write<- more_plot_data_all %>% mutate(on= str_replace(on,'\n', ' '))

  shap_file <-paste0(here("data","results","shap.analysis/"), dataset, "_median_shap_value.tsv")
  
  # Write the file
  write_tsv(more_plot_data_all_to_write, shap_file)

}
  

for (dataset in c(
  "Alldata",
  "EO-CRC",
  "LO-CRC"
)) {
  get_shap_median(dataset = dataset, label_case = "1", model_types_to_evaluate = c("RF"))
}
