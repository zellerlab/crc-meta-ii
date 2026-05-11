# Try from /g/scb2/zeller/karcher/mambaforge/envs/siamcat
library(mlr3)
library(SIAMCAT)
library(tidyverse)
library(here)

# For 3 types of models, this runs for around 10 minutes

model_type <- 'RF'

for (dataset in c("EO-CRC",'LO-CRC','Alldata', 'AD')) {

    dataset_name_path_map <- list(
      'EO-CRC'= here('data','results','Training.eo.lo.crc.rf.models.Rdata'),
      'LO-CRC'=  here('data','results','Training.eo.lo.crc.rf.models.Rdata'),
      "Alldata" = here('data','results','Training.unified.crc.model.Rdata'),
      'AD' = here('data','results','Training.ad.ctr.rf.model.Rdata')
    )

    dataset_name_model_object_map <- list(
      "Alldata" = 'models.all.rf',
      'EO-CRC'= 'models.eo.rf',
      'LO-CRC'='models.lo.rf',
      'AD' ='models.ad'
    )   

    extractFold <- function(siamcatObject, what_fold = NULL) {
        if (what_fold == "training") {
            sc_folds <- siamcatObject@data_split$training.folds
        } else if (what_fold == "testing") {
            sc_folds <- siamcatObject@data_split$test.folds
        } else {
            stop("what_fold must be either 'training' or 'testing'")
        }
        all_folds <- map2(sc_folds, 1:length(sc_folds), \(resamp, resamp_id) {
            folds <- map2(resamp, 1:length(resamp), \(fold, fold_id) {
                fold <- data.frame(sampleID = fold)
                fold$fold_id <- fold_id
                return(fold)
            })
            folds <- do.call('rbind', folds)
            folds$resamp_id <- resamp_id
            return(folds)
        })
        all_folds <- do.call('rbind', all_folds)
        return(all_folds)
    }

    load(dataset_name_path_map[[dataset]])
    assign("data_raw", get(dataset_name_model_object_map[[dataset]]), envir = .GlobalEnv)
    meta <- as.data.frame(data_raw@label$label)
    colnames(meta)[1] <- 'Condition'
    tmp <- as.data.frame(data_raw@label$info)
    colnames(tmp)[1] <- "code"
    tmp$trans <- rownames(tmp)
    meta <- meta %>%
        rownames_to_column('sampleID') %>%
        inner_join(tmp, by = c("Condition" = "code")) %>%
        mutate(Condition = trans) %>%
        select(-trans)
    norm_feat <- data_raw@norm_feat$norm.feat %>% t() %>% as.data.frame()
    models_list <- data_raw@model_list$models
    n <- names(models_list)
    models_list <- map(models_list, \(x) x$model)
    names(models_list) <- n

    fold_ids <- extractFold(data_raw, what_fold = "testing")
    fold_ids_training <- extractFold(data_raw, what_fold = "training")

    for (foldIndex in 1:data_raw@data_split$num.folds) {
        for (repeatIndex in 1:data_raw@data_split$num.resample) {
            print(repeatIndex)
            test_data <- fold_ids %>%
                filter(resamp_id == repeatIndex) %>%
                filter(fold_id == foldIndex) %>%
                select(-fold_id, -resamp_id) %>%
                inner_join(
                    norm_feat %>%
                        rownames_to_column('sampleID'),
                    by = 'sampleID'
                ) %>%
                inner_join(meta, by = 'sampleID') %>%
                relocate(sampleID, Condition)

            train_data <- fold_ids_training %>%
                filter(resamp_id == repeatIndex) %>%
                filter(fold_id == foldIndex) %>%
                select(-fold_id, -resamp_id) %>%
                inner_join(
                    norm_feat %>%
                        rownames_to_column('sampleID'),
                    by = 'sampleID'
                ) %>%
                inner_join(meta, by = 'sampleID') %>%
                relocate(sampleID, Condition)
            
            write_tsv(
              test_data,
              file.path(
                here("data", "results","shap.analysis","fold_info"),
                paste0("test_data_fold", foldIndex, "__repeat_", repeatIndex, "__", dataset, ".tsv")
              )
            )
            
            # Save training data
            write_tsv(
              train_data,
              file.path(
                here("data", "results","shap.analysis","fold_info"),
                paste0("training_data_fold", foldIndex, "__repeat_", repeatIndex, "__", dataset, ".tsv")
              )
            )
            
            # Save trained model
            write_rds(
              models_list[[paste0("cv_fold", foldIndex, "_rep", repeatIndex)]],
              file.path(
                here("data","results","shap.analysis", "models"),
                paste0("model__fold_id_", foldIndex, "__repeat_", repeatIndex, "__", dataset, "__", model_type, ".rds")
              )
            )

                 }
        }
}

