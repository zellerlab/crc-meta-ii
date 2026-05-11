######################
# Extended data figure 3
######################

source(here('src','utils.R'))
load(here('data','results', 'lmm.tables.16S.WGS.Rdata'))

lmm.table.motu.filtered <- lmm.table.motu %>%
  filter(P.adj < 0.05) %>%
  dplyr::rename(LME = Effect.size) %>%
  arrange(LME)

# Select top 10 positive and 10 negative species
top_positive_species <- lmm.table.motu.filtered %>%
  arrange(desc(LME)) %>%
  slice(1:10) %>%
  pull(species)

top_negative_species <- lmm.table.motu.filtered %>%
  arrange(LME) %>%
  slice(1:10) %>%
  pull(species)

# Combine top species
top_species <- c(top_negative_species, rev(top_positive_species))

# Filter for top species
filtered_table <- lmm.table.motu.filtered %>%
  filter(species %in% top_species) %>%
  mutate(species = factor(species, levels = top_species))

# Pivot longer to get other datasets' effect sizes
long_data <- filtered_table %>%
  pivot_longer(
    cols = starts_with('Effect.size_'),
    names_to = 'Dataset',
    values_to = 'Effect_size'
  ) %>%
  filter(!is.na(Effect_size)) %>%
  mutate(Dataset= str_remove(Dataset,'Effect.size_'))


family_palette <- colorRampPalette(brewer.pal(11, 'Spectral'))(13)

# Generate forest plot 

forest_plot_wgs <- ggplot() +
  geom_tile(
    data = filtered_table,
    aes(x = -3.2, y = species, fill = family),
    width = 0.2, height = 0.8, alpha = 0.6
  ) +
  geom_point(
    data = filtered_table,
    aes(x = LME, y = species),
    size = 5, color = '#C44600', alpha = 1, shape = 18
  ) +
  geom_point(
    data = long_data,
    aes(x = Effect_size, y = species, shape = Dataset),
    size = 3, color = '#C44600', alpha = 0.5, show.legend = TRUE
  ) +
  geom_line(
    data = filtered_table %>%
      select(lowCI = `conf.int_2.5 %`, highCI = `conf.int_97.5 %`, species) %>%
      pivot_longer(-c(species)),
    aes(x = value, y = as.numeric(species), group = species),
    color = 'black'
  ) +
  scale_fill_manual(values = family_palette, name = 'Family') +
  scale_shape_manual(
    values = c(0, 1, 3, 4, 7, 8, 9, 10, 12, 11, 13, 14, 16, 17, 18,
               25, 42, 35, 94, 95, 62, 43)
  ) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  coord_cartesian(xlim = c(-3, 3), clip = 'off') +
  geom_vline(xintercept = 0, linetype = 'dotted', size = 1) +
  xlab('Effect Size') +
  ylab('Species') +
  theme_paper +
  theme(
    legend.position = 'right',
    legend.box.spacing = unit(1.5, 'lines'),
    plot.margin = margin(t=10,r=10,b=10,l=20), # Adjust margin for family bars
    axis.text.y = element_text(margin=margin(r=20)) # Space for family bars
  )

# Save the plot
ggsave(plot = forest_plot_wgs, filename = here('figures','extended.data.figure3','Extended.Data.Figure3a.pdf'),width = 16, height = 14)


######################
# Extended data figure 3b

load(here('data','results','Training.wgs.rf.models.Rdata'))

models <- list(model.rf.wgs.motus, model.rf.wgs.genus)
labels <- c('CRC mOTU level CV:', 'CRC genus level CV:')
trained_on <- list(NULL, NULL)
colours <- c('darkorchid' ,'gray66')

roc_species_genus_comp <- plot_roc_siamcat_models(models, labels, colours, trained_on, alpha=0.5)

ggsave(plot = roc_species_genus_comp, filename = here('figures','extended.data.figure3','Extended.Data.Figure3b.pdf'),width = 7, height = 7)

######################
# Extended data figure 3c

# Load FOBT/FIT metadata
FOBT_meta <- read.table(here('data','Metadata.all.samples.balanced.tsv'), sep = '\t', header = T) %>%
  mutate(FIT= ifelse(FIT > 10, 'Positive', 'Negative')) %>% 
  filter(FOBT=='Positive'|FOBT=='Negative' | FIT=='Positive'|FIT=='Negative') %>% 
  filter(Condition=='CRC'| Condition=='CTR')

unique(FOBT_meta$Cohort)

# Load LOSO models to extract the number of samples on training and test sets

model_path <- here('data','results','scv.loso','crc.loso.test')

selected_models <- c('BaxterNT_2016', 'HanniganGD_2018', 'YangJ_2019', 'ZellerG_2014')

rdata_files <- list.files(model_path, pattern = '.Rdata$', full.names = TRUE)

selected_files <- rdata_files[grepl(paste(selected_models, collapse = '|'), rdata_files)]

loso.eval <- list()

for (file in selected_files) {
  load(file) 
  model_name <- gsub('.*Model\\.(.*?)\\.CRC.LOSO.*', '\\1', basename(file))  # Extract the model name 
  loso.eval[[model_name]] <- loso.eval.crc[[model_name]]
}

# Calculate FPR and TPR for FOBT/FIT status
result_list<- list()

for(i in 1:length(loso.eval)){
  
  metadata <- FOBT_meta %>%  filter(Cohort==names(loso.eval)[i])
  head(metadata)
  
  if (all(metadata$Cohort == 'BaxterNT_2016')) {  
    
    metadata <- metadata %>% data.frame() %>% 
      mutate(FIT = case_when( FIT == 'Positive' ~ 1, FIT =='Negative' ~ 0, TRUE ~ as.numeric(FIT)))
    
    predictor <- metadata$FIT
    
  } else {
    
    metadata <- metadata %>%data.frame() %>%  
      mutate(FOBT = case_when( FOBT == 'Positive' ~ 1,  FOBT == 'Negative' ~ 0,TRUE ~ as.numeric(FOBT))) 
    
    predictor <- metadata$FOBT
    
  }
  
  true_labels <- metadata$Condition
  
  true_labels <- ifelse(true_labels == 'CRC', 1, 0)
  
  valid_indices <- !is.na(predictor) & !is.na(true_labels)
  predictor <- predictor[valid_indices]
  true_labels <- true_labels[valid_indices]
  
  # Calculate TP, FP, TN, FN
  TP <- sum(predictor == 1 & true_labels == 1, na.rm = TRUE) 
  FP <- sum(predictor == 1 & true_labels == 0, na.rm = TRUE) 
  TN <- sum(predictor == 0 & true_labels == 0, na.rm = TRUE)  
  FN <- sum(predictor == 0 & true_labels == 1, na.rm = TRUE)
  
  TPR <- TP / (TP + FN)  
  FPR <- FP / (FP + TN)  
  
  # Extract cohort name
  cohort_name <- unique(metadata$Cohort)
  
  # Extract model name
  model_name <- names(loso.eval)[i]
  
  # Store results in a structured list
  result_list[[i]] <- list(
    Cohort = cohort_name,
    Model = model_name,
    TPR = TPR,
    FPR = FPR
  )
}

result_df <- do.call(rbind, lapply(result_list, as.data.frame)) 

#Load trained siamcat object for loso models (to get training sample size)

model_path <- here('data','results','scv.loso','crc.loso.train')

rdata_files <- list.files(model_path, pattern = '.Rdata$', full.names = TRUE)

selected_files <- rdata_files[grepl(paste(selected_models, collapse = '|'), rdata_files)]

loso.training <- list()

for (file in selected_files) {
  data_env <- new.env()
  load(file, envir = data_env)  
  
  model_name <- gsub('.*Model\\.(.*?)\\.CRC.LOSO.*', '\\1', basename(file))
  
  obj_name <- ls(data_env)[1]  
  loso.training[[model_name]] <- get(obj_name, envir = data_env)
}

print(loso.training)

# Load also unified CRC model
load(here('data','results','Training.unified.crc.model.Rdata'))

models <- list(loso.eval$BaxterNT_2016, loso.eval$HanniganGD_2018, loso.eval$YangJ_2019, loso.eval$ZellerG_2014, models.all.rf )
labels <- c('BaxterNT_2016', 'HanniganGD_2018', 'YangJ_2019','ZellerG_2014', 'Universal CRC classifier')
trained_on <- list(loso.training[[1]]$BaxterNT_2016_LOSO,loso.training[[2]]$HanniganGD_2018_LOSO,loso.training[[3]]$YangJ_2019_LOSO,loso.training[[4]]$ZellerG_2014_LOSO, NULL)
colours <- c('#77AADD',  '#FFAABB', '#99DDFF', '#44BB99', 'black')
linetypes <- c('solid', 'solid', 'solid', 'solid', '11')
FOBT_tpr <- result_df

plot<-plot_roc_siamcat_models_with_fobt(models = models , labels = labels, linetypes = linetypes,colours = colours, FOBT_tpr =FOBT_tpr ,trained_on = trained_on )


ggsave(plot = plot, filename = here('figures','extended.data.figure3','Extended.Data.Figure3c.pdf'),width = 7, height = 7)

######################
# Extended data figure 3d:  Heatmap of study to study transfer(SST) models  

# Load SST models 
load(here('data','results','scv.loso','crc.scv.sst','Models.SCV.SST.Rdata'))

SCV.CRC<- models.scv
SST.eval.CRC<- holdout.evaluation.crc

# Get AUC's and number of samples
SST_metrics <- get_siamcat_holdout_metrics(siamcat_list_cv = SCV.CRC, siamcat_list_holdout = SST.eval.CRC)

SST_metrics$Train<- factor(SST_metrics$Train, levels = rev(unique(SST_metrics$Train)))

SST_metrics$Test<- factor(SST_metrics$Test, levels = (unique(SST_metrics$Train)))

# Plot heatmap
sst_heatmap <- ggplot(SST_metrics, aes(x = Test, y = Train, fill = AUC)) +
  geom_tile(color = 'grey90') +
  geom_tile(data = SST_metrics %>% filter(Type == 'CV'),
            aes(x = Test, y = Train),
            color = 'white', size = 1.2, fill = NA) +
  geom_text(aes(label = round(AUC, 2)), size = 4, color = 'white') +
  scale_fill_viridis_c(name = 'AUC', option = 'viridis', direction = -1) +
  labs(title = 'Study to study transfer evaluation',
       x = 'Test Dataset', y = 'Train Dataset') +
  theme_paper +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())


ggsave(sst_heatmap, filename = here('figures','extended.data.figure3','Extended.Data.Figure3d.pdf'), height = 15, width = 15 )

























