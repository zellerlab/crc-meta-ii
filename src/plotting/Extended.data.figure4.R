######################
# Extended data figure 4
######################

source(here('src','utils.R'))
params <- yaml::read_yaml(here("src", "parameters.yml"))
plotting <- params$plotting

######################
# Extended data figure 4a

meta.ad <-  read_tsv(here('data','Metadata.all.samples.tsv')) %>% 
  filter(Condition=='smallAD' | Condition=='AdvAD' | Condition=='AD'| Condition=='CTR') %>% 
  mutate(Condition= gsub('AdvAD','AD',Condition)) %>% 
  mutate(Condition= gsub('smallAD','AD',Condition))

all.data <- read.table(here('data','Relab.all.samples.tsv'),sep='\t',check.names = F)  %>%
  rownames_to_column('genus') %>%
  filter(genus!='unassigned') %>%
  column_to_rownames('genus')

all.data <- all.data[which(rowSums(all.data > 0) / ncol(all.data) > 0.1),]

all.data.ad <- all.data[,meta.ad$Sample_ID]


lmm.table.ad.vs.ctr <- run_lmem(
  data_df = all.data.ad,
  meta_df  = meta.ad, column_name='Condition', ref_group='CTR', feature_column_name = 'Taxon')


write.table(lmm.table.ad.vs.ctr, file= here('data','results' ,'lmm.tables.ad.ctr.tsv'), sep='\t', quote = F,row.names = F)

condition_colors <- params$plotting$condition_colors

volcano_ad_ctr <- plot_volcano(
  plot_df = lmm.table.ad.vs.ctr %>%
    select(Taxon, P.val, P.adj, Effect.size, pr.shift, pr.AD, pr.CTR, n.CTR, n.AD),
  group_case = 'AD',
  group_control = 'CTR',
  feature_column_name = 'Taxon',
  min_segment_length = 0.05, nudge_y = 0.05, max.overlaps = 5,
  color_vector = c(AD = condition_colors$AD, CTR = condition_colors$CTR, 'n.s.' = 'white') # Custom colors
) +
  xlab('Adenoma enrichment effect size')+
  theme(
    legend.box = 'horizontal',
    legend.spacing.y = unit(0.1, 'cm'),
    legend.position = c(0.35, .99),
    legend.key.size = unit(0.75, 'lines'),
    legend.justification = c(0, 1)
  )

ggsave(volcano_ad_ctr,file= here('figures','extended.data.figure4','Extended.Data.Figure4a.pdf'), height = 5, width = 5)

######################
# Extended data figure 4b

# Load all LOSO models from a single directory
rdata_dir <- here('data','results','scv.loso','ad.loso.test')

flist <- list.files(rdata_dir, pattern ='.Rdata', full.names = TRUE)

# Parse cohort and AD group
model_info <- str_match(basename(flist), "^Model\\.(.*?)\\.(AD|advAD|smallAD)\\.LOSO")[, 2:3]
cohort_names <- model_info[, 1]
ad_type <- model_info[, 2]

# Load all models into a named list
all.eval <- setNames(
  lapply(flist, function(x) {
    env <- new.env()
    load(x, envir = env)
    as.list(env)
  }),
  paste(cohort_names, ad_type, sep = "_")
)

# Helper function to extract model objects by AD type
get_models_by_type <- function(ad_type_label, loso_key) {
  relevant_names <- grep(paste0("_", ad_type_label, "$"), names(all.eval), value = TRUE)
  model_list <- list()
  for (n in relevant_names) {
    cohort <- str_replace(n, paste0("_", ad_type_label), "")
    model <- all.eval[[n]][[loso_key]][[cohort]]
    model_list[[cohort]] <- model
  }
  return(model_list)
}

# Extract models
models_ad     <- get_models_by_type("AD", "loso.eval.ad")
models_AdvAD  <- get_models_by_type("advAD", "loso.eval.adv.ad")
models_smallAD<- get_models_by_type("smallAD", "loso.eval.small.ad")

# Function to extract AUC and sample size
extract_model_metrics <- function(models, model_names) {
  if (length(models) != length(model_names)) {
    stop("Mismatch: Number of models and names are not equal.")
  }
  data.frame(
    Model = model_names,
    AUC = sapply(models, function(model) {
      if (!is.null(model@eval_data$auroc)) model@eval_data$auroc else NA
    }),
    Test_Samples = sapply(models, function(model) {
      if (!is.null(model@phyloseq@sam_data)) nrow(model@phyloseq@sam_data) else NA
    })
  )
}

# Extract metrics
model_names_ad     <- names(models_ad)
model_names_adv    <- names(models_AdvAD)
model_names_small  <- names(models_smallAD)

AD_metrics     <- extract_model_metrics(models_ad, model_names_ad)
AdvAD_metrics  <- extract_model_metrics(models_AdvAD, model_names_adv)
smallAD_metrics<- extract_model_metrics(models_smallAD, model_names_small)

# Merge wide format table
AD_loso_all_auc_df <- AD_metrics %>%
  left_join(AdvAD_metrics, by = "Model", suffix = c(".AD", ".AdvAD")) %>%
  left_join(smallAD_metrics, by = "Model") %>%
  rename(AUC.SmallAD = AUC, Test_Samples.SmallAD = Test_Samples)

# Save to file
# write.table(AD_loso_all_auc_df, here('data' ,'results', 'LOSO.AD.all.AUC.table.tsv', quote = FALSE, sep = '\t', row.names = FALSE)

# Reshape to long format for plotting
AD_loso_all_auc_long <- AD_loso_all_auc_df %>%
  select(Model, AUC.AD, AUC.AdvAD, AUC.SmallAD) %>%
  pivot_longer(cols = starts_with("AUC."),
               names_to = "AD_group", values_to = "AUC") %>%
  mutate(AD_group = recode(AD_group,
                           "AUC.AD" = "AD",
                           "AUC.AdvAD" = "AdvAD",
                           "AUC.SmallAD" = "SmallAD"))

AD_loso_all_samples_long <- AD_loso_all_auc_df %>%
  select(Model, Test_Samples.AD, Test_Samples.AdvAD, Test_Samples.SmallAD) %>%
  pivot_longer(cols = starts_with("Test_Samples."),
               names_to = "AD_group", values_to = "Test_Samples") %>%
  mutate(AD_group = recode(AD_group,
                           "Test_Samples.AD" = "AD",
                           "Test_Samples.AdvAD" = "AdvAD",
                           "Test_Samples.SmallAD" = "SmallAD"))

# Merge long format AUC and sample count
AD_loso_all_long <- left_join(AD_loso_all_auc_long, AD_loso_all_samples_long,
                              by = c("Model", "AD_group")) %>% drop_na(AUC)




# #count number of models generated for each dataset, if there are two these are same models so remove AD model
# the reason of this is I combined small ad and Adv ad samples as AD and trained a model 
model_auc_counts <- AD_loso_all_long %>% 
  group_by(Model) %>%
  summarize(n_auc = sum(!is.na(AUC)))

df_filtered <- AD_loso_all_long %>%
  left_join(model_auc_counts, by = "Model") %>%
  filter(!(n_auc == 2 & AD_group == "AD")) %>%
  select(-n_auc)  

colours_ad <- unlist(plotting$color_mapping_ad)

# Create the plot
auc_all_ad_plot<-ggplot(df_filtered, aes(x = Model, y = AUC, shape = AD_group, color=Model)) +
  geom_point(size = 5, alpha= 0.7) + 
  geom_hline(yintercept = 0.5, linetype='dashed', color='grey',size=0.8)+
  labs(title = "LOSO AD evaluations (trained on CRC)",
       x = "Model",
       y = "AUC",
       color = "Test",
       shape= 'Group') + 
  scale_color_manual(values=colours_ad) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


auc_all_ad_plot


ggsave(auc_all_ad_plot, file=here('figures','figure4','Figure4b.pdf'),width = 7, height = 7)


######################
# Extended data figure 4c

load(here('data','results','Training.ad.ctr.rf.model.Rdata'))

models <- list(models.ad)
labels <- c("AD Classifier CV:")
trained_on <- list(NULL)
colours <- c( "#FFBA08")

ad_auc_plot<- plot_roc_siamcat_models(models, labels, colours, trained_on,alpha=0.8)

ggsave(ad_auc_plot, file=here('figures','extended.data.figure4','Extended.Data.Figure4c.pdf'),width = 6, height = 6)


######################
# Extended data figure 4d

load(here('data','results','Training.unified.crc.model.Rdata'))

# calculate prevalance of top 4 genera for each stage

top_genera <- c("Fusobacterium", "Peptostreptococcus", "Parvimonas", "Porphyromonas")

# Get metadata and prediction matrix
all.meta <- models.all.rf@phyloseq@sam_data %>% data.frame()
all.meta$sampleID <- rownames(all.meta)

CvpredMatrix <- models.all.rf@pred_matrix %>%
  as.data.frame() %>%
  rownames_to_column('sampleID') %>%
  pivot_longer(-sampleID, names_to = "Fold", values_to = "Pred") %>%
  group_by(sampleID) %>%
  summarize(medianPredictionProb = median(Pred), .groups = "drop") %>%
  left_join(all.meta, by = "sampleID")

# Compute threshold (90% specificity) using ROC on CRC vs CTR
roc_curve <- roc(
  response = CvpredMatrix$Condition == "CRC",
  predictor = CvpredMatrix$medianPredictionProb)
threshold <- roc_curve$thresholds[which(roc_curve$specificities >= 0.9)[1]]

meta.ad <- read_tsv(here("data", "Metadata.all.samples.tsv")) %>%
  filter(Condition %in% c("smallAD", "AdvAD", "AD"))
meta.yachida <- meta.ad %>% filter(Cohort == "YachidaS_2019")

all.data <- read.table(here("data", "Relab.all.samples.tsv"), sep = '\t', check.names = FALSE) %>%
  rownames_to_column("genus") %>%
  filter(genus != "unassigned") %>%
  column_to_rownames("genus")

# Get 90th percentile threshold from CTR samples
ctr_data <- models.all.rf@phyloseq@otu_table %>%
  as.data.frame() %>%
  rownames_to_column("genus") %>%
  pivot_longer(-genus, names_to = "sampleID", values_to = "Relab") %>%
  filter(genus %in% top_genera) %>%
  left_join(all.meta, by = c("sampleID")) %>%
  filter(Condition == "CTR")

prevalence_threshold <- ctr_data %>%
  group_by(genus) %>%
  summarize(threshold = quantile(Relab, 0.9), .groups = "drop")


# Get CRC stage-wise prevalence for Yachida
crc_prevalence_data_stages_yachida <- models.all.rf@phyloseq@otu_table %>%
  as.data.frame() %>%
  rownames_to_column("genus") %>%
  pivot_longer(-genus, names_to = "sampleID", values_to = "Relab") %>%
  filter(genus %in% top_genera) %>%
  left_join(all.meta, by = "sampleID") %>%
  filter(Cohort == "YachidaS_2019") %>%
  filter(!(Condition == "CTR" & Stage == 0)) %>%
  left_join(prevalence_threshold, by = "genus") %>%
  group_by(Stage, genus) %>%
  summarize(Prevalence = mean(Relab >= threshold) * 100, .groups = "drop") %>%
  mutate(Stage=as.character(Stage))

# Get adenoma stage-wise prevalence for Yachida
ad_prevalence_data_stages_yachida <- all.data[, meta.yachida$Sample_ID] %>%
  rownames_to_column("genus") %>%
  pivot_longer(-genus, names_to = "SampleID", values_to = "Relab") %>%
  filter(genus %in% top_genera) %>%
  left_join(meta.ad %>% select(Stage = Condition, Sample_ID), by = c("SampleID"='Sample_ID' )) %>%
  left_join(prevalence_threshold, by = "genus") %>%
  group_by(Stage, genus) %>%
  summarize(Prevalence = mean(Relab >= threshold) * 100, .groups = "drop") 

# Combine prevalence data
prevalance_all_stages_yachida <- bind_rows(crc_prevalence_data_stages_yachida, ad_prevalence_data_stages_yachida) %>%
  drop_na() %>%
  mutate(Stage = case_when(
    Stage == 0 ~ "Stage 0",
    Stage == 1 ~ "Stage I",
    Stage == 2 ~ "Stage II",
    Stage == 3 ~ "Stage III",
    Stage == 4 ~ "Stage IV",
    TRUE ~ as.character(Stage)
  )) %>%
  mutate(Stage = factor(Stage, levels = c("smallAD", "Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV"))) %>%
  mutate(genus = factor(genus, levels = rev(top_genera)))


prev_heatmap_stages_yachida <- ggplot(prevalance_all_stages_yachida, aes(x = genus, y = Stage, fill = Prevalence)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", Prevalence)), vjust = 1, color = "black") +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 100)) +
  theme_minimal() +
  labs(x = "Tumor Stage", y = "Genus", fill = "Positivity (%)") +
  theme_paper +
  theme(
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 12),
    axis.ticks.x = element_line(),
    axis.text.x.bottom = element_text(size = 12, angle = 90),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.ticks.y = element_blank()
  ) +
  scale_y_discrete(limits = rev(levels(prevalance_all_stages_yachida$Stage)))

# Get TPR per stage (CRC)
df.plot.yachida <- CvpredMatrix %>%
  filter(Cohort == "YachidaS_2019") %>%
  mutate(Stage = case_when(
    Stage == 0 ~ "Stage 0",
    Stage == 1 ~ "Stage I",
    Stage == 2 ~ "Stage II",
    Stage == 3 ~ "Stage III",
    Stage == 4 ~ "Stage IV",
    TRUE ~ as.character(Stage)
  )) %>%
  mutate(Stage = if_else(Condition == "CTR" & Stage == "Stage 0", NA_character_, Stage)) %>%
  group_by(Stage) %>%
  summarize(tp = sum(medianPredictionProb > threshold), n = n(), .groups = "drop") %>%
  mutate(tpr = tp / n, StageLabel = paste0(Stage, "\n (N:", n, ")")) %>%
  drop_na(Stage)


# Get TPR per stage (Adenomas)
AD.predictions <- read_tsv(here("data", "results", "Adenoma.CRC.microbiome.signature.scores.tsv"))

Ad.df.plot.yachida <- AD.predictions %>%
  left_join(meta.ad %>% select(Sample_ID, Cohort,Condition), by = c("sampleID" = "Sample_ID")) %>%
  filter(Cohort == "YachidaS_2019") %>%
  rename(Stage = Condition) %>%
  group_by(Stage) %>%
  summarize(tp = sum(medianPredictionProb > threshold), n = n(), .groups = "drop") %>%
  mutate(tpr = tp / n, StageLabel = paste0(Stage, "\n (N:", n, ")")) %>%
  drop_na(Stage)

# Combine TPR data
df.plot.yachida <- bind_rows(df.plot.yachida, Ad.df.plot.yachida) %>%
  mutate(Stage = factor(Stage, levels = c("smallAD", "Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV")),
         StageLabel = factor(StageLabel, levels = rev(c(
           "smallAD\n (N:67)", "Stage 0\n (N:69)", "Stage I\n (N:73)",
           "Stage II\n (N:35)", "Stage III\n (N:51)", "Stage IV\n (N:20)"
         ))))

# Plot TPR barplot
g.y <- ggplot(df.plot.yachida, aes(x = tpr, y = StageLabel, fill = Stage)) +
  geom_bar(stat = "identity", color = "black") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1)) +
  xlab("") +
  ylab("True positive rate") +
  scale_fill_manual(values = c("#FFBA08", "#F48C06", "#E85D04", "#DC2F02", "#D00000", "#9D0208", "#6A040F"), guide = FALSE) +
  theme_paper +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    panel.grid.major.y = element_line(colour = "lightgrey")
  )

# Combine both panels
final_plot_yachida <- plot_grid(
  g.y, prev_heatmap_stages_yachida,
  align = "h",
  axis = "l",
  rel_heights = c(1, 1),
  rel_widths = c(0.6, 1),
  ncol = 2
)

# Save figure
ggsave(final_plot_yachida, filename = here("figures", "extended.data.figure4", "Extended.Data.Figure4d.pdf"),width = 6, height = 6)









