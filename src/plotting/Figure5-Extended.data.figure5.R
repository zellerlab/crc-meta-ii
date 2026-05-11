######################
# Figure 5 and Extended data figure 5: Correlation of fiber intake and CRC microbiome signature score in different cohorts
######################

params <- yaml::read_yaml(here("src", "parameters.yml"))
plotting <- params$plotting
source(here('src','utils.R'))

######################
# Figure 5a: Correlation of fiber intake (YachidaS_2019) and CRC microbiome signature score

# Load universal CRC model
load(here('data','results','Training.unified.crc.model.Rdata'))

CVModelPredictions <- pred_matrix(models.all.rf) 

cvPredMatrix <- CVModelPredictions %>%
  as.data.frame() %>%
  rownames_to_column('sampleID') %>%
  as_tibble() %>%
  pivot_longer(-sampleID) %>%
  group_by(sampleID) %>%
  summarize(medianPredictionProb = median(value)) %>%
  left_join(models.all.rf@label$label %>% as.data.frame() %>% dplyr::rename(label = '.') %>% rownames_to_column('sampleID'), by = 'sampleID') %>%
  mutate(label = ifelse(label == '1', "CRC", "Control")) %>%
  mutate(type = 'reference (CV)') %>% 
  mutate(sampleID=as.character(sampleID))

# Load metadata to get fiber intake info for YachidaS_2019
all.meta<- read_tsv(here('data','Metadata.all.samples.balanced.tsv')) %>% filter(Condition=='CRC'| Condition=='CTR') %>% 
  as.data.frame()  %>% rename(sampleID=Sample_ID) %>% select(Cohort, sampleID, Condition)

cvPredMatrix <- cvPredMatrix %>% left_join(all.meta, by = "sampleID") %>% filter(Cohort=='YachidaS_2019')

Yachida_fiber<- read_xlsx(here('data','fiber_data','ijc34398-sup-0003-tables3.xlsx'),sheet = 1,skip = 1)  %>%
  rename('Sample_Name'= 'Patient ID') %>% 
  mutate(Sample_Name=as.character(Sample_Name)) %>% 
  left_join(read.table(file=here('data','fiber_data','SraRunTable_Yachida_2019.txt'), sep = ',',header=T) %>% select(Run,Sample_Name)  %>%
              mutate(Sample_Name=as.character(Sample_Name))
            , by='Sample_Name') %>% dplyr::rename(sampleID= Run) %>% select(sampleID, Fiber) 

Yachida_fiber<- Yachida_fiber %>% 
  left_join(cvPredMatrix %>% 
              select(sampleID, medianPredictionProb,Condition), by='sampleID') %>% 
  drop_na(medianPredictionProb) %>% 
  filter(Fiber <100)     #seems like one samples with really high Fiber intake (131) so it seems like typo 

# Calculate Spearman correlations for each group
cor_all <- cor.test(~ medianPredictionProb + Fiber, data = Yachida_fiber, method = "spearman")
cor_crc <- cor.test(~ medianPredictionProb + Fiber, data = filter(Yachida_fiber, Condition == "CRC"), method = "spearman")
cor_ctr <- cor.test(~ medianPredictionProb + Fiber, data = filter(Yachida_fiber, Condition == "CTR"), method = "spearman")

# Format results
label_all <- paste0("rho = ", round(cor_all$estimate, 2), ", p = ", signif(cor_all$p.value, 3))
label_crc <- paste0("rho = ", round(cor_crc$estimate, 2), ", p = ", signif(cor_crc$p.value, 3))
label_ctr <- paste0("rho = ", round(cor_ctr$estimate, 2), ", p = ", signif(cor_ctr$p.value, 3))

# Define label positions
x_pos <- max(Yachida_fiber$Fiber, na.rm = TRUE) * 0.6  # fixed x
y_max <- max(Yachida_fiber$medianPredictionProb, na.rm = TRUE)
spacing <- 0.05  

# Y positions for stacked text
y_pos_all <- y_max
y_pos_crc <- y_max - spacing
y_pos_ctr <- y_max - 2 * spacing

condition_colors <- plotting$condition_colors

fiber_corr_plot <- ggplot(Yachida_fiber, aes(x = Fiber, y = medianPredictionProb)) +
  geom_point(
    aes(fill = Condition),
    shape = 21,
    size = 4,
    color = "black",
    alpha = 0.7
  ) +
  geom_smooth(aes(color = Condition), method = "lm", se = FALSE) +
  geom_smooth(aes(group = 1), method = "lm", color = "black", se = FALSE) +
  scale_fill_manual(values = condition_colors) +
  scale_color_manual(values = condition_colors) +
  annotate("text", x = x_pos, y = y_pos_all, label = label_all, hjust = 0, color = "black", size = 4) +
  annotate("text", x = x_pos, y = y_pos_crc, label = label_crc, hjust = 0, color = condition_colors$CRC, size = 4) +
  annotate("text", x = x_pos, y = y_pos_ctr, label = label_ctr, hjust = 0, color = condition_colors$CTR, size = 4) +
  xlim(0, 50) +
  ylim(0.1, 0.95) +
  labs(
    title = "YachidaS_2019",
    x = "Fiber intake (estimated g/day)",
    y = "CRC microbiome siganture score",
    fill = "Group",     # fill replaces color in legend title
    color = "Group"
  ) +
  theme_paper +
  theme(
    legend.position = c(0.98, 0.02),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
  )

ggsave(plot=fiber_corr_plot, here('figures','figure5','Figure5a.pdf'), width = 6, height = 6)

######################
# Figure 5b: Correlation of fiber intake (SowahSA_2022) and CRC microbiome signature score

# Load all scores for fiber studies
Fiber_scores <- read_tsv(here('data','results','CRC_scores_fiber_datasets.tsv'))

Sowah_meta <- read.table(here('data','fiber_data','SraRunTable_Sowah.csv'), sep=',', header = T) %>%
  select(Run,Sample.Name, Sample_name) %>%
  mutate(Time_points=  str_split_fixed(Sample_name,'_',4)[,2],
         Batch=  str_split_fixed(Sample_name,'_',4)[,3],
         Patient_ID = str_split_fixed(Sample_name,'_',4)[,1],
         Identifier = str_split_fixed(Sample_name,'_',4)[,4]) 


Sowah_scores <- Fiber_scores %>% filter(label == 'SowahSA_2022') %>% 
  left_join(readRDS(here('data', 'fiber_data','SowahSA_2022_res_IDTaxa.rds')) %>%
              as.data.frame() %>%
              rownames_to_column('taxon') %>%
              pivot_longer(-taxon) %>%
              mutate(genus = str_split_fixed(taxon, "\\|", n = 7)[, 6],
                     genus = str_replace(genus, "g__", "")) %>%
              dplyr::rename(sampleID = name) %>%
              mutate(genus = ifelse(genus == "", "unassigned", genus)) %>%
              filter(genus != 'unassigned') %>%
              group_by(sampleID, genus) %>%
              summarise(value = sum(value), .groups = "drop") %>%
              select(-genus) %>% 
              group_by(sampleID) %>%
              mutate(value = sum(value)) %>% 
              ungroup() %>% unique(), by='sampleID') %>% 
  left_join(Sowah_meta, by=c('sampleID'='Run')) 

# Load also sample groups used in this study
Groups <- read_xlsx(here('data', 'fiber_data', '13073_2022_1030_MOESM4_ESM.xlsx'), sheet = 1) %>% 
  select(Patient_ID=Sample_id, Treatment_group) %>% 
  mutate(Patient_ID= str_split_fixed(Patient_ID,'_',2)[,1] ) %>% 
  distinct()

Sowah_scores <- Sowah_scores %>% 
  left_join(Groups, by='Patient_ID') %>% 
  mutate(Time_points= factor(Time_points,levels= c('T0','T1','T2','T3'))) %>% 
  drop_na(Treatment_group)

Sowah_scores <- Sowah_scores %>%
  group_by(Patient_ID, Time_points) %>%
  slice_max(value, n = 1, with_ties = FALSE) %>%  # one sample per timepoint
  ungroup()

# Load dietary metadata
diet_long<-read.table(here('data','fiber_data', 'Diet all long.csv'), sep=',',header = T)

diet_labels<- read_csv(here('data','fiber_data','Diet labels all.csv'))

label_map <- setNames(diet_labels$Description, diet_labels$Label)

# Replace column names if they are in label_map
colnames(diet_long) <- ifelse(
  colnames(diet_long) %in% names(label_map),
  label_map[colnames(diet_long)],
  colnames(diet_long)  # keep original if no match
)

# First  take the baseline of samples and correlate with fiber intake 
Sowah_all_diet<- Sowah_scores %>%
  mutate(Patient_ID=as.integer(Patient_ID)) %>%
  select(medianPredictionProb, Time_points,Patient_ID, Treatment_group,Batch) %>% 
  left_join(diet_long, by=c('Patient_ID'='HELENA_ID', 'Time_points'='Timepoint'))

# Select baseline samples 
Sowah_all_diet_T0 <- Sowah_all_diet %>% filter(Time_points=='T0')

plot_data <- Sowah_all_diet_T0 %>% 
  select(Patient_ID, medianPredictionProb, Total_fiber) %>%
  rename("Fiber_total"='Total_fiber') %>%
  pivot_longer(cols = starts_with("Fiber_"), names_to = "Fiber_Type", values_to = "Fiber_Value")

corr_stats <- plot_data %>%
  group_by(Fiber_Type) %>%
  summarise(
    rho = cor(medianPredictionProb, Fiber_Value, method = "spearman", use = "complete.obs"),
    p_value = cor.test(medianPredictionProb, Fiber_Value, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    subtitle = paste0("rho = ", round(rho, 2),
                      ",p = ", format.pval(p_value, digits = 2, eps = .001))
  )

# Plot the correlation for total fiber 
Sowah_baseline_cor <- ggplot(plot_data, aes(x = Fiber_Value, y = medianPredictionProb)) +
  geom_point(
    shape = 21,             
    size = 4,              
    fill = "grey55",         
    color = "black",        
    alpha = 0.7             
  ) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate("text", x = 35, y = 0.9, label = corr_stats$subtitle, hjust = 0, color = "black", size = 4) +
  xlim(0, 50) +
  ylim(0.1, 0.95) +
  labs(
    x = "Fiber intake (estimated g/day)",
    y = "CRC microbiome signature score",
    title = "SowahSA_2021"
  ) +
  theme_presentation()

ggsave(plot=Sowah_baseline_cor, here('figures','figure5','Figure5b.pdf'), width = 6, height = 6)

######################
# Extended Data Figure 5a: Other diet variables correlated with CRC microbiome signature score at baseline in SowahSA_2021

# Select baseline (T0) samples 
Sowah_all_diet_T0 <- Sowah_all_diet %>% filter(Time_points == 'T0')

# Put all diet vars into long format
diet_vars <- setdiff(names(Sowah_all_diet_T0),
                     c("medianPredictionProb","Time_points","Patient_ID","Treatment_group","Batch"))

plot_data <- Sowah_all_diet_T0 %>%
  select(Patient_ID, medianPredictionProb, all_of(diet_vars)) %>%
  pivot_longer(cols = all_of(diet_vars),
               names_to = "DietVar", values_to = "DietValue")

# Spearman correlations per diet variable
corr_stats <- plot_data %>%
  group_by(DietVar) %>%
  summarise(
    n = sum(!is.na(medianPredictionProb) & !is.na(DietValue)),
    rho = ifelse(n >= 3,
                 cor(medianPredictionProb, DietValue,
                     method = "spearman", use = "complete.obs"),
                 NA_real_),
    p_value = ifelse(n >= 3,
                     suppressWarnings(cor.test(medianPredictionProb, DietValue,
                                               method = "spearman")$p.value),
                     NA_real_),
    .groups = "drop"
  ) %>%
  filter(!is.na(p_value)) %>%
  mutate(
    signif = p_value <= 0.05,
    subtitle = paste0("rho = ", round(rho, 2),
                      ", p = ", format.pval(p_value, digits = 2, eps = .001))
  )

# Filter to significant variables and other diet vars of interest 
sig_vars <- corr_stats %>%
  filter(signif | DietVar=='Fat'| DietVar=='Menu_comp_mostly_animal_prod'|
         DietVar=='Sausages_and_processed_meat'| DietVar=='Cholesterol') %>%
  arrange(p_value) %>%    # sorted by significance
  filter(!str_detect(DietVar, 'fiber|Fiber'))

# Inspect which variables are significant
sig_vars %>% select(DietVar, rho, p_value)

# Plot function 
make_sowah_plot <- function(dvar) {
  df <- plot_data %>% filter(DietVar == dvar)
  lab <- (sig_vars %>% filter(DietVar == dvar) %>% pull(subtitle))[1]
  
  # annotation positions based on observed ranges
  xr <- range(df$DietValue, na.rm = TRUE)
  yr <- range(df$medianPredictionProb, na.rm = TRUE)
  x_pos <- xr[1] + 0.6 * diff(xr)
  y_pos <- yr[2]
  
  ggplot(df, aes(x = DietValue, y = medianPredictionProb)) +
    geom_point(shape = 21, size = 4, fill = "grey55", color = "black", alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    annotate("text", x = x_pos, y = y_pos, label = lab,
             hjust = 0, color = "black", size = 4) +
    labs(
      x = paste0(dvar),
      y = NULL,  # no per-plot y label
      title = paste0("SowahSA_2021\n", dvar)
    ) +
    theme_presentation() +
    theme(axis.title.y = element_blank())
}

plots_list <- purrr::map(sig_vars$DietVar, make_sowah_plot) %>%
  rlang::set_names(sig_vars$DietVar)

# Arrange grid of plots (4 columns)
sig_grid <- cowplot::plot_grid(
  plotlist = plots_list,
  ncol = 4,          
  labels = NULL,
  align = "hv"
)

# Shared y-axis label
shared_y <- cowplot::ggdraw() +
  cowplot::draw_label(
    "CRC microbiome signature score",
    angle = 90, vjust = 0.5, fontface = "bold"
  )

# Combine shared y-label and the grid of plots
sig_grid_with_ylabel <- cowplot::plot_grid(
  shared_y, sig_grid,
  ncol = 2,
  rel_widths = c(0.06, 1)
)

# Save the combined plot
ggsave(
  here("figures","Extended.Data.Figure5","Extended.Data.Figure5a.pdf"),
  sig_grid_with_ylabel,
  width = 18, height = 16, dpi = 300
)

######################
# Figure 5c: Habitual diet and CRC microbiome signature score in the DeFilippis_2016 cohort

# Load metadata of the study
DeFilippis_meta<- read.table(here('data','fiber_data','SraRunTable-DeFilippis_2016.csv'),sep=',', header = T) %>% 
  filter(env_biome=='feces') %>% 
  select(Run, Label, Sample.Name) %>% 
  mutate(host_diet = str_split_fixed(Label, " subject", 2)[,1],
         host_diet = str_split_fixed(host_diet, "from ",2) [,2]) %>% 
  mutate(host_diet_general = ifelse(host_diet=='omnivorous','omnivorous','vegetarian/vegan'))

# Load CRC microbiome signature scores and merge with metadata
DeFilippis_scores <- Fiber_scores %>% filter(label == 'DeFilippis_2016') %>%
  left_join(DeFilippis_meta %>% select(Run,host_diet_general,Sample.Name), by = c("sampleID" = "Run")) %>%
  mutate(
    host_diet_general = factor(host_diet_general, levels = c("omnivorous", "vegetarian/vegan"))) %>% 
  mutate(
    x = as.numeric(host_diet_general),  # ensure order
    x_jittered = x + runif(n(), -0.1, 0.1)
  )

# Calculate p value using Wilcoxon test between the two diet groups
pval_df <- DeFilippis_scores %>%
  wilcox_test(medianPredictionProb ~ host_diet_general) %>%
  adjust_pvalue(method = "BH") %>%        
  add_xy_position(x = "host_diet_general") %>%
  mutate(
    label = paste0("p = ", signif(p, 3)),  # now works
    y.position = y.position - 0.02,
    x = (as.numeric(factor(group1, 
                           levels = levels(DeFilippis_scores$host_diet_general))) +
           as.numeric(factor(group2, 
                             levels = levels(DeFilippis_scores$host_diet_general)))) / 2
  )

# Get group sizes
group.n_DeFilippis <- DeFilippis_scores %>%
  count(host_diet_general)

# Plot
DeFlippis_habitual<- ggplot() +
  geom_boxplot(
    data = DeFilippis_scores,
    aes(x = host_diet_general, y = medianPredictionProb),
    width = 0.5, outlier.shape = NA
  ) +
  geom_jitter(
    data = DeFilippis_scores,
    aes(x = host_diet_general, y = medianPredictionProb),
    size = 4, shape = 21, fill = "grey80", color = "black",
    alpha = 0.5, width = 0.3
  ) +
  geom_text(
    data = group.n_DeFilippis,
    aes(x = host_diet_general,
        y = max(DeFilippis_scores$medianPredictionProb) * 1.05,
        label = paste0("(N = ", n, ")")),
    inherit.aes = FALSE,
    vjust = 0
  ) +
  geom_text(
    data = pval_df,
    aes(x = x, y = y.position, label = label),
    inherit.aes = FALSE,
    vjust = 0, size = 4
  ) +
  scale_x_discrete(drop = TRUE) +
  labs(y = "CRC microbiome signature score") +
  theme_presentation() +
  theme(
    legend.position = 'top',
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank()
  ) +
  coord_cartesian(ylim = c(NA, max(DeFilippis_scores$medianPredictionProb) * 1.15))

ggsave(plot=DeFlippis_habitual, here('figures','figure5','Figure5c.pdf'),width = 3, height = 6, dpi = 300)

######################
# Figure 5d: CRC microbiome signature score changes with fiber enriched Mediterranean diet (BarberC_2021)

# Load metadata of the study
Barber_meta<-read.table(here('data','fiber_data','BarberC2021.tsv'), sep='\t', header = T) %>% 
  select(c(sampleID=run_accession, sample_title)) %>% 
  mutate(time_points= str_split_fixed(sample_title,' ',3)[,3]) %>%
  mutate(individual= str_split_fixed(sample_title,' ',3)[,2]) 

# Load CRC microbiome signature scores and merge with metadata
Barber_scores <- Fiber_scores %>% filter(label == 'BarberC_2021') %>%
  left_join(Barber_meta, by='sampleID') %>%
  mutate(individual = as.factor(individual),
         time_points = factor(time_points, levels = c('WD', 'FMD'))) 

Barber_scores_paired <- Barber_scores %>%
  group_by(individual) %>%
  filter(n_distinct(time_points) == 2) %>%
  ungroup()

Barber_wide_scores <- Barber_scores_paired %>%
  select(individual, time_points, medianPredictionProb) %>%
  pivot_wider(names_from = time_points, values_from = medianPredictionProb)

# Perform Wilcoxon signed-rank test for paired samples
wilcox_res <- wilcox.test(
  Barber_wide_scores$WD,
  Barber_wide_scores$FMD,
  paired = TRUE,
  exact = FALSE
)

# Prepare data for stat_pvalue_manual
pval_df <- data.frame(
  group1 = "WD",
  group2 = "FMD",
  label = paste0("p = ", round(signif(wilcox_res$p.value, 3),3)),
  y.position = max(Barber_scores_paired$medianPredictionProb, na.rm = TRUE) * 1.05
)

Barber_scores_paired <- Barber_scores_paired %>%
  mutate(
    x = as.numeric(factor(time_points, levels = c("WD", "FMD"))),  # ensure order
    x_jittered = x + runif(n(), -0.1, 0.1)
  )

# Plot
Barber_scores_plot <- ggplot() +
  geom_boxplot(data = Barber_scores_paired, aes(x = time_points, y = medianPredictionProb), 
               width = 0.5, outlier.shape = NA) +
  geom_point(data = Barber_scores_paired,
             aes(x = x_jittered, y = medianPredictionProb),
             size = 4,
             shape = 21,
             fill = "grey",
             color = "black",
             alpha = 0.7) +
  stat_pvalue_manual(
    pval_df,
    label = "label",  
    tip.length = 0.01,
    size = 4
  ) +
  labs(y = "CRC microbiome signature score") +
  theme_presentation() +
  ylim(c(0,0.65)) +
  theme(legend.position = 'top', 
        axis.title.x = element_blank())

ggsave(plot=Barber_scores_plot, here('figures','figure5','Figure5d.pdf'),width = 3, height = 6, dpi = 300)

######################
# Extended Data Figure 5d : Changes in specific genera with fiber enriched Mediterranean diet (BarberC_2021)

selected_genera <- c("CAG-41","Agathobaculum", "Lachnospira")

mOTUS.3.1.metadata<- read_tsv(here('data','mOTUs3.1.0.genome_metadata_edited.tsv')) 
# this is the metadata file for mOTUs3.1.0, which contains taxonomic annotation for the mOTUs. 
# We need this to get the genus level annotation for the mOTUs in BarberC_2021 dataset, which is profiled with mOTUs3.1.0 profiler and thus has mOTUs as features.

Barber_relab <- readRDS(here('data','fiber_data','Barber_2021_res_mOTUs.rds')) %>%
  as.data.frame() %>%
  rownames_to_column('taxon') %>%
  pivot_longer(-taxon) %>%
  group_by(name) %>%
  mutate(mOTUs = str_split_fixed(taxon, '\\|', n = 8)[, 8]) %>%
  mutate(mOTUs= ifelse(taxon=='unassigned','unassigned',mOTUs)) %>% 
  dplyr::rename(sampleID = name) %>%
  group_by(sampleID, mOTUs) %>%
  filter(mOTUs != '') %>%
  pivot_wider(id_cols = mOTUs, names_from = sampleID, values_from = value) %>%
  as.data.frame %>%
  pivot_longer(-mOTUs) %>% 
  dplyr::rename('sampleID'= 'name', 'count'='value') %>% 
  left_join(mOTUS.3.1.metadata, by= c('mOTUs' = 'GENOME_SOURCE')) %>% 
  filter(DOMAIN!='Unclassified') %>% 
  dplyr::rename('genus'= 'GTDB_genus') %>% 
  mutate(genus=ifelse(mOTUs=='unassigned', 'unassigned',genus)) %>% 
  filter(mOTUs!='meta_mOTU_v31_12968') %>% # these specific motus don't have any entry in mOTUs metadata, so we can't assign them to any genus, so we just exclude this motu from the analysis.
  select(-DOMAIN) %>% 
  group_by(sampleID, genus) %>% 
  dplyr::summarise(count=sum(count)) %>% 
  pivot_wider(id_cols = genus , names_from = sampleID, values_from = count, values_fill = 0)  %>% 
  filter(genus!='unassigned') %>% 
  column_to_rownames('genus')  %>%
  mutate(across(everything(), ~ . / sum(.))) 

Barber_relab_selected <- Barber_relab %>%
  select(Barber_scores_paired$sampleID) %>%
  rownames_to_column('genus') %>%
  pivot_longer(-genus, names_to = "SampleID", values_to = "Relab") %>%
  filter(genus %in% selected_genera) %>%
  left_join(Barber_scores_paired %>% select(time_points, sampleID, individual), 
            by = c('SampleID' = 'sampleID')) %>%
  mutate(logRelab = log10(Relab + 1e-4)) %>% 
  mutate(genus= factor(genus,levels=c('CAG-41','Agathobaculum','Lachnospira')))  %>% 
  select(-SampleID)

df_wide <- Barber_relab_selected %>%
  pivot_wider(
    id_cols = c(individual, genus),
    names_from = time_points,
    values_from = logRelab
  )

# Now, for each genus, perform the paired Wilcoxon test
wilcox_results <- df_wide %>%
  group_by(genus) %>%
  summarise(
    p_value = wilcox.test(FMD, WD, paired = TRUE, exact = FALSE)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))


df_log_annot <- Barber_relab_selected %>%
  left_join(wilcox_results, by = "genus")

annotation_df <- df_log_annot %>%
  group_by(genus) %>%
  summarise(
    y_pos = max(logRelab, na.rm = TRUE) + 0.2,
    p_label = paste0("p = ", signif(first(p_adj), 2))
  )

light_colors <- c("#d4e6f1", "#f9e79f", "#fadbd8")

# Plot the abundance changes with boxplots and jittered points, faceted by genus, and annotated with p-values
abundance_plot <- ggplot(df_log_annot, aes(x = time_points, y = logRelab)) +
  geom_boxplot(outlier.shape = NA, color = "black", fill = NA, size = 0.3) + 
  geom_jitter(aes(fill = genus), width = 0.1, size = 3, alpha = 0.8, shape = 21, color = 'black') +
  geom_text(data = annotation_df, aes(x = 1.5, y = y_pos, label = p_label),
            inherit.aes = FALSE, size = 3, hjust = 0.5, fontface = "bold") +
  facet_wrap(~genus, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = light_colors) +
  labs(y = "log10(relative abundance + 1e-4)") +
  theme_paper +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    axis.title.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "none"
  )

ggsave(plot=abundance_plot, here('figures','extended.data.figure5','Extended.Data.Figure5d.pdf'),width = 4, height = 6, dpi = 300)

######################
# Figure 5f: CRC microbiome signature score changes with longitudinal fiber enriched  diet (Delannoy-Bruno_2021a)

Delannoy_meta<-read.table(here('data','fiber_data','DelannoyBruno2021.csv'), sep=',',header = T) %>%
  filter(host_scientific_name !='Mus musculus') %>% #exclude mouse samples 
  filter(Assay.Type !='AMPLICON') %>%
  select(Run, Sample_name, Diet, host_scientific_name) 

Delannoy_scores<- Fiber_scores %>% filter(label == 'Delannoy-Bruno_2021') %>%
 left_join(Delannoy_meta, by=c('sampleID'='Run'))  %>%
  mutate(
    Group = ifelse(str_detect(Sample_name, 'HS'), 'HS', 'PF'),
    Diet = ifelse(str_detect(Sample_name, 'HS'), '2Fib4Fib', Diet),
  ) %>% 
  mutate(individual=case_when(str_detect(Sample_name, 'HS2') ~str_split_fixed(Sample_name, '_',4)[,3],
                              str_detect(Sample_name, 'PF') ~str_split_fixed(Sample_name, '_',3)[,2])) %>% 
  mutate(individual=str_remove(individual,'S.')) %>% 
  mutate(day=case_when(str_detect(Sample_name, 'HS2') ~str_split_fixed(Sample_name, '_',4)[,4],
                       str_detect(Sample_name, 'PF') ~str_split_fixed(Sample_name, '_',3)[,3])) %>% 
  mutate(day= str_remove(day, 'Day')) %>%
  mutate(day = factor(day, levels = sort(as.numeric(unique(day))))) %>% 
  mutate(individual = as.character(as.numeric(individual))) %>%
  mutate(individual = factor(individual, levels = sort(as.numeric(unique(individual)))))

Delannoy_2021_scores_PF <- Delannoy_scores %>% 
  filter(Group=='PF') %>%
  mutate(fiber_amount = case_when(str_detect(Diet,'HiSF-LoFV') ~ 0, 
                                  str_detect(Diet, '1 snack') ~ 1,
                                  str_detect(Diet, '2 snack') ~ 2,
                                  str_detect(Diet, '3 snack') ~ 3)) %>% 
  arrange(individual, day) %>%
  group_by(individual, Diet) %>%
  mutate(time_points = case_when(
    str_detect(Sample_name, 'PF') ~ paste0(str_replace_all(Diet, "\\s+", "_"), "_", row_number()),
    TRUE ~ as.character(day)  # fallback for non-PF samples
  )) %>%
  ungroup() %>% 
  mutate(Diet_general = str_remove(Diet, " \\(1 snack a day\\)| \\(2 snacks a day\\)| \\(3 snacks a day\\)")) %>% 
  mutate(Diet_general = str_replace(Diet_general, " \\(pre-intervention\\)", "")) %>%
  mutate(Diet_general = str_replace(Diet_general, " \\(post-intervention\\)", "")) %>% 
  filter(Diet != "Free diet" ) %>% 
  mutate(time= case_when(time_points == 'HiSF-LoFV_(pre-intervention)_1' ~ 'TP-1', 
                         time_points == 'HiSF-LoFV_(pre-intervention)_2' ~ 'TP-2',
                         time_points == 'HiSF-LoFV_(pre-intervention)_3' ~ 'TP-3',
                         time_points == 'Pea_fiber_(1_snack_a_day)_1' ~ 'TP-4',
                         time_points == 'Pea_fiber_(1_snack_a_day)_2' ~ 'TP-5',
                         time_points == 'Pea_fiber_(2_snacks_a_day)_1' ~ 'TP-6',
                         time_points == 'Pea_fiber_(2_snacks_a_day)_2' ~ 'TP-7',
                         time_points == 'Pea_fiber_(3_snacks_a_day)_1' ~ 'TP-8',
                         time_points == 'Pea_fiber_(3_snacks_a_day)_2' ~ 'TP-9',
                         time_points == 'Pea_fiber_(3_snacks_a_day)_3' ~ 'TP-10',
                         time_points == 'Pea_fiber_(3_snacks_a_day)_4' ~ 'TP-11',
                         time_points == 'Pea_fiber_(3_snacks_a_day)_5' ~ 'TP-12',
                         time_points == 'HiSF-LoFV_(post-intervention)_1' ~ 'TP-13',
                         time_points == 'HiSF-LoFV_(post-intervention)_2' ~ 'TP-14',
                         time_points == 'HiSF-LoFV_(post-intervention)_3' ~ 'TP-15',
  )) %>% 
  mutate(time = factor(time, levels=c('TP-1', 'TP-2', 'TP-3', 'TP-4','TP-5' ,'TP-6' ,'TP-7', 'TP-8' ,'TP-9', 'TP-10', 'TP-11' ,'TP-12', 'TP-13', 'TP-14', 'TP-15')))

# Get ordered time points
tp_levels <- levels(Delannoy_2021_scores_PF$time)

# Define positions for rects
stage_rects <- data.frame(
  xmin = seq_along(tp_levels) - 0.5,
  xmax = seq_along(tp_levels) + 0.5,
  stage = tp_levels
)

# Assign colors
stage_rects <- stage_rects %>% 
  mutate(fill= case_when( 
    stage == 'TP-1' | stage =='TP-2' | stage =='TP-3' | stage =='TP-13' | stage =='TP-14' | stage =='TP-15' ~ "grey90",
    stage == 'TP-4' | stage =='TP-5'  ~  "#d9f0d3", # one snack a day
    stage == 'TP-6' | stage =='TP-7'   ~ "#a6dba0", # 2 snacks a day
    stage == 'TP-8' | stage =='TP-9' |  stage =='TP-10' | stage =='TP-11' | stage =='TP-12' ~ "#5aae61", # 3 snacks a day
    TRUE ~ NA_character_
  ))


# Generate separate linear mixed effect models for both fiber amount and diet type as fixed effects 
model_fiber_snack <- lmer(medianPredictionProb ~ fiber_amount + (1 | individual) + (1| time_points) , data = Delannoy_2021_scores_PF)
coefs <- summary(model_fiber_snack)$coefficients
pval <- coefs["fiber_amount", "Pr(>|t|)"]

model_diet <- lmer(medianPredictionProb ~ Diet_general + (1 | individual) + (1| time_points) , data = Delannoy_2021_scores_PF)
coefs_diet <- summary(model_diet)$coefficients
pval_diet <- coefs_diet["Diet_generalPea fiber", "Pr(>|t|)"]

# Format p-values for display
label_fiber <- paste0("p-val (snack amount as a fixed effect) = ", signif(pval, 2))
label_diet <- paste0("p-val (diet type as a fixed effect)= ", signif(pval_diet, 2))

Delannoy_2021_scores_PF <- Delannoy_2021_scores_PF %>%
  mutate(
    Diet_block = case_when(Diet_general=='Pea fiber' ~ 2, 
                           Diet == 'HiSF-LoFV (pre-intervention)'  ~ 1,
                           Diet == 'HiSF-LoFV (post-intervention)'  ~ 3,)
  )

bottom_rects <- Delannoy_2021_scores_PF %>%
  group_by(Diet_block, Diet_general) %>%
  summarise(
    xmin = min(as.numeric(time)-0.5),
    xmax = max(as.numeric(time)+0.5),
    .groups = "drop"
  ) %>%
  mutate(
    fill = rep(c("#8C8C8C", "#ADC2A9", "#8C8C8C"), length.out = n()),
    label = Diet_general,
    ymin = -0.03,
    ymax = -0.005 
  )

p_Delannoy_PF <- ggplot(Delannoy_2021_scores_PF, aes(x = time, y = medianPredictionProb)) +
  geom_rect(data = stage_rects, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 0.3, show.legend = FALSE) +
  geom_point(
    shape = 21,             
    size = 4,               
    fill = "grey",         
    color = "black",       
    alpha = 0.7             
  ) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "blue", size = 1) +
  geom_rect(data = bottom_rects, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            show.legend = FALSE, color = "black", alpha = 0.8) +
  geom_text(data = bottom_rects, inherit.aes = FALSE,
            aes(x = (xmin + xmax)/2, y = -0.035, label = label),
            size = 3.5, vjust = 1) +
  theme_paper +
  scale_fill_identity() +
  scale_x_discrete(limits = tp_levels)+
  labs(x = 'Time points', y = 'CRC microbiome signature score') +
  theme(legend.position = 'none') +
  annotate("text", x = 10.5, y = Inf, vjust = 2, hjust = 0, label = label_fiber, parse = FALSE, size = 4) +
  annotate("text", x = 10.5, y = Inf, vjust = 4, hjust = 0, label = label_diet, parse = FALSE, size = 4)

fiber_dots <- Delannoy_2021_scores_PF %>% 
  select(time, fiber_amount) %>%
  unique() %>%
  mutate(
    time = factor(time, levels = tp_levels),
    fiber_amount = as.integer(fiber_amount)
  ) %>%
  uncount(weights = fiber_amount) %>%
  group_by(time) %>%
  mutate(y = row_number()) %>%
  ungroup()

p0_fiber_dots <- ggplot(fiber_dots, aes(x = time, y = y)) +
  geom_point(shape = 21, size = 2.5, fill = "#4A7C59", color = "black") +
  scale_x_discrete(limits = tp_levels) +
  #scale_y_continuous(expand = c(0.1, 0.1)) +
  labs(y = "Fiber snack") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(), 
    axis.title.x = element_blank(), 
    panel.grid = element_blank(),
    plot.margin = margin(t = -10, b = 5)
  )

p0_combined <- p_Delannoy_PF / p0_fiber_dots + plot_layout(heights = c(5, 1.3))

ggsave(plot=p0_combined, here('figures','figure5','Figure5f.pdf'),width =9, height = 9, dpi = 300)

######################
# Figure 5g: CRC microbiome signature score changes with longitudinal fiber enriched  diet (Delannoy-Bruno_2021b)

Delannoy_2021_scores_HS <- Delannoy_scores %>% 
  filter(Group=='HS') %>%
  mutate(Diet_day= case_when(day==1~ paste0('Free_diet day:',day), 
                             day==9 | day==10  |day==11 ~ paste0('HiSF-LoFV_(pre-intervention day:', day),
                             day==23 | day==24 |day==25 ~ paste0('2_snack_a_day day:', day),
                             day==33 | day==34 |day==35 ~ paste0('HiSF-LoFV_(pre-intervention day:', day),
                             day==47 | day==48 |day==49 ~ paste0('4_snack_a_day day:', day)
  )) %>% drop_na(Diet_day) %>%
  filter(Diet_day !='Free_diet day:1') %>% 
  mutate(Diet_day = factor(Diet_day, levels=c("HiSF-LoFV_(pre-intervention day:9" , "HiSF-LoFV_(pre-intervention day:10"    ,
                                              "HiSF-LoFV_(pre-intervention day:11",  "2_snack_a_day day:23"  ,"2_snack_a_day day:24","2_snack_a_day day:25"   ,
                                              "HiSF-LoFV_(pre-intervention day:33","HiSF-LoFV_(pre-intervention day:34" , "HiSF-LoFV_(pre-intervention day:35",
                                              "4_snack_a_day day:47",  "4_snack_a_day day:48" ,"4_snack_a_day day:49" )))  %>%
  mutate(fiber_amount = case_when(str_detect(Diet_day,'HiSF-LoFV') ~ 0,
                                  str_detect(Diet_day, '2_snack') ~ 2,
                                  str_detect(Diet_day, '4_snack') ~ 4)) %>% 
  mutate(time= case_when(Diet_day == 'HiSF-LoFV_(pre-intervention day:9' ~ 'TP-1', 
                         Diet_day == 'HiSF-LoFV_(pre-intervention day:10' ~ 'TP-2',
                         Diet_day == 'HiSF-LoFV_(pre-intervention day:11' ~ 'TP-3',
                         Diet_day == '2_snack_a_day day:23' ~ 'TP-4',
                         Diet_day == '2_snack_a_day day:24' ~ 'TP-5',
                         Diet_day == '2_snack_a_day day:25' ~ 'TP-6',
                         Diet_day == 'HiSF-LoFV_(pre-intervention day:33' ~ 'TP-7',
                         Diet_day == 'HiSF-LoFV_(pre-intervention day:34' ~ 'TP-8',
                         Diet_day == 'HiSF-LoFV_(pre-intervention day:35' ~ 'TP-9',
                         Diet_day == '4_snack_a_day day:47' ~ 'TP-10',
                         Diet_day == '4_snack_a_day day:48' ~ 'TP-11',
                         Diet_day == '4_snack_a_day day:49' ~ 'TP-12'
  )) %>% 
  mutate(time = factor(time, levels=c('TP-1', 'TP-2', 'TP-3', 'TP-4','TP-5' ,'TP-6' ,'TP-7', 'TP-8' ,'TP-9', 'TP-10', 'TP-11' ,'TP-12')))


# Generate seperate linear mixed effect models for both fiber amount and diet type as fixed effects
model_fiber_snack_HS <- lmer(medianPredictionProb ~ fiber_amount + (1 | individual) + (1| day) , data = Delannoy_2021_scores_HS)
coefs_hs <- summary(model_fiber_snack_HS)$coefficients
pval_hs <- coefs_hs["fiber_amount", "Pr(>|t|)"]

Delannoy_2021_scores_HS <-Delannoy_2021_scores_HS %>% mutate(Diet_general =ifelse(str_detect(Diet_day, 'HiSF-LoFV'), 'HiSF-LoFV','Fiber'))

model_diet_HS <- lmer(medianPredictionProb ~ Diet_general + (1 | individual) + (1| day) , data = Delannoy_2021_scores_HS)
coefs_diet_HS <- summary(model_diet_HS)$coefficients
pval_diet_hs <- coefs_diet_HS["Diet_generalHiSF-LoFV", "Pr(>|t|)"]

# Format p-values for display
label_fiber2 <- paste0("p-val (snack amount as a fixed effect) = ", signif(pval_hs, 2))
label_diet2 <- paste0("p-val (diet type as a fixed effect)= ", signif(pval_diet_hs, 2))

Delannoy_2021_scores_HS <- Delannoy_2021_scores_HS %>%
  mutate(
    Diet_block = case_when(Diet_day=='HiSF-LoFV_(pre-intervention day:9' | Diet_day=='HiSF-LoFV_(pre-intervention day:10' |Diet_day=='HiSF-LoFV_(pre-intervention day:11'  ~ 1, 
                           fiber_amount == 2 ~ 2,
                           Diet_day=='HiSF-LoFV_(pre-intervention day:33' | Diet_day=='HiSF-LoFV_(pre-intervention day:34' |Diet_day=='HiSF-LoFV_(pre-intervention day:35'  ~ 3, 
                           fiber_amount == 4 ~ 4)
  )

# Get ordered time points
tp_levels2 <- levels(Delannoy_2021_scores_HS$time)

stage_rects2 <- data.frame(
  xmin = seq_along(tp_levels2) - 0.5,
  xmax = seq_along(tp_levels2) + 0.5,
  stage = tp_levels2
)

# Assign colors
stage_rects2 <- stage_rects2 %>% 
  mutate(fill= case_when( 
    stage == 'TP-1' | stage =='TP-2' | stage =='TP-3' ~ "grey90",
    stage == 'TP-4' | stage =='TP-5' | stage =='TP-6' ~ "#d9f0d3",
    stage == 'TP-7' | stage =='TP-8' | stage =='TP-9' ~ "grey90",
    stage == 'TP-10' | stage =='TP-11' | stage =='TP-12' ~ "#a6dba0",
    TRUE ~ NA_character_
  ))

bottom_rects2 <- Delannoy_2021_scores_HS %>%
  group_by(Diet_block, Diet_general) %>%
  summarise(
    xmin = min(as.numeric(time)-0.5),
    xmax = max(as.numeric(time)+0.5),
    .groups = "drop"
  ) %>%
  mutate(
    fill = rep(c("#8C8C8C", "#ADC2A9", "#8C8C8C","#ADC2A9"), length.out = n()),
    label = Diet_general,
    ymin = -0.03, 
    ymax = -0.005  
  )

p_Delannoy_HS <- ggplot(Delannoy_2021_scores_HS, aes(x = time, y = medianPredictionProb)) +
  geom_rect(data = stage_rects2, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 0.3, show.legend = FALSE) +
  
  geom_point(
    shape = 21,             
    size = 4,               
    fill = "grey",         
    color = "black",       
    alpha = 0.7             
  ) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "blue", size = 1) +
  geom_rect(data = bottom_rects2, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            show.legend = FALSE, color = "black", alpha = 0.8) +
  geom_text(data = bottom_rects2, inherit.aes = FALSE,
            aes(x = (xmin + xmax)/2, y = -0.035, label = label),
            size = 3.5, vjust = 1) +
  theme_paper +
  scale_x_discrete(limits = tp_levels2)+
  scale_fill_identity() +
  labs(x = 'Time points', y = 'CRC microbiome signature score') +
  theme(legend.position = 'none')+
  annotate("text", x = 3.5, y = 0.1, vjust = 2, hjust = 0, label = label_fiber2, parse = FALSE, size = 4) +
  annotate("text", x = 3.5, y = 0.15, vjust = 4, hjust = 0, label = label_diet2, parse = FALSE, size = 4)

Delannoy_2021_scores_HS <- Delannoy_2021_scores_HS %>%
  mutate(time = factor(time, levels = tp_levels2))

fiber_dots_df <- Delannoy_2021_scores_HS %>%
  select(time, fiber_amount) %>%
  unique() %>%
  mutate(
    time = factor(time, levels = tp_levels2),
    fiber_amount = as.integer(fiber_amount)
  ) %>%
  uncount(weights = fiber_amount) %>%
  group_by(time) %>%
  mutate(y = row_number()) %>%
  ungroup()

p_fiber_dots <- ggplot(fiber_dots_df, aes(x = time, y = y)) +
  geom_point(shape = 21, size = 2.5, fill = "#4A7C59", color = "black") +
  scale_x_discrete(limits = tp_levels2) +
  labs(y = "Fiber snack") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(), 
    axis.title.x = element_blank(), 
    panel.grid = element_blank(),
    plot.margin = margin(t = -10, b = 5)
  )

p0_combined2 <- p_Delannoy_HS / p_fiber_dots + plot_layout(heights = c(5, 1.3))

ggsave(plot=p0_combined2, here('figures','figure5','Figure5g.pdf'),width =9, height = 9, dpi = 300)

##################
# Extended Data Figure 5e: CRC microbiome signature score changes with longitudinal fiber enriched  diet (HealeyG_2018)

# Load metadata of the study
Healey_meta <- read.table(here("data", "fiber_data", "HealeyG_2018_metadata.txt"), sep='\t', header = T)  %>% 
  mutate(Run=str_split_fixed(sampleid, '_',2)[,2])

Healey_sra <- read.table(here("data", "fiber_data", "SraRunTable-Healey.csv"), sep=',', header = T)  %>% select(Run,  Host_Diet)

Healey_meta<- Healey_meta %>% left_join(Healey_sra, by='Run')

# Load the CRC microbiome signature scores of the samples in this study, and join with the metadata to get the time point annotation for each sample. 
Healey_score <- Fiber_scores %>% filter(label == 'HealeyG_2018') %>%
  right_join(
    Healey_meta %>%
      select(subject_id, Run, treatment, gender, time_days, number, timepoint, timepoint_id, timepoint_numeric,Host_Diet),
    by = c("sampleID" = "Run")
  ) %>%
  mutate(
    subject_id  = as.factor(subject_id),
    timepoint_id = factor(timepoint_id, levels = c("before_1","after_2","before_3","after_4"))
  ) %>%
  group_by(subject_id) %>%
  mutate(
    Intervention_start = case_when(
      any(timepoint_id == "before_1" & treatment == "control") ~ "Control",
      any(timepoint_id == "before_1" & treatment == "fiber")   ~ "Fiber",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(
    timepoint_annotated = case_when(
      Intervention_start == "Control" & timepoint_id == "before_1" ~ "before placebo",
      Intervention_start == "Control" & timepoint_id == "after_2"  ~ "after placebo",
      Intervention_start == "Control" & timepoint_id == "before_3" ~ "before fiber",
      Intervention_start == "Control" & timepoint_id == "after_4"  ~ "after fiber",
      Intervention_start == "Fiber"   & timepoint_id == "before_1" ~ "before fiber",
      Intervention_start == "Fiber"   & timepoint_id == "after_2"  ~ "after fiber",
      Intervention_start == "Fiber"   & timepoint_id == "before_3" ~ "before placebo",
      Intervention_start == "Fiber"   & timepoint_id == "after_4"  ~ "after placebo",
      TRUE ~ NA_character_
    ),
    tp_key = case_when(
      Intervention_start == "Control" & timepoint_annotated == "before placebo" ~ "1_before placebo",
      Intervention_start == "Control" & timepoint_annotated == "after placebo"  ~ "2_after placebo",
      Intervention_start == "Control" & timepoint_annotated == "before fiber"   ~ "3_before fiber",
      Intervention_start == "Control" & timepoint_annotated == "after fiber"    ~ "4_after fiber",
      Intervention_start == "Fiber"   & timepoint_annotated == "before fiber"   ~ "1_before fiber",
      Intervention_start == "Fiber"   & timepoint_annotated == "after fiber"    ~ "2_after fiber",
      Intervention_start == "Fiber"   & timepoint_annotated == "before placebo" ~ "3_before placebo",
      Intervention_start == "Fiber"   & timepoint_annotated == "after placebo"  ~ "4_after placebo",
      TRUE ~ NA_character_
    ),
    Intervention_start = factor(Intervention_start, levels = c("Control","Fiber"))
  ) %>%
  mutate(
    tp_key = factor(
      tp_key,
      levels = c(
        # control order
        "1_before placebo","2_after placebo","3_before fiber","4_after fiber",
        # fiber order
        "1_before fiber","2_after fiber","3_before placebo","4_after placebo"
      )
    )
  )


Healey_score <- Healey_score %>%
  mutate(tp_num = as.integer(sub("^([0-9]+).*", "\\1", as.character(tp_key))))

# Compare before and after fiber intake 
only_before_after <- Healey_score %>% 
  filter(timepoint_annotated %in% c("before fiber", "after fiber")) %>% 
  mutate(timepoint_annotated= factor(timepoint_annotated, levels=c('before fiber','after fiber'))) %>%
  select(value=medianPredictionProb, timepoint_annotated, Host_Diet, subject_id) 

# Perform Wilcoxon signed-rank test 
stat.test_fiber <- only_before_after %>%
  wilcox_test(value ~ timepoint_annotated) %>%
  adjust_pvalue(method = "BH") %>%
  add_xy_position(x = "timepoint_annotated") %>%
  mutate(
    p.adj.label = paste0("p = ", signif(p.adj, 3)),
    y.position = y.position ,
    x = (as.numeric(factor(group1, levels = levels(factor(only_before_after$timepoint_annotated)))) +
           as.numeric(factor(group2, levels = levels(factor(only_before_after$timepoint_annotated))))) / 2
  )

# Get group sizes
group.n_fiber <- only_before_after %>%
  dplyr::count(timepoint_annotated) %>%
  select(n) %>% 
  unique() 

# Plot
Healey_fiber_figure <- ggplot() +
  geom_boxplot(
    data = only_before_after,
    aes(x = timepoint_annotated, y = value),
    width = 0.5, outlier.shape = NA
  ) +
  geom_jitter(
    data = only_before_after,
    aes(x = timepoint_annotated, y = value),
    size = 3.5, shape = 21, fill = "grey80", color = "black", alpha = 0.8, width = 0.2
  ) +
  geom_text(
    data = stat.test_fiber,
    aes(x = x, y = y.position + 0.04, label = p.adj.label),
    vjust = 0
  ) +  
  geom_text(
      data = group.n_fiber,
      aes(x = 1.5, y = max(only_before_after$value) * 1.05,
          label = paste0("(N = ", n, ')')),
      inherit.aes = FALSE,
      vjust = 0
    ) +
  labs(y = "CRC microbiome signature score", x = "Fiber intake") +
  theme_paper +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  )

# Save the figure
ggsave(plot=Healey_fiber_figure, here('figures','figure5','Figure5e.pdf'),width = 2, height = 4, dpi = 300)

##################
# Extended Data Figure 5c: CRC microbiome signature score changes with longitudinal fiber enriched  diet (HealeyG_2018) - only placebo arm

# Extract only the before and after placebo samples 
only_before_after_placebo <- Healey_score %>% 
  filter(timepoint_annotated %in% c("before placebo", "after placebo")) %>% 
  mutate(timepoint_annotated= factor(timepoint_annotated, levels=c('before placebo','after placebo'))) %>%
  select(value=medianPredictionProb, timepoint_annotated, Host_Diet, subject_id) 

# Perform Wilcoxon signed-rank test
stat.test_placebo <- only_before_after_placebo %>%
  wilcox_test(value ~ timepoint_annotated) %>%
  adjust_pvalue(method = "BH") %>%
  add_xy_position(x = "timepoint_annotated") %>%
  mutate(
    p.adj.label = paste0("p = ", signif(p.adj, 3)),
    y.position = y.position ,
    x = (as.numeric(factor(group1, levels = levels(factor(only_before_after_placebo$timepoint_annotated)))) +
           as.numeric(factor(group2, levels = levels(factor(only_before_after_placebo$timepoint_annotated))))) / 2
  )

# Get group sizes
group.n_placebo<- only_before_after_placebo %>%
  dplyr::count(timepoint_annotated) %>%
  select(n) %>% 
  unique() 

# Plot
Healey_placebo_figure <- ggplot() +
  geom_boxplot(
    data = only_before_after_placebo,
    aes(x = timepoint_annotated, y = value),
    width = 0.5, outlier.shape = NA
  ) +
  geom_jitter(
    data = only_before_after_placebo,
    aes(x = timepoint_annotated, y = value),
    size = 3.5, shape = 21, fill = "grey80", color = "black", alpha = 0.8, width = 0.2
  ) +
  geom_text(
    data = stat.test_placebo,
    aes(x = x, y = y.position + 0.04, label = p.adj.label),
    vjust = 0
  ) +  
  geom_text(
    data = group.n_placebo,
    aes(x = 1.5, y = max(only_before_after_placebo$value) * 1.4,
        label = paste0("(N = ", n, ')')),
    inherit.aes = FALSE,
    vjust = 0
  ) +
  labs(y = "CRC microbiome signature score", x = "Fiber intake") +
  theme_paper +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  )

# Save the figure
ggsave(plot=Healey_placebo_figure, here('figures','extended.data.figure5','Extended.data.Figure5c.pdf'),width = 3.5, height = 6, dpi = 300)

##################
# Data figure 5b: CRC microbiome signature score changes with longitudinal fiber enriched  diet (HealeyG_2018) - only baseline comparison between high and low fiber intake

Healey_score_BL <- Healey_score %>%
  filter(timepoint_numeric == 1) %>%
  select(medianPredictionProb, Host_Diet) %>%
  pivot_longer(-Host_Diet, names_to = "measure", values_to = "value") %>% 
  mutate(Host_Diet = ifelse(Host_Diet=='High dietary fibre intake', 'HFD','LFD')) %>% 
  mutate(Host_Diet= factor(Host_Diet, levels=c('LFD', 'HFD')))

# Perform Wilcoxon rank-sum test
stat.test <- Healey_score_BL %>% 
  wilcox_test(value ~ Host_Diet) %>%
  adjust_pvalue(method = "BH") %>%
  add_xy_position(x = "Host_Diet") %>%
  dplyr::mutate(p.adj.label = paste0("p = ", signif(p.adj, 3))) %>% 
  mutate(y.position= y.position-0.02) %>%
  mutate(
    x = (as.numeric(factor(group1, levels = levels(Healey_score_BL$Host_Diet))) +
           as.numeric(factor(group2, levels = levels(Healey_score_BL$Host_Diet)))) / 2
  )

# Get group sizes
group.n <- Healey_score_BL %>%
  dplyr::count(Host_Diet)

# Plot
Healey_BL_figure <-ggplot() +
  geom_boxplot(
    data = Healey_score_BL,
    aes(x = Host_Diet, y = value),
    width = 0.5, outlier.shape = NA
  ) +
  geom_jitter(
    data = Healey_score_BL,
    aes(x = Host_Diet, y = value),
    size = 3.5, shape = 21, fill = "grey80", color = "black", alpha = 0.8, width = 0.2
  ) +
  geom_text(
    data = group.n,
    aes(x = Host_Diet, y = max(Healey_score_BL$value) * 1.05,
        label = paste0("(N = ", n,')')),
    inherit.aes = FALSE,
    vjust = 0
  ) +
  geom_text(
    data = stat.test,
    aes(x = x, y = y.position, label = p.adj.label),
    vjust = 0
  ) +  
  labs(y = "CRC microbiome signature score", x = "Host Diet") +
  theme_paper +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  )
 
# Save the figure
ggsave(plot=Healey_BL_figure, here('figures','extended.data.figure5','Extended.data.Figure5b.pdf'),width = 3, height = 6, dpi = 300)

##################
# Extended Data Figure 5e: CRC microbiome signature score changes with longitudinal fiber enriched  diet (Delannoy-Bruno_2022)

Delannoy_2022_scores<- Fiber_scores %>% filter(label == 'Delannoy-Bruno_2022')

# Load the metadata of the study
meta<-read.table(here('data','fiber_data','Delannoy-Bruno_2022.csv'), sep=',',header = T) %>%
  filter(Assay.Type !='AMPLICON') %>% # only keep the WGS samples
  select(Run, Sample_name, Diet, host_scientific_name) 

Delannoy_2022_scores<- Delannoy_2022_scores %>% right_join(meta, by=c('sampleID'='Run')) 

Delannoy_2022_scores_wgs <- Delannoy_2022_scores %>%
  mutate(individual= str_split_fixed(Sample_name, '_',4)[,1]) %>% 
  mutate(time_points= str_split_fixed(Sample_name, '_',4)[,2]) %>% 
  mutate(Group= str_split_fixed(Sample_name, '_',4)[,3])  %>% 
  mutate(Group= str_remove(Group, 'Study1')) %>% 
  mutate(Group= str_remove(Group, 'Study2')) %>%  
  mutate(time_points= str_remove(time_points, "Week")) %>% 
  mutate(time_points = factor(time_points, levels = sort(as.numeric(unique(time_points))))) 

Delannoy_2022_scores_wgs_grouped_PF <- Delannoy_2022_scores_wgs %>% 
  filter(Group=='PeaFiber') %>%
  arrange(individual, time_points) %>%
  group_by(individual, Diet) %>%
  mutate(diet_time_points = case_when(
    str_detect(Sample_name, 'PeaFiber') ~ paste0(str_replace_all(Diet, "\\s+", "_"), "_", row_number()),
    TRUE ~ as.character(time_points)  
  )) %>%
  ungroup() %>% 
  mutate(diet_time_points= factor(diet_time_points, levels=c("Free_diet_1" ,  "Free_diet_2"  ,
                                                             "Pea_fiber_(1_snack_a_day)_1" , "Pea_fiber_(2_snacks_a_day)_1",
                                                             "Pea_fiber_(3_snacks_a_day)_1", "Pea_fiber_(3_snacks_a_day)_2" ,'Free_diet_3'))) %>% 
  mutate(Diet_general= ifelse(str_detect(diet_time_points, 'Free'), 'Free_diet', 'Pea fiber' )) %>% 
  mutate(fiber_amount = case_when(str_detect(diet_time_points, 'Free') ~ 0,
                                  str_detect(diet_time_points, '1') ~ 1,
                                  str_detect(diet_time_points, '2') ~ 2,
                                  str_detect(diet_time_points, '3') ~ 3
  )) %>% 
  mutate(Twin_id= str_split_fixed(individual,"-",2)[,1])

# Generate separate linear mixed effect models for both fiber amount and diet type as fixed effects 
model1_snack <- lmer(medianPredictionProb ~ fiber_amount + (1 | individual) + (1| time_points) , data = Delannoy_2022_scores_wgs_grouped_PF)
coefs <- summary(model1_snack)$coefficients
pval <- coefs["fiber_amount", "Pr(>|t|)"]

model1_diet <- lmer(medianPredictionProb ~ Diet_general + (1 | individual) + (1| time_points) , data = Delannoy_2022_scores_wgs_grouped_PF)
coefs_diet <- summary(model1_diet)$coefficients
pval1_diet <- coefs_diet["Diet_generalPea fiber", "Pr(>|t|)"]

# generate linear mixed effect model also for the twin pairing as a random effect to see if the significance of the fixed effect of fiber amount is affected by the twin pairing.
model1_snack_twin <- lmer(medianPredictionProb ~ fiber_amount + (1 | individual) + (1| time_points) + (1 | Twin_id) , data = Delannoy_2022_scores_wgs_grouped_PF)
coefs2 <- summary(model1_snack_twin)$coefficients
pval2 <- coefs2["fiber_amount", "Pr(>|t|)"]

# Get ordered time points
tp_levels <- levels(Delannoy_2022_scores_wgs_grouped_PF$diet_time_points)

# Define positions for rects (assuming each level is one unit wide)
stage_rects <- data.frame(
  xmin = seq_along(tp_levels) - 0.5,
  xmax = seq_along(tp_levels) + 0.5,
  stage = tp_levels
)

# Assign colors
stage_rects$fill_color <- case_when(
  grepl("Free_diet*", stage_rects$stage) ~ "grey90",
  grepl("Pea_fiber_\\(1_snack_a_day\\)", stage_rects$stage) ~ "#d9f0d3",
  grepl("Pea_fiber_\\(2_snacks_a_day\\)", stage_rects$stage) ~ "#a6dba0",
  grepl("Pea_fiber_\\(3_snacks_a_day\\)", stage_rects$stage) ~ "#5aae61",
  TRUE ~ NA_character_
)

# Format p-values for display
label_fiber <- paste0("p-val (snack amount as a fixed effect) = ", signif(pval, 2))
label_diet <- paste0("p-val (diet type as a fixed effect)= ", signif(pval1_diet, 2))

df <- Delannoy_2022_scores_wgs_grouped_PF %>%
  arrange(individual, as.integer(diet_time_points))

pos <- position_jitter(width = 0.2, height = 0, seed = 42)

# Plot 
p_Delannoy_PF<- ggplot(Delannoy_2022_scores_wgs_grouped_PF, aes(x = diet_time_points, y = medianPredictionProb)) +

  geom_rect(data = stage_rects, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_color),
            alpha = 0.3, show.legend = FALSE) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "blue", size = 1) +
  theme_paper +
  scale_fill_identity() +
  labs(x = 'Time points', y = 'Colorectal cancer microbiome signature score') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90))+
  annotate("text", x = 1.5, y = Inf, vjust = 2, hjust = 0, label = label_fiber, parse = FALSE, size = 4) +
  annotate("text", x = 1.5, y = Inf, vjust = 4, hjust = 0, label = label_diet, parse = FALSE, size = 4)

# Orange fiber  intervention 

Delannoy_2022_scores_wgs_grouped_OF <- Delannoy_2022_scores_wgs %>% 
  filter(Group=='OrangeFiber') %>% # only keep the orange fiber intervention samples
  arrange(individual, time_points) %>%
  group_by(individual, Diet) %>%
  mutate(diet_time_points = case_when(
    str_detect(Sample_name, 'OrangeFiber') ~ paste0(str_replace_all(Diet, "\\s+", "_"), "_", row_number()),
    TRUE ~ as.character(time_points)  
  )) %>%
  ungroup() %>% 
  mutate(diet_time_points= factor(diet_time_points, levels=c("Free_diet_1", "Free_diet_2"  , "Orange_fiber_(3_snacks_a_day)_1" , "Orange_fiber_(3_snacks_a_day)_2" ))) %>%
  mutate(Diet_general= ifelse(str_detect(diet_time_points, 'Free'), 'Free_diet', 'Orange fiber' )) %>% 
  mutate(fiber_amount = case_when(str_detect(diet_time_points, 'Free') ~ 0,
                                  str_detect(diet_time_points, '3') ~ 3
  )) %>% 
  mutate(Twin_id= str_split_fixed(individual,"-",2)[,1])

# Generate separate linear mixed effect models for both fiber amount and diet type as fixed effects 
model2_snack <- lmer(medianPredictionProb ~ fiber_amount + (1 | individual) + (1| time_points) , data = Delannoy_2022_scores_wgs_grouped_OF)
pval <- summary(model2_snack)$coefficients["fiber_amount", "Pr(>|t|)"]

model2_diet <- lmer(medianPredictionProb ~ Diet_general + (1 | individual) + (1| time_points) , data = Delannoy_2022_scores_wgs_grouped_OF)
pval2_diet <- summary(model2_diet)$coefficients["Diet_generalOrange fiber", "Pr(>|t|)"]

label_fiber2 <- paste0("p-val (snack amount as a fixed effect) = ", signif(pval, 2))
label_diet2 <- paste0("p-val (diet type as a fixed effect)= ", signif(pval2_diet, 2))

# Get ordered time points
tp_levels <- levels(Delannoy_2022_scores_wgs_grouped_OF$diet_time_points)

stage_rects <- data.frame(
  xmin = seq_along(tp_levels) - 0.5,
  xmax = seq_along(tp_levels) + 0.5,
  stage = tp_levels
)

# Assign colors
stage_rects$fill_color <- case_when(
  grepl("Free_diet*", stage_rects$stage) ~ "grey90",
  grepl("Orange_fiber_\\(1_snack_a_day\\)", stage_rects$stage) ~ "#d9f0d3",
  grepl("Orange_fiber_\\(3_snacks_a_day\\)", stage_rects$stage) ~ "#5aae61",
  TRUE ~ NA_character_
)

df_of <- Delannoy_2022_scores_wgs_grouped_OF %>%
  arrange(individual, as.integer(diet_time_points))

pos_of <- position_jitter(width = 0.1, height = 0, seed = 42)

p_Delannoy_OF <- ggplot(df_of, aes(x = diet_time_points, y = medianPredictionProb)) +
  geom_rect(data = stage_rects, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_color),
            alpha = 0.3, show.legend = FALSE) +
  geom_path(aes(group = individual, color = Twin_id),
            alpha = 0.4, position = pos_of) +
  geom_point(aes(color = Twin_id), size = 3, alpha = 0.6, position = pos_of) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "blue", size = 1) +
  theme_paper +
  scale_fill_identity() +
  labs(x = "Time points", y = "Colorectal cancer microbiome signature score") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90)) +
  annotate("text", x = 1.5,
           y = max(df_of$medianPredictionProb, na.rm = TRUE) * 1.05,
           hjust = 0, label = label_fiber2, size = 4) +
  annotate("text", x = 1.5,
           y = max(df_of$medianPredictionProb, na.rm = TRUE) * 1.10,
           hjust = 0, label = label_diet2, size = 4)

# combine the two figures for PF and OF and save
p_combined <- p_Delannoy_PF / p_Delannoy_OF + plot_layout(heights = c(1, 1))
ggsave(plot=p_combined, here('figures','extended.data.figure5','Extended.Data.Figure5ef.pdf'),width = 15, height = 9, dpi = 300)