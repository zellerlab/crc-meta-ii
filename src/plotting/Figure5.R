######################
# Figure 5: The CRC microbiome signature score is inversely associated with dietary fibre intake and can be modulated by dietary intervention
######################

source(here('src','utils.R'))
source(here('src','cancerness.utils.R'))

######################
# Figure 5a:Correlation of fiber intake (YachidaS_2019) and CRC microbiome signature score

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
# Figure 5b:Correlation of fiber intake (SowahSA_2022) and CRC microbiome signature score

Sowah_relab <-readRDS(here('data','fiber_data','SowahSA_2022_res_IDTaxa.rds')) %>% 
  as.data.frame() %>%
  rownames_to_column('taxon') %>%
  pivot_longer(-taxon) %>%
  mutate(genus = str_split_fixed(taxon, "\\|", n = 7)[, 6],
         genus = str_replace(genus, "g__", "")) %>%
  dplyr::rename(sampleID = name) %>%
  mutate(genus = ifelse(genus == "", "unassigned", genus)) %>% 
  group_by(sampleID,genus) %>%
  summarise(value=sum(value)) %>% 
  pivot_wider(names_from = sampleID,values_from = value) %>%
  column_to_rownames('genus') %>%
  mutate(across(everything(), ~ . / sum(.)))


# Load metadata
Sowah_meta <- read.table(here('data','fiber_data','SraRunTable_Sowah.csv'), sep=',', header = T) %>%
  select(Run,Sample.Name, Sample_name) %>%
  mutate(Time_points=  str_split_fixed(Sample_name,'_',4)[,2],
         Batch=  str_split_fixed(Sample_name,'_',4)[,3],
         Patient_ID = str_split_fixed(Sample_name,'_',4)[,1],
         Identifier = str_split_fixed(Sample_name,'_',4)[,4]) 

# assign CRC scores to the Sowah samples
dummyControls <- data.frame(dummySample = rpois(n = dim(Sowah_relab)[1], lambda = 10))


rownames(dummyControls) <- rownames(Sowah_relab)
testProfile <- cbind(Sowah_relab, dummyControls)

tmp <- lapply(colnames(testProfile), function(x) {
  tmp <- data.frame(dummyColumnName = rep(0, sum(!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile))))
  rownames(tmp) <- norm_params(models.all.rf)$retained.feat[!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile)]
  colnames(tmp)[1] <- x
  return(tmp)
})

testProfile <- rbind(testProfile, as.data.frame(tmp, check.names = F))
dummyMeta <- data.frame(label = c(rep("CRC", dim(testProfile)[2] - 1 ), "Control"))
rownames(dummyMeta) <- colnames(testProfile)

Sowah_score_plot <- get_cancerness(siam = models.all.rf, pro = testProfile, meta = dummyMeta, evaluated.label='SowahSA_2022')

Sowah_scores <- Sowah_score_plot[[2]] %>%  
  filter(type=='evaluated') %>%
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


#load also sample groups used in this study

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


diet_long<-read.table(here('data','fiber_data', 'Diet all long.csv'), sep=',',header = T)

diet_labels<- read.table(here('data','fiber_data','Diet labels all.csv'), sep=',',header = T) 

label_map <- setNames(diet_labels$Description, diet_labels$Label)

# Replace column names if they are in label_map
colnames(diet_long) <- ifelse(
  colnames(diet_long) %in% names(label_map),
  label_map[colnames(diet_long)],
  colnames(diet_long)  # keep original if no match
)

# First  take the baseline of samples and correlate with fiber intake 
Sowah_all_diet<- Sowah_scores%>%
  mutate(Patient_ID=as.integer(Patient_ID)) %>%
  select(medianPredictionProb, Time_points,Patient_ID, Treatment_group,Batch) %>% 
  left_join(diet_long, by=c('Patient_ID'='HELENA_ID', 'Time_points'='Timepoint'))

# select baseline samples 
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
  theme_paper

ggsave(plot=Sowah_baseline_cor, here('figures','figure5','Figure5b.pdf'), width = 6, height = 6)

######################
# Figure 5c: CRC microbiome signature score changes with fiber enriched Mediterranean diet (BarberC_2021)

mOTUS.3.1.metadata<- read_tsv(here('data','mOTUs3.1.0.genome_metadata_edited.tsv'))

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
  filter(mOTUs!='meta_mOTU_v31_12968') %>% # these specific motus don't have any entry in mOTUs3.1.0.genome_metadata.tsv 
  select(-DOMAIN) %>% 
  group_by(sampleID, genus) %>% 
  dplyr::summarise(count=sum(count)) %>% 
  pivot_wider(id_cols = genus , names_from = sampleID, values_from = count, values_fill = 0)  %>% 
  filter(genus!='unassigned') %>% 
  column_to_rownames('genus')  %>%
  mutate(across(everything(), ~ . / sum(.))) 


# Calculate CRC scores

dummyControls <- data.frame(dummySample = rpois(n = dim(Barber_relab)[1], lambda = 10))

rownames(dummyControls) <- rownames(Barber_relab)
testProfile <- cbind(Barber_relab, dummyControls)

tmp <- lapply(colnames(testProfile), function(x) {
  tmp <- data.frame(dummyColumnName = rep(0, sum(!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile))))
  rownames(tmp) <- norm_params(models.all.rf)$retained.feat[!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile)]
  colnames(tmp)[1] <- x
  return(tmp)
})

testProfile <- rbind(testProfile, as.data.frame(tmp, check.names = F))
dummyMeta <- data.frame(label = c(rep("CRC", dim(testProfile)[2] - 1 ), "Control"))
rownames(dummyMeta) <- colnames(testProfile)


Barber_score_plot <- get_cancerness(siam = models.all.rf, pro = testProfile, meta = dummyMeta, evaluated.label='BarberC_2021')

Barber_scores<- Barber_score_plot[[2]] %>% filter(type=='evaluated')

# Load metadaata
Barber_meta<-read.table(here('data','fiber_data','BarberC2021.tsv'), sep='\t', header = T) %>% 
  select(c(sampleID=run_accession, sample_title)) %>% 
  mutate(time_points= str_split_fixed(sample_title,' ',3)[,3]) %>%
  mutate(individual= str_split_fixed(sample_title,' ',3)[,2]) 


Barber_scores <-Barber_scores %>% left_join(Barber_meta, by='sampleID') %>%
  mutate(individual = as.factor(individual),
         time_points = factor(time_points, levels = c('WD', 'FMD'))) 

Barber_scores_paired <- Barber_scores %>%
  group_by(individual) %>%
  filter(n_distinct(time_points) == 2) %>%
  ungroup()


Barber_wide_scores <- Barber_scores_paired %>%
  select(individual, time_points, medianPredictionProb) %>%
  pivot_wider(names_from = time_points, values_from = medianPredictionProb)

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

Barber_scores_plot <- ggplot() +
  geom_boxplot(data = Barber_scores_paired, aes(x = time_points, y = medianPredictionProb), 
               width = 0.5, outlier.shape = NA) +
  geom_line(data = Barber_scores_paired, aes(x = x_jittered, y = medianPredictionProb, group = individual),
            color = "grey", alpha = 0.5) +
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
  theme_paper +
  ylim(c(0,0.65)) +
  theme(legend.position = 'top', 
        axis.title.x = element_blank())

ggsave(plot=Barber_scores_plot, here('figures','figure5','Figure5c.pdf'),width = 3, height = 6, dpi = 300)


######################
# Figure 5d: Assess fiber intake changes in CRC-depleted genera (Lachnospira, Agathobaculum, CAG-41) 

selected_genera <- c("CAG-41","Agathobaculum", "Lachnospira")

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

ggsave(plot=abundance_plot, here('figures','figure5','Figure5d.pdf'),width = 4, height = 6, dpi = 300)

######################
# Figure 5e: CRC microbiome signature score changes with longitudinal fiber enriched  diet (Delannoy-Bruno_2021a)

Delannoy_relab <- readRDS(here('data','fiber_data','Delannoy-Bruno_2021_res_mOTUs.rds')) %>%
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
  filter(mOTUs!='meta_mOTU_v31_12968') %>% # these specific motus don't have any entry in mOTUs3.1.0.genome_metadata.tsv 
  select(-DOMAIN) %>% 
  group_by(sampleID, genus) %>% 
  dplyr::summarise(count=sum(count)) %>% 
  pivot_wider(id_cols = genus , names_from = sampleID, values_from = count, values_fill = 0)  %>% 
  filter(genus!='unassigned') %>% 
  column_to_rownames('genus')  %>%
  mutate(across(everything(), ~ . / sum(.))) 


# Calculate CRC scores
dummyControls <- data.frame(dummySample = rpois(n = dim(Delannoy_relab)[1], lambda = 10))
rownames(dummyControls) <- rownames(Delannoy_relab)
testProfile <- cbind(Delannoy_relab, dummyControls)

tmp <- lapply(colnames(testProfile), function(x) {
  tmp <- data.frame(dummyColumnName = rep(0, sum(!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile))))
  rownames(tmp) <- norm_params(models.all.rf)$retained.feat[!norm_params(models.all.rf)$retained.feat %in% rownames(testProfile)]
  colnames(tmp)[1] <- x
  return(tmp)
})


testProfile <- rbind(testProfile, as.data.frame(tmp, check.names = F))
dummyMeta <- data.frame(label = c(rep("CRC", dim(testProfile)[2] - 1 ), "Control"))
rownames(dummyMeta) <- colnames(testProfile)

Delannoy_score_plot <- get_cancerness(siam = models.all.rf, pro = testProfile, meta = dummyMeta, evaluated.label='Delannoy-Bruno_2021')

Delannoy_scores<- Delannoy_score_plot[[2]] %>% 
  filter(type=='evaluated')

Delannoy_meta<-read.table(here('data','fiber_data','DelannoyBruno2021.csv'), sep=',',header = T) %>%
  filter(host_scientific_name !='Mus musculus') %>%
  filter(Assay.Type !='AMPLICON') %>%
  select(Run, Sample_name, Diet, host_scientific_name) 


Delannoy_scores<- Delannoy_scores %>% left_join(Delannoy_meta, by=c('sampleID'='Run'))  %>%
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

# Define positions for rects (assuming each level is one unit wide)
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


# generate linear mixed effect model 
model_fiber_snack <- lmer(medianPredictionProb ~ fiber_amount + (1 | individual) + (1| time_points) , data = Delannoy_2021_scores_PF)
summary(model_fiber_snack)

coefs <- summary(model_fiber_snack)$coefficients
pval <- coefs["fiber_amount", "Pr(>|t|)"]
print(pval)


model_diet <- lmer(medianPredictionProb ~ Diet_general + (1 | individual) + (1| time_points) , data = Delannoy_2021_scores_PF)
summary(model_diet)


coefs_diet <- summary(model_diet)$coefficients
pval_diet <- coefs_diet["Diet_generalPea fiber", "Pr(>|t|)"]
print(pval_diet)


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
    ymin = -0.03,  # below data range
    ymax = -0.005  # thin band
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
  
  # Diet phase labels
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

ggsave(plot=p0_combined, here('figures','figure5','Figure5e.pdf'),width =9, height = 9, dpi = 300)


######################
# Figure 5f: CRC microbiome signature score changes with longitudinal fiber enriched  diet (Delannoy-Bruno_2021b)


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


# generate linear mixed effect model 
model_fiber_snack_HS <- lmer(medianPredictionProb ~ fiber_amount + (1 | individual) + (1| day) , data = Delannoy_2021_scores_HS)
summary(model_fiber_snack_HS)


coefs_hs <- summary(model_fiber_snack_HS)$coefficients
pval_hs <- coefs_hs["fiber_amount", "Pr(>|t|)"]
print(pval_hs)


Delannoy_2021_scores_HS <-Delannoy_2021_scores_HS %>% mutate(Diet_general =ifelse(str_detect(Diet_day, 'HiSF-LoFV'), 'HiSF-LoFV','Fiber'))

model_diet_HS <- lmer(medianPredictionProb ~ Diet_general + (1 | individual) + (1| day) , data = Delannoy_2021_scores_HS)
summary(model_diet_HS)


coefs_diet_HS <- summary(model_diet_HS)$coefficients
pval_diet_hs <- coefs_diet_HS["Diet_generalHiSF-LoFV", "Pr(>|t|)"]
print(pval_diet_hs)


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

# Define positions for rects (assuming each level is one unit wide)
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
  
  # Diet phase labels
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


ggsave(plot=p0_combined2, here('figures','figure5','Figure5f.pdf'),width =9, height = 9, dpi = 300)













