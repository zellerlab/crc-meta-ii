######################
# Figure 4
######################

source(here('src','utils.R'))
source(here('src','cancerness.utils.R'))

######################
# Figure 4a: CRC microbiome scores of adenoma samples 

set.seed(200)
meta.ad <-  read_tsv(here('data','Metadata.all.samples.tsv')) %>% 
  filter(Condition=='smallAD' | Condition=='AdvAD' | Condition=='AD') %>% 
  mutate(Condition_general= gsub('AdvAD','AD',Condition)) %>% 
  mutate(Condition_general= gsub('smallAD','AD',Condition))


all.data <- read.table(here('data','Relab.all.samples.tsv'),sep='\t',check.names = F)  %>%
  rownames_to_column('genus') %>%
  filter(genus!='unassigned') %>%
  column_to_rownames('genus')


ad.data<- all.data[,meta.ad$Sample_ID]


load(here('data', 'results', 'Training.unified.crc.model.Rdata'))


dummyControls <- data.frame(dummySample = rpois(n = dim(ad.data)[1], lambda = 10))

rownames(dummyControls) <- rownames(ad.data)
ad.data <- cbind(ad.data, dummyControls)


# SIAMCAT fails if not all training features are found at testing time. so add them.
tmp <- lapply(colnames(ad.data), function(x) {
  tmp <- data.frame(dummyColumnName = rep(0, sum(!norm_params(models.all.rf)$retained.feat %in% rownames(ad.data))))
  rownames(tmp) <- norm_params(models.all.rf)$retained.feat[!norm_params(models.all.rf)$retained.feat %in% rownames(ad.data)]
  colnames(tmp)[1] <- x
  return(tmp)
})

ad.data <- rbind(ad.data, as.data.frame(tmp, check.names = F))
dummyMeta_ad <- data.frame(label = c(rep('CRC', dim(ad.data)[2] - 1 ), 'Control'))
rownames(dummyMeta_ad) <- colnames(ad.data)


all_features<-unique(norm_params(models.all.rf)$retained.feat)

ad.data <- ad.data %>% 
  rownames_to_column('rn') %>%
  complete(rn = union(rn, all_features[which(!all_features%in%rownames(ad.data))])) %>% 
  column_to_rownames('rn') %>%
  replace(is.na(.), 0)


ad.crc.score <- get_cancerness(siam = models.all.rf, pro = ad.data, meta = dummyMeta_ad ,evaluated.label='Adenoma', evaluated.color = '#FFBA08')


write.table(ad.crc.score[[2]], file=here('data', 'results', 'Adenoma.CRC.microbiome.signature.scores.tsv'), sep='\t',quote = F,row.names = F)

ggsave(ad.crc.score[[1]], file= here('figures','figure4','Figure4a.pdf'),width = 7, height = 7)


######################
# Figure 4b: Stage evaluation with positivity rate

all.meta<- (models.all.rf@phyloseq@sam_data) 

all.meta$sampleID <- rownames(all.meta)

CvpredMatrix <- models.all.rf@pred_matrix %>%
  as.data.frame() %>%
  rownames_to_column('sampleID') %>%
  as_tibble() %>%
  pivot_longer(-sampleID) %>%
  group_by(sampleID) %>%
  summarize(medianPredictionProb = median(value)) %>%
  mutate(type = 'evaluated') %>%
  mutate(alpha = 1) %>% 
  left_join(all.meta %>% as.data.frame(), by = 'sampleID') 



roc.curve <- roc(cases=CvpredMatrix %>% 
                   filter(Condition=='CRC') %>% 
                   pull(medianPredictionProb),
                 controls=CvpredMatrix %>% 
                   filter(Condition=='CTR') %>% 
                   pull(medianPredictionProb))

threshold <- roc.curve$thresholds[which(roc.curve$specificities >= 0.9)[1]]

df.plot <- CvpredMatrix %>% 
  mutate(Stage = case_when(
    Stage == 0 ~ 'Stage 0', 
    Stage == 1 ~ 'Stage I', 
    Stage == 2 ~ 'Stage II',
    Stage == 3 ~ 'Stage III',
    Stage == 4 ~ 'Stage IV',
    TRUE ~ as.character(Stage)
  )) %>% 
  mutate(Stage= case_when(Condition=='CTR' & Stage=='Stage 0' ~ NA, 
                          TRUE ~ Stage)) %>% 
  group_by(Stage) %>% 
  summarize(tp=sum(medianPredictionProb > threshold), n=n()) %>% 
  mutate(tpr=tp/n) %>% 
  mutate(StageLabel= paste0(as.character(Stage), '\n (N:' , as.character(n),')' )) %>% 
  drop_na(Stage)

# Load AD predictions 

ad.predictions <- ad.crc.score[[2]] %>% 
  filter(type=='evaluated')  %>% 
  left_join(meta.ad %>% select(Sample_ID, Condition), by=c('sampleID'='Sample_ID')) %>% 
  filter(Condition=='smallAD'| Condition=='AdvAD') 


ad.roc <- roc(cases= ad.predictions %>% 
                      filter(label=='Adenoma') %>% 
                      pull(medianPredictionProb),
                    controls=CvpredMatrix %>% 
                      filter(Condition=='CTR') %>% 
                      pull(medianPredictionProb)) 

threshold <- ad.roc$thresholds[which(ad.roc$specificities >= 0.9)[1]]

ad.df.plot <- ad.predictions %>% 
  rename(Stage = Condition) %>% 
  mutate(Stage = as.factor(Stage)) %>% 
  group_by(Stage) %>% 
  summarize(tp=sum(medianPredictionProb > threshold), n=n()) %>% 
  mutate(tpr=tp/n) %>% 
  mutate(StageLabel= paste0(as.character(Stage), '\n (N:' , as.character(n),')' )) %>% 
  drop_na(Stage)

df.plot<- df.plot %>%  rbind(ad.df.plot) %>%  
  mutate(Stage = factor(Stage, levels = c(
    'smallAD',
    'AdvAD',
    'Stage 0',
    'Stage I',
    'Stage II',
    'Stage III',
    'Stage IV'
  )))



# calculate prevalance of top 4 genera for each stage

top_genera <- c('Fusobacterium', 'Peptostreptococcus', 'Parvimonas', 'Porphyromonas')

# Calculate genus prevalence in CTR Samples to determine  positivity threshold
ctr_data <- models.all.rf@phyloseq@otu_table %>%
  as.data.frame() %>%
  rownames_to_column('genus') %>%
  pivot_longer(-genus, names_to = 'SampleID', values_to = 'Relab') %>%
  filter(genus %in% top_genera) %>%
  left_join(models.all.rf@phyloseq@sam_data %>% data.frame() %>%  rownames_to_column('SampleID'), by = 'SampleID', multiple = 'any') %>%
  filter(Condition == 'CTR')

# Calculate a threshold for positivity (e.g., 90 th percentile of relative abundance in CTR)
prevalence_threshold <- ctr_data %>%
  group_by(genus) %>%
  summarize(threshold = quantile(Relab, probs = 0.9), .groups = 'drop')


crc_prevalence_data_stages <- models.all.rf@phyloseq@otu_table %>%
  as.data.frame() %>%
  rownames_to_column('genus') %>%
  pivot_longer(-genus, names_to = 'SampleID', values_to = 'Relab') %>%
  filter(genus %in% top_genera) %>%
  left_join(models.all.rf@phyloseq@sam_data %>% data.frame() %>%  rownames_to_column('SampleID'), by = 'SampleID', multiple = 'any') %>%
  filter(!(Condition=='CTR' & Stage==0)) %>%
  left_join(prevalence_threshold, by = 'genus') %>% # Join the threshold
  group_by(Stage, genus) %>%
  summarize(Prevalence = mean(Relab >= threshold) * 100, .groups = 'drop')



all.data <- read.table(here('data','Relab.all.samples.tsv'),sep='\t',check.names = F)  %>%
  rownames_to_column('genus') %>%
  filter(genus!='unassigned') %>%
  column_to_rownames('genus')


ad_prevalence_data_stages<- all.data[,meta.ad$Sample_ID] %>% rownames_to_column('genus') %>%
  pivot_longer(-genus, names_to = 'SampleID', values_to = 'Relab') %>%
  filter(genus %in% top_genera) %>%
  left_join(meta.ad %>%  select(Stage=Condition, Sample_ID) ,  by = c('SampleID'='Sample_ID')) %>%
  left_join(prevalence_threshold, by = 'genus') %>% # Join the threshold
  group_by(Stage, genus) %>%
  summarize(Prevalence = mean(Relab >= threshold) * 100, .groups = 'drop') %>% # Apply threshold
  filter(Stage!='AD')





prevalance_all_stages<- rbind(crc_prevalence_data_stages, ad_prevalence_data_stages) %>% 
  drop_na() %>% 
  mutate(Stage= case_when(Stage==0 ~ 'Stage 0',
                          Stage==1 ~ 'Stage I',
                          Stage==2 ~'Stage II',
                          Stage ==3 ~ 'Stage III',
                          Stage == 4 ~ 'Stage IV',
                          TRUE ~Stage)) %>%
  left_join(df.plot %>% select(Stage, StageLabel), by='Stage') %>%
  mutate(StageLabel = factor(StageLabel, levels = c(
    'smallAD\n (N:590)',
    'AdvAD\n (N:766)',
    'Stage 0\n (N:228)',
    'Stage I\n (N:418)',
    'Stage II\n (N:329)',
    'Stage III\n (N:453)',
    'Stage IV\n (N:185)'
  ))) %>%
  mutate(genus=factor(genus, levels=rev(c('Peptostreptococcus','Fusobacterium','Parvimonas','Porphyromonas'))))


prev_heatmap_stages<- ggplot(prevalance_all_stages, aes(x = genus, y = StageLabel, fill = Prevalence)) +
  geom_tile(color = 'white') +  # White grid lines for separation
  geom_text(aes(label = sprintf('%.1f%%', Prevalence)), vjust = 1,color='black') +
  scale_fill_gradient(low = 'white', high = 'black',limits = c(0, 100)) +  # Heatmap colors from low to high prevalence
  theme_minimal() +
  labs(
    x = 'Tumor Stage',
    y = 'Genus',
    fill = 'Positivity (%)') +
  theme_paper +
  theme(
    #axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = 'right',
    axis.text.x = element_text(size = 12),
    axis.ticks.x = element_line(),
    axis.ticks.x.top = element_line(),
    axis.ticks.x.bottom = element_blank(),
    axis.text.x.top = element_text(size = 12, angle = 90),
    axis.text.x.bottom = element_blank(),
    axis.line.x.top = element_line(),
    axis.line.x.bottom = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.ticks.y = element_blank()) +
  scale_x_discrete(position = 'top') +
  scale_y_discrete(limits = rev(levels(prevalance_all_stages$StageLabel)))


g.a <- df.plot %>%
  drop_na(Stage) %>%
  mutate(StageLabel = factor(StageLabel, levels = rev(c(
    'smallAD\n (N:590)',
    'AdvAD\n (N:766)',
    'Stage 0\n (N:228)',
    'Stage I\n (N:418)',
    'Stage II\n (N:329)',
    'Stage III\n (N:453)',
    'Stage IV\n (N:185)'
  )))) %>%
  ggplot(aes(x = tpr , y = StageLabel , fill = Stage)) +
  geom_bar(stat = 'identity', col = 'black') + 
  theme_classic() + 
  scale_x_continuous(limits = c(0, 1)) +
  xlab('') +  
  ylab('True positive rate') +  # Label y-axis
  scale_fill_manual(values = c('#FFBA08','#F48C06','#E85D04','#DC2F02','#D00000','#9D0208','#6A040F'), guide = FALSE) +  # Manual fill colors
  theme_paper + 
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 14),
    axis.text.y= element_text(size = 14),
    panel.grid.major.y = element_line(colour = 'lightgrey')
  ) 

final_plot <- plot_grid(
  prev_heatmap_stages, g.a,
  align = 'v',
  axis = 'l', 
  rel_heights  = c(0.3, 1),
  rel_widths = c(1,1),
  nrow = 2)

ggsave(final_plot, filename = here('figures','figure4','Figure4b.pdf'), width = 7 ,height = 12)

######################
# Figure 4c: Variance explained by tumor stage, location, geographic location, age and T2D status

