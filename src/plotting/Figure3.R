######################
# Figure 3
######################

source(here('src','utils.R'))
params <- yaml::read_yaml(here('src', 'parameters.yml'))
plotting <- params$plotting

######################
# Figure 3a: World map showing sample distribution 

options(scipen = 999)

world <- map_data('world')

meta<- read.table(here( 'data', 'Metadata.all.samples.balanced.tsv'), header = T, sep = '\t') %>% select(Country, Cohort) 


country.count<- meta %>%  group_by(Country) %>% unique() %>% 
  summarise(Number_of_Cohorts = n())


country_names <- c(
  ARG = 'Argentina', AUS = 'Australia', CAN = 'Canada', CHL = 'Chile', CHN = 'China',
  FRA = 'France', GER = 'Germany', IND = 'India', ITA = 'Italy', JPN = 'Japan',
  LUX = 'Luxembourg', TWN = 'Taiwan', UK = 'UK', USA = 'USA', VNM = 'Vietnam'
)

sample_counts_per_country <- meta %>%
  dplyr::count(Country)


sample_counts_per_country <- sample_counts_per_country %>%
  mutate(
    region = case_when(
      Country %in% names(country_names) ~ country_names[Country],
      TRUE ~ Country  # Keep original name if not in mapping
    )
  )

# Merge world map data with country counts
world_data <- left_join(world, sample_counts_per_country, by = 'region')  %>% filter(region != 'Antarctica') %>% 
  mutate(Cohort_avail=ifelse(is.na(Country), 'No', 'Yes')) %>% 
  left_join(country.count, by='Country')

# First, calculate centroids for each country
country_centroids <- world_data %>%
  group_by(region) %>%
  summarise(
    long = mean(long, na.rm = TRUE),
    lat = mean(lat, na.rm = TRUE),
    Country = first(Country),
    Number_of_Cohorts = first(Number_of_Cohorts),
    n = first(n)
  ) %>%
  filter(!is.na(Country))

# Create a function to generate offset coordinates for labels
generate_offset <- function(x, y, angle, distance) {
  data.frame(
    x = x + cos(angle) * distance,
    y = y + sin(angle) * distance
  )
}

# Generate offset coordinates for each country
set.seed(123)  # for reproducibility
country_centroids <- country_centroids %>%
  rowwise() %>%
  mutate(
    angle = runif(1, 0, 2*pi),
    distance = runif(1, 10, 40),
    offset = list(generate_offset(long, lat, angle, distance))
  ) %>%
  unnest(offset)

# Define the countries you want to adjust
countries_to_adjust <- c('FRA', 'GER', 'ITA', 'CHN', 'TWN','UK')

# Adjust the x and y coordinates for specific countries
country_centroids <- country_centroids %>%
  mutate(
    x = ifelse(Country %in% countries_to_adjust, x + 3, x),
    y = ifelse(Country %in% countries_to_adjust, y + 1, y)
  )

# Now plot with the adjusted coordinates
worldplot <- ggplot() +
  geom_polygon(data = world_data, aes(x = long, y = lat, group = group, fill = Cohort_avail), 
               color = 'black', linewidth = 0.2) +
  geom_segment(data = country_centroids, 
               aes(x = long, y = lat, xend = x, yend = y),
               color = 'gray50', linewidth = 0.5) +
  geom_label_repel(data = country_centroids,
                            aes(x = x, y = y, 
                                label = paste0(Country, ': ', Number_of_Cohorts, 
                                               ' ', ifelse(Number_of_Cohorts == 1, 'study', 'studies'), 
                                               ', N = ', n)),
                            size = 3.5, fill = 'white', color = 'black', fontface = 'bold',
                            label.padding = unit(0.4, 'lines'),
                            max.overlaps = Inf,
                            segment.color = 'gray50')+
  coord_fixed(1.3) +
  scale_fill_manual(values = c(No = 'gray', Yes = '#4A0E4E'),
                    name = 'Countries',
                    labels = c(No = 'Not included', Yes = 'Included')) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = 'bottom'
  )


ggsave(plot = worldplot, file= here('figures','figure3','Figure3a.pdf') , height = 8,width = 10)

######################
# Figure 3b: AUROC for unified CRC model

# Load unified model

load(here('data','results','Training.unified.crc.model.Rdata'))

models <- list(models.all.rf)
labels <- c('Classifier cross validated on CRC')
colours <- c( 'black' )

crc_model_auc_plot<- plot_roc_siamcat_models(models, labels, colours, alpha=0.8)

ggsave(plot = crc_model_auc_plot, file= here('figures','figure3','Figure3b.pdf') ,width = 6, height = 6)

######################
# Figure 3c: Performance of SCV and LOSO models across individual 16S and WGS datasets

# Load CRC LOSO models
path.CRC.LOSO.test<- list.files(path = here('data','results','scv.loso','crc.loso.test') ,
                                     pattern = 'Rdata',full.names = T)

LOSO.evaluated.CRC <- list()

# Loop through each RData file, load it, and extract the relevant list element
for (path in path.CRC.LOSO.test) {
  load(path)
  list_name <- names(loso.eval.crc)[1] 
  LOSO.evaluated.CRC[[list_name]] <- loso.eval.crc[[list_name]]
}


# Load SCV and SST models 
load(here('data','results','scv.loso','crc.scv.sst','Models.SCV.SST.Rdata'))

SCV.CRC<- models.scv
SST.eval.CRC<- holdout.evaluation.crc

# Get AUC's and number of samples
SCV_metrics <- get_siamcat_holdout_metrics(siamcat_list_cv = SCV.CRC, siamcat_list_holdout = SST.eval.CRC)

LOSO_CRC_metrics<- get_siamcat_loso_metrics(LOSO.evaluated.CRC)


# Combine SCV and LOSO dataframes
df <- rbind(SCV_metrics, LOSO_CRC_metrics)

# Get the assay type for hold-out test datasets from metadata
assay <- read.table(here('data', 'Metadata.all.samples.balanced.tsv'), sep = '\t', header = T) %>%
  select(Cohort, Assay) %>%
  unique()

# Merge and process the data
df <- df %>%
  mutate(
    Type = factor(Type, levels = c('CV', 'SST', 'LOSO'), ordered = TRUE),
    AUC = as.numeric(AUC)
  ) %>%
  left_join(assay, by = c('Test' = 'Cohort')) %>%
  mutate(Assay = factor(Assay, levels = c('WGS', '16S'))) %>%
  mutate(Test = factor(Test, levels = unique(df %>% arrange(Assay) %>% pull(Test)))) %>%
  arrange(Assay)


all.meta<- read.table(here('data','Metadata.all.samples.balanced.tsv'), sep = '\t', header = T)  %>%
  filter(Condition=='CRC'| Condition=='CTR') %>% as.data.frame() 

# Load seqdepth and genus richness 


seqdepth<- read.table(here('data','Sequencing_depth.tsv'),sep = '\t',header = T)

seqdepth.mean <- seqdepth %>% 
  select(Sample_ID, Assay,Seq_depth, Cohort) %>% 
  mutate(Cohort =factor(Cohort, levels=c(unique(df$Test)))) %>% 
  mutate(Assay= factor(Assay, levels=c('WGS','16S'))) %>% 
  group_by(Cohort) %>% 
  summarise(mean_seqdepth= mean(Seq_depth))


# Summarize the data for plotting
df_summary <- df %>%
  group_by(Type, Test) %>%
  summarise(mean_AUC = mean(AUC),
            sd_AUC = sd(AUC),
            .groups = 'drop')

df_average <- df_summary %>%
  group_by(Type) %>%
  summarise(std_AUC = sd(mean_AUC),
            mean_AUC = mean(mean_AUC),
            .groups = 'drop')


# Precompute ymin and ymax for the CV and LOSO rectangles
cv_ymin <- subset(df_average, Type == 'CV')$mean_AUC - subset(df_average, Type == 'CV')$std_AUC
cv_ymax <- subset(df_average, Type == 'CV')$mean_AUC + subset(df_average, Type == 'CV')$std_AUC

loso_ymin <- subset(df_average, Type == 'LOSO')$mean_AUC - subset(df_average, Type == 'LOSO')$std_AUC
loso_ymax <- subset(df_average, Type == 'LOSO')$mean_AUC + subset(df_average, Type == 'LOSO')$std_AUC


# Arrange Tests based on Assay type
df <- df %>% 
  mutate(Test = factor(Test, levels = unique(Test[order(Assay)])))

# Define color gradients for each assay
color_gradients <- list(
  'WGS' = c('white', '#9DC4C1', '#01665E'),
  '16S' = c('white', '#83652B', '#543005')
)

get_limits <- function(df, assay_type) {
  df %>%
    filter(Assay == assay_type) %>%
    pull(mean_libsize) %>%
    range()
}

theme_paper <- ggembl::theme_presentation() +
  theme(
    axis.title = element_text(face = 'bold', size = 12), # Bold axis titles
    panel.border = element_rect(fill=NA, colour='black', size=1.5),
    axis.text = element_text(face = 'bold', size = 12),
  )

# Combine Test and Number_of_Samples for x-axis labels
df <- df %>%
  mutate(Test_Label = case_when(Type=='CV' ~ paste0(Test, ' (N = ', Number_of_Samples, ')'), 
                                TRUE ~ NA )) %>% 
  mutate(Test_Label=factor(Test_Label, levels=c(
    'FengQ_2015 (N = 92)', 'GaoR_2021 (N = 32)', 
    'HanniganGD_2018 (N = 32)', 'LiuNN_2022 (N = 158)', 
    'ThomasAM_2019_2 (N = 56)', 'ThomasAM_2019_1 (N = 48)', 
    'VogtmannE_2016 (N = 104)', 'WirbelJ_2019 (N = 120)', 
    'XuJ_2022 (N = 60)', 'YachidaS_2019 (N = 498)', 
    'YangJ_2019 (N = 190)', 'YangY_2021c (N = 198)', 
    'YuJ_2017 (N = 108)', 'ZellerG_2014 (N = 104)', 
    'BaxterNT_2016 (N = 240)', 'DengX_2018 (N = 34)', 
    'DuX_2022 (N = 54)', 'KonishiY_2022 (N = 874)', 
    'OkumuraS_2021 (N = 1026)', 'TernesD_2022 (N = 92)', 
    'UchinoY_2021 (N = 102)', 'YangY_2021a (N = 702)', 
    'YangY_2021b (N = 246)', 'YangTW_2019 (N = 74)', 
    'YoungC2_2021 (N = 868)', 'YoungC1_2021 (N = 80)', 
    'ZhouZ_2022 (N = 36)'
  )))



write.table(df %>% select(-Test_Label), file = here('data','results','ALl.predictions.LOSO.SST.SCV.models.tsv'), sep = '\t' ,quote = F )

# Create the plot
Auc.loso.cv.scv.plot <- ggplot(df) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = cv_ymin, ymax = cv_ymax,
           fill = 'white',color= '#0072B2', alpha = 0.5) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = loso_ymin, ymax = loso_ymax,
           fill = 'white', color= '#882255',alpha = 0.5) +
  geom_hline(data = df_average %>% filter(Type == 'CV'), aes(yintercept = mean_AUC),
             linetype = 'dashed', linewidth = 1.3, color = '#0072B2',alpha=0.5) +
  geom_hline(data = df_average %>% filter(Type == 'LOSO'), aes(yintercept = mean_AUC),
             linetype = 'dashed', linewidth = 1.3, color = '#882255',alpha=0.5) +
  geom_boxplot(data = df %>% filter(Type == 'SST'), aes(x = Test, y = AUC),
               fill = 'white', color = 'black', trim = TRUE, size = 0.75) +
  geom_point(data = df %>% filter(Type == 'LOSO' | Type == 'CV'), aes(x = Test, y = AUC, fill = Type, shape = Type),
             size = 3.5, position = position_dodge2(width = 0.7, preserve = 'total')) +
  scale_shape_manual(name = 'ML Approach', values = c(21, 25)) +
  scale_fill_manual(name = 'ML Approach', values = c(
    'CV' = '#0072B2',
    'LOSO' = '#882255'
  )) +
  ggnewscale::new_scale_fill() + # add new scale for assay
  geom_tile(data = df %>% select(Test, Assay) %>% unique(),
            aes(x = Test, y = 0.39, fill = Assay, width = 1, height = 0.02),
            inherit.aes = FALSE, color = 'black',size=0.5) +
  scale_fill_manual(name = 'Assay', values = c('WGS' = '#287D76', '16S' = '#775721')) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0.38, 1)) +
  scale_x_discrete(labels = unique(na.omit(df$Test_Label)))+  # Use combined labels for x-axis
  theme_paper +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    title = element_text(size = 16),
    legend.background = element_blank(),
    legend.position = 'top',
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  ) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', linewidth = 1.3, color = 'gray') 

# Print the plot
print(Auc.loso.cv.scv.plot)         

# Create the box plot with individual sample points
seqdepth<- seqdepth %>% left_join(df %>% select(Cohort= Train, Test_Label) %>% drop_na(Test_Label), by='Cohort')

seq_depth_plot<-ggplot(seqdepth, aes(x = Test_Label, y = log10(Seq_depth))) +
  geom_boxplot(outlier.shape = NA, color = 'black') +  
  geom_jitter(width = 0.2, alpha = 0.2, aes(color = Assay)) + 
  labs(y = 'log(seq depth)') +
  facet_wrap(~ Assay, scales = 'free',) + 
  scale_color_manual(name = 'Assay', values = c('WGS' = '#287D76', '16S' = '#775721')) +  
  theme_paper +  
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    title = element_text(size = 16),
    strip.text.x = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = 'grey90'),
    panel.grid.minor.y = element_line(color = 'grey90')
  )+
  guides(color='none')



# Adjust Auc.loso.cv.scv.plot to ensure consistent x-axis
Auc.loso.cv.scv.plot <- Auc.loso.cv.scv.plot +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.x = element_blank()  
  )

# Adjust seq_depth_plot to match x-axis alignment
seq_depth_plot <- seq_depth_plot +
  theme(
    axis.text.x = element_text(size = 12, angle = 60, hjust = 1), 
    axis.ticks.x = element_line(),  
    axis.title.x = element_blank()  
  )


diversity_richness_df<- read.table(file = '/g/scb/zeller/pekel/meta_analysis/data/results/Diversity.richness.CRC.CTR.groups.tsv',sep = '\t', header = T) %>% 
  select(-Test_Label)


diversity_richness_df<- diversity_richness_df %>% 
  left_join(df %>% select(Test, Test_Label) %>% drop_na(), 
            by=c('Cohort'='Test')) %>% 
  mutate(Assay=factor(Assay, levels=c('WGS','16S')))


richness<-ggplot(diversity_richness_df, aes(x = Test_Label, y = Genus.richness)) +
  geom_boxplot(outlier.shape = NA, color = 'black') + 
  geom_jitter(width = 0.2, alpha = 0.2, aes(color = Assay)) +  
  scale_color_manual(name = 'Assay', values = c('WGS' = '#287D76', '16S' = '#775721')) +  
  geom_hline(yintercept = 0, linetype= 'dashed', color='grey') +
  theme_paper +  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x= element_blank(),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    title = element_text(size = 16),
    strip.text.x = element_blank(),
    #legend.position = 'left',
    legend.background = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = 'grey90'),
    panel.grid.minor.y = element_line(color = 'grey90')
  )+
  guides(color='none')



# Combine plots with vertical alignment on the x-axis
final_plot <- plot_grid(
  Auc.loso.cv.scv.plot,
  richness,
  seq_depth_plot,
  align = 'v',   # Align vertically to share x-axis
  axis = 'lr',   # Align along left and right axes
  rel_heights = c(0.6,0.15, 0.4),  # Adjust relative heights if needed
  nrow = 3      # Stack plots vertically
)

ggsave(plot = final_plot, file= here('figures','figure3','Figure3c.pdf') , height = 10,width = 10)


######################
# Figure 3d: 

shap_values<- read_tsv(here("data","results","shap.analysis", "Alldata_median_shap_value.tsv")) %>% 
  mutate(feature= as.factor(feature))


#calculate the percantage of the mean absolute shap value too see how much contribute each feature to the unified model

perc_mean_shap<-shap_values %>%
  select(feature, shap_value) %>%
  group_by(feature) %>%
  summarise(mean_abs_shap = mean(abs(shap_value)), .groups = "drop") %>%
  mutate(percentage = mean_abs_shap / sum(mean_abs_shap) * 100)

l <- levels(shap_values$feature) 

shap_values$feature <- factor(shap_values$feature, levels = rev(l))

feature_order <- shap_values %>%
  group_by(feature) %>%
  summarize(
    mean_abs_shap = mean(abs(shap_value), na.rm = TRUE),
    spearman_sign = unique(spearman_sign),
    .groups = "drop"
  ) %>%
  arrange(
    desc(spearman_sign),  # Positive correlations first
    ifelse(spearman_sign == 1, -mean_abs_shap, mean_abs_shap)
  ) %>%
  mutate(feature_ordered = factor(feature, levels = feature))


top_features <- feature_order %>%
  arrange(desc(mean_abs_shap)) %>%
  slice(1:10) %>%
  pull(feature)


shap_values <- shap_values %>%
  filter(feature %in% top_features) %>%
  left_join(feature_order %>% select(feature, feature_ordered), by = "feature") %>%
  mutate(
    feature_ordered = factor(feature, levels = levels(feature_order$feature_ordered))
  )

plot <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(data = shap_values, aes(x = shap_value, y = feature_ordered, color = feature_value),
             position = position_jitter(height = 0.45, width = 0), alpha = 0.5, size = 1) +
  geom_point(data = shap_values %>%
               group_by(feature) %>%
               summarize(n = mean(abs(shap_value)) * spearman_sign) %>%
               distinct(), aes(x = n, y = feature), shape = 18, color = 'black', size = 1.5, inherit.aes = F) +
  geom_point(data = shap_values %>%
               group_by(feature) %>%
               summarize(n = mean(abs(shap_value)) * spearman_sign) %>%
               distinct(), aes(x = n, y = feature), shape = 5, color = 'black', size = 1.5, inherit.aes = F) +
  theme_presentation() +
  coord_flip() +
  scale_color_gradientn(
    colors = c("dodgerblue3", "white", "#C44600"),
    limits = c(-3, 3),
    oob = scales::squish
  ) +
  scale_shape_manual(values = c("CRC" = 16, "CTR" = 1)) +
  xlab("SHAP value") +
  ylab("Genus") +
  theme(
    axis.text.x = element_text(size=6, angle=30, hjust=1), # Smaller x-axis text
    axis.text.y = element_text(size=6), # Smaller y-axis text
    axis.title.x = element_blank(), # Smaller x-axis title
    axis.title.y = element_text(size=7), # Smaller y-axis title
    axis.ticks.length = unit(0.5, "mm"), # Make tick marks shorter
    legend.position = c(0.95, 0.95), # Top-right corner inside the plot
    legend.justification = c("right", "top"), # Align legend to top-right
    legend.direction = "horizontal", # Arrange legend items horizontally
    legend.text = element_text(size=6), # Reduce legend text size
    legend.title = element_text(size=7), # Reduce legend title size
    legend.key.size = unit(4, "mm") # Reduce legend key size
  )

plot_up_data <- shap_values %>%
  select(feature_ordered, spearman_sign) %>%
  distinct() %>%
  left_join(perc_mean_shap %>% select(feature, percentage), by = c('feature_ordered' = 'feature')) %>%
  mutate(perc_signed = percentage * spearman_sign)


plot_up <- ggplot(plot_up_data, aes(x = perc_signed, y = feature_ordered)) +
  geom_col(fill = 'grey66') +
  geom_text(
    aes(
      label = paste0(round(percentage, 1), '%'),
      x = perc_signed + ifelse(perc_signed >= 0, 1, -1),  # small nudge
      hjust = ifelse(perc_signed >= 0, 0, 1)  # align text away from bar
    ),
    fontface = 'bold',
    color = 'black',
    size = 2
  )+

  theme_presentation() +
  coord_flip() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 6), 
    axis.title.y = element_text(size = 7)
  ) +
  xlab("Relative contribution")

# Save the combined plot with adjusted sizes
ggsave(plot=plot_up + plot + plot_layout(heights=c(1,0.8), guides='keep'), file=here('figures','figure3','Figure3d.pdf'), dpi=300, width=2.4, height=4.2, units="in")
