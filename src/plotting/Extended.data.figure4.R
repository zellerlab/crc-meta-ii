######################
# Extended data figure 4
######################

source(here('src','utils.R'))
params <- yaml::read_yaml(here("src", "parameters.yml"))
plotting <- params$plotting
condition_colors <- params$plotting$condition_colors

######################
# Extended data figure 4a: Volcano plot comparing CRC vs CTR

all.meta<- read_tsv(here('data','Metadata.all.samples.tsv')) %>% 
  filter(Condition=='CRC' | Condition=='CTR') 

all.data <- read.table(here('data','Relab.all.samples.tsv'),sep='\t',check.names = F)  %>%
  rownames_to_column('genus') %>%
  filter(genus!='unassigned') %>%
  column_to_rownames('genus')

all.data <- all.data[which(rowSums(all.data > 0) / ncol(all.data) > 0.1),]

all.data <- all.data[,all.meta$Sample_ID]

lmm.table.general.crc <- run_lmem(
  data_df = all.data,
  meta_df  = all.meta, column_name='Condition', ref_group='CTR', feature_column_name = 'Taxon')

write.table(lmm.table.general.crc, file= here('data','results' ,'lmm.tables.ctr.crc.allsamples.tsv'), sep='\t', quote = F,row.names = F)


volcano_crc_ctr <- plot_volcano(
  plot_df = lmm.table.general.crc %>%
    select(Taxon, P.val, P.adj, Effect.size, pr.shift, pr.CTR, pr.CRC, n.CRC, n.CTR),
  group_case = 'CRC',
  group_control = 'CTR',
  feature_column_name = 'Taxon',
  min_segment_length = 0.05, nudge_y = 0.05, max.overlaps = 5,
  color_vector = c( condition_colors$CRC, condition_colors$CTR,  'n.s.' = 'white') # Custom colors
) +
  xlab('Effect size')+
  theme(
    legend.box = 'horizontal',
    legend.spacing.y = unit(0.1, 'cm'),
    legend.position = c(0, .99),
    legend.key.size = unit(0.75, 'lines'),
    legend.justification = c(0, 1)
  )

ggsave(volcano_crc_ctr,file= here('figures','extended.data.figure4','Extended.Data.Figure4a.pdf'), height = 5, width = 5)

######################
# Extended data figure 4b: Volcano plot comparing AD vs CRC

meta.ad.crc <-  read_tsv(here('data','Metadata.all.samples.tsv')) %>% 
  filter(Condition=='smallAD' | Condition=='AdvAD' | Condition=='AD'| Condition=='CRC') %>% 
  mutate(Condition= gsub('AdvAD','AD',Condition)) %>% 
  mutate(Condition= gsub('smallAD','AD',Condition))

all.data.ad.crc <- all.data[,meta.ad.crc$Sample_ID]

lmm.table.ad.vs.crc <- run_lmem(
  data_df = all.data.ad.crc,
  meta_df  = meta.ad.crc, column_name='Condition', ref_group='CRC', feature_column_name = 'Taxon')

write.table(lmm.table.ad.vs.crc, file= here('data','results' ,'lmm.tables.ad.crc.tsv'), sep='\t', quote = F,row.names = F)
lmm.table.ad.vs.crc<- read_tsv( here('data','results' ,'lmm.tables.ad.crc.tsv'))


volcano_ad_crc <- plot_volcano(
  plot_df = lmm.table.ad.vs.crc %>%
    select(Taxon, P.val, P.adj, Effect.size, pr.shift, pr.AD, pr.CRC, n.CRC, n.AD),
  group_case = 'CRC',
  group_control = 'AD',
  feature_column_name = 'Taxon',
  min_segment_length = 0.05, nudge_y = 0.05, max.overlaps = 5,
  color_vector = c( condition_colors$CRC, condition_colors$AD,  'n.s.' = 'white') # Custom colors
) +
  xlab('Effect size')+
  theme(
    legend.box = 'horizontal',
    legend.spacing.y = unit(0.1, 'cm'),
    legend.position = c(0, .99),
    legend.key.size = unit(0.75, 'lines'),
    legend.justification = c(0, 1)
  )

ggsave(volcano_ad_crc,file= here('figures','extended.data.figure4','Extended.Data.Figure4b.pdf'), height = 5, width = 5)

######################
# Extended data figure 4c: Volcano plot comparing AD vs CTR

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

volcano_ad_ctr <- plot_volcano(
  plot_df = lmm.table.ad.vs.ctr %>%
    select(Taxon, P.val, P.adj, Effect.size, pr.shift, pr.AD, pr.CTR, n.CTR, n.AD),
  group_case = 'AD',
  group_control = 'CTR',
  feature_column_name = 'Taxon',
  min_segment_length = 0.05, nudge_y = 0.05, max.overlaps = 5,
  color_vector = c(AD = condition_colors$AD, CTR = condition_colors$CTR, 'n.s.' = 'white') # Custom colors
) +
  xlab('Effect size')+
  theme(
    legend.box = 'horizontal',
    legend.spacing.y = unit(0.1, 'cm'),
    legend.position = c(0.35, .99),
    legend.key.size = unit(0.75, 'lines'),
    legend.justification = c(0, 1)
  )

  ggsave(volcano_ad_ctr,file= here('figures','extended.data.figure4','Extended.Data.Figure4c.pdf'), height = 5, width = 5)

######################
# Extended data figure 4d: Scatter plot comparing effect sizes in CRC vs CTR and AD vs CTR

scatter.data <- lmm.table.general.crc %>% 
  select(Taxon, Effect.size,P.adj) %>% 
  left_join(lmm.table.ad.vs.ctr %>%select(Taxon, Effect.size,P.adj), by='Taxon' ,suffix = c('.CRC','.AD')) 


plot_comparison_scatter_ad <- function(
    data, x_col, y_col, x_label, y_label,
    p_col_x, p_col_y, feature_column_name = "Bacteria",
    effect_thresh = 0, padj_thresh = 0.05, theme_obj = theme_paper
) {
  
  x <- data %>%
    dplyr::mutate(
      sig_x = !is.na(.data[[p_col_x]]) & (.data[[p_col_x]] < padj_thresh),
      sig_y = !is.na(.data[[p_col_y]]) & (.data[[p_col_y]] < padj_thresh),
      up_x  = !is.na(.data[[x_col]])   & (.data[[x_col]]  >  effect_thresh),
      dn_x  = !is.na(.data[[x_col]])   & (.data[[x_col]]  < -effect_thresh),
      up_y  = !is.na(.data[[y_col]])   & (.data[[y_col]]  >  effect_thresh),
      dn_y  = !is.na(.data[[y_col]])   & (.data[[y_col]]  < -effect_thresh),
      enriched_in = dplyr::case_when(
        (sig_x & up_x) & (sig_y & up_y) ~ "both in CRC and AD",
        (sig_x & dn_x) & (sig_y & dn_y) ~ "CTR",
        TRUE                             ~ "n.s."
      ),
      enriched_in = factor(enriched_in, levels = c("both in CRC and AD", "CTR", "n.s.")),
      label = case_when(enriched_in!="n.s." ~.data[[feature_column_name]],
                        Taxon %in% c("Fusobacterium","Parvimonas","Peptostreptococcus","Porphyromonas","Gemella") ~ .data[[feature_column_name]],
                        TRUE ~ ""),
      font = case_when(enriched_in!="n.s." ~ "bold.italic",
                        Taxon %in% c("Fusobacterium","Parvimonas","Peptostreptococcus","Porphyromonas","Gemella") ~ "bold.italic",
                        TRUE ~ "")) %>%
    dplyr::select(-sig_x, -sig_y, -up_x, -dn_x, -up_y, -dn_y)
  
  axis_limits <- range(c(data[[x_col]], data[[y_col]]), na.rm = TRUE)
  padding <- diff(axis_limits) * 0.05   # 5% margin
  axis_limits <- c(axis_limits[1] - padding, axis_limits[2] + padding)
  
  ggplot(x, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(aes(fill = enriched_in), shape = 21, color = "black", alpha = 0.75, size = 3) +
    geom_hline(yintercept = 0, color = "grey", lty = "dashed", lwd = 0.3) +
    geom_vline(xintercept = 0, color = "grey", lty = "dashed", lwd = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "grey", lty = "dashed", lwd = 0.3) +
    ggrepel::geom_text_repel(
      aes(label = label, fontface = font),
      color = "black", segment.color = "black",
      segment.size = 0.2, segment.ncp = 3, max.overlaps = 25,
      min.segment.length = 0.2, nudge_x = 0.05, nudge_y = 0.05,
      size = 3, seed = 123
    ) +
    scale_fill_manual(values = c("#E19448",  "dodgerblue3", "white")) +
    theme_obj +
    theme(
      legend.box = "horizontal",
      legend.spacing.y = unit(0.1, "cm"),
      legend.position = c(.01, .99),
      legend.key.size = unit(0.75, "lines"),
      legend.justification = c(0, 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    labs(fill = paste0("p.adj < ", padj_thresh), x = x_label, y = y_label) +
    xlim(axis_limits) + ylim(axis_limits) +
    coord_equal()
}


plot.comp.effect.ad.crc<- plot_comparison_scatter_ad(data=scatter.data, 
                                                      x_col='Effect.size.CRC', y_col='Effect.size.AD',
                                                      p_col_x = 'P.adj.CRC', p_col_y = 'P.adj.AD',
                                                      x_label= 'Effect size (CTR vs CRC)' , y_label = 'Effect size (CTR vs AD)',
                                                      feature_column_name = 'Taxon')


ggsave(plot.comp.effect.ad.crc, file= here('figures','extended.data.figure4','Extended.Data.Figure4d.pdf'), height = 5, width = 5)


######################
# Extended data figure 4e: CRC microbiome signature scores in AD samples

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

# Load the universal CRC trained model
load(here('data', 'results', 'Training.unified.crc.model.Rdata'))

dummyControls <- data.frame(dummySample = rpois(n = dim(ad.data)[1], lambda = 10))

rownames(dummyControls) <- rownames(ad.data)
ad.data <- cbind(ad.data, dummyControls)


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

# Save the plot
ggsave(ad.crc.score[[1]], file= here('figures','extended.data.figure4','Extended.Data.Figure4e.pdf'),width = 7, height = 7)

######################
# Extended data figure 4h: Shap value calculation for AD model 

shap_values<- read_tsv(here("data","results","shap.analysis", "AD_median_shap_value.tsv")) %>% 
  mutate(feature= as.factor(feature))

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
  slice(1:25) %>%
  pull(feature)

shap_values <- shap_values %>%
  filter(feature %in% top_features) %>%
  left_join(feature_order %>% select(feature, feature_ordered), by = "feature") %>%
  mutate(
    feature_ordered = factor(feature, levels = levels(feature_order$feature_ordered))
  )

# Plotting
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
    colors = c("dodgerblue3", "white", "#FEE391"),
    limits = c(-3, 3),
    oob = scales::squish
  ) +
  scale_shape_manual(values = c("AD" = 16, "CTR" = 1)) +
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


plot_up <- ggplot(plot_up_data, aes(x = perc_signed / 100, y = feature_ordered)) +
  geom_col(fill = 'grey66') +
  geom_text(
    aes(
      label = paste0(round(percentage, 1), '%'),
      x = perc_signed + ifelse(perc_signed >= 0, 1, -1), 
      hjust = ifelse(perc_signed >= 0, 0, 1)  
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
ggsave(plot, file=here('figures','extended.data.figure4','Extended.Data.Figure4h.pdf'), dpi=300, width=8, height=6)
