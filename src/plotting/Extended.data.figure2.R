######################
# Extended data figure 2
######################

source(here('src','utils.R'))
params <- yaml::read_yaml(here("src", "parameters.yml"))
plotting <- params$plotting

# Load LMM tables
load(here('data','results', 'lmm.tables.eo.lo.Rdata'))


######################
# Extended data figure 2a

volcano_lo <- plot_volcano(
  plot_df = lmm.table.lo %>%
    select(Taxon, P.val, P.adj, Effect.size, pr.shift, `pr.LO-CRC`, `pr.LO-CTR`, `n.LO-CTR`, `n.LO-CRC`),
  fdr_thresh = 0.05,
  group_case = 'LO-CRC',
  group_control = 'LO-CTR',
  feature_column_name = 'Taxon',
  max.overlaps = 20,
  nudge_y = 10,
  color_vector = c('LO-CRC' = plotting$condition_colors$`LO-CRC`, 'LO-CTR'= plotting$condition_colors$`LO-CTR`, "n.s." = "white")) +
  scale_x_continuous(limits=c(-0.5,0.9), breaks = c(-0.5, -0.25, 0, 0.25, 0.5,0.75))+
  xlab('Late-onset enrichment effect size')

ggsave(volcano_lo,file= here('figures','extended.data.figure2','Extended.Data.Figure2a.pdf'), height = 5, width = 5)

######################
# Extended data figure 2b

volcano_eo <- plot_volcano(
  plot_df = lmm.table.eo %>%
    select(Taxon, P.val, P.adj, Effect.size, pr.shift, `pr.EO-CRC`, `pr.EO-CTR`, `n.EO-CTR`, `n.EO-CRC`),
  fdr_thresh = 0.05,
  group_case = 'EO-CRC',
  group_control = 'EO-CTR',
  feature_column_name = 'Taxon',
  max.overlaps = 20,
  color_vector = c('EO-CRC' = plotting$condition_colors$`EO-CRC`, 'EO-CTR' = plotting$condition_colors$`EO-CTR`, "n.s." = "white")) +
  scale_x_continuous(limits=c(-0.5,0.9), breaks = c(-0.5, -0.25, 0, 0.25, 0.5,0.75))+
  xlab("Early-onset enrichment effect size")


ggsave(volcano_eo, file=here('figures','extended.data.figure2','Extended.Data.Figure2b.pdf'), height = 5, width = 5)

######################
# Extended data figure 2c

lmm.table.eo.lo.crc <- lmm.table.eo.lo.crc %>%  #some dataset does not have either EO-CRC or LO-CRC samples
  select(where(~ !all(is.na(.))))

volcano_eolo_crc<- plot_volcano(
  plot_df = lmm.table.eo.lo.crc %>% 
    select(Taxon, P.val, P.adj, Effect.size, pr.shift, 'pr.LO-CRC', 'pr.EO-CRC', 'n.LO-CRC', 'n.EO-CRC'),
  group_case = 'EO-CRC',
  group_control = 'LO-CRC',
  feature_column_name = 'Taxon',
  max.overlaps = 10,
  min_segment_length = 0.05,
  nudge_y = 2,
  color_vector = c( 'EO-CRC' = plotting$condition_colors$`EO-CRC`, 'LO-CRC' = plotting$condition_colors$`LO-CRC`,"n.s." = "white") ) +
  xlab('Enrichment effect size')+
  theme(
    legend.box = "horizontal",
    legend.spacing.y = unit(0.1, "cm"),
    legend.spacing.x = unit(0.03, "cm"),
    legend.position = c(0.44, .99),
    legend.key.size = unit(0.75, "lines"),
    legend.justification = c(0, 1)
  )+
  scale_x_continuous(limits = c(-0.5,0.5),breaks = c(-0.5, -0.25, 0, 0.25, 0.5))


ggsave(volcano_eolo_crc, file=here('figures','extended.data.figure2','Extended.Data.Figure2c.pdf'), height = 5, width = 5)

######################
# Extended data figure 2d

lmm.table.eo.lo.ctr <- lmm.table.eo.lo.ctr %>%  #some dataset does not have either EO-CRC or LO-CRC samples
  select(where(~ !all(is.na(.))))

volcano_eolo_ctr<- plot_volcano(
  plot_df = lmm.table.eo.lo.ctr %>% 
    select(Taxon, P.val, P.adj, Effect.size, pr.shift, 'pr.LO-CTR', 'pr.EO-CTR', 'n.LO-CTR', 'n.EO-CTR'),
  group_case = 'LO-CTR',
  group_control = 'EO-CTR',
  feature_column_name = 'Taxon',
  max.overlaps = 10,
  min_segment_length = 0.05, 
  nudge_y = 0.1,
  color_vector = c('LO-CTR' = plotting$condition_colors$`LO-CTR`, 'EO-CTR' = plotting$condition_colors$`EO-CTR`, "n.s." = "white") # Custom colors
) +
  xlab('Enrichment effect size')+
  theme(
    legend.box = "horizontal",
    legend.spacing.y = unit(0.1, "cm"),
    legend.spacing.x = unit(0.03, "cm"),
    legend.position = c(0.02, .99),
    legend.key.size = unit(0.75, "lines"),
    legend.justification = c(0, 1)
  )+
  scale_x_continuous(limits = c(-0.25,0.25),breaks = c(-0.25,-0.1, 0, 0.1,0.25))

ggsave(volcano_eolo_ctr, file=here('figures','extended.data.figure2','Extended.Data.Figure2d.pdf'), height = 5, width = 5)

######################
# Extended data figure 2e

load(here('data','results','Training.eo.lo.crc.rf.models.Rdata'))

models <- list(models.lo.rf, siamcat.test.evaluated.lo.holdout.rf)
labels <- c("Classifier cross validated on LO-CRC", 'Classifier trained on EO-CRC and tested on LO-CRC"')
trained_on <- list(NULL, models.eo.rf)
colours <- c( "darkred" , 'black')

lo_models_auc_plot<- plot_roc_siamcat_models(models, labels, colours, trained_on,alpha=0.8)

ggsave(lo_models_auc_plot, file=here('figures','extended.data.figure2','Extended.Data.Figure2e.pdf'), height = 7, width = 7)


