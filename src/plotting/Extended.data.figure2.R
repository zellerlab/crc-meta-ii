######################
# Extended data figure 2
######################

source(here('src', 'utils.R'))
params <- yaml::read_yaml(here('src', 'parameters.yml'))
plotting <- params$plotting

# Load LMM tables
load(here('data', 'results', 'lmm.tables.eo.lo.Rdata'))

######################
# Extended data figure 2a

volcano_lo <- plot_volcano(
  plot_df = lmm.table.lo %>%
    select(Bacteria, P.val, P.adj, Effect.size, pr.shift, `pr.LO-CRC`, `pr.LO-CTR`, `n.LO-CTR`, `n.LO-CRC`),
  fdr_thresh = 0.05,
  group_case = 'LO-CRC',
  group_control = 'LO-CTR',
  feature_column_name = 'Bacteria',
  max.overlaps = 20,
  nudge_y = 10,
  color_vector = c('LO-CRC' = plotting$condition_colors$`LO-CRC`, 'LO-CTR' = plotting$condition_colors$`LO-CTR`, 'n.s.' = 'white')) +
  scale_x_continuous(limits = c(-0.5, 0.9), breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  xlab('Late-onset enrichment effect size')

ggsave(volcano_lo, file = here('figures', 'extended.data.figure2', 'Extended.Data.Figure2a.pdf'), height = 5, width = 5)

######################
# Extended data figure 2b

volcano_eo <- plot_volcano(
  plot_df = lmm.table.eo %>%
    select(Bacteria, P.val, P.adj, Effect.size, pr.shift, `pr.EO-CRC`, `pr.EO-CTR`, `n.EO-CTR`, `n.EO-CRC`),
  fdr_thresh = 0.05,
  group_case = 'EO-CRC',
  group_control = 'EO-CTR',
  feature_column_name = 'Bacteria',
  max.overlaps = 20,
  color_vector = c('EO-CRC' = plotting$condition_colors$`EO-CRC`, 'EO-CTR' = plotting$condition_colors$`EO-CTR`, 'n.s.' = 'white')) +
  scale_x_continuous(limits = c(-0.5, 0.9), breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  xlab('Early-onset enrichment effect size')

ggsave(volcano_eo, file = here('figures', 'extended.data.figure2', 'Extended.Data.Figure2b.pdf'), height = 5, width = 5)

######################
# Extended data figure 2c

# Remove columns where all values are NA (dataset completeness)
lmm.table.eo.lo.ctr <- lmm.table.eo.lo.ctr %>%
  select(where(~ !all(is.na(.))))

volcano_eolo_ctr <- plot_volcano(
  plot_df = lmm.table.eo.lo.ctr %>% 
    select(Bacteria, P.val, P.adj, Effect.size, pr.shift, 'pr.LO-CTR', 'pr.EO-CTR', 'n.LO-CTR', 'n.EO-CTR'),
  group_case = 'LO-CTR',
  group_control = 'EO-CTR',
  feature_column_name = 'Bacteria',
  max.overlaps = 10,
  min_segment_length = 0.05, 
  nudge_y = 0.1,
  color_vector = c('LO-CTR' = plotting$condition_colors$`LO-CTR`, 'EO-CTR' = plotting$condition_colors$`EO-CTR`, 'n.s.' = 'white')) +
  xlab('Enrichment effect size') +
  theme(
    legend.box = 'horizontal',
    legend.spacing.y = unit(0.1, 'cm'),
    legend.spacing.x = unit(0.03, 'cm'),
    legend.position = c(0.02, 0.99),
    legend.key.size = unit(0.75, 'lines'),
    legend.justification = c(0, 1)
  ) +
  scale_x_continuous(limits = c(-0.25, 0.25), breaks = c(-0.25, -0.1, 0, 0.1, 0.25))

ggsave(volcano_eolo_ctr, file = here('figures', 'extended.data.figure2', 'Extended.Data.Figure2c.pdf'), height = 5, width = 5)

######################
# Extended data figure 2e

lmm.table.eo.lo.crc <- lmm.table.eo.lo.crc %>%
  select(where(~ !all(is.na(.))))

volcano_eolo_crc <- plot_volcano(
  plot_df = lmm.table.eo.lo.crc %>% 
    select(Bacteria, P.val, P.adj, Effect.size, pr.shift, 'pr.LO-CRC', 'pr.EO-CRC', 'n.LO-CRC', 'n.EO-CRC'),
  group_case = 'EO-CRC',
  group_control = 'LO-CRC',
  feature_column_name = 'Bacteria',
  max.overlaps = 10,
  min_segment_length = 0.05,
  nudge_y = 2,
  color_vector = c('EO-CRC' = plotting$condition_colors$`EO-CRC`, 'LO-CRC' = plotting$condition_colors$`LO-CRC`, 'n.s.' = 'white')) +
  xlab('Enrichment effect size') +
  theme(
    legend.box = 'horizontal',
    legend.spacing.y = unit(0.1, 'cm'),
    legend.spacing.x = unit(0.03, 'cm'),
    legend.position = c(0.44, 0.99),
    legend.key.size = unit(0.75, 'lines'),
    legend.justification = c(0, 1)
  ) +
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, -0.25, 0, 0.25, 0.5))

ggsave(volcano_eolo_crc, file = here('figures', 'extended.data.figure2', 'Extended.Data.Figure2e.pdf'), height = 5, width = 5)


######################
# Extended data figure 2f

load(here('data', 'results', 'Training.eo.lo.crc.rf.models.updated.Rdata'))

models <- list(models.lo.rf, siamcat.test.evaluated.lo.holdout.rf)
labels <- c('Classifier cross validated on LO-CRC', 'Classifier trained on EO-CRC and tested on LO-CRC')
trained_on <- list(NULL, models.eo.rf)
colours <- c('darkred', 'black')

lo_models_auc_plot <- plot_roc_siamcat_models(models, labels, colours, trained_on, alpha = 0.8)

ggsave(lo_models_auc_plot, file = here('figures', 'extended.data.figure2', 'Extended.Data.Figure2f.pdf'), height = 7, width = 7)

######################
# Extended data figure 2g

# Load LMM tables for EO with different cutoffs
load(here('data','results', 'lmm.tables.eo.diff.cutoff.Rdata')) 

alpha_fdr <- 0.05
library(GGally)

# Build a comparison table EO vs LO for each cutoff
compare_list <- imap(lmm.eo.list, function(tab_eo, cutoff) {
  left_join(
    lmm.table.lo %>%
      select(Bacteria, Effect.size, P.adj),
    tab_eo %>%
      select(Bacteria, Effect.size, P.adj),
    by = "Bacteria",
    suffix = c(".LO", ".EO")
  ) %>%
    mutate(
      cutoff = as.numeric(cutoff)) %>% 
    mutate(
      enriched_in = case_when(
        # Both tests show significant enrichment
        (P.adj.LO < 0.05 & Effect.size.LO > 0) & 
          (P.adj.EO <0.05 & Effect.size.EO > 0) ~ 'both enriched',
        
        # Both tests show significant depletion
        (P.adj.LO < 0.05 & Effect.size.LO < 0) & 
          (P.adj.EO <0.05 & Effect.size.EO < 0) ~ 'both depleted',
        
        # Both tests non-significant
        P.adj.EO >= 0.05 & P.adj.LO > 0.05 ~ 'both n.s.',
        
        # disagreements or partial significance
        TRUE ~ "one sig"
      )
    )
})

compare_all <- bind_rows(compare_list)

target_taxa <- c(
  "Fusobacterium","Peptostreptococcus","Parvimonas","Porphyromonas",
  "Gemella","Campylobacter_A","Hungatella"
)

taxa_pair <- c(
  "Peptostreptococcus" = "#86DF9B",
  "Fusobacterium"      = "#FFB9C7",
  "Parvimonas"         = "#75D9FF",
  "Porphyromonas"      = "#CCD276",
  "Gemella"            = "#FAC18D",
  "Campylobacter_A"    = "#53E2D2",
  "Hungatella"         = "#CDC5FF",
  "Other"              = "grey66"   # neutral for all non-listed taxa
)

lo_effects <- compare_all %>%
  select(Bacteria, Effect.size.LO) %>%
  distinct() %>%
  rename(LO = Effect.size.LO)

eo_effects_wide <- compare_all %>%
  select(Bacteria, cutoff, Effect.size.EO) %>%
  pivot_wider(
    names_from  = cutoff,
    values_from = Effect.size.EO,
    names_glue  = "EO (< {cutoff})"
  )

mat_lo_eo <- lo_effects %>%
  inner_join(eo_effects_wide, by = "Bacteria") %>%
  as.data.frame() %>% 
  select(c(Bacteria, LO , `EO (< 40)` ,`EO (< 45)`,  `EO (< 50)`,`EO (< 55)`))

plot_df <- mat_lo_eo %>%
  mutate(taxa = if_else(Bacteria %in% target_taxa, Bacteria, "Other")) %>%
  select(taxa, LO, `EO (< 40)`, `EO (< 45)`, `EO (< 50)`, `EO (< 55)`)

plot_df$taxa <- factor(plot_df$taxa, levels = c(names(taxa_pair)[names(taxa_pair)!="Other"], "Other"))

df <- mat_lo_eo %>% select(-Bacteria) 
vars <- names(df)

cor_df <- compare_all %>%
  group_by(cutoff) %>%
  summarise(
    corr = cor(Effect.size.LO, Effect.size.EO,
               method = "pearson", use = "complete.obs"),
    r2   = corr^2
  )

jac_df <- compare_all %>%
  group_by(cutoff) %>%
  summarise(
    numerator   = sum(enriched_in %in% c("both enriched", "both depleted")),
    denominator = n() - sum(enriched_in == "both n.s."),
    jaccard     = numerator / denominator
  )

df_label <- cor_df %>%
  inner_join(jac_df, by = "cutoff")

df_label

safe_cor <- function(data, mapping, ..., method = "pearson") {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  r <- if (length(x) >= 2) suppressWarnings(cor(x, y, use = "complete.obs", method = method)) else NA_real_
  lab <- if (is.finite(r)) sprintf("Corr : %.2f", r) else "r = NA"
  
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = lab, size = 4) +
    theme_void()
}

# Make sure factor levels & palette are consistent
plot_df$taxa <- factor(plot_df$taxa, levels = c(setdiff(names(taxa_pair), "Other"), "Other"))
pal <- taxa_pair[levels(plot_df$taxa)]

# Main plot
p_pairs <- ggpairs(
  data    = plot_df,
  columns = 2:6,
  mapping = aes(fill = taxa, color = taxa),
  lower   = list(
    continuous = wrap(
      "points",
      shape = 21, size = 3, alpha = 0.5,
      colour = "black",
      position = position_jitter(width = 0.1, height = 0.1)
    )
  ),
  upper = list(continuous = safe_cor),     # custom correlation panel
  diag  = list(continuous = "densityDiag"),
  title = "Comparison of EO (different cut-offs) signatures to LO",
  legend = 1
) +
  scale_fill_manual(values = pal, breaks = levels(plot_df$taxa), name = "Genus") +
  scale_color_manual(values = pal, guide = "none") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

p_pairs

ggsave(p_pairs, file= here('figures','extended.data.figure2','Extended.Data.Figure2g.pdf'), width = 10, height = 10)

######################
# Extended data figure 2h
load(here("data", "results",'Testing.eo.crc.by.age.cutoffs.Rdata'))    

models = c(eval_eo_tested[['40']],eval_eo_tested[['45']], eval_eo_tested[['50']],eval_eo_tested[['55']])
labels <- c("Classifier trained on LO and tested on \n EO-CRC (< 40)",
            "EO-CRC (< 45)",
            "EO-CRC (< 50)",
            "EO-CRC (< 55)")

trained_on <- list(models.lo.rf,models.lo.rf, models.lo.rf,  models.lo.rf)
colours <- c("#FFD08A", "#FFA500", "black", "#CC7000")

cross_tech_diff_eo_cutoff <- plot_roc_siamcat_models(models, labels, colours, trained_on, alpha=0.6)
ggsave(cross_tech_diff_eo_cutoff, file=here('figures','extended.data.figure2','Extended.Data.Figure2h.pdf'), height=7, width=7)

######################
# Extended data figure 2i: Shannon diversity and richness across EO/LO samples

all.meta<- read.table(here('data','Metadata.all.samples.tsv'), sep = '\t', header = T)  %>%
  filter(Condition=='CRC'| Condition=='CTR') %>% as.data.frame() 

diversity_richness_df <- read_tsv(here('data','results','Diversity.richness.tsv')) %>%
  left_join(all.meta %>% select(Sample_ID, Group), by = "Sample_ID") %>%
  drop_na(Group) %>%
  mutate(Cohort = as.factor(Cohort),
         Group = as.factor(Group))

# Blocked Wilcoxon p-values for Genus richness
wilcox_by_condition <- diversity_richness_df %>%
  group_by(Condition) %>%
  summarise(
    p_value = tryCatch({
      wt <- coin::wilcox_test(
        `Genus richness` ~ Group | Cohort,
        data = cur_data(),
        distribution = "approximate"
      )
      as.numeric(pvalue(wt))
    }, error = function(e) NA_real_)
  )

blocked_pval_richness <- diversity_richness_df %>%
  group_by(Condition) %>%
  summarise(groups_present = paste(sort(unique(Condition)), collapse = " vs ")) %>%
  left_join(wilcox_by_condition, by = "Condition")

p_crc <- blocked_pval_richness %>% filter(Condition == "CRC") %>% pull(p_value)
p_ctr <- blocked_pval_richness %>% filter(Condition == "CTR") %>% pull(p_value)

# Filter cohorts that have all 4 groups
required_groups <- c("EO-CRC", "LO-CRC", "EO-CTR", "LO-CTR")

div_filt <- diversity_richness_df %>%
  group_by(Cohort) %>%
  filter(all(required_groups %in% unique(Group))) %>%
  ungroup()

# Build plot_df for Genus richness (keeps Sample_ID)
plot_df <- bind_rows(
  div_filt %>% mutate(Cohort_facet = "All cohorts"),
  div_filt %>% mutate(Cohort_facet = as.character(Cohort))
) %>%
  mutate(alpha_val = ifelse(Cohort_facet == 'All cohorts', 'dark', 'light'),
         Group = factor(Group, levels = c("EO-CRC", "LO-CRC", "EO-CTR", "LO-CTR")))

upper_lim <- max(plot_df$`Genus richness`, na.rm = TRUE) * 1.3

seg_all <- tibble::tribble(
  ~Cohort_facet,   ~x, ~xend, ~y,                     ~yend,
  "All cohorts",    1,    1,  upper_lim * 0.80,  upper_lim * 0.85,
  "All cohorts",    2,    2,  upper_lim * 0.80,  upper_lim * 0.85,
  "All cohorts",    1,    2,  upper_lim * 0.85,  upper_lim * 0.85,
  "All cohorts",    3,    3,  upper_lim * 0.80,  upper_lim * 0.85,
  "All cohorts",    4,    4,  upper_lim * 0.80,  upper_lim * 0.85,
  "All cohorts",    3,    4,  upper_lim * 0.85,  upper_lim * 0.85
)

lab_all <- tibble::tribble(
  ~Cohort_facet,   ~x,  ~y,                    ~label,
  "All cohorts",  1.5, upper_lim * 0.9, paste0("Blocked \n p = ", signif(p_crc, 3)),
  "All cohorts",  3.5, upper_lim * 0.9, paste0("Blocked \n p = ", signif(p_ctr, 3))
)

# Blocked Wilcoxon p-values for Shannon diversity
wilcox_by_condition_sd <- diversity_richness_df %>%
  group_by(Condition) %>%
  summarise(
    p_value = tryCatch({
      wt <- coin::wilcox_test(
        `Shannon diversity` ~ Group | Cohort,
        data = cur_data(),
        distribution = "approximate"
      )
      as.numeric(pvalue(wt))
    }, error = function(e) NA_real_)
  )

blocked_pval_shannon <- diversity_richness_df %>%
  group_by(Condition) %>%
  summarise(groups_present = paste(sort(unique(Group)), collapse = " vs ")) %>%
  left_join(wilcox_by_condition_sd, by = "Condition")

p_crc_sd <- blocked_pval_shannon %>% filter(Condition == "CRC") %>% pull(p_value)
p_ctr_sd <- blocked_pval_shannon %>% filter(Condition == "CTR") %>% pull(p_value)

# Build plot_df_sd for Shannon diversity
plot_df_sd <- bind_rows(
  div_filt %>% mutate(Cohort_facet = as.character(Cohort)),
  div_filt %>% mutate(Cohort_facet = "All cohorts")
) %>%
  mutate(alpha_val = ifelse(Cohort_facet == 'All cohorts', 'dark', 'light'),
         Group = factor(Group, levels = c("EO-CRC", "LO-CRC", "EO-CTR", "LO-CTR")))

upper_lim_sd <- max(plot_df_sd$`Shannon diversity`, na.rm = TRUE) * 1.3

seg_all_sd <- tibble::tribble(
  ~Cohort_facet,   ~x, ~xend, ~y,                     ~yend,
  "All cohorts",    1,    1,  upper_lim_sd * 0.70,  upper_lim_sd * 0.75,
  "All cohorts",    2,    2,  upper_lim_sd * 0.70,  upper_lim_sd * 0.75,
  "All cohorts",    1,    2,  upper_lim_sd * 0.75,  upper_lim_sd * 0.75,
  "All cohorts",    3,    3,  upper_lim_sd * 0.70,  upper_lim_sd * 0.75,
  "All cohorts",    4,    4,  upper_lim_sd * 0.70,  upper_lim_sd * 0.75,
  "All cohorts",    3,    4,  upper_lim_sd * 0.75,  upper_lim_sd * 0.75
)

lab_all_sd <- tibble::tribble(
  ~Cohort_facet,   ~x,  ~y,                     ~label,
  "All cohorts",  1.5, upper_lim_sd * 0.80, paste0("Blocked \n p = ", signif(p_crc_sd, 3)),
  "All cohorts",  3.5, upper_lim_sd * 0.80, paste0("Blocked \n p = ", signif(p_ctr_sd, 3))
)

# Prepare combined "All cohorts" data for both metrics
only_all_cohorts <- plot_df %>%
  filter(Cohort_facet == 'All cohorts') %>%
  left_join(
    plot_df_sd %>% filter(Cohort_facet == 'All cohorts') %>% select(Sample_ID, `Shannon diversity`),
    by = 'Sample_ID', suffix = c('', '')
  ) %>%
  select(-c(Cohort_facet, alpha_val, Cohort, Assay, Condition, Group)) %>%
  pivot_longer(-Sample_ID, names_to = 'Metric') %>%
  left_join(all.meta %>% select(Sample_ID, Group), by = 'Sample_ID') %>%
  mutate(Group = factor(Group, levels = c("EO-CRC", "LO-CRC", "EO-CTR", "LO-CTR")))

only_all_data_plot <- ggplot(
  only_all_cohorts,
  aes(x = Group, y = value)
) +
  geom_boxplot(aes(fill = Group), color = "black", width = 0.7) +
  scale_fill_manual(values = plotting$condition_colors) +
  geom_segment(
    data = seg_all %>% mutate(Metric = 'Genus richness') %>% bind_rows(seg_all_sd %>% mutate(Metric = 'Shannon diversity')),
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.2,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = lab_all %>% mutate(Metric = 'Genus richness') %>% bind_rows(lab_all_sd %>% mutate(Metric = 'Shannon diversity')),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 5
  ) +
  facet_wrap(~ Metric, ncol = 2, scales = 'free_y') +
  theme_paper +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(
  only_all_data_plot,
  file = here("figures", "extended.data.figure2", "Extended.Data.Figure2i.pdf"),
  height = 8, width = 8
)

