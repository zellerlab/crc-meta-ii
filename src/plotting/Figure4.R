######################
# Figure 4
######################

source(here('src','utils.R'))
params <- yaml::read_yaml(here("src", "parameters.yml"))
plotting <- params$plotting

######################
# Figure 4a: Volcano plot for CRC tissue vs adjacent normal tissue comparison

# Read in the LMM tables for tissue comparisons
lmm.tissue <- read_tsv(here('data','results','lmm.tables.tissue.tsv')) 

volcano_tissue <- plot_volcano(
  plot_df = lmm.tissue %>% select(Taxon, P.val, P.adj, Effect.size,`n.Adj. non-tumor`,`n.Primary tumor`,`pr.Adj. non-tumor`,`pr.Primary tumor`),
  fdr_thresh = 0.05, 
  group_case = 'Primary tumor', 
  group_control = 'Adj. non-tumor',
  min_segment_length=0.2,max.overlaps = 10,
  color_vector = c(plotting$condition_colors$CRC, plotting$condition_colors$CTR,'white'),
  feature_column_name = 'Taxon') +
  xlab('Enrichment effect size (CRC tissue)')

ggsave(volcano_tissue, file=here('figures','figure4','Figure4a.pdf'), height = 5, width = 5)

######################
# Figure 4b: Scatter plot comparing effect sizes in tissue and fecal samples

# Load Effect size CRC vs CTR in fecal samples        
load(here('data','results','lmm.table.crc.ctr.Rdata'))

# Combine LMM tables for scatter plot
lmm.table.combined<- full_join(lmm.tissue %>% select(Taxon, P.val, P.adj, Effect.size) ,
                               lmm.table.general.crc %>% select(Taxon, P.val, P.adj, Effect.size), 
                               by='Taxon', suffix = c('.tissue','.fecal')) 

# Function to create scatter plot comparing effect sizes in tissue and fecal samples
plot_comparison_scatter_tissue <- function(
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
        (sig_x & up_x) & (sig_y & up_y) ~ "CRC",
        (sig_x & dn_x) & (sig_y & dn_y) ~ "CTR",
        TRUE                             ~ "n.s."
      ),
      enriched_in = factor(enriched_in, levels = c("CRC", "CTR", "n.s.")),
      label = dplyr::if_else(enriched_in!="n.s.", .data[[feature_column_name]], ""),
      font  = dplyr::if_else(enriched_in!="n.s.", "bold.italic", "")
    ) %>%
    dplyr::select(-sig_x, -sig_y, -up_x, -dn_x, -up_y, -dn_y)
  
  axis_limits <- range(c(data[[x_col]], data[[y_col]]), na.rm = TRUE)
  padding <- diff(axis_limits) * 0.05   # 5% margin
  axis_limits <- c(axis_limits[1] - padding, axis_limits[2] + padding)
  
  ggplot(x, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(aes(fill = enriched_in), shape = 21, color = "black", alpha = 0.75, size = 3) +
    geom_hline(yintercept = 0, color = "grey", lty = "solid", lwd = 0.3) +
    geom_vline(xintercept = 0, color = "grey", lty = "solid", lwd = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "grey", lty = "solid", lwd = 0.3) +
    ggrepel::geom_text_repel(
      aes(label = label, fontface = font),
      color = "black", segment.color = "black",
      segment.size = 0.2, segment.ncp = 3, max.overlaps = 25,
      min.segment.length = 0.2, nudge_x = 0.05, nudge_y = 0.05,
      size = 3, seed = 123
    ) +
    scale_fill_manual(values = c("CRC" = "#C44600", "CTR" = "dodgerblue3", "n.s." = "white")) +
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

# Create scatter plot comparing effect sizes in tissue and fecal samples
scatter_plot<- plot_comparison_scatter_tissue(data = lmm.table.combined,x_col = 'Effect.size.tissue',y_col = 'Effect.size.fecal',p_col_x = "P.adj.tissue",
                            p_col_y = "P.adj.fecal", feature_column_name = 'Taxon',x_label = 'Effect size tissue',y_label = 'Effect size fecal')

ggsave(scatter_plot, file=here('figures','figure4','Figure4b.pdf'), height = 5, width = 5)

######################
# Figure 4c: Mean relative abundance of key CRC-associated genera in fecal and tissue samples

selected_genera <- factor(c( "Peptostreptococcus", "Fusobacterium" , "Parvimonas", 'Porphyromonas' , "Gemella","Hungatella","Campylobacter_A" ),
  levels=rev(c( "Fusobacterium", "Gemella" , 'Porphyromonas',"Peptostreptococcus",  "Campylobacter_A" ,"Parvimonas","Hungatella"))) 

# Load relative abundance data for tissue samples
count_rel_mat2 <- count_rel_mat[which(rowSums(count_rel_mat > 0) / ncol(count_rel_mat) > 0.1),] # 10 % prev threshold

meta_df2 <- meta_df %>% 
  filter(sample_type=='Adj. non-tumor' | sample_type=='Primary tumor')  %>% rename(Cohort='study_name') 

count_rel_mat2 <- count_rel_mat2[,meta_df2$Sample_ID]
count_rel_mat2 <- count_rel_mat2 %>% 
  as.data.frame() %>%
  rownames_to_column('Taxon') %>% 
  filter(Taxon!='unassigned') %>% 
  column_to_rownames('Taxon')

abun_df <- count_rel_mat2 %>% 
  as_tibble(rownames = "bac") %>%
  gather(-bac, key = "Sample_ID", value = "rel")  %>% 
  mutate(is_prev = rel > 1e-4) %>% 
  mutate(l10_rel = log10(rel + 1e-5)) %>%
  left_join(meta_df %>% dplyr::select(Sample_ID, sample_type)) %>% 
  group_by(bac, sample_type) %>% 
  summarise(
    log_mean_rel = mean(l10_rel, na.rm = TRUE),
    mean_rel = mean(rel, na.rm = TRUE),
    sd_rel = sd(rel, na.rm = TRUE),
    prev = sum(is_prev) / length(unique(Sample_ID))
  ) %>% 
  mutate(in_core = bac %in% selected_genera) %>% 
  ungroup() %>% 
  filter(bac != "unassigned") %>% 
  filter(in_core==TRUE)

# Load relative abundance data for fecal samples 
fecal.meta<- read.table(here('data','Metadata.all.samples.tsv'), sep = '\t', header = T)  %>%
  filter(Condition=='CRC'| Condition=='CTR') %>% as.data.frame() 

fecal.data <- read.table(here('data','Relab.all.samples.tsv'),sep='\t',check.names = F)  %>% 
  rownames_to_column('genus') %>% 
  filter(genus!='unassigned') %>% 
  column_to_rownames('genus') 

fecal.data<-fecal.data[,fecal.meta$Sample_ID]

abun_df_fecal<- fecal.data %>% 
  as_tibble(rownames = "bac") %>%
  gather(-bac, key = "Sample_ID", value = "rel")  %>% 
  mutate(is_prev = rel > 1e-4) %>% 
  mutate(l10_rel = log10(rel + 1e-5)) %>%
  left_join(fecal.meta %>% dplyr::select(Sample_ID, Cohort, Condition)) %>% 
  group_by(bac, Condition) %>% 
  summarise(
    log_mean_rel = mean(l10_rel, na.rm = TRUE),
    mean_rel = mean(rel, na.rm = TRUE),
    sd_rel = sd(rel, na.rm = TRUE),
    prev = sum(is_prev) / length(unique(Sample_ID))
  ) %>% 
  mutate(in_core = bac %in% selected_genera) %>% 
  ungroup() %>% 
  filter(bac != "unassigned") %>% 
  filter(in_core==TRUE) 

# Prepare group-wise means for the alluvial plot: 
# One row per group (fecal healthy, fecal CRC, Adj. non-tumor, Primary tumor) per taxon, with the mean relative abundance for that group and taxon.

abun_df_tissue <- abun_df %>%
  filter(sample_type %in% c("Adj. non-tumor", "Primary tumor")) %>%
  transmute(bac, group = sample_type, mean_rel)

# fecal means -> rename groups to match the figure
abun_df_fecal2 <- abun_df_fecal %>%
  transmute(
    bac,
    group = recode(Condition, CTR = "Fecal healthy", CRC = "Fecal CRC"),
    mean_rel
  )

means_all <- bind_rows(abun_df_fecal2, abun_df_tissue) %>%
  filter(bac %in% selected_genera)

group_order <- c("Fecal healthy", "Fecal CRC", "Adj. non-tumor", "Primary tumor")

n_fecal <- fecal.meta %>%
  mutate(group = recode(Condition, CTR = "Fecal healthy", CRC = "Fecal CRC")) %>%
  count(group, name = "N")

n_tissue <- meta_df %>%
  filter(sample_type %in% c("Adj. non-tumor", "Primary tumor")) %>%
  transmute(group = sample_type) %>%
  count(group, name = "N")

n_per_group <- bind_rows(n_fecal, n_tissue)

plot_df <- means_all %>%
  left_join(n_per_group, by = "group") %>%
  mutate(
    group_lbl = paste0(as.character(group), "\n(N=", N, ")")
  ) %>% 
  mutate(bac=factor(bac,levels=(c( "Fusobacterium" , "Gemella" , 'Porphyromonas',"Peptostreptococcus", "Campylobacter_A" ,"Parvimonas","Hungatella"))))

stack_height <- plot_df %>%
  group_by(group_lbl) %>%
  summarise(total = sum(mean_rel), .groups = "drop")

group_order <- c("Fecal healthy\n(N=3263)", "Fecal CRC\n(N=3516)", "Adj. non-tumor\n(N=449)", "Primary tumor\n(N=457)")

# Prepare data for alluvial plot
plot_df_alluvial <- plot_df %>%
  select(group_lbl, bac, mean_rel) %>%
  mutate(
    group_lbl = factor(as.character(group_lbl), levels = group_order),
  ) %>%
  group_by(group_lbl, bac) %>%
  summarise(mean_rel = sum(mean_rel, na.rm = TRUE), .groups = "drop") %>%
  tidyr::complete(group_lbl, bac, fill = list(mean_rel = 0)) %>%   # fill missing combos with 0
  arrange(bac, group_lbl)

# Total stack height per group for the black outline boxes
stack_height <- plot_df_alluvial %>%
  group_by(group_lbl) %>%
  summarise(total = sum(mean_rel), .groups = "drop")

base_cols <- c("#ED90A4", "#D8A06A" ,"#ABB150" ,"#62BE79", "#00C1B2", "#48B8DE", "#ACA2EC")
pale_cols <- lighten(base_cols, amount = 0.4)   # amount: 0 (no change) → 1 (white)

col_map <- setNames(pale_cols[seq_along(selected_genera)], selected_genera)

# Plotting
p_alluvial <- ggplot(
  plot_df_alluvial,
  aes(x = group_lbl, stratum = bac, alluvium = bac, y = mean_rel, fill = bac)) + 
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray")+
  geom_stratum() +
  {
    lab_df <- plot_df_alluvial %>%
      group_by(group_lbl) %>%
      arrange(desc(mean_rel), .by_group = TRUE) %>%
      mutate(ypos = cumsum(mean_rel) - mean_rel/2) %>%
      ungroup() %>%
      filter(grepl("^Primary tumor tissue", group_lbl),
             bac %in% c("Fusobacterium", "Streptococcus"))
    geom_text(data = lab_df, aes(y = ypos, label = bac),
              fontface = "italic", size = 4.3)
  } +
  scale_y_continuous( "Mean relative abundance", labels = scales::label_number(accuracy = 0.01)) +
  scale_x_discrete(NULL, expand = expansion(mult = c(0.06, 0.06))) +
  scale_fill_manual(values = col_map) +
  guides(fill = guide_legend(title = NULL, ncol = 1)) +
  theme_paper +
  theme(
    axis.text.x = element_text(lineheight = 1.1),
    legend.position = "right",
    legend.key.size = unit(6, "pt")
  )

main <- p_alluvial + theme_paper

ggsave(main, file=here('figures','figure4','Figure4c.pdf'), width = 10,height = 7)

######################
# Figure 4d-g: FC of top taxa across CRC stages and tumor locations in fecal and tissue samples
library(MuMIn)

taxa_levels  <- c("Fusobacterium", "Porphyromonas", "Peptostreptococcus", "Parvimonas")
stage_levels <- c("I", "II", "III", "IV")
loc_levels   <- c("proximal", "distal", "rectum")  
pseudo       <- 1e-5

# SIAMCAT association parameters
probs.fc <- seq(0.1, 0.9, 0.05)
log.n0   <- 1e-5

# Colors (lightened)
taxa_cols <- c("#ED90A4", "#62BE79", "#48B8DE", "#ABB150")
taxa_cols <- lighten(taxa_cols, amount = 0.4)

# Map stage / localization to numeric for regression
stage_to_num <- function(x) {
  dplyr::case_when(
    x == "0"   ~ 0,
    x == "I"   ~ 1,
    x == "II"  ~ 2,
    x == "III" ~ 3,
    x == "IV"  ~ 4,
    TRUE       ~ NA_real_
  )
}

loc_to_num <- function(x) {
  dplyr::case_when(
    x == "proximal" ~ 1,
    x == "distal"   ~ 2,
    x == "rectum"   ~ 3,
    TRUE            ~ NA_real_
  )
}

# Generic LMM trend per taxon: response ~ 1 + x_var + (1|group_var)
fit_lmm_trend <- function(df, response, group_var, x_var, taxon) {
  df <- df %>% tidyr::drop_na(.data[[x_var]])
  
  if (length(unique(df[[x_var]])) < 2) {
    return(tibble(
      taxa = taxon,
      n    = nrow(df),
      slope = NA_real_,
      intercept = NA_real_,
      R2_marginal_m0    = NA_real_,
      R2_conditional_m0 = NA_real_,
      R2_marginal_m1    = NA_real_,
      R2_conditional_m1 = NA_real_,
      LRT_chisq         = NA_real_,
      LRT_df            = NA_real_,
      LRT_p             = NA_real_
    ))
  }
  
  form0 <- as.formula(
    paste(response, "~ 1 + (1 |", group_var, ")")
  )
  form1 <- as.formula(
    paste(response, "~ 1 +", x_var, "+ (1 |", group_var, ")")
  )
  
  m0 <- lmer(form0, data = df, REML = FALSE)
  m1 <- lmer(form1, data = df, REML = FALSE)
  
  aa    <- anova(m0, m1)
  r2_m0 <- suppressWarnings(r.squaredGLMM(m0))
  r2_m1 <- suppressWarnings(r.squaredGLMM(m1))
  
  slope     <- unname(fixef(m1)[x_var])
  intercept <- unname(fixef(m1)["(Intercept)"])
  
  tibble(
    taxa = taxon,
    n    = nrow(df),
    slope = slope,
    intercept = intercept,
    R2_marginal_m0    = as.numeric(r2_m0[1]),
    R2_conditional_m0 = as.numeric(r2_m0[2]),
    R2_marginal_m1    = as.numeric(r2_m1[1]),
    R2_conditional_m1 = as.numeric(r2_m1[2]),
    LRT_chisq = aa$Chisq[2],
    LRT_df    = aa$`Chi Df`[2],
    LRT_p     = aa$`Pr(>Chisq)`[2]
  )
}

make_trend_labels <- function(trend_tbl) {
  trend_tbl %>%
    mutate(
      taxa      = factor(taxa, levels = taxa_levels, ordered = TRUE),
      slope_fmt = ifelse(is.na(slope), "NA", sprintf("%.3f", slope)),
      r2_fmt    = ifelse(is.na(R2_conditional_m1), "NA", sprintf("%.3f", R2_conditional_m1)),
      label     = paste0("slope = ", slope_fmt, "\nR\u00B2 = ", r2_fmt)
    )
}

# Read tissue metadata & counts ----
meta_raw <- read_tsv(here('data', 'tissue_data',"Amplicon_meta-analysis_meta_clean_2501007_filtered.tsv")) %>%
  mutate(
    Localization = anatomic_location_prox_dist,
    Localization = ifelse(
      anatomic_location_prox_dist == "distal" & anatomic_location_segment == "rectum",
      anatomic_location_segment,
      anatomic_location_prox_dist
    )
  )

data_raw <- read_tsv(here('data', 'tissue_data',"Amplicon_tissue_meta-analysis_data_251007.tsv")) %>%
  tibble::column_to_rownames("genus") %>%
  mutate(across(everything(), ~ .x / sum(.x)))   # relab per sample

# Studies with stage info
studies_tissue <- meta_raw %>%
  filter(!is.na(stage)) %>%
  pull(study_name) %>%
  unique()

# Patients with both primary & adjacent samples
patients_tissue <- meta_raw %>%
  filter(study_name %in% studies_tissue) %>%
  group_by(patient_id) %>%
  filter(n_distinct(sample_type) == 2) %>%
  pull(patient_id)

# Stage info from primary tumor per patient
primary_stage <- meta_raw %>%
  filter(study_name %in% studies_tissue) %>%
  mutate(
    sample_type = str_trim(sample_type),
    stage       = na_if(str_trim(stage), "")
  ) %>%
  filter(sample_type == "Primary tumor", !is.na(stage)) %>%
  group_by(patient_id) %>%
  summarise(stage_primary = first(stage), .groups = "drop")

# Cleaned tissue meta with filled stage in Adj. non-tumor
meta_tissue <- meta_raw %>%
  filter(study_name %in% studies_tissue,
         patient_id %in% patients_tissue) %>%
  select(Sample_ID, study_name, patient_id, disease_status,
         sample_type, stage, specimen_type) %>%
  mutate(
    sample_type = str_trim(sample_type),
    stage       = na_if(str_trim(stage), "")
  ) %>%
  left_join(primary_stage, by = "patient_id") %>%
  mutate(
    stage = if_else(
      sample_type == "Adj. non-tumor" & is.na(stage),
      stage_primary,
      stage
    )
  ) %>%
  select(-stage_primary)

# Subset taxa and samples
data_tissue <- data_raw[, meta_tissue$Sample_ID] %>%
  tibble::rownames_to_column("taxa") %>%
  filter(taxa %in% taxa_levels) %>%
  tibble::column_to_rownames("taxa")

# Long + paired
data_long_tissue <- data_tissue %>%
  tibble::rownames_to_column("taxa") %>%
  pivot_longer(-taxa, names_to = "Sample_ID", values_to = "abund") %>%
  left_join(meta_tissue, by = "Sample_ID") %>%
  group_by(patient_id, taxa) %>%
  filter(n_distinct(sample_type) == 2) %>%
  ungroup()

# TISSUE: FC (tumor tissue vs adjacent non-tumor) across CRC stages
paired_fc <- data_long_tissue %>%
  select(taxa, study_name, patient_id, abund, sample_type, stage) %>%
  filter(taxa %in% taxa_levels,
         !is.na(stage),
         stage != "0") %>%
  pivot_wider(names_from = sample_type, values_from = abund) %>%
  mutate(
    fold_change = (`Primary tumor` + pseudo) / (`Adj. non-tumor` + pseudo),
    log10_fc    = log10(fold_change)
  )

# Per-study mean log10 FC per stage & taxa
summary_taxa_stage_study <- paired_fc %>%
  group_by(taxa, study_name, stage) %>%
  summarise(
    n_pairs         = sum(!is.na(`Primary tumor`) & !is.na(`Adj. non-tumor`)),
    mean_log10_fc   = mean(log10_fc, na.rm = TRUE),
    median_log10_fc = median(log10_fc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    stage     = factor(stage, levels = stage_levels, ordered = TRUE),
    stage_num = stage_to_num(as.character(stage)),
    taxa      = factor(taxa, levels = taxa_levels, ordered = TRUE)
  )

# Weighted mean per stage (across studies)
wm_stage_taxa <- summary_taxa_stage_study %>%
  group_by(taxa, stage, stage_num) %>%
  summarise(
    wmean_log10_fc = weighted.mean(mean_log10_fc, w = n_pairs, na.rm = TRUE),
    .groups = "drop"
  )

# LMM slopes per taxon (tissue, stage)
fc_trend_table_tissue <- purrr::map_dfr(
  taxa_levels,
  ~ fit_lmm_trend(
    df        = summary_taxa_stage_study %>% filter(taxa == .x),
    response  = "mean_log10_fc",
    group_var = "study_name",
    x_var     = "stage_num",
    taxon     = .x
  )
) %>%
  mutate(taxa = factor(taxa, levels = taxa_levels, ordered = TRUE))

# Regression segments for plotting
fc_trend_segments_tissue <- fc_trend_table_tissue %>%
  mutate(
    taxa    = factor(taxa, levels = taxa_levels, ordered = TRUE),
    x_start = 1,
    x_end   = 4,
    y_start = intercept + slope * x_start,
    y_end   = intercept + slope * x_end
  )

# Labels (slope + R2)
tissue_labels <- make_trend_labels(fc_trend_table_tissue)

# Tissue stage plot
p_study_tissue <- ggplot(summary_taxa_stage_study,
                          aes(x = stage_num, y = mean_log10_fc)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_jitter(aes(fill = taxa),
              color = "black", size = 2.7, alpha = 0.4,
              shape = 21, stroke = 0.8, width = 0.1) +
  geom_segment(
    data = fc_trend_segments_tissue,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
    color = "black", linewidth = 0.7, inherit.aes = FALSE
  ) +
  geom_point(
    data = wm_stage_taxa,
    aes(x = stage_num, y = wmean_log10_fc, fill = taxa),
    inherit.aes = FALSE,
    shape = 23, size = 3, color = "black", stroke = 1
  ) +
  geom_text(
    data = tissue_labels,
    aes(label = label),
    x = 1.05, y = 2,
    hjust = 0, vjust = 1,
    size = 3.5,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ taxa, ncol = 4, scales = "fixed", drop = FALSE) +
  scale_fill_manual(values = taxa_cols) +
  scale_x_continuous(
    breaks = 0:4,
    labels = c("0", "I", "II", "III", "IV")
  ) +
  labs(
    x = "Stage",
    y = "Mean log10 fold change per study (Tumor vs Adj.NT)",
    shape = "Study",
    title = "Tissue"
  ) +
  ylim(-1.5, 2) +
  theme(
    axis.text.x           = element_text(angle = 45, hjust = 1),
    panel.grid.major.x    = element_blank(),
    strip.text            = element_text(face = "bold"),
    axis.title.x          = element_text(margin = margin(t = 8)),
    axis.title.y          = element_text(margin = margin(r = 8)),
    legend.key.height     = unit(0.9, "lines"),
    legend.key.width      = unit(1.2, "lines")
  ) +
  theme_presentation()

# FECAL: gFC (CRC vs CTR) across stages in fecal samples

# Load Fecal metadata and relative abundance data, subset to taxa of interest and samples with stage info

fecal_meta <- read_tsv(here('data',"Metadata.all.samples.260304.tsv")) %>%
  filter(Condition %in% c("CRC", "CTR")) %>%
  mutate(
    Stage = if_else(Condition == "CTR", NA_character_, as.character(Stage)),
    Stage = case_when(
      Stage %in% c("0")        ~ "0",
      Stage %in% c("1", "I")   ~ "I",
      Stage %in% c("2", "II")  ~ "II",
      Stage %in% c("3", "III") ~ "III",
      Stage %in% c("4", "IV")  ~ "IV",
      TRUE ~ Stage
    )
  )

# Cohorts with stage info
studies_fecal <- fecal_meta %>%
  filter(!is.na(Stage)) %>%
  pull(Cohort) %>%
  unique()

fecal_meta_stage <- fecal_meta %>%
  filter(Cohort %in% studies_fecal) %>%
  filter(!(Condition == "CRC" & is.na(Stage))) %>%
  select(Sample_ID, Cohort, Condition, Stage) %>%
  filter(!(Condition == "CRC" & Stage == "0")) %>%
  as.data.frame()


tbl_counts_fecal <- fecal_meta_stage %>%
  group_by(Cohort, Stage) %>%
  summarise(n = n(), .groups = "drop") %>%
  drop_na(Stage) %>%
  tidyr::pivot_wider(names_from = Stage, values_from = n)

tbl_counts_fecal

#  Fecal data (4 taxa)
fecal_data <- read.table(here('data',"Relab.all.samples.tsv"),
                         sep = "\t", check.names = FALSE, header = TRUE) %>%
  tibble::rownames_to_column("taxa") %>%
  filter(taxa %in% taxa_levels) %>%
  tibble::column_to_rownames("taxa")

fecal_data <- fecal_data[, intersect(colnames(fecal_data), fecal_meta_stage$Sample_ID), drop = FALSE]

#  SIAMCAT helper: one Cohort × Stage 
run_assoc_stage_cohort <- function(cohort, stage) {
  meta_sub <- fecal_meta_stage %>%
    filter(
      Cohort == cohort,
      Condition %in% c("CRC", "CTR"),
      (Condition == "CTR") | (Condition == "CRC" & Stage == stage)
    )
  
  n_pos <- sum(meta_sub$Condition == "CRC")
  n_neg <- sum(meta_sub$Condition == "CTR")
  
  if (n_pos == 0 || n_neg == 0) {
    return(tibble(
      taxa   = rownames(fecal_data),
      Cohort = cohort,
      Stage  = stage,
      n_pos  = n_pos,
      n_neg  = n_neg,
      gFC    = NA_real_,
      p.val  = NA_real_,
      p.adj  = NA_real_
    ))
  }
  
  sel  <- meta_sub$Sample_ID
  feat <- fecal_data[, sel, drop = FALSE]
  
  lab <- create.label(
    meta    = meta_sub %>% tibble::column_to_rownames("Sample_ID"),
    label   = "Condition",
    case    = "CRC",
    control = "CTR"
  )
  
  sc <- siamcat(
    feat = feat,
    meta = meta_sub %>% tibble::column_to_rownames("Sample_ID"),
    label = lab
  )
  
  sc <- check.associations(
    sc,
    test         = "wilcoxon",
    feature.type = "original",
    probs.fc     = probs.fc,
    log.n0       = log.n0,
    verbose      = 0
  )
  
  associations(sc) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("taxa") %>%
    transmute(
      taxa,
      Cohort = cohort,
      Stage  = stage,
      n_pos  = n_pos,
      n_neg  = n_neg,
      gFC    = fc,
      p.val,
      p.adj
    )
}

#  Run SIAMCAT across Cohort × Stage 
cohorts <- sort(unique(fecal_meta_stage$Cohort))

gfc_stage_cohort_tbl <- purrr::map_dfr(
  cohorts,
  \(coh) purrr::map_dfr(stage_levels, \(stg) run_assoc_stage_cohort(coh, stg))
) %>%
  arrange(taxa, Cohort, Stage) %>%
  mutate(
    taxa      = factor(taxa, levels = taxa_levels, ordered = TRUE),
    Stage     = factor(Stage, levels = stage_levels, ordered = TRUE),
    stage_num = stage_to_num(as.character(Stage))
  )

# Weighted average across cohorts per stage
gfc_stage_weighted <- gfc_stage_cohort_tbl %>%
  group_by(taxa, Stage) %>%
  summarise(
    total_n_pos       = sum(n_pos, na.rm = TRUE),
    total_n_neg       = sum(n_neg, na.rm = TRUE),
    weighted_mean_gFC = weighted.mean(gFC, w = n_pos, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    taxa      = factor(taxa, levels = taxa_levels, ordered = TRUE),
    Stage     = factor(Stage, levels = stage_levels, ordered = TRUE),
    stage_num = stage_to_num(as.character(Stage))
  )

# LMM slopes per taxon (fecal, stage)
fc_trend_table_fecal <- purrr::map_dfr(
  taxa_levels,
  ~ fit_lmm_trend(
    df        = gfc_stage_cohort_tbl %>% filter(taxa == .x),
    response  = "gFC",
    group_var = "Cohort",
    x_var     = "stage_num",
    taxon     = .x
  )
) %>%
  mutate(taxa = factor(taxa, levels = taxa_levels, ordered = TRUE))

fc_trend_segments_fecal <- fc_trend_table_fecal %>%
  mutate(
    taxa    = factor(taxa, levels = taxa_levels, ordered = TRUE),
    x_start = 1,
    x_end   = 4,
    y_start = intercept + slope * x_start,
    y_end   = intercept + slope * x_end
  )

fecal_labels <- make_trend_labels(fc_trend_table_fecal)

# Fecal stage plot
p_study_fecal <- ggplot(gfc_stage_cohort_tbl,
                         aes(x = stage_num, y = gFC)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_segment(
    data = fc_trend_segments_fecal,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
    color = "black", linewidth = 0.7, inherit.aes = FALSE
  ) +
  geom_jitter(aes(fill = taxa),
              color = "black", size = 2.7, alpha = 0.4,
              shape = 21, stroke = 0.8, width = 0.1) +
  geom_point(
    data = gfc_stage_weighted,
    aes(x = stage_num, y = weighted_mean_gFC, fill = taxa),
    inherit.aes = FALSE,
    shape = 23, size = 3, color = "black", stroke = 1
  ) +
  geom_text(
    data = fecal_labels,
    aes(label = label),
    x = 1.05, y = 2,
    hjust = 0, vjust = 1,
    size = 3.5,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ taxa, ncol = 4, scales = "fixed", drop = FALSE) +
  scale_fill_manual(values = taxa_cols, limits = taxa_levels, breaks = taxa_levels) +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("I", "II", "III", "IV"),
    limits = c(1, 4)
  ) +
  labs(
    x = "Stage",
    y = "gFC (CRC vs CTR) per study",
    shape = "Study",
    title = "Fecal"
  ) +
  ylim(-1.5, 2) +
  theme(
    axis.text.x           = element_text(angle = 45, hjust = 1),
    panel.grid.major.x    = element_blank(),
    strip.text            = element_text(face = "bold"),
    axis.title.x          = element_text(margin = margin(t = 8)),
    axis.title.y          = element_text(margin = margin(r = 8)),
    legend.key.height     = unit(0.9, "lines"),
    legend.key.width      = unit(1.2, "lines")
  ) +
  theme_presentation()

# Combine stage plots & save

p_together <- p_study_fecal + p_study_tissue + plot_layout(guides = "collect")

ggsave(
  p_together,
  file  = here("figures","figure4","Figure4d-e.pdf"),
  width = 14,
  height = 6
)

# TISSUE: FC (tumor tissue vs adjacent non-tumor) across localization (proximal, distal, rectum)

# Fill localization per patient
primary_loc <- meta_raw %>%
  filter(study_name %in% studies_tissue) %>%
  mutate(sample_type = str_trim(sample_type)) %>%
  filter(sample_type == "Primary tumor", !is.na(Localization)) %>%
  group_by(patient_id) %>%
  summarise(Localization_primary = first(Localization), .groups = "drop")

meta_tissue_loc <- meta_raw %>%
  filter(study_name %in% studies_tissue,
         patient_id %in% patients_tissue) %>%
  select(Sample_ID, study_name, patient_id, disease_status,
         sample_type, stage, Localization, specimen_type) %>%
  mutate(sample_type = str_trim(sample_type)) %>%
  left_join(primary_loc, by = "patient_id") %>%
  mutate(
    Localization = if_else(
      sample_type == "Adj. non-tumor" & is.na(Localization),
      Localization_primary,
      Localization
    )
  ) %>%
  select(-Localization_primary)

# Long + paired (localization)
data_long_tissue_loc <- data_tissue %>%
  tibble::rownames_to_column("taxa") %>%
  pivot_longer(-taxa, names_to = "Sample_ID", values_to = "abund") %>%
  left_join(meta_tissue_loc, by = "Sample_ID") %>%
  group_by(patient_id, taxa) %>%
  filter(n_distinct(sample_type) == 2) %>%
  ungroup()

paired_fc_loc <- data_long_tissue_loc %>%
  select(taxa, study_name, patient_id, abund, sample_type, Localization) %>%
  filter(Localization %in% loc_levels, !is.na(Localization)) %>%
  pivot_wider(names_from = sample_type, values_from = abund) %>%
  mutate(
    fold_change = (`Primary tumor` + pseudo) / (`Adj. non-tumor` + pseudo),
    log10_fc    = log10(fold_change)
  )

summary_taxa_loc_study <- paired_fc_loc %>%
  group_by(taxa, study_name, Localization) %>%
  summarise(
    n_pairs         = sum(!is.na(`Primary tumor`) & !is.na(`Adj. non-tumor`)),
    mean_log10_fc   = mean(log10_fc, na.rm = TRUE),
    median_log10_fc = median(log10_fc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    taxa         = factor(taxa, levels = taxa_levels, ordered = TRUE),
    Localization = factor(Localization, levels = loc_levels, ordered = TRUE),
    loc_num      = loc_to_num(as.character(Localization))
  )

wm_loc_taxa <- summary_taxa_loc_study %>%
  group_by(taxa, Localization, loc_num) %>%
  summarise(
    wmean_log10_fc = weighted.mean(mean_log10_fc, w = n_pairs, na.rm = TRUE),
    .groups = "drop"
  )

fc_trend_loc_tissue <- purrr::map_dfr(
  taxa_levels,
  ~ fit_lmm_trend(
    df        = summary_taxa_loc_study %>% filter(taxa == .x),
    response  = "mean_log10_fc",
    group_var = "study_name",
    x_var     = "loc_num",
    taxon     = .x
  )
) %>%
  mutate(taxa = factor(taxa, levels = taxa_levels, ordered = TRUE))

fc_trend_segments_loc_tissue <- fc_trend_loc_tissue %>%
  mutate(
    taxa    = factor(taxa, levels = taxa_levels, ordered = TRUE),
    x_start = 1,
    x_end   = 3,
    y_start = intercept + slope * x_start,
    y_end   = intercept + slope * x_end
  )

tissue_loc_labels <- make_trend_labels(fc_trend_loc_tissue)

p_tissue_loc <- ggplot(summary_taxa_loc_study,
                       aes(x = loc_num, y = mean_log10_fc)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_jitter(aes(fill = taxa),
              color = "black", size = 2.7, alpha = 0.4,
              shape = 21, stroke = 0.8, width = 0.15) +
  geom_segment(
    data = fc_trend_segments_loc_tissue,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
    color = "black", linewidth = 0.7, inherit.aes = FALSE
  ) +
  geom_point(
    data = wm_loc_taxa,
    aes(x = loc_num, y = wmean_log10_fc, fill = taxa),
    inherit.aes = FALSE,
    shape = 23, size = 3, color = "black", stroke = 1
  ) +
  geom_text(
    data = tissue_loc_labels,
    aes(label = label),
    x = 1.05, y = 2,
    hjust = 0, vjust = 1,
    size = 3.5,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ taxa, ncol = 4, scales = "fixed", drop = FALSE) +
  scale_fill_manual(values = taxa_cols) +
  scale_x_continuous(
    breaks = 1:3,
    labels = loc_levels,
    limits = c(1, 3)
  ) +
  labs(
    x = "Localization",
    y = "Mean log10 fold change per study (Tumor vs Adj.NT)",
    shape = "Study",
    title = "Tissue"
  ) +
  ylim(-1, 2) +
  theme_presentation() +
  theme(
    axis.text.x        = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    strip.text         = element_text(face = "bold"),
    axis.title.x       = element_text(margin = margin(t = 8)),
    axis.title.y       = element_text(margin = margin(r = 8)),
    legend.key.height  = unit(0.9, "lines"),
    legend.key.width   = unit(1.2, "lines")
  )

# FECAL: gFC (CRC vs CTR) across localization (proximal, distal, rectum)
fecal_meta_loc <- read.table(here("data","Metadata.all.samples.tsv"),
                             sep = "\t", header = TRUE) %>%
  filter(Condition %in% c("CRC", "CTR")) %>%
  select(Sample_ID, Cohort, Condition, Localization)

# Cohorts with localization info
studies_fecal_loc <- fecal_meta_loc %>%
  filter(!is.na(Localization)) %>%
  pull(Cohort) %>%
  unique()

fecal_meta_loc <- fecal_meta_loc %>%
  filter(Cohort %in% studies_fecal_loc) %>%
  filter(!(Condition == "CRC" & is.na(Localization))) %>%
  mutate(
    Localization = case_when(
      Localization == "Left_Distal"    ~ "distal",
      Localization == "Right_Proximal" ~ "proximal",
      Localization == "Rectum"         ~ "rectum",
      TRUE                             ~ Localization
    ),
    Localization = if_else(
      Localization %in% loc_levels,
      Localization,
      NA_character_
    )
  ) 

# Sample counts per Localization
tbl_counts_fecal_loc <- fecal_meta_loc %>%
  group_by(Cohort, Localization) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Localization, values_from = n)


fecal_data_loc <- read.table(here("data","Relab.all.samples.tsv"),
                             sep = "\t", check.names = FALSE, header = TRUE) %>%
  tibble::rownames_to_column("taxa") %>%
  filter(taxa %in% taxa_levels) %>%
  tibble::column_to_rownames("taxa")

fecal_data_loc <- fecal_data_loc[, intersect(colnames(fecal_data_loc), fecal_meta_loc$Sample_ID), drop = FALSE]

# SIAMCAT helper: one Cohort × Localization 
run_assoc_loc_cohort <- function(cohort, loc) {
  meta_sub <- fecal_meta_loc %>%
    filter(
      Cohort == cohort,
      Condition %in% c("CRC", "CTR"),
      (Condition == "CTR") | (Condition == "CRC" & Localization == loc)
    )
  
  n_pos <- sum(meta_sub$Condition == "CRC")
  n_neg <- sum(meta_sub$Condition == "CTR")
  
  if (n_pos == 0 || n_neg == 0) {
    return(tibble(
      taxa         = rownames(fecal_data_loc),
      Cohort       = cohort,
      Localization = loc,
      n_pos        = n_pos,
      n_neg        = n_neg,
      gFC          = NA_real_,
      p.val        = NA_real_,
      p.adj        = NA_real_
    ))
  }
  
  sel  <- meta_sub$Sample_ID
  feat <- fecal_data_loc[, sel, drop = FALSE]
  
  lab <- create.label(
    meta    = meta_sub %>% tibble::column_to_rownames("Sample_ID"),
    label   = "Condition",
    case    = "CRC",
    control = "CTR"
  )
  
  sc <- siamcat(
    feat  = feat,
    meta  = meta_sub %>% tibble::column_to_rownames("Sample_ID"),
    label = lab
  )
  
  sc <- check.associations(
    sc,
    test         = "wilcoxon",
    feature.type = "original",
    probs.fc     = probs.fc,
    log.n0       = log.n0,
    verbose      = 0
  )
  
  associations(sc) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("taxa") %>%
    transmute(
      taxa,
      Cohort       = cohort,
      Localization = loc,
      n_pos        = n_pos,
      n_neg        = n_neg,
      gFC          = fc,
      p.val,
      p.adj
    )
}

cohorts_loc <- sort(unique(fecal_meta_loc$Cohort))

gfc_loc_cohort_tbl <- purrr::map_dfr(
  cohorts_loc,
  \(coh) purrr::map_dfr(loc_levels, \(loc) run_assoc_loc_cohort(coh, loc))
) %>%
  arrange(taxa, Cohort, Localization) %>%
  mutate(
    taxa         = factor(taxa, levels = taxa_levels, ordered = TRUE),
    Localization = factor(Localization, levels = loc_levels, ordered = TRUE),
    loc_num      = loc_to_num(as.character(Localization))
  )

gfc_loc_weighted <- gfc_loc_cohort_tbl %>%
  group_by(taxa, Localization) %>%
  summarise(
    total_n_pos       = sum(n_pos, na.rm = TRUE),
    total_n_neg       = sum(n_neg, na.rm = TRUE),
    weighted_mean_gFC = weighted.mean(gFC, w = n_pos, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    taxa         = factor(taxa, levels = taxa_levels, ordered = TRUE),
    Localization = factor(Localization, levels = loc_levels, ordered = TRUE),
    loc_num      = loc_to_num(as.character(Localization))
  )

fc_trend_loc_fecal <- purrr::map_dfr(
  taxa_levels,
  ~ fit_lmm_trend(
    df        = gfc_loc_cohort_tbl %>% filter(taxa == .x),
    response  = "gFC",
    group_var = "Cohort",
    x_var     = "loc_num",
    taxon     = .x
  )
) %>%
  mutate(taxa = factor(taxa, levels = taxa_levels, ordered = TRUE))

fc_trend_segments_loc_fecal <- fc_trend_loc_fecal %>%
  mutate(
    taxa    = factor(taxa, levels = taxa_levels, ordered = TRUE),
    x_start = 1,
    x_end   = 3,
    y_start = intercept + slope * x_start,
    y_end   = intercept + slope * x_end
  )

fecal_loc_labels <- make_trend_labels(fc_trend_loc_fecal)


p_fecal_loc <- ggplot(gfc_loc_cohort_tbl,
                      aes(x = loc_num, y = gFC)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_segment(
    data = fc_trend_segments_loc_fecal,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
    color = "black", linewidth = 0.7, inherit.aes = FALSE
  ) +
  geom_jitter(aes(fill = taxa),
              color = "black", size = 2.7, alpha = 0.4,
              shape = 21, stroke = 0.8, width = 0.15) +
  geom_point(
    data = gfc_loc_weighted,
    aes(x = loc_num, y = weighted_mean_gFC, fill = taxa),
    inherit.aes = FALSE,
    shape = 23, size = 3, color = "black", stroke = 1
  ) +
  geom_text(
    data = fecal_loc_labels,
    aes(label = label),
    x = 1.05, y = 2,
    hjust = 0, vjust = 1,
    size = 3.5,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ taxa, ncol = 4, scales = "fixed", drop = FALSE) +
  scale_fill_manual(values = taxa_cols, limits = taxa_levels, breaks = taxa_levels) +
  scale_x_continuous(
    breaks = 1:3,
    labels = loc_levels,
    limits = c(1, 3)
  ) +
  labs(
    x = "Localization",
    y = "gFC (CRC vs CTR) per study",
    shape = "Study",
    title = "Fecal"
  ) +
  ylim(-1, 2) +
  theme_presentation() +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    strip.text         = element_text(face = "bold"),
    axis.title.x       = element_text(margin = margin(t = 8)),
    axis.title.y       = element_text(margin = margin(r = 8)),
    legend.key.height  = unit(0.9, "lines"),
    legend.key.width   = unit(1.2, "lines")
  )

# Combine localization plots & save
p_tissue_loc <- p_tissue_loc + labs(title = "Tissue")
p_fecal_loc  <- p_fecal_loc  + labs(title = "Fecal")

p_loc_together <- p_tissue_loc + p_fecal_loc + plot_layout(guides = "collect")

ggsave(
  p_loc_together,
  file  = here('figures','figure4',"Figure4f-g.pdf"),
  width = 14,
  height = 6
)
