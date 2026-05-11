######################
# Figure 6: Functional meta-analysis of CRC and Virulence Factors (VFs)
######################

# Load dependencies and parameters
source(here('src','utils.R'))
params <- yaml::read_yaml(here('src', 'parameters.yml'))
plotting <- params$plotting

######################
# Figure 6a: CRC enrichments of VFs in WGS data

# Load linear mixed model results
load(here('data','results','lmm.tables.functional.profiles.Rdata'))

# Load VF annotation and category data
vf.meta <- read.table(here('data','functional_data','VFGs_all.14Feb25_VS_GMGC.selected.table.tsv'),
                      sep = '\t', header = TRUE, quote = "", fill = TRUE, comment.char = "")

vfdb.category <- readxl::read_xls(here('data','functional_data','VFs.xlsx'), skip = 1)

# Join VF annotations and clean VF names
lmm.table.vf <- lmm.table.vf %>%
  left_join(vf.meta %>% select(VF_old_names, VFgene_protein), by = c("VF" = "VFgene_protein"), multiple = "any") %>%
  left_join(vf.meta %>% select(VF_Name = VF, VFgene_protein), by = c("VF" = "VFgene_protein"), multiple = "any") %>%
  mutate(
    VF_Name = sub(".*\\|\\s*", "", VF_Name),                       
    VF_Name = gsub("\\s*\\([^)]*\\)", "", VF_Name)              
  ) %>%
  left_join(vfdb.category %>% select(VFcategory, Bacteria, VF_Name, VF_FullName),
            by = "VF_Name", multiple = "any")

# Filter significant VFs (P.adj < 0.05) with positive effect sizes
lmm.table.vf.plot <- lmm.table.vf %>%
  filter(Effect.size > 0, P.adj < 0.05) %>%
  arrange(P.adj)

# Fill missing VFcategory values manually
lmm.table.vf.plot <- lmm.table.vf.plot %>%
  mutate(VFcategory = case_when(
    VF %in% c('fadA', 'fap2', 'fplA', 'radD', '', 'bcfc', 'bcfb') ~ 'Adherence',
    VF %in% c('fn1079', 'fn1079 (fn-dps)', 'fss2') ~ 'Stress survival',
    VF == 'lytC22' ~ 'Immune modulation',
    VF %in% c('mend', 'pduc', 'pduD', 'baiA', 'baiK', 'baiH') ~ 'Nutritional/Metabolic factor',
    VF %in% c('icsP/sopA', 'aha_3493') ~ 'Exoenzyme',
    TRUE ~ VFcategory
  ))

# Create effect size dataframe for plotting
effect.size.df <- lmm.table.vf.plot %>%
  select(VF, contains("Effect"), `conf.int_2.5 %`, `conf.int_97.5 %`) %>%
  rename_with(~ str_remove(.x, "Effect_")) %>%
  pivot_longer(
    cols = -c(VF, `conf.int_2.5 %`, `conf.int_97.5 %`),
    names_to = "Dataset", values_to = "Effect.size"
  ) %>%
  mutate(Dataset = ifelse(Dataset == "Effect.size", "LME", Dataset),
         Dataset = str_remove(Dataset, "Effect_size_"))

# Order VFs based on LME effect size
VF_order <- effect.size.df %>%
  filter(Dataset == "LME") %>%
  arrange(Effect.size) %>%
  pull(VF)

effect.size.df <- effect.size.df %>%
  mutate(VF = factor(VF, levels = VF_order))

# Main effect size plot
effect.size <- ggplot() +
  geom_point(data = effect.size.df %>% filter(Dataset == "LME"),
             aes(x = Effect.size, y = VF), size = 2, color = "black") +
  geom_col(width = 0.8) +
  geom_jitter(data = effect.size.df %>% filter(Dataset != "LME"),
              aes(x = Effect.size, y = VF), size = 1, color = "grey") +
  geom_errorbar(data = effect.size.df %>% filter(Dataset == "LME"),
                aes(xmin = `conf.int_2.5 %`, xmax = `conf.int_97.5 %`, y = VF),
                color = "black", size = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  theme_paper +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    legend.box = "vertical"
  ) +
  scale_x_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75))

# VFDB2 availability plot
VFDB2_availability <- lmm.table.vf.plot %>%
  left_join(vf.meta %>% select(VF_old_names, pos), by = "VF_old_names", multiple = "any") %>%
  mutate(VFDB2_avail = ifelse(is.na(pos), "No", "Yes")) %>%
  select(VF, VFcategory, VFDB2_avail) %>%
  mutate(
    VF = factor(VF, levels = rev(unique(lmm.table.vf.plot$VF))),
    VFcategory = factor(VFcategory, levels = unique(VFcategory))
  )

vfdb2_plot <- ggplot(VFDB2_availability, aes(x = 1, y = VF, fill = VFDB2_avail)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = c("Yes" = "black", "No" = "white")) +
  scale_y_discrete(limits = rev(VFDB2_availability$VF)) +
  labs(x = "VFDB2 Availability", fill = "Available") +
  theme_paper +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_text(margin = margin(t = 15, b = 0, l = 90)),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

# Significance and VF labels
category_bar_data <- lmm.table.vf.plot %>%
  mutate(
    VF = factor(VF, levels = rev(unique(VF))),
    VFcategory = factor(VFcategory, levels = unique(VFcategory)),
    VF_label = paste0(VF, " (", format(P.adj, scientific = TRUE, digits = 2), ")")
  ) %>%
  select(VF, VFcategory, VF_label)

# Custom category colors
custom_colors <- c(
  "Adherence" = "#ebac23", 
  "Immune modulation" = "#008cf9", 
  "Others" = "#bdbdbd",  
  "Exoenzyme" = "#5954d6", 
  "Exotoxin" = "#006e00",
  "Stress survival" = "#878500",
  "Nutritional/Metabolic factor" = "#00a76c",
  "Effector delivery system" = "#ff9287",
  "Biofilm" = "#b80058",
  "Biofilm formation" = "#d163e6"
)

# VF category bar plot
category_bar_plot <- ggplot(category_bar_data, aes(x = 1, y = VF, fill = VFcategory)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = custom_colors) +
  scale_y_discrete(labels = rev(category_bar_data$VF_label)) +
  labs(x = "VF Category", fill = "Category") +
  theme_paper +
  theme(
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.position = "left",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

# Combine all three plots into one figure
final_plot <- cowplot::plot_grid(
  category_bar_plot, vfdb2_plot, effect.size,
  align = "h", axis = "tb", rel_widths = c(0.85, 0.1, 1), ncol = 3
)

# Save final figure
ggsave(final_plot, file = here('figure', 'figure6', 'Figure6a.pdf'), height = 14, width = 10)

######################
# Figure 6c: CRC Enrichment in Fusobacterium subspecies across continents

load(here('data','results','lmm.tables.16S.WGS.Rdata'))

# Extract country information from metadata
country= read.table(here('data','Metadata.all.samples.tsv'),sep='\t',header=T) %>% filter(Assay=='WGS') %>% select(Country, Cohort) %>% unique() 

all.meta<- read.table(here('data','Metadata.all.samples.tsv'), sep = '\t', header = T)  %>%
  filter(Condition=='CRC'| Condition=='CTR') %>% as.data.frame() 

all.data.motus <- read.table(here('data','Raw.counts.wgs.motus.all.samples.tsv'),
                             sep='\t', header = TRUE, check.names = FALSE)

all.meta.wgs <- all.meta %>% filter(Assay == 'WGS')

# keep only WGS samples present in metadata
all.data.motus <- all.data.motus[, colnames(all.data.motus) %in% all.meta.wgs$Sample_ID]

# convert to relative abundances
all.data.motus <- all.data.motus %>%
  rownames_to_column('motus') %>%
  filter(rowSums(select(., where(is.numeric))) != 0) %>%
  column_to_rownames('motus')

all.data.motus <- all.data.motus %>% mutate(across(everything(), ~ . / sum(.)))

all.meta.wgs <- all.meta.wgs %>%
  select(Sample_ID, Cohort, Country, Condition) %>%
  mutate(
    patient_id = Sample_ID,
    Continent = case_when(
      Country %in% c("AUS","GER","FRA","ITA") ~ "Europe",
      Country %in% c("CAN","USA")             ~ "North America",
      Country %in% c("CHN","JPN")             ~ "Asia",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Continent))


all.data.motus<- normalize_data(all.data.motus)

# Keep fuso-only mOTUs in the matrix
all.data.motus <- all.data.motus %>%
  rownames_to_column("motu") %>%
  dplyr::filter(str_detect(motu, "11278|14406|05384|05385|13245|14121|01000|00999|01001|01005|01002|01003|01004|04102|01837")) %>%
  column_to_rownames("motu")

fuso_levels <- c(
  "Fusobacterium mortiferum [ref_mOTU_v31_11278]",
  "Fusobacterium sp. [meta_mOTU_v31_14406]",
  "Fusobacterium ulcerans [ref_mOTU_v31_05384]",
  "Fusobacterium varium [ref_mOTU_v31_05385]",
  "Fusobacterium sp. [meta_mOTU_v31_13245]",
  "Fusobacterium sp. [meta_mOTU_v31_14121]",
  "Fusobacterium sp. oral taxon 370 [ref_mOTU_v31_01000]",
  "Fusobacterium periodonticum [ref_mOTU_v31_00999]",
  "Fusobacterium nucleatum subsp. animalis [ref_mOTU_v31_01001]",
  "Fusobacterium nucleatum [ref_mOTU_v31_01005]",
  "Fusobacterium nucleatum subsp. vincentii [ref_mOTU_v31_01002]",
  "Fusobacterium nucleatum subsp. nucleatum [ref_mOTU_v31_01003]",
  "Fusobacterium hwasookii/nucleatum [ref_mOTU_v31_01004]",
  "Fusobacterium equinum/gonidiaformans [ref_mOTU_v31_04102]",
  "Fusobacterium necrophorum [ref_mOTU_v31_01837]"
)

continents  <- c("North America","Europe","Asia")
continent_cols  <- c("Asia"="red","Europe"="#1f77b4","North America"="#2ca02c")
continent_order <- c("North America","Europe","Asia")

run_lmem <- function(data_df, meta_df,ref_group,column_name, feature_column_name = "Bacteria" ) {
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lme4)
  
  
  # Initialize result containers
  features <- p.val <- NULL
  effect.size <- list()
  conf.int <- list()
  pr.shift <- list()
  sample.counts <- list()
  
  # Check if `Cohort` column exists
  if (!"Cohort" %in% colnames(meta_df)) {
    stop("The `meta_df` dataset must contain a `Cohort` column. Please ensure your metadata is correctly formatted.")
  }
  
  # Check if the group column has more than two levels
  unique_groups <- unique(meta_df[[column_name]])
  if (length(unique_groups) > 2) {
    stop(paste(
      "The column", column_name, "has more than two levels:",
      paste(unique_groups, collapse = ", "),
      "\nPlease reduce the levels to two before running the function."
    ))
  }
  
  # Reshape and join metadata
  data <- data_df %>%
    tibble::rownames_to_column(feature_column_name) %>%
    pivot_longer(-all_of(feature_column_name), names_to = "Sample_ID", values_to = "value") %>%
    left_join(meta_df %>% select(all_of(column_name), Cohort, Sample_ID), by = "Sample_ID")
  
  Features <- rownames(data_df)
  
  # Iterate over each feature
  for (i in seq_along(Features)) {
    feature_name <- Features[i]
    tmp <- data %>% filter(.data[[feature_column_name]] == feature_name)
    
    # Skip if subset is empty or group has fewer than 2 levels
    if (nrow(tmp) == 0 || length(unique(tmp[[column_name]])) < 2) {
      next
    }
    
    # Set reference group
    tmp[[column_name]] <- as.factor(tmp[[column_name]])
    tmp[[column_name]] <- relevel(tmp[[column_name]], ref = ref_group)
    
    # Extract non-reference group
    non_ref_group <- setdiff(levels(tmp[[column_name]]), ref_group)
    
    # Count samples in each group
    n.group1 <- sum(tmp[[column_name]] == ref_group)
    n.group2 <- sum(tmp[[column_name]] == non_ref_group)
    sample.counts[[i]] <- tibble(n.group1 = n.group1, n.group2 = n.group2)
    
    # Per-cohort estimates
    datasetWiseEstimates <- list()
    uq <- unique(tmp$Cohort)
    for (dsIndex in seq_along(uq)) {
      dataset <- tmp %>% filter(Cohort == uq[dsIndex])
      if (length(unique(dataset[[column_name]])) < 2) {
        datasetWiseEstimates[[length(datasetWiseEstimates) + 1]] <- NA
      } else {
        # Build formula safely
        form <- as.formula(paste("value ~", column_name))
        datasetWiseEstimates[[length(datasetWiseEstimates) + 1]] <- coefficients(summary(lm(form, data = dataset)))[2, 1]
      }
    }
    names(datasetWiseEstimates) <- uq
    
    # General linear mixed-effects model
    lmem_form <- as.formula(paste("value ~", column_name, "+ (1 | Cohort)"))
    lmem <- lmer(lmem_form, data = tmp)
    features[i] <- unique(tmp[[feature_column_name]])
    p.val[i] <- coefficients(summary(lmem))[2, "Pr(>|t|)"]
    effect.size[[i]] <- c(coefficients(summary(lmem))[2, "Estimate"], unlist(datasetWiseEstimates))
    conf.int[[i]] <- tryCatch({
      confint(lmem, level = 0.95)[4, ]
    }, error = function(e) {
      c(NA, NA)  # Fallback if confidence intervals cannot be computed
    })
    
    # Prevalence shift
    x.pos <- tmp %>% filter(.data[[column_name]] == non_ref_group)
    x.neg <- tmp %>% filter(.data[[column_name]] == ref_group)
    pr.n <- mean(x.neg$value >= 1e-04)
    pr.p <- mean(x.pos$value >= 1e-04)
    pr.shift[[i]] <- c(pr.p - pr.n, pr.n, pr.p)
  }
  
  # Combine results into a table
  result_table <- tibble(
    !!feature_column_name := features,
    P.val = p.val,
    Effect.size = effect.size,
    conf.int = conf.int,
    pr.shift = pr.shift,
    sample.counts = sample.counts
  )
  
  # Expand nested results
  expanded_table <- result_table %>%
    unnest_wider(c(Effect.size, conf.int, pr.shift, sample.counts), names_sep = "_") %>%
    dplyr::rename(
      Effect.size = Effect.size_1,
      pr.shift = pr.shift_1,
      !!paste0("pr.", ref_group) := pr.shift_2,
      !!paste0("pr.", non_ref_group) := pr.shift_3,
      !!paste0("n.", ref_group) := sample.counts_n.group1,
      !!paste0("n.", non_ref_group) := sample.counts_n.group2
    ) %>%
    mutate(P.adj = p.adjust(P.val, method = "BH"))
  
  # Return the final table
  return(expanded_table)
}


# container
continent_results <- list()

for (CT in c("North America","Europe","Asia")) {
  meta_sub <- all.meta.wgs %>% filter(Continent == CT)
  keep_samples <- intersect(colnames(all.data.motus), meta_sub$Sample_ID)
  
  if (length(keep_samples) < 8) {
    message("Skipping ", CT, ": too few samples")
    next
  }
  
  data_sub<- all.data.motus[, keep_samples, drop = FALSE]
  
  res<- run_lmem(
    data_df = data_sub,
    meta_df  = meta_sub, ref_group='CTR',column_name='Condition')
  
  if (nrow(res) > 0) {
    res$Continent <- CT
    continent_results[[CT]] <- res
  }
}

# bind into one table
continent_results <- dplyr::bind_rows(continent_results)

write_tsv(continent_results %>% rename(Taxa= Bacteria), file=here('data','results', 'lmm.tables.fuso.geography.tsv') )

species_lookup <- lmm.table.motu %>%
  dplyr::distinct(Bacteria, species) %>%
  dplyr::mutate(species = stringr::str_trim(species))

# Prep continent results for plotting
fuso_continent_lmm2 <- continent_results %>%
  left_join(species_lookup, by = "Bacteria") %>%
  mutate(
    species   = coalesce(species, Bacteria) |> str_trim(),
    species   = factor(species, levels = rev(fuso_levels)),
    Continent = factor(Continent, levels = continent_order)
  ) %>%
  select(Bacteria, species, Continent, Effect.size,`conf.int_2.5 %`, `conf.int_97.5 %`) %>%
  rename(value = Effect.size) 


fuso_dataset_pts2 <- continent_results %>%
  # add species names
  left_join(species_lookup, by = "Bacteria") %>%
  mutate(species = coalesce(species, Bacteria) |> str_trim()) %>%
  filter(species %in% fuso_levels) %>%
  select(species, starts_with("Effect.size_")) %>%
  pivot_longer(
    cols = starts_with("Effect.size_"),
    names_to = "Cohort", values_to = "value"
  ) %>%
  mutate(Cohort = sub("^Effect.size_", "", Cohort)) %>%
  left_join(country, by = "Cohort", multiple = "any") %>%
  mutate(
    Continent = case_when(
      Country %in% c("AUS","GER","FRA","ITA") ~ "Europe",
      Country %in% c("CAN","USA")             ~ "North America",
      Country %in% c("CHN","JPN")             ~ "Asia",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Continent)) %>%
  mutate(
    species   = factor(species, levels = rev(fuso_levels)),
    Continent = factor(Continent, levels = continent_order)
  ) 

# keep only finite rows for plotting
fuso_dataset_pts2_plot <- fuso_dataset_pts2 %>%
  dplyr::filter(is.finite(value))

fuso_continent_lmm2_plot <- fuso_continent_lmm2 %>%
  dplyr::filter(is.finite(value), is.finite(`conf.int_2.5 %`), is.finite( `conf.int_97.5 %`))

lanes_per_species <- length(continent_order)   # 3
lane_step   <- 8                               # vertical distance between lanes
lane_margin <- 0.4                             # inner padding for the background band
lane_idx    <- setNames(seq_along(continent_order), continent_order)

mk_y <- function(df, lvl) {
  df %>%
    mutate(
      species    = factor(species, levels = lvl),
      species_id = as.integer(species),
      Continent  = factor(Continent, levels = continent_order),
      lane_id    = lane_idx[as.character(Continent)],
      species_center = ((species_id - 1) * lanes_per_species + (lanes_per_species + 1) / 2) * lane_step,
      y = species_center + (lane_id - (lanes_per_species + 1) / 2) * lane_step
    )
}

# Tidy continent-level table & rename CI columns -> conf.low/conf.high 
fuso_continent_lmm2 <- fuso_continent_lmm2 %>%
  rename(conf.low = `conf.int_2.5 %`, conf.high = `conf.int_97.5 %`) %>%
  mutate(
    species   = factor(species, levels = rev(fuso_levels)),
    Continent = factor(Continent, levels = continent_order)
  )

fuso_dataset_pts2 <- fuso_dataset_pts2 %>%
  mutate(
    species   = factor(species, levels = rev(fuso_levels)),
    Continent = factor(Continent, levels = continent_order)
  )

lvl <- rev(fuso_levels)
fuso_continent_lmm2_plot <- fuso_continent_lmm2 %>%
  filter(is.finite(value)) %>%
  mk_y(lvl) %>%
  filter(is.finite(conf.low), is.finite(conf.high))

fuso_dataset_pts2_plot <- fuso_dataset_pts2 %>%
  filter(is.finite(value)) %>%
  mk_y(lvl)

# Species-level alternating background (covers all lanes for each species)
stripe_half <- 1.5 * lane_step - lane_margin
bg_species <- tibble::tibble(
  species        = lvl,
  species_id     = seq_along(lvl),
  species_center = ((species_id - 1) * lanes_per_species + (lanes_per_species + 1) / 2) * lane_step,
  ymin           = species_center - stripe_half,
  ymax           = species_center + stripe_half,
  fillcol        = rep(c("white", "grey95"), length.out = length(lvl))
)

# Plot
forest_plot_fuso_datasets_plus_continent <- ggplot() +
  geom_rect(
    data = bg_species,
    aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, fill = fillcol),
    inherit.aes = FALSE, alpha = 1
  ) +
  scale_fill_identity() +
  ggnewscale::new_scale_fill() +
  geom_jitter(
    data = fuso_dataset_pts2_plot,
    aes(x = value, y = y, fill = Continent),
    size = 2.8, color = "black", alpha = 0.6, shape = 21,
    width = 0, height = 0
  ) +
  geom_errorbarh(
    data = fuso_continent_lmm2_plot,
    aes(y = y, xmin = conf.low, xmax = conf.high),
    height = 0, linewidth = 1.2, color = "black", alpha = 0.9
  ) +
  geom_point(
    data = fuso_continent_lmm2_plot,
    aes(x = value, y = y, fill = Continent),
    shape = 23, size = 4.5, stroke = 1.1, color = "black"
  ) +
  scale_fill_manual(values = continent_cols, name = "Continent") +
  scale_y_continuous(
    breaks = bg_species$species_center,
    labels = bg_species$species,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(oob = scales::oob_keep) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.7) +
  xlab("Effect Size") + ylab("Species") +
  theme_paper +
  theme(
    axis.title.y     = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    axis.text.y      = element_text(face = "bold"),
    panel.grid       = element_blank()
  )


ggsave(
  forest_plot_fuso_datasets_plus_continent,
  file = here('figure','figure6','Figure6c.pdf'),
  width = 15, height = 12
)

#  Per-continent LMM + subgroup Q-test via meta package 

library(meta)
cohort_counts <- all.meta.wgs %>%
  group_by(Continent) %>%
  summarise(
    n_cohorts = n_distinct(Cohort),
    n_samples = n(),
    .groups   = "drop"
  )

message("Cohorts per continent:")
print(cohort_counts)

# Which continents are "thin" (< 3 cohorts)
thin_continents <- cohort_counts %>%
  filter(n_cohorts < 3) %>%
  pull(Continent)

if (length(thin_continents) > 0) {
  message("following continents have < 3 cohorts and will be ",
          "treated with caution:\n  ", paste(thin_continents, collapse = ", "))
}


# This is how I generated 3 lmm for 3 continents
continent_results <- list()

for (CT in continent_order) {
  meta_sub     <- all.meta.wgs %>% filter(Continent == CT)
  keep_samples <- intersect(colnames(all.data.motus), meta_sub$Sample_ID)
  
  if (length(keep_samples) < 8) {
    message("Skip ", CT, ": too few samples")
    next
  }
  
  data_sub <- all.data.motus[, keep_samples, drop = FALSE]
  
  res <- run_lmem(
    data_df            = data_sub,
    meta_df            = meta_sub,
    ref_group          = "CTR",
    column_name        = "Condition"
  )
  
  if (!is.null(res) && nrow(res) > 0) {
    n_coh        <- cohort_counts$n_cohorts[cohort_counts$Continent == CT]
    res$Continent <- CT
    res$n_cohorts <- n_coh
    
    # Include the flag estimates from thin subgroups
    res$reliability <- ifelse(
      n_coh < 3,
      "low — only 2 cohorts, tau^2 unreliable, CI wide",
      ifelse(n_coh < 5, "moderate — 3-4 cohorts", "good")
    )
    
    continent_results[[CT]] <- res
  }
}

continent_results <- dplyr::bind_rows(continent_results) %>% 
  filter(Bacteria!='ref_mOTU_v31_01000')


motus <- unique(continent_results$Bacteria)

use_common_tau <- length(thin_continents) > 0
message("\ntau.common = ", use_common_tau,
        if (use_common_tau) " (forced because >= 1 continent has only 2 cohorts)" else "")

subgroup_results <- lapply(setNames(motus, motus), function(m) {
  
  df <- continent_results %>%
    filter(Bacteria == m) %>%
    filter(
      is.finite(Effect.size),
      is.finite(`conf.int_2.5 %`),
      is.finite(`conf.int_97.5 %`)
    ) %>%
    mutate(
      # Reconstruct SE from 95% CI
      SE = (`conf.int_97.5 %` - `conf.int_2.5 %`) / (2 * 1.96)
    ) %>%
    filter(is.finite(SE), SE > 0)
  
  # Need at least 2 continents with valid estimates
  if (nrow(df) < 2) return(NULL)
  
  tryCatch({
    m_sub <- metagen(
      TE         = Effect.size,
      seTE       = SE,
      data       = df,
      studlab    = Continent,
      subgroup   = Continent,
      sm         = "MD",
      common     = FALSE,
      random     = TRUE,
      tau.common = use_common_tau
    )
    
    print(m_sub)
    
    tibble(
      motu             = m,
      n_continents     = nrow(df),
      # Between-subgroup Q-test — this IS the subgroup heterogeneity test
      Q_between        = m_sub$Q.b.random,
      df_between       = nrow(df) - 1,
      p_subgroup       = m_sub$pval.Q.b.random,
      I2_overall       = round(m_sub$I2 * 100, 1),
      tau_common_used  = use_common_tau,
      thin_subgroup    = paste(thin_continents, collapse = "; "),
      # Per-continent effects for easy reading
      effect_NorthAm   = df$Effect.size[df$Continent == "North America"] %||% NA_real_,
      effect_Europe    = df$Effect.size[df$Continent == "Europe"]        %||% NA_real_,
      effect_Asia      = df$Effect.size[df$Continent == "Asia"]          %||% NA_real_,
      se_NorthAm       = df$SE[df$Continent == "North America"]          %||% NA_real_,
      se_Europe        = df$SE[df$Continent == "Europe"]                 %||% NA_real_,
      se_Asia          = df$SE[df$Continent == "Asia"]                   %||% NA_real_
    )
  }, error = function(e) {
    message("metagen failed for motu: ", m, " — ", conditionMessage(e))
    NULL
  })
})

# Helper for safely extracting single values (base R doesn't have %||% for length-0)
`%||%` <- function(x, y) if (length(x) == 0 || all(is.na(x))) y else x

subgroup_results <- bind_rows(subgroup_results) %>%
  mutate(
    p_subgroup_adj = p.adjust(p_subgroup, method = "BH"),
    sig = case_when(
      p_subgroup_adj < 0.001 ~ "***",
      p_subgroup_adj < 0.01  ~ "**",
      p_subgroup_adj < 0.05  ~ "*",
      TRUE                   ~ "ns"
    )
  ) %>%
  left_join(
    lmm.table.motu %>%
      distinct(Bacteria, species) %>%
      rename(motu = Bacteria) %>%
      mutate(species = str_trim(species)),
    by = "motu"
  ) %>%
  mutate(species = coalesce(species, motu)) %>%
  relocate(species, .after = motu)

write_tsv(
  subgroup_results,
  here("data", "results", "subgroup_meta_Q_test_continents.tsv")
)

# Summary

cat("\n── Subgroup Q-test results (are continent effects heterogeneous?) ──\n")
cat("   tau.common =", use_common_tau, "\n")
cat("   Thin subgroup(s):", paste(thin_continents, collapse = ", "), "\n\n")

