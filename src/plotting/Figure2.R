######################
# Figure 2
######################

library(here)

source(here('src','utils.R'))
params <- yaml::read_yaml(here("src", "parameters.yml"))
plotting <- params$plotting

# Load LMM tables for 16S and WGS
load(here('data','results', 'lmm.table.crc.ctr.Rdata'))
load(here('data','results', 'lmem.tables.eo.lo.Rdata'))

######################
# Figure 2a: Forest plot for EO/LO-CRC  vs EO/LO-CTR

# Import genus-level data and metadata

all.meta <- read.table(here('data','Metadata.all.samples.tsv'), sep = '\t', header = TRUE) %>%
  filter(Condition %in% c('CRC', 'CTR')) %>%
  as.data.frame()

all.data <- read.table(here('data','Relab.all.samples.tsv'), sep = '\t', check.names = FALSE) %>%
  rownames_to_column('genus') %>%
  filter(genus != 'unassigned') %>%
  column_to_rownames('genus')

# Filter for genera present in >10% of samples
all.data <- all.data[rowSums(all.data > 0) / ncol(all.data) > 0.1, ]

# EO vs LO comparisons within CRC and CTR

# Subset EO samples
all.meta.eo <- all.meta %>% filter(Group %in% c('EO-CRC', 'EO-CTR'))
all.data.eo <- all.data[, all.meta.eo$Sample_ID]

# Subset LO samples
all.meta.lo <- all.meta %>% filter(Group %in% c('LO-CRC', 'LO-CTR'))
all.data.lo <- all.data[, all.meta.lo$Sample_ID]

lmm.table.eo <- lmm.table.eo %>%  #some dataset does not have either EO-CRC or EO-CTR samples
  select(where(~ !all(is.na(.))))

# Prepare Early Onset data
lmm.table.eo <- lmm.table.eo %>%
  filter(P.adj < 0.05) %>%
  dplyr::rename(LME = Effect.size) %>%
  arrange(LME) %>%
  rename_with(~ gsub("Effect.size_", "", .), starts_with("Effect.size")) %>%
  mutate(Group = "EO") %>%
  rename_with(~ gsub("EO-", "", .x))

# Prepare Late Onset data
lmm.table.lo <- lmm.table.lo %>%
  filter(P.adj < 0.05) %>%
  dplyr::rename(LME = Effect.size) %>%
  arrange(LME) %>%
  rename_with(~ gsub("Effect.size_", "", .), starts_with("Effect.size")) %>%
  mutate(Group = "LO") %>%
  rename_with(~ gsub("LO-", "", .x))

# Merge datasets
datasets <- unique(c(all.meta.eo$Cohort, all.meta.lo$Cohort))

merged <- bind_rows(lmm.table.eo, lmm.table.lo)

# Order by absolute LME values
merged_ordered <- merged %>%
  arrange(desc((LME))) %>%
  filter(Group %in% c("EO", "LO"))%>%
  select(-c(pr.shift, pr.CRC, pr.CTR , P.adj))

# Select bacteria present in both groups
common_bacteria <- merged_ordered %>%
  group_by(Bacteria) %>%
  filter(n() == 2) %>%
  pivot_longer(c(LME, datasets)) %>%
  dplyr::rename(group = name, `model effect size` = value)

# Add size and alpha for points
common_bacteria <- common_bacteria %>%
  mutate(
    size = ifelse(group == "LME", 2, 1.5),
    alpha = ifelse(group == "LME", 1, 0.8)
  ) %>%
  dplyr::rename(Dataset = group)

# Remove NA values
common_bacteria <- common_bacteria %>%
  filter(!is.na(`model effect size`), !is.na(Bacteria))


# Select top 10 positive and top 10 negative bacteria based on LME values
top_positive_bacteria <- common_bacteria %>%
  filter(Dataset == "LME") %>%
  group_by(Bacteria) %>%
  summarize(max_LME = max(`model effect size`, na.rm = TRUE)) %>%
  arrange(desc(max_LME)) %>%
  slice(1:10) %>% # Top 10 positive
  pull(Bacteria)

top_negative_bacteria <- common_bacteria %>%
  filter(Dataset == "LME") %>%
  group_by(Bacteria) %>%
  summarize(min_LME = min(`model effect size`, na.rm = TRUE)) %>%
  arrange(min_LME) %>% # Top 10 negative
  slice(1:10) %>%
  pull(Bacteria)

# Top_negative_bacteria
# Combine the top positive and negative bacteria
top_bacteria <- c( top_negative_bacteria,rev(top_positive_bacteria))

# Filter the common_bacteria dataset to include only these bacteria
common_bacteria_filtered <- common_bacteria %>%
  filter(Bacteria %in% top_bacteria) %>%
  mutate(Bacteria = factor(Bacteria, levels=top_bacteria))

shapes <- data.frame(
  shape = c(0:6),
  x = 0:6 %/% 5,
  y = -(0:6 %% 5)
)

# Generate forest plot 
forest_plot <- ggplot() +
  geom_point(data = common_bacteria_filtered %>%
               filter(Dataset != "LME") %>%
               filter(Group == "EO"),
             aes(x = `model effect size`, y = Bacteria, color = Group), alpha = 0, size = 0, show.legend = FALSE) +
  geom_point(data = common_bacteria_filtered %>%
               filter(Dataset != "LME") %>%
               filter(Group == "EO"),
             aes(x = `model effect size`, y = as.numeric(Bacteria) + 0.2,
                 shape = Dataset, color = Group), size = 2, show.legend = TRUE) +
  geom_point(data = common_bacteria_filtered %>%
               filter(Dataset == "LME") %>%
               filter(Group == "EO") %>%
               arrange(desc(`model effect size`)),
             aes(x = `model effect size`, y = as.numeric(Bacteria) + 0.2,
                 shape = Dataset, color = Group), size = 4.5, show.legend = FALSE) +
  geom_point(data = common_bacteria_filtered %>%
               filter(Dataset != "LME") %>%
               filter(Group == "LO"),
             aes(x = `model effect size`, y = as.numeric(Bacteria) - 0.2,
                 shape = Dataset, color = Group), size = 2, show.legend = FALSE) +
  geom_point(data = common_bacteria_filtered %>%
               filter(Dataset == "LME") %>%
               filter(Group == "LO") %>%
               arrange(desc(`model effect size`)),
             aes(x = `model effect size`, y = as.numeric(Bacteria) - 0.2,
                 shape = Dataset, color = Group), size = 4.5, show.legend = FALSE) +
  geom_line(data = common_bacteria_filtered %>%
              filter(Dataset == "LME") %>%
              filter(Group == "EO") %>%
              select(lowCI = `conf.int_2.5 %`, highCI = `conf.int_97.5 %`, Bacteria, Group) %>%
              pivot_longer(-c(Bacteria, Group)),
            aes(x = value, y = as.numeric(Bacteria) + 0.2, group = interaction(Bacteria, Group)),
            color = 'black') +
  geom_line(data = common_bacteria_filtered %>%
              filter(Dataset == "LME") %>%
              filter(Group == "LO") %>%
              select(lowCI = `conf.int_2.5 %`, highCI = `conf.int_97.5 %`, Bacteria, Group) %>%
              pivot_longer(-c(Bacteria, Group)),
            aes(x = value, y = as.numeric(Bacteria) - 0.2, group = interaction(Bacteria, Group)),
            color = 'black') +
  scale_shape_manual(
    values = c(15, 0, 1, 3, 4, 7, 8, 9, 10, 12, 11,13,14,16,17,18,
               25,42,35,94,95,62,43),
    breaks = c("LME", datasets)) +
  scale_color_manual(values = c("EO" = "darkorange", "LO" = "darkred")) +
  scale_alpha_continuous(guide = 'none') +
  ylab('Genus') +
  xlab('Model Effect Size') +
  geom_vline(xintercept = 0, linetype = 'dotted', size = 1) +
  theme_paper +
  theme(
    legend.box = "vertical",
    legend.spacing.y = unit(0.1, "cm"),
    legend.position = c(.75, .55),
    legend.key.size = unit(0.75, "lines"),
    legend.justification = c(0, 1)) +
  guides(
    color = guide_legend(title = "Group", order=1,override.aes = list(size = 5,shape = 15)), # Ensure Group legend is included
    shape = guide_legend(override.aes = list(size = 5), ncol=1,title = "Datasets",order = 2))


ggsave(plot = forest_plot, file= here('figures','figure2','Figure2a.pdf') , width = 7.5,  height = 10)

######################
# Figure 2b: Scatter plot

load(here('data','results', 'lmm.tables.eo.lo.Rdata'))

lmm.table.all<- left_join(lmm.table.lo %>% select(Bacteria,P.val, P.adj, Effect.size, pr.shift,'pr.LO-CRC','pr.LO-CTR','n.LO-CTR', 'n.LO-CRC') ,
                           lmm.table.eo %>% select(Bacteria,P.val, P.adj, Effect.size, pr.shift,'pr.EO-CRC','pr.EO-CTR', 'n.EO-CTR', 'n.EO-CRC'),
                           by='Bacteria', suffix = c('.LO','.EO'))


corr <- cor(
  lmm.table.all$Effect.size.EO,
  lmm.table.all$Effect.size.LO,
  method = "pearson",
  use = "complete.obs"
)

# Format the label
corr_label <- paste0("Pearson r = ", round(corr, 2))

plot_comparison <- function(data, x_col, y_col, x_label, y_label) {
  # Prepare data with color coding based on criteria
  data <- data %>%
    mutate(
      point_color = case_when(
        Effect.size.LO > 0.1 & Effect.size.LO > 0.1 & P.adj.EO < 0.05 & P.adj.LO < 0.05 ~ "dodgerblue3",
        Effect.size.EO < -0.1 & Effect.size.LO < -0.1 & P.adj.EO < 0.05 & P.adj.LO < 0.05 ~ "dodgerblue3",
        TRUE ~ "white" # Remaining points are white
      ),
      enriched_in =
        case_when(
          P.adj.EO < 0.05 & Effect.size.EO > 0.1 |  P.adj.LO < 0.05 & Effect.size.LO > 0.1 ~ 'CRC',
          P.adj.EO < 0.05 & Effect.size.EO < -0.1 |  P.adj.LO < 0.05 & Effect.size.LO < -0.1  ~ 'CTR',
          TRUE ~ "n.s."
        ),
      label = ifelse(P.adj.EO < 0.05 | P.adj.LO < 0.05, Bacteria, ""), # Highlight significant points
      font = ifelse(P.adj.EO < 0.05 | P.adj.LO < 0.05, "bold.italic", "italic") # Font style for labels
    )
  
  head(data)
  axis_limits <- range(c(data[[x_col]], data[[y_col]]), na.rm = TRUE)
  
  ggplot(data, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(
      aes(
        size = 3, # You can adjust size dynamically by linking it to a column
        fill = enriched_in
      ),
      shape = 21, color = "black", alpha = 0.75
    ) +
    geom_hline(yintercept = 0, color = "grey", lty = "solid", lwd = 0.3) +
    geom_vline(xintercept = 0, color = "grey", lty = "solid", lwd = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "grey", lty = "solid", lwd = 0.3) +
    ggrepel::geom_text_repel(
      aes(label = label, fontface = font),
      color = "black",
      segment.color = "black",
      segment.size = 0.3,         # Make label lines thinner
      segment.ncp = 3,            # Number of control points (smooth lines)
      max.overlaps = 20,         # Ensure all significant labels appear
      min.segment.length = 0.2,   # Reduce the minimum length of connecting lines
      nudge_x = 0.05,             # Slightly adjust text placement to prevent overlapping
      nudge_y = 0.05,
      size = 3,
      seed = 123
    ) +
    scale_fill_manual(values = c('#C44600', 'dodgerblue3', 'white')) +
    theme_paper +
    theme(
      legend.box = "horizontal",
      legend.spacing.y = unit(0.1, "cm"),
      legend.position = c(.01, .99),
      legend.key.size = unit(0.75, "lines"),
      legend.justification = c(0, 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    labs(
      fill = paste0("P adj. < 0.05"),
      x = x_label,
      y = y_label
    ) +
    guides(
      size = "none",
      fill = guide_legend(order = 2, override.aes = list(size = 3.5))
    ) +
    scale_size(range = c(2, 5)) + 
    coord_equal() +
    xlab(x_label) +
    ylab(y_label)
  }

scatter <-plot_comparison(
  data = lmm.table.all,
  x_col = "Effect.size.EO",
  y_col = "Effect.size.LO",
  x_label = "Effect Size (EO)",
  y_label = "Effect Size (LO)")+
  scale_x_continuous(limits = c(-0.5,1),breaks = c(-0.5, -0.25, 0, 0.25, 0.5,0.75)) +
  scale_y_continuous(limits=c(-0.5,1),breaks = c(-0.5, -0.25, 0, 0.25, 0.5,0.75))
 

jaccard_table <- lmm.table.all  %>%
  mutate(
    enriched_in = case_when(
      # Both tests show significant enrichment
      (P.adj.EO < 0.05 & Effect.size.EO > 0) & 
        (P.adj.LO < 0.05 & Effect.size.LO > 0) ~ 'both enriched',
      
      # Both tests show significant depletion
      (P.adj.EO < 0.05 & Effect.size.EO < 0) & 
        (P.adj.LO < 0.05 & Effect.size.LO < 0) ~ 'both depleted',
      
      # Both tests non-significant
      P.adj.EO >= 0.05 & P.adj.LO >= 0.05 ~ 'both n.s.',
      
      # All other cases (disagreements or partial significance)
      TRUE ~ "one sig"
    ),
    label = ifelse(P.adj.EO < 0.05 | P.adj.LO < 0.05, Bacteria, ""),
    font = ifelse(P.adj.EO < 0.05 | P.adj.LO < 0.05, "bold.italic", "italic")
  ) %>% 
  select(Bacteria, P.adj.EO, P.adj.LO, Effect.size.LO, Effect.size.EO, enriched_in)

jaccard_index <- jaccard_table %>% 
  summarise(
    numerator = sum(enriched_in %in% c("both enriched", "both depleted")),
    denominator = n() - sum(enriched_in == "both n.s."),  # Exclude doubly non-sig features
    jaccard = numerator / denominator
  ) %>%
  pull(jaccard)

# Output the result
print(jaccard_index)

ggsave(scatter,file=here('figures','figure2','Figure2b.pdf'), height=6, width=6)

######################
# Figure 2d: AUROC for EO-CRC model and, evaluation of EO samples using LO CRC/CTR trained model  

# Load classifiers 

load(here('data','results','Training.eo.lo.crc.rf.models.Rdata'))
models <- list(models.eo.rf, siamcat.test.evaluated.eo.holdout.rf)
labels <- c("Classifier cross validated on EO-CRC","Classifier trained on LO-CRC and tested on EO-CRC")
trained_on <- list(NULL, models.lo.rf)
colours <- c( "darkorange", "black" )

eo_models_auc_plot <- plot_roc_siamcat_models(models, labels, colours, trained_on, alpha=0.7)


ggsave(eo_models_auc_plot, file=here('figures','figure2','Figure2d.pdf'),width = 7, height = 7)

# Figure 2e: Mean absolute SHAP value comparison of EO- and LO-CRC models

# Load EO-CRC and LO-CRC models shap median values

eo.model <- read.table(paste0(here('data','results','shap.analysis' ,'EO-CRC_median_shap_value.tsv')), sep='\t', header = TRUE) %>%
  select(feature, shap_value, spearman_sign) %>%
  group_by(feature) %>%
  summarise(mean_abs_shap = mean(abs(shap_value)) * unique(spearman_sign))


lo.model<- read.table(paste0(here('data','results','shap.analysis' ,'LO-CRC_median_shap_value.tsv')), sep='\t',header = T) %>%
  select(feature, shap_value, spearman_sign) %>%
  group_by(feature) %>%
  summarise(mean_abs_shap = mean(abs(shap_value)) * unique(spearman_sign))

merged_data<-NULL

# Merge the datasets by feature

merged_data <- eo.model %>%
  left_join(lo.model, by = "feature", suffix = c("_eo", "_lo")) %>%
  mutate(
    color = case_when(
      mean_abs_shap_eo > 0 & mean_abs_shap_lo > 0 ~ "Positive in both models",
      mean_abs_shap_eo < 0 & mean_abs_shap_lo < 0 ~ "Negative in both models",
      TRUE ~ "Other"
    ),
    abs_mean_abs_shap_eo = abs(mean_abs_shap_eo),
    abs_mean_abs_shap_lo = abs(mean_abs_shap_lo),
    mean_abs_shap_combined = (abs_mean_abs_shap_eo + abs_mean_abs_shap_lo) / 2 # Calculate mean of absolute SHAP values
  ) %>%
  mutate(color = factor(color, levels = c("Positive in both models", "Negative in both models", "Other")))

# Select top 15 features based on mean absolute SHAP values
top_features <- merged_data %>%
  arrange(desc(mean_abs_shap_combined)) %>% # Rank by combined mean SHAP values
  slice(1:15) %>%
  pull(feature)

# Add a column to indicate whether a feature should be labeled
merged_data <- merged_data %>%
  mutate(label = ifelse(feature %in% top_features, feature, NA))

# Create the scatter plot
scatter_plot <- ggplot(merged_data, aes(x = mean_abs_shap_eo, y = mean_abs_shap_lo)) +
  geom_point(
    aes(fill = color),
    shape = 21, 
    color = "black", 
    size = 3, 
    alpha = 0.75
  ) +
  geom_hline(yintercept = 0, color = "grey", linetype = "solid", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey", linetype = "solid", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "solid", linewidth = 0.3) +
  ggrepel::geom_text_repel(
    aes(label = label),
    fontface = "bold",
    color = "black",
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0.1,
    nudge_x = 0.005,
    nudge_y = 0.005,
    size = 3,
    seed = 123
  ) +
  scale_fill_manual(
    values = c(
      "Positive in both models" = '#C44600', 
      "Negative in both models" = 'dodgerblue3', 
      "Other" = 'grey'
    ),
    name = "Spearman sign"
  ) +
  theme_paper +
  theme(
    legend.box = "horizontal",
    legend.spacing.y = unit(0.1, "cm"),
    legend.position = c(.01, .99),
    legend.key.size = unit(0.75, "lines"),
    legend.justification = c(0, 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(linewidth = 0.3, color = "black")
  ) +
  labs(
    x = "Mean absolute SHAP values of EO-CRC model",
    y = "Mean absolute SHAP values of LO-CRC model"
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3.5))) +
  coord_equal()

print(scatter_plot)


ggsave(scatter_plot,file=  here('figures','figure2','Figure2e.pdf'), height=6, width=6)







