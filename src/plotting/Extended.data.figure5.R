######################
# Extended data figure 5
######################

source(here('src','utils.R'))
load(here('data','results','lmm.tables.functional.profiles.Rdata'))
######################
# Extended data figure 5a

gmm_data_long <- lmm.GMMs %>%
  filter(P.adj < 0.01) %>%  # Filter significant results
  pivot_longer(
    cols = starts_with("Effect.size"),
    names_to = "Dataset",
    values_to = "Effect.size"
  ) %>%
  mutate(
    Dataset = ifelse(Dataset == "Effect.size", "LME", gsub("Effect.size_", "", Dataset)) 
  )


# Reorder features (Description) by LME
description_order <- gmm_data_long %>%
  filter(Dataset == "LME") %>%
  arrange(desc(abs(Effect.size))) %>%
  slice(1:30) %>% 
  arrange(desc((Effect.size))) %>%
  pull(Description)


gmm_data_long <- gmm_data_long %>%
  mutate(Description = factor(Description, levels = rev(description_order))) %>% drop_na(Description)

# Get unique datasets
unique_datasets <- unique(gmm_data_long$Dataset)

# Generate a list of shapes, ensuring enough shapes for all datasets
shapes <- c(15, 0, 1, 3, 4, 7, 8, 9, 10, 12, 11, 2, 5, 6, 13, 14, 18, 16, 17, 19, 20, 21)

# Map datasets to shapes
if (length(unique_datasets) > length(shapes)) {
  stop("Not enough shapes for all datasets. Please add more shape codes.")
}

dataset_shapes <- setNames(shapes[1:length(unique_datasets)], unique_datasets)

# Plot the forest plot
forest_plot_gmms <- ggplot() +
  # Plot LME as larger points
  geom_point(
    data = gmm_data_long %>% filter(Dataset == "LME"),
    aes(x = Effect.size, y = Description, shape = Dataset),
    size = 4, color = '#C44600'
  ) +
  # Plot other datasets as smaller points
  geom_point(
    data = gmm_data_long %>% filter(Dataset != "LME"),
    aes(x = Effect.size, y = Description, shape = Dataset),
    size = 2, color = '#C44600'
  ) +
  # Add error bars for LME (if applicable)
  geom_errorbarh(
    data = gmm_data_long %>%
      filter(Dataset == "LME"),
    aes(xmin = `conf.int_2.5 %`, xmax = `conf.int_97.5 %`, y = Description),
    height = 0.2, color = "black"
  ) +
  # Add vertical line for zero effect
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 1) +
  scale_alpha_continuous(guide = 'none') +
  scale_shape_manual(
    values = dataset_shapes  # Dynamically generated shapes
  ) +
  theme_minimal() +
  xlab("Model Effect Size") +
  ylab("GMMs") +
  theme_paper +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    legend.box = "vertical"
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 4), title = "Cohort")
  )

ggsave(forest_plot_gmms,file= here('figures','extended.data.figure5','Extended.Data.Figure5a.pdf'), width = 9,  height = 10)

######################
# Extended data figure 5b

# # Generate a forest plot

# # Reshape the data for all datasets including LME
KEGG.pathways.long <- lmm.KEGG.pathways %>%
  filter(P.adj < 0.05) %>%  # Filter significant results
  pivot_longer(
    cols = starts_with("Effect.size"),
    names_to = "Dataset",
    values_to = "Effect.size"
  ) %>%
  mutate(
    Dataset = ifelse(Dataset == "Effect.size", "LME", gsub("Effect.size_", "", Dataset))  # Rename LME for clarity
  ) %>% drop_na('KEGG_pathways')


# Reorder features (Description) by LME
kegg_pathway_order <- KEGG.pathways.long %>%
  filter(Dataset == "LME") %>%
  arrange(desc(abs(Effect.size))) %>%
  slice(1:30) %>%
  arrange(desc( (Effect.size))) %>%
  pull(KEGG_pathways)

# Get unique datasets
unique_datasets <- unique(KEGG.pathways.long$Dataset)

# Generate a list of shapes, ensuring enough shapes for all datasets
shapes <- c(15, 0, 1, 3, 4, 7, 8, 9, 10, 12, 11, 2, 5, 6, 13, 14, 18, 16, 17, 19, 20, 21)

dataset_shapes <- setNames(shapes[1:length(unique_datasets)], unique_datasets)

# Select first 15 and last 15 module names
selected_pathways <- kegg_pathway_order


# Filter the dataset to include only the selected modules
KEGG.pathway.long.filtered <- KEGG.pathways.long %>%
  filter(KEGG_pathways %in% selected_pathways) %>%
  mutate(KEGG_pathways = factor(KEGG_pathways, levels = rev(selected_pathways)))


forest_plot_ko_pathway <- ggplot() +
  geom_point(
    data = KEGG.pathway.long.filtered %>% filter(Dataset == "LME"),
    aes(x = Effect.size, y = KEGG_pathways, shape = Dataset),
    size = 4, color = '#C44600'
  ) +
  geom_point(
    data = KEGG.pathway.long.filtered %>% filter(Dataset != "LME"),
    aes(x = Effect.size, y = KEGG_pathways, shape = Dataset),
    size = 2, color = '#C44600'
  ) +
  geom_errorbarh(
    data = KEGG.pathway.long.filtered %>% filter(Dataset == "LME"),
    aes(xmin = `conf.int_2.5 %`, xmax = `conf.int_97.5 %`, y = KEGG_pathways),
    height = 0.2, color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 1) +
  scale_alpha_continuous(guide = 'none') +
  scale_shape_manual(values = dataset_shapes) +
  theme_minimal() +
  xlab("Model Effect Size") +
  ylab("KEGG Pathways") +
  theme_paper +
  theme(
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    legend.box = "vertical"
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 4), title = "Cohort")
  )


ggsave(forest_plot_ko_pathway,file= here('figures','extended.data.figure5','Extended.Data.Figure5b.pdf'), width = 14,  height = 10)


######################
# Extended data figure 5c

# Reshape the data for all datasets including LME
KEGG.modules.long <- lmm.KEGG.modules %>%
  filter(P.adj < 0.05) %>%  # Filter significant results
  pivot_longer(
    cols = starts_with("Effect.size"),
    names_to = "Dataset",
    values_to = "Effect.size"
  ) %>%
  mutate(
    Dataset = ifelse(Dataset == "Effect.size", "LME", gsub("Effect.size_", "", Dataset))  # Rename LME for clarity
  ) %>% drop_na('KEGG_modules')


# Reorder features (Description) by LME
kegg_module_order <- KEGG.modules.long %>%
  filter(Dataset == "LME") %>%
  arrange(desc(abs(Effect.size))) %>%
  slice(1:30) %>% 
  arrange(desc((Effect.size))) %>%
  pull(KEGG_modules)

# Get unique datasets
unique_datasets <- unique(KEGG.modules.long$Dataset)

# Generate a list of shapes, ensuring enough shapes for all datasets
shapes <- c(15, 0, 1, 3, 4, 7, 8, 9, 10, 12, 11, 2, 5, 6, 13, 14, 18, 16, 17, 19, 20, 21)

dataset_shapes <- setNames(shapes[1:length(unique_datasets)], unique_datasets)

selected_modules <-kegg_module_order

# Filter the dataset to include only the selected modules
KEGG.module.long.filtered <- KEGG.modules.long %>%
  filter(KEGG_modules %in% selected_modules) %>%
  mutate(KEGG_modules = factor(KEGG_modules, levels = rev(selected_modules)))


forest_plot_ko_module <- ggplot() +
  geom_point(
    data = KEGG.module.long.filtered %>% filter(Dataset == "LME"),
    aes(x = Effect.size, y = KEGG_modules, shape = Dataset),
    size = 4, color = '#C44600'
  ) +
  geom_point(
    data = KEGG.module.long.filtered %>% filter(Dataset != "LME"),
    aes(x = Effect.size, y = KEGG_modules, shape = Dataset),
    size = 2, color = '#C44600'
  ) +
  geom_errorbarh(
    data = KEGG.module.long.filtered %>% filter(Dataset == "LME"),
    aes(xmin = `conf.int_2.5 %`, xmax = `conf.int_97.5 %`, y = KEGG_modules),
    height = 0.2, color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 1) +
  scale_alpha_continuous(guide = 'none') +
  scale_shape_manual(values = dataset_shapes) +
  theme_minimal() +
  xlab("Model Effect Size") +
  ylab("KEGG modules") +
  theme_paper +
  theme(
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    legend.box = "vertical"
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 4), title = "Cohort")
  )


ggsave(forest_plot_ko_module,file= here('figures','extended.data.figure5','Extended.Data.Figure5c.pdf'), height = 5, width = 5)


























