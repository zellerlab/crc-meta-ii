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

vfdb.category <- readxl::read_xls(here('data','functional_data','VFs.xls'), skip = 1)

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
