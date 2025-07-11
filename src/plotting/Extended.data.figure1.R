######################
# Extended data figure 1 
######################

# Generate Forest plot for  Assay types
# Function to clean and prepare LME data

load(here('data','results','lmm.tables.16S.WGS.Rdata'))

params <- yaml::read_yaml(here("src", "parameters.yml"))
plotting <- params$plotting


theme_paper <- ggembl::theme_presentation() +
  theme(
    axis.title = element_text(face = "bold", size = 12), # Bold axis titles
    panel.border = element_rect(fill=NA, colour='black', size=1.5),
    axis.text = element_text(face = "bold", size = 12),
  )


prepare_data <- function(df, assay) {
  df %>%
    rename(LME = Effect.size) %>%
    rename_with(~ gsub("Effect.size_", "", .), starts_with("Effect.size")) %>%
    arrange(LME) %>%
    mutate(Group = assay)
}

# Apply function to both datasets
merged <- bind_rows(
  prepare_data(lmm.table.wgs, "WGS"),
  prepare_data(lmm.table.16S, "16S")
) %>%
  arrange(desc(LME)) %>%
  select(-c(pr.shift, pr.CRC, pr.CTR, P.adj))

# Select top bacteria based on LME
top_bacteria <- merged %>%
  filter(Group=='WGS') %>% 
  slice(c(1:10, (n() - 9):n())) %>%
  pull(Taxon)

bacteria_filtered <- merged %>% filter(Taxon %in% top_bacteria)

# Function to extract effect sizes and confidence intervals
extract_effects <- function(df, assay) {
  df %>%
    filter(Group == assay) %>%
    select(-c(P.val, n.CTR, n.CRC, Group)) %>%
    pivot_longer(cols = -Taxon, names_to = "Cohort", values_to = "model effect size") %>%
    mutate(Assay = assay) %>%
    drop_na(`model effect size`)
}

# Merge WGS and 16S effect sizes (including CI)
combined_effects <- bind_rows(
  extract_effects(bacteria_filtered, "WGS"),
  extract_effects(bacteria_filtered, "16S")
)

# Order bacteria by LME (from 16S)
bacteria_order <- bacteria_filtered %>%
  filter(Group == "WGS") %>%
  distinct(Taxon, LME) %>%
  arrange(desc(LME)) %>%
  pull(Taxon)

# Apply ordering and Y-axis offset
combined_effects <- combined_effects %>%
  filter(Cohort!='conf.int_97.5 %' &Cohort!='conf.int_2.5 %') %>% 
  mutate(Taxon = factor(Taxon, levels = rev(bacteria_order)),
         Y_position = as.numeric(Taxon) + ifelse(Assay == "WGS", 0.2, -0.2))

# Get unique bacteria labels (to display only once)
bacteria_labels <- combined_effects %>%
  distinct(Taxon, .keep_all = TRUE) %>%
  select(Taxon, Y_position)

# Extract confidence intervals for plotting lines
ci_lines <- bacteria_filtered %>%
  select(Taxon, Assay = Group, lowCI = `conf.int_2.5 %`, highCI = `conf.int_97.5 %`) %>%
  mutate(Y_position = as.numeric(factor(Taxon, levels = rev(bacteria_order))) + 
           ifelse(Assay == "WGS", 0.2, -0.2)) %>%
  pivot_longer(cols = c(lowCI, highCI), names_to = "CI_type", values_to = "value")

ci_wide <- ci_lines %>%
  pivot_wider(names_from = CI_type, values_from = value) %>%
  mutate(
    Taxon = as.character(Taxon),
    Assay = as.character(Assay)
  )

# Define WGS and 16S cohort lists
cohorts_wgs <- c("FengQ_2015", "GaoR_2021", "HanniganGD_2018", "LiuNN_2022", 
                 "ThomasAM_2019_2", "ThomasAM_2019_1", "VogtmannE_2016", "WirbelJ_2019", 
                 "XuJ_2022", "YachidaS_2019", "YangJ_2019", "YangY_2021c", "YuJ_2017", "ZellerG_2014")

cohorts_16s <- c("YangY_2021a", "YangY_2021b", "BaxterNT_2016", "OkumuraS_2021", "KonishiY_2022", 
                 "DuX_2022", "TernesD_2022", "YoungC1_2021", "YoungC2_2021", "UchinoY_2021", 
                 "DengX_2018", "ZhouZ_2022", "YangTW_2019")

# Ensure Cohorts are factors and ordered by Assay type 
combined_effects <- combined_effects %>%
  mutate(Cohort = factor(Cohort, levels = c("LME", cohorts_wgs, cohorts_16s)),  # LME is first
         Assay = factor(Assay, levels = c("WGS", "16S"))) 

# Define shared shapes for WGS and 16S cohorts
shared_shapes <- c(16, 17, 18, 19, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# Assign distinct shapes for Cohorts, ensuring LME = 15
cohort_shapes <- c(
  "LME" = 15, 
  setNames(shared_shapes[1:length(cohorts_wgs)], cohorts_wgs), 
  setNames(shared_shapes[1:length(cohorts_16s)], cohorts_16s)  
)

assay_colors <- plotting$assay_colors

# Generate forest plot
forest_plot <- ggplot() +
  geom_segment(data = ci_wide,
               aes(x = lowCI, xend = highCI, y = Y_position, yend = Y_position, color = Assay),
               alpha = 0.8, size = 1)+
  geom_point(data = combined_effects %>% filter(Cohort != "LME"), 
             aes(x = `model effect size`, y = Y_position, color = Assay, shape = Cohort),
             alpha = 0.5, size = 2, show.legend = TRUE) +  
  geom_point(data = combined_effects %>% filter(Cohort == "LME"), 
             aes(x = `model effect size`, y = Y_position, color = Assay, shape = Cohort),
             alpha = 0.8, size = 4.5, show.legend = TRUE) + 
  labs(
    x = "Model Effect Size",
    y = "Genus") +
  
  scale_y_continuous(breaks = bacteria_labels$Y_position, labels = bacteria_labels$Taxon) +
  scale_color_manual(values = assay_colors) +
  scale_shape_manual(values = cohort_shapes) +
  guides(
    shape = guide_legend(order = 1, title = "Cohorts", override.aes = list(size = 4)), 
    color = guide_legend(order = 2, title = "Assay", override.aes = list(size = 5)),  
    shape = guide_legend(override.aes = list(shape = cohort_shapes, size = 4)) 
  ) +
  geom_vline(xintercept = 0, linetype = 'dotted', size = 1) +
  theme_paper +
  theme(
    legend.box = "vertical",
    legend.spacing.y = unit(0.1, "cm"),
    legend.position = c(.01, .99),
    legend.key.size = unit(0.75, "lines"),
    legend.justification = c(0, 1)
  )

ggsave( plot = forest_plot,  
        filename =here('figures','extended.data.figure1','Extended.Data.Figure1d.pdf')  ,width = 9,  height = 12)


