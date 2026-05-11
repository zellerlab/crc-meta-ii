######################
# PCoA analysis
######################
# This script performs Principal Coordinates Analysis (PCoA) on genus-level microbiome data
# using Bray-Curtis and Canberra distances, grouped by dataset, disease condition, and hypervariable regions of 16S data.

source(here('src','utils.R'))
params <- yaml::read_yaml(here("src", "parameters.yml"))
plotting <- params$plotting

# Load data
all.meta<- read.table(here('data','Metadata.all.samples.tsv'), sep = '\t', header = T)  %>%
  filter(Condition=='CRC'| Condition=='CTR') %>% as.data.frame() 

all.data <- read.table(here('data','Relab.all.samples.tsv'),sep='\t',check.names = F)  %>% 
  rownames_to_column('genus') %>% 
  filter(genus!='unassigned') %>% 
  column_to_rownames('genus') 

# Sanity check
all.data <- all.data[,colnames(all.data) %in% (all.meta$Sample_ID)]
stopifnot(all(sort(all.meta$Sample_ID) == sort(colnames(all.data))))

# Load plotting parameters
condition_colors <- unlist(plotting$condition_colors)
brown_cohorts <- unlist(plotting$brown_cohorts)
cyan_cohorts <- unlist(plotting$cyan_cohorts)
brown_palette <- unlist(plotting$brown_palette)
cyan_palette <- unlist(plotting$cyan_palette)

# Count samples per cohort
cohort_colors <- c(
  setNames(cyan_palette[1:length(cyan_cohorts)], cyan_cohorts),
  setNames(brown_palette[1:length(brown_cohorts)], brown_cohorts)
)

cohort_counts <- all.meta %>%
  group_by(Cohort) %>%
  summarise(Sample_Count = n())

all.meta_with_counts <- all.meta %>%
  left_join(cohort_counts, by = "Cohort")

# Format cohort labels
cohort_labels <- cohort_counts %>% 
  mutate(Label = paste0(Cohort, " (N=", Sample_Count, ")")) %>%
  left_join(all.meta %>% select(Cohort,Assay) %>%  unique(), by='Cohort') %>% 
  arrange(desc(Assay == "WGS"), Assay) %>% 
  pull(Label)


names(cohort_labels) <- cohort_counts %>% 
  mutate(Label = paste0(Cohort, " (N=", Sample_Count, ")")) %>%
  left_join(all.meta %>% select(Cohort,Assay) %>%  unique(), by='Cohort') %>% arrange(desc(Assay == "WGS"), Assay)  %>% pull(Cohort)

color_combinations <- list(
  list(cyan = cyan_palette, brown = brown_palette)
)

cohort_colors <- c(setNames(color_combinations[[1]]$cyan[1:14], cyan_cohorts),
                   setNames(color_combinations[[1]]$brown[1:13], brown_cohorts))

######################
# Figure 1a: PCoA by Cohort using Bray-Curtis distance matrix

pt_pcoa_study <- plot_PCoA(
  meta_df = all.meta_with_counts,
  mat = all.data,
  method = "bray",
  grouping = 'Cohort',
  transformed = 'no'
)

current_plot_bray <- pt_pcoa_study +
  scale_color_manual(values = alpha(cohort_colors, 0.9), breaks = names(cohort_labels), labels = cohort_labels) +
  scale_fill_manual(values = alpha(cohort_colors, 0.9), breaks = names(cohort_labels), labels = cohort_labels) +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2))


ggsave(current_plot_bray, file= here('figures','figure1','Figure1a.pdf'), height = 8,width = 10)

######################
# Extended Data Figure 1a: PCoA by Cohort using Canberra distance matrix

pt_pcoa_study_canberra <- plot_PCoA(
  meta_df = all.meta_with_counts,
  mat = all.data,
  method = "canberra",
  grouping = 'Cohort',
  transformed = 'no'
)

current_plot_canberra <- pt_pcoa_study_canberra +
  scale_color_manual(values = alpha(cohort_colors, 0.9), breaks = names(cohort_labels), labels = cohort_labels) +
  scale_fill_manual(values = alpha(cohort_colors, 0.9), breaks = names(cohort_labels), labels = cohort_labels) +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2))

ggsave(current_plot_canberra, file=  here('figures','extended.data.figure1','Extended.Data.Figure1a.pdf'), height = 8,width = 10)

######################
# Extended Data Figure 1b: PCoA by Disease Condition using Bray-Curtis distance matrix

pt_pcoa_disease <- plot_PCoA(
  meta_df = all.meta_with_counts,
  mat = all.data,
  method = "bray",
  grouping = 'Condition',
  transformed = 'no',
  point_alpha= 0.6
)

pt_pcoa_disease_bray<- pt_pcoa_disease +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors) +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2))

ggsave(pt_pcoa_disease_bray, file=here('figures','extended.data.figure1','Extended.Data.Figure1b.pdf'), height = 8,width = 10)

######################
# Extended Data Figure 1c: PCoA by Disease Condition using Canberra distance matrix

pt_pcoa_disease_can <- plot_PCoA(
  meta_df = all.meta_with_counts,
  mat = all.data,
  method = "canberra",
  grouping = 'Condition',
  transformed = 'no',
  point_alpha= 0.6
)

pt_pcoa_disease_can<- pt_pcoa_disease_can +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors) +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2))

ggsave(pt_pcoa_disease_can, file= here('figures','extended.data.figure1','Extended.Data.Figure1c.pdf'), height = 8,width = 10)

######################
# Extended Data Figure 1e: PCoA by Hypervariable Regions using Canberra distance matrix (only 16S data)

all.meta.16s <- all.meta %>%  filter(Assay=='16S')
all.data.16s <- all.data[,all.meta.16s$Sample_ID]

# using canberra distance matrix
pt_pcoa_hypervariable <- plot_PCoA(
  meta_df = all.meta.16s,
  mat = all.data.16s, method = "canberra", grouping = 'Hypervariable_region',transformed='no')


pt_pcoa_hypervariable <- pt_pcoa_hypervariable +
  scale_fill_manual(values=c('#C9CBA3' ,'#FFE1A8','#E26D5C', '#723D46')) +
  coord_fixed(ratio = 1)

ggsave(pt_pcoa_hypervariable, file=here('figures','extended.data.figure1','Extended.Data.Figure1e.pdf'), height = 6,width = 8)

######################
# Extended Data Figure 1f: PCoA by Hypervariable Regions using Bray-Curtis distance matrix (only 16S data)

pt_pcoa_hypervariable_bray <- plot_PCoA(
  meta_df = all.meta.16s,
  mat = all.data.16s, method = "bray", grouping = 'Hypervariable_region',transformed='no')

pt_pcoa_hypervariable_bray <- pt_pcoa_hypervariable_bray +
  scale_fill_manual(values=c('#C9CBA3' ,'#FFE1A8','#E26D5C', '#723D46')) +
  coord_fixed(ratio = 1)

ggsave(pt_pcoa_hypervariable_bray, file= here('figures','extended.data.figure1','Extended.Data.Figure1f.pdf'), height = 6,width = 8)

######################
# Extended Data Figure 2d: PCoA by disease condition for EO-CRC and LO-CRC samples only

all.meta.eo.lo.crc <- all.meta %>%  filter(Group== "EO-CRC" | Group=="LO-CRC")
all.data.eo.lo.crc <- all.data[,all.meta.eo.lo.crc$Sample_ID]


pt_pcoa_eo_lo_crc <- plot_PCoA(
  meta_df = all.meta.eo.lo.crc,
  mat = all.data.eo.lo.crc, method = "bray", grouping = 'Group',transformed='no')

ggsave(pt_pcoa_eo_lo_crc + scale_fill_manual(values = plotting$condition_colors),
       file= here('figures','extended.data.figure2','Extended.Data.Figure2d.pdf'), height = 6,width = 8)
