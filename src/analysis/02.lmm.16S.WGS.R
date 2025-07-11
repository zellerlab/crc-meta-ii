######################
# Linear Mixed Effect Model for sequencing approaches
######################
# This script performs linear mixed effects modeling (LMM) to compare microbial taxonomic differences
# between CRC and CTR samples, stratified by sequencing approaches (16S vs WGS).
# It also includes LMM analysis at the mOTU level (for WGS data only), with taxonomic annotation.
# Output tables are saved as an .Rdata object for downstream analysis and visualization.

# Load libraries and setup

source(here('src','utils.R'))


# Import genus level data and metadata

all.meta<- read.table(here('data','Metadata.all.samples.tsv'), sep = '\t', header = T)  %>%
  filter(Condition=='CRC'| Condition=='CTR') %>% as.data.frame() 


all.data <- read.table(here('data','Relab.all.samples.tsv'),sep='\t',check.names = F)  %>% 
  rownames_to_column('genus') %>% 
  filter(genus!='unassigned') %>% 
  column_to_rownames('genus') 

all.data<-  all.data[which(rowSums(all.data > 0) / ncol(all.data) > 0.1),]

# Subset by sequencing assay type (16S and WGS)

all.meta.16s <- all.meta %>%  filter(Assay=='16S')
all.data.16s<- all.data[,all.meta.16s$Sample_ID] 

all.meta.wgs <- all.meta %>%  filter(Assay=='WGS')  
all.data.wgs<- all.data[,all.meta.wgs$Sample_ID] 

lmm.table.wgs <- run_lmem(
data_df = all.data.wgs,
meta_df  = all.meta.wgs,ref_group='CTR',column_name='Condition')

lmm.table.16S <- run_lmem(
data_df = all.data.16s,
meta_df  = all.meta.16s, ref_group='CTR',column_name='Condition')

# Also run on mOTU level

# Load raw mOTU count data (WGS only)

all.data.motus <- read.table(here('data','Raw.counts.wgs.motus.all.samples.tsv'), sep='\t', header = T, check.names = F)

# Keep only WGS samples that are present in metadata

all.data.motus<-all.data.motus[colnames(all.data.motus)%in% all.meta.wgs$Sample_ID]

# Convert relative abundance table

all.data.motus<- all.data.motus %>%
  rownames_to_column('motus') %>%
  filter(rowSums(select(., where(is.numeric))) != 0)%>%
  column_to_rownames('motus')

all.data.motus<- all.data.motus %>% mutate(across(everything(), ~ . / sum(.)))

# Run LMM on normalized mOTU data to compare CRC vs CTR

lmm.table.motu <- run_lmem( data_df = all.data.motus,  meta_df  = all.meta, column_name='Condition', ref_group='CTR')

# Load mOTU taxonomy annotations

taxtable<-read_tsv(here('data','mOTUs_taxonomy_table.tsv')) %>% 
  mutate(species = str_split_fixed(taxon, "\\|", n = 8)[, 7],
         species = str_replace(species, "s__", ""),
         species = str_trim(species, side = "left"),  
         mOTUs = str_split_fixed(taxon, "\\|", n = 8)[, 8],
         mOTUs = str_trim(mOTUs, side = "left")  ,
         family = str_split_fixed(taxon, "\\|", n = 8)[, 5],
         family = str_replace(family, "f__", ""),
         family = str_trim(family, side = "left")
  ) %>% 
  mutate(
    species = str_replace(species, "^\\[([^\\]]+)\\]", "\\1")
  ) %>% 
  mutate(
    species =
      paste0(str_trim(str_replace(species, "\\[.*", "")), " [", mOTUs, "]")
  ) %>% select(-taxon) %>%
  mutate(species=str_replace(species,'species incertae sedis','sp.'))


lmm.table.motu<-lmm.table.motu %>% left_join(taxtable, by=c(Bacteria='mOTUs')) %>%
  mutate(species = str_trim(species))

# Save LMM results for WGS genus-level, 16S genus-level, and WGS mOTU-level to file
save(lmm.table.wgs, lmm.table.16S, lmm.table.motu, file = here('data','results', 'lmm.tables.16S.WGS.Rdata'))


