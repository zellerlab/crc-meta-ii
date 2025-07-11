######################
# Linear Mixed Effect Model for functional profiles
######################
# This script performs linear mixed effects modeling (LMM) to compare functional microbiome features
# (virulence factors, KEGG orthologs, KEGG modules, KEGG pathways, and gut metabolic modules)
# between colorectal cancer (CRC) and control (CTR) samples in WGS datasets.


source(here('src','utils.R'))

######################

# Read and subset metadata for WGS CRC/CTR samples
all.meta.wgs<- read.table(here('data','Metadata.all.samples.tsv'), sep = '\t', header = T)  %>%
  filter(Condition=='CRC'| Condition=='CTR') %>% as.data.frame()  %>% 
  filter(Assay=='WGS') 

# Load relative abundance of virulence factors
VFs_relab<- read.table(here('data','functional_data','VFGs.relab.wgs.samples.tsv'), sep = '\t', header = T, check.names = F)  

# Run LMM for virulence factor profiles
lmm.table.vf <- run_lmem_functional(
  data_df = VFs_relab,
  meta_df  = all.meta.wgs.vf %>% select(-Sample_ID) %>% rename(Sample_ID=Sample_ID_VFs), 
  column_name='Condition', ref_group='CTR',feature_column_name = 'VF')


# KEGG Orthologs: Prepare KO matrix for KEGG Mapper & LMM
KO_normalized_matrix <- read.table(here('data', 'functional_data','All.KEGG.ko.combined.scaled.together.tsv'), sep = '\t', header = T) 
colnames(KO_normalized_matrix)<- gsub("\\.", "-", colnames(KO_normalized_matrix))

# Subset and filter KOs
KO_normalized_matrix <- KO_normalized_matrix %>%  select(any_of(all.meta.wgs$Sample_ID_KOs))
KO_normalized_matrix<-  KO_normalized_matrix[which(rowSums(KO_normalized_matrix > 0) / ncol(KO_normalized_matrix) > 0.1),]

# Run LMM for KEGG orthologs
lmm.KEGG.ko <- run_lmem_functional(
  data_df = KO_normalized_matrix,
  meta_df  = all.meta.wgs %>% select(-Sample_ID) %>% rename(Sample_ID=Sample_ID_KOs),
  column_name='Condition', ref_group='CTR',feature_column_name = 'KO')

# Generate a table for KEGG Mapper visualization (significant red/blue, nonsignificant grey)
KEGG_mapper_table<- lmm.KEGG.ko %>%
  filter(P.adj< 0.05) %>% filter(Effect.size> 0 ) %>% pull(Bacteria) %>%  as.data.frame() %>% mutate(color= "red") %>%
  rbind(lmm.KEGG.ko %>%
          filter(P.adj< 0.05) %>% filter(Effect.size< 0) %>% pull(Bacteria) %>%  as.data.frame() %>% mutate(color= "dodgerblue")) %>%
  rbind(lmm.KEGG.ko %>%
          filter(P.adj > 0.05) %>% pull(Bacteria) %>%  as.data.frame() %>% mutate(color= "grey")) %>%
  write.table(file = here('data','results','KEGG_mapper_table.tsv'), row.names = F, quote = F,col.names = F)

# Annotate LMM results with KEGG mapper color
lmm.KEGG.ko <- lmm.KEGG.ko %>% inner_join(KEGG_mapper_table %>% rename(c('KO'='.','KEGG_mapper_color' ='color')), by='KO')

# KEGG Modules: Prepare KEGG module matrix & LMM

KO_modules<- read.table(here('data','functional_data','All.KEGG.modules.combined.scaled.together.tsv'), sep='\t', check.names = F)
KO_modules<-  KO_modules[which(rowSums(KO_modules > 0) / ncol(KO_modules) > 0.1),]
colnames(KO_modules)<- gsub("\\.", "-", colnames(KO_modules))
KO_modules <- KO_modules %>%  select(any_of(all.meta.wgs$Sample_ID_KOs))

# Load module descriptions
Kegg_modules_0.6 <- read.table(here('data','functional_data','kegg_modules_filtered_0_6.csv'), sep = ',', header=T) %>% 
  select(module_ko=num , description) %>% unique() %>% 
  mutate(description= str_replace(description, "^\\s+", ""))

# Aggregate values by module description
KO_modules <- KO_modules %>% rownames_to_column('module') %>% pivot_longer(-module) %>% 
  left_join(Kegg_modules_0.6, by = c(module='module_ko')) %>% 
  filter(module %in%  Kegg_modules_0.6$module_ko) %>%
  select(-module) %>% 
  group_by(description, name) %>% 
  summarise(value=sum(value)) %>%
  pivot_wider(names_from = name, values_from = value)

KO_modules <- KO_modules %>% column_to_rownames('description')

# Run LMM for KEGG modules
lmm.KEGG.modules <- run_lmem_functional(
  data_df = KO_modules,feature_column_name = 'KEGG_modules',
  meta_df  = all.meta.wgs %>% select(-Sample_ID) %>% rename(Sample_ID=Sample_ID_KOs) , column_name='Condition', ref_group='CTR')


# KEGG Pathways: Prepare KEGG pathways matrix & LMM

KO_pathways<- read.table(here('data', 'functional_data','All.KEGG.pathways.combined.scaled.together.tsv'), sep='\t', check.names = F)

Kegg_pathways_0.5 <- read.table(here('data', 'functional_data','kegg_pathways_filtered_0_5.csv'), sep = ',') %>% 
  select(path_ko=V1 , description=V4) %>%
  mutate(path_ko= str_remove(path_ko, 'path:')) %>% unique()

# Aggregate values by pathway description
KO_pathways <- KO_pathways %>% rownames_to_column('ko') %>% pivot_longer(-ko) %>% 
  left_join(Kegg_pathways_0.5, by = c(ko='path_ko')) %>% 
  filter(ko %in%  Kegg_pathways_0.5$path_ko) %>%
  select(-ko) %>% 
  group_by(description, name) %>% 
  summarise(value=sum(value)) %>%
  pivot_wider(names_from = name, values_from = value)


colnames(KO_pathways)<- gsub("\\.", "-", colnames(KO_pathways))
KO_pathways <- KO_pathways %>%  select(any_of(all.meta.wgs$Sample_ID_KOs)) %>% column_to_rownames('description')
KO_pathways<-  KO_pathways[which(rowSums(KO_pathways > 0) / ncol(KO_pathways) > 0.1),]

# Run LMM for KEGG pathways
lmm.KEGG.pathways <- run_lmem_functional(
  data_df = KO_pathways,
  feature_column_name = 'KEGG_pathways',
  meta_df  = all.meta.wgs %>% select(-Sample_ID) %>% rename(Sample_ID=Sample_ID_KOs) ,
  column_name='Condition', ref_group='CTR')

# Gut Metabolic Modules (GMMs): LMM analysis
GMMs_relab<- readRDS(here('data','functional_data','merged_GMMs.RDS')) %>% select(-Description) %>% column_to_rownames('Module')
colnames(GMMs_relab) <- gsub("\\.", "-", colnames(GMMs_relab))
GMMs_relab <- GMMs_relab %>%  select(any_of(all.meta.wgs$Sample_ID_KOs))

# Run LMM for GMMs
lmm.GMMs <- run_lmem_functional(
  data_df = GMMs_relab,
  meta_df  = all.meta.wgs %>% select(-Sample_ID) %>% rename(Sample_ID=Sample_ID_KOs) , column_name='Condition', ref_group='CTR',feature_column_name = 'Module')

lmm.GMMs<- lmm.GMMs %>% 
  left_join(readRDS(here('data','functional_data','merged_GMMs.RDS')) %>% select(Description, Module), by='Module')


# Save lmm tables for all functional profiles
save(lmm.table.vf,lmm.KEGG.ko,lmm.KEGG.modules, lmm.KEGG.pathways,lmm.GMMs, file=here('data','results','lmm.tables.functional.profiles.Rdata'))
     


