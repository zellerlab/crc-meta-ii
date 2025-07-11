# Load utility functions
source(here('src','utils.R'))

# Load and filter data/metadata to include only CRC and Control samples

all.meta<- read.table(here('data','Metadata.all.samples.tsv'), sep = '\t', header = T)  %>%
  filter(Condition=='CRC'| Condition=='CTR') %>% as.data.frame() 

all.data <- read.table(here('data','Relab.all.samples.tsv'),sep='\t',check.names = F)  %>% 
  rownames_to_column('genus') %>% 
  filter(genus!='unassigned') %>% 
  column_to_rownames('genus') 

all.data <- all.data[,all.meta$Sample_ID]

# Load raw count data

all.data.raw <- read.table(here('data','Raw.counts.all.samples.tsv'),sep='\t',check.names = F)  %>% 
  rownames_to_column('genus') %>% 
  filter(genus!='unassigned') %>% 
  column_to_rownames('genus')

all.data.raw<- all.data.raw[,all.meta$Sample_ID]

# Rarefaction: normalize read depth to 1000 reads per sample
# Calculate Shannon diversity and genus richness on rarefied data

set.seed(1)
# Perform rarefaction and calculate diversity and richness

count_rar <- as.matrix(t(vegan::rrarefy(t(round(all.data.raw[,all.meta$Sample_ID],0)),1000)))

count_rar_rel <- prop.table(count_rar,2)
div <- vegan::diversity(t(count_rar_rel), index = "shannon") %>% enframe(name = "Sample_ID",value = "Shannon diversity")
rich <- colSums(count_rar_rel > 0) %>% enframe(name = "Sample_ID",value = "Genus richness")

# Combine diversity and richness data

diversity_richness_df <- full_join(div,rich,by = "Sample_ID")

diversity_richness_df<- diversity_richness_df %>% left_join(all.meta %>% select(Cohort,Sample_ID, Assay,Condition), by='Sample_ID') %>% 
  mutate(Cohort=as.factor(Cohort)) 


# Save the final combined data table

write.table(diversity_richness_df,file = here('data','results','Diversity.richness.tsv'),sep = '\t',row.names = F)






















