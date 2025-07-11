######################
# Linear Mixed Effect Model for EO/LO CRC and CTR comparisons
######################
# This script performs linear mixed effects modeling (LMM) to compare microbial taxonomic differences
# between early-onset (EO) and late-onset (LO) colorectal cancer (CRC) and control (CTR) samples.
# Stratified comparisons include: CRC vs CTR (overall), EO-CRC vs EO-CTR, LO-CRC vs LO-CTR,
# EO-CRC vs LO-CRC, and EO-CTR vs LO-CTR

# Load libraries and setup

source(here('src','utils.R'))
params <- yaml::read_yaml(here("src", "parameters.yml"))
plotting <- params$plotting


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


# General CRC vs CTR comparison (all ages)

lmm.table.general.crc <- run_lmem(
  data_df = all.data[, all.meta$Sample_ID],
  meta_df = all.meta,
  column_name = 'Condition',
  ref_group = 'CTR'
)

save(lmm.table.general.crc, file = here('data', 'results', 'lmm.table.crc.ctr.Rdata'))


# EO vs LO comparisons within CRC and CTR

# Subset EO samples
all.meta.eo <- all.meta %>% filter(Group %in% c('EO-CRC', 'EO-CTR'))
all.data.eo <- all.data[, all.meta.eo$Sample_ID]

# Subset LO samples
all.meta.lo <- all.meta %>% filter(Group %in% c('LO-CRC', 'LO-CTR'))
all.data.lo <- all.data[, all.meta.lo$Sample_ID]

# EO-CRC vs EO-CTR comparison
lmm.table.eo <- run_lmem(
  data_df = all.data.eo,
  meta_df  = all.meta.eo, column_name='Group', ref_group='EO-CTR')

lmm.table.eo <- lmm.table.eo %>%  #some dataset does not have either EO-CRC or EO-CTR samples
  select(where(~ !all(is.na(.))))

# LO-CRC vs LO-CTR comparison
lmm.table.lo <- run_lmem(
  data_df = all.data.lo,
  meta_df  = all.meta.lo, column_name='Group', ref_group='LO-CTR')


# EO vs LO within CRC
all.meta.eo.lo.crc <- all.meta %>% filter(Group %in% c('EO-CRC', 'LO-CRC'))
all.data.eo.lo.crc <- all.data[, all.meta.eo.lo.crc$Sample_ID]

lmm.table.eo.lo.crc <- run_lmem(
  data_df = all.data.eo.lo.crc,
  meta_df  = all.meta.eo.lo.crc,
  column_name = 'Group',
  ref_group = 'LO-CRC') %>%
  select(where(~ !all(is.na(.))))  


# EO vs LO within CTR
all.meta.eo.lo.ctr <- all.meta %>% filter(Group %in% c('EO-CTR', 'LO-CTR'))
all.data.eo.lo.ctr <- all.data[, all.meta.eo.lo.ctr$Sample_ID]

lmm.table.eo.lo.ctr <- run_lmem(
  data_df = all.data.eo.lo.ctr,
  meta_df  = all.meta.eo.lo.ctr,
  column_name = 'Group',
  ref_group = 'LO-CTR') %>%
  select(where(~ !all(is.na(.)))) 


# Save EO/LO comparison results
save(lmm.table.lo, lmm.table.eo, lmm.table.eo.lo.crc, lmm.table.eo.lo.ctr,
     file = here('data', 'results', 'lmm.tables.eo.lo.Rdata'))
