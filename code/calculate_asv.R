#### 1. Set-up ####

# Load libraries
library(phyloseq)
library(tidyverse)

# Read in data
ps <- readRDS("data/Manistee_MiSeq_data/asv_data/phyloseq.rds")

# Source functions
source("code/functions.R")

#### 2. Format ASV table ####

# Get ASV table from phyloseq object
ps.ASV.tbl <- as.data.frame(otu_table(ps)) %>% 
  as_tibble(., rownames = "SAMPLE.ID")

# Get sample data from phyloseq object
ps.sample.tbl <- as.data.frame(sample_data(ps)) %>% 
  as_tibble(.)

# Get taxonomy data
ps.tax.tbl <- as.data.frame(tax_table(ps)) %>% 
  as_tibble(., rownames = "ASV") %>% 
  rename_all(toupper) %>% 
  pivot_longer(-ASV, names_to = "RANKING", values_to = "TAXON") %>% 
  mutate(TAXON = sapply(TAXON, tax_clean, simplify = TRUE)) %>% 
  pivot_wider(id_cols = ASV, names_from = "RANKING", values_from = TAXON) %>% 
  mutate(SPECIES = mapply(s_unclass, KINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES, SIMPLIFY = TRUE), 
         GENUS = mapply(g_unclass, KINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS, SIMPLIFY = TRUE), 
         FAMILY = mapply(f_unclass, KINGDOM, PHYLUM, CLASS, ORDER, FAMILY, SIMPLIFY = TRUE), 
         ORDER = mapply(o_unclass, KINGDOM, PHYLUM, CLASS, ORDER, SIMPLIFY = TRUE), 
         CLASS = mapply(c_unclass, KINGDOM, PHYLUM, CLASS, SIMPLIFY = TRUE), 
         PHYLUM = mapply(p_unclass, KINGDOM, PHYLUM, SIMPLIFY = TRUE))

# Format
ASV.all.tbl <- ps.ASV.tbl %>% 
  pivot_longer(-SAMPLE.ID, names_to = "ASV", values_to = "COUNT") %>% 
  # Add sample data
  inner_join(ps.sample.tbl, ., by = "SAMPLE.ID") %>% 
  # Remove mock communities and negative controls
  filter(TYPE == "environmental") %>% 
  # Combine sequences by substrate (i.e., roots or soil)
  group_by(SUBSTRATE, PLOT, STAND, ASV) %>% 
  summarise(COUNT = sum(COUNT)) %>% 
  ungroup() %>% 
  # Remove any ASVs with no sequences (after mock and negative controls removed)
  group_by(ASV) %>% 
  mutate(ASV.TOTAL = sum(COUNT)) %>% 
  ungroup() %>% 
  filter(ASV.TOTAL > 0) %>% 
  select(-ASV.TOTAL)

#### 3. Save and clean up data ####

# Unload phyloseq
detach("package:phyloseq", unload = TRUE)

# Save ASV table
saveRDS(ps.tax.tbl, file = "data/working/ps.tax.tbl.RData")
saveRDS(ASV.all.tbl, file = "data/working/ASV.all.tbl.Rdata")

# Remove unnecessary data
rm(ps, ps.ASV.tbl, ps.sample.tbl)
