#### 1. Set-up ####

# Load libraries
library(tidyverse)

# Read in data
ASV.all.tbl <- readRDS("data/working/ASV.all.tbl.RData")
ps.tax.tbl <- readRDS("data/working/ps.tax.tbl.RData")

# Source functions
source("code/functions.R")

#### 2. Format taxonomy data

# Update tax table
ps.tax.modified.tbl <- ps.tax.tbl %>% 
  mutate(FAMILY = ifelse(GENUS == "Cryptococcus" & FAMILY == "Cryptococcaceae", "Tremellaceae", FAMILY)) %>% 
  mutate(GENUS = ifelse(is.na(GENUS) & CLASS == "Lecanoromycetes", "Lecanoromycetes incertae sedis", GENUS)) %>% 
  mutate(GENUS = ifelse(is.na(GENUS) & ORDER == "Helotiales", "Helotiales incertae sedis", GENUS)) %>% 
  mutate(GENUS = ifelse(is.na(GENUS) & ORDER == "Hypocreales", "Hypocreales incertae sedis", GENUS)) %>% 
  mutate(GENUS = ifelse(is.na(GENUS) & PHYLUM == "Rozellomycota","Rozellomycota incertae sedis", GENUS)) %>% 
  mutate(GENUS = ifelse(GENUS == "unclassified GS35" & ORDER == "unclassified GS35","unclassified GS35 (c)", GENUS)) %>% 
  mutate(GENUS = ifelse(GENUS == "unclassified GS35" & ORDER == "GS35","unclassified GS35 (o)", GENUS)) %>% 
  mutate(GENUS = ifelse(GENUS == "unclassified GS18" & ORDER == "unclassified GS18","unclassified GS18 (c)", GENUS)) %>% 
  mutate(GENUS = ifelse(GENUS == "unclassified GS18" & ORDER == "GS18","unclassified GS18 (c)", GENUS))

#### 3. Calculate genus-level abundances ####

# Abundances
genus.all.tbl <- ASV.all.tbl %>% 
  inner_join(ps.tax.modified.tbl, ., by = "ASV") %>% 
  # Calculate by plot
  group_by(GENUS, SUBSTRATE, STAND, PLOT) %>% 
  summarize(COUNT = sum(COUNT)) %>% 
  ungroup() %>% 
  select(GENUS, SUBSTRATE, STAND, PLOT, COUNT)

# Total sequences by substrate after outliers removed
genus.total.seqs.no.outliers.tbl <- genus.all.tbl %>% 
  filter(PLOT != "Plot_6" & PLOT != "Plot_19" & PLOT != "Plot_43" & PLOT != "Plot_60") %>% 
  group_by(SUBSTRATE) %>% 
  summarize(SUBSTRATE.TOTAL.COUNT = sum(COUNT))

#### 5. Calculate percentage of total sequences by substrate ####

# Percentage of total seqs, soil
soil.genus.perc.tbl <- genus.all.tbl %>% 
  filter(PLOT != "Plot_6" & PLOT != "Plot_19" & PLOT != "Plot_43" & PLOT != "Plot_60") %>% 
  filter(SUBSTRATE == "soil") %>% 
  group_by(GENUS) %>% 
  summarise(GENUS.COUNT = sum(COUNT)) %>% 
  ungroup() %>% 
  mutate(TOTAL.COUNT = sum(GENUS.COUNT)) %>% 
  mutate(GENUS.PERC = 100 * (GENUS.COUNT / TOTAL.COUNT)) %>% 
  select(GENUS, GENUS.PERC)

# Percentage of total seqs, roots
roots.genus.perc.tbl <- genus.all.tbl %>% 
  filter(PLOT != "Plot_6" & PLOT != "Plot_19" & PLOT != "Plot_43" & PLOT != "Plot_60") %>% 
  filter(SUBSTRATE == "decomp.roots") %>% 
  group_by(GENUS) %>% 
  summarise(GENUS.COUNT = sum(COUNT)) %>% 
  ungroup() %>% 
  mutate(TOTAL.COUNT = sum(GENUS.COUNT)) %>% 
  mutate(GENUS.PERC = 100 * (GENUS.COUNT / TOTAL.COUNT)) %>% 
  select(GENUS, GENUS.PERC)

#### 5. Hellinger transform and removed unclassified genera ####

# Transform soil genera abundances
soil.genus.no.outliers.tmp.tbl <- genus.all.tbl %>% 
  # Drop plots that were SOM outliers
  filter(PLOT != "Plot_6" & PLOT != "Plot_19" & PLOT != "Plot_43" & PLOT != "Plot_60") %>% 
  filter(SUBSTRATE == "soil") %>% 
  # Drop genera that were not present
  group_by(GENUS) %>% 
  mutate(GENUS.TOTAL = sum(COUNT)) %>% 
  ungroup() %>% 
  filter(GENUS.TOTAL > 0) %>% 
  # Proportion
  group_by(PLOT) %>% 
  mutate(PLOT.COUNT = sum(COUNT)) %>% 
  ungroup() %>% 
  group_by(GENUS, PLOT) %>% 
  mutate(PROP = sum(COUNT) / PLOT.COUNT) %>% 
  ungroup() %>% 
  # Hellinger transform abundances for nmds and titan
  mutate(HELLINGER = sqrt(PROP))

# Remove sequences not classified at the genus level
soil.genus.tbl <- soil.genus.no.outliers.tmp.tbl %>% 
  separate(col = GENUS, into = c("UNCLASSIFIED", "UNCLASSIFIED.2"), remove = FALSE, sep = " ") %>% 
  filter(UNCLASSIFIED != "unclassified" & is.na(UNCLASSIFIED.2)) %>% 
  select(-UNCLASSIFIED, -UNCLASSIFIED.2)

# Transform roots genera
roots.genus.no.outliers.tmp.tbl <- genus.all.tbl %>% 
  # Drop plots that were SOM outliers
  filter(PLOT != "Plot_6" & PLOT != "Plot_19" & PLOT != "Plot_43" & PLOT != "Plot_60") %>% 
  filter(SUBSTRATE == "decomp.roots") %>% 
  # Drop genera that were not present
  group_by(GENUS) %>% 
  mutate(GENUS.TOTAL = sum(COUNT)) %>% 
  ungroup() %>% 
  filter(GENUS.TOTAL > 0) %>% 
  # Proportion
  group_by(PLOT) %>% 
  mutate(PLOT.COUNT = sum(COUNT)) %>% 
  ungroup() %>% 
  group_by(GENUS, PLOT) %>% 
  mutate(PROP = sum(COUNT) / PLOT.COUNT) %>% 
  ungroup() %>% 
  # Hellinger transform abundances for nmds and titan
  mutate(HELLINGER = sqrt(PROP))

# Remove sequences not classified at the genus level
roots.genus.tbl <- roots.genus.no.outliers.tmp.tbl %>% 
  separate(col = GENUS, into = c("UNCLASSIFIED", "UNCLASSIFIED.2"), remove = FALSE, sep = " ") %>% 
  filter(UNCLASSIFIED != "unclassified" & is.na(UNCLASSIFIED.2)) %>% 
  select(-UNCLASSIFIED, -UNCLASSIFIED.2)

#### 6. Trim to retain only abundant genera ####

# Soil, retain only abundant genera
soil.genus.trim.tbl <- soil.genus.tbl %>% 
  mutate(PA = ifelse(COUNT > 0, 1, 0)) %>% 
  group_by(GENUS) %>% 
  inner_join(., soil.genus.perc.tbl, by = "GENUS") %>%
  mutate(OCCURRENCE = sum(PA)) %>% 
  ungroup() %>% 
  # Trim out genera that are present in fewer than 5 plots and < 0.01% of reads
  filter(GENUS.PERC >= 0.1 & OCCURRENCE >= 5)

# Roots, retain only abundant genera
roots.genus.trim.tbl <- roots.genus.tbl %>% 
  mutate(PA = ifelse(COUNT > 0, 1, 0)) %>% 
  group_by(GENUS) %>% 
  inner_join(., roots.genus.perc.tbl, by = "GENUS") %>%
  mutate(OCCURRENCE = sum(PA)) %>% 
  ungroup() %>% 
  # Trim out genera that are present in fewer than 5 plots and < 0.01% of reads
  filter(GENUS.PERC >= 0.1 & OCCURRENCE >= 5)
