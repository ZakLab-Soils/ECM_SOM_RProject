#### 1. Set-up ####

# Load libraries
library(tidyverse)
library(data.table)

#### 2. Read in and format data ####

# Read in sample information
py.gcms.info.tbl <- as_tibble(fread("data/Manistee_env_data/pyGCMS/py_gcms_pre_decay_root_som_sample_info.txt", sep = "\t", header = TRUE))
# Read in data
py.gcms.som.tbl <- as_tibble(fread("data/Manistee_env_data/pyGCMS/py_gcms_pre_decay_root_som_data.txt", sep = "\t", header = TRUE)) %>% 
  # Update column names
  rename("COMPOUND" = Compound, "TYPE" = Type, "SOURCE" = Source) %>% 
  # Update lignin and phenols
  mutate(SOURCE = ifelse(SOURCE == "Lignin+TMAH", "Lignin", SOURCE), 
         SOURCE = ifelse(SOURCE == "Phenol+TMAH", "Phenol", SOURCE)) %>% 
  # Long format
  pivot_longer(-c(COMPOUND, TYPE, SOURCE), names_to = "SAMPLE.ID", values_to = "PROPORTION") %>% 
  group_by(SAMPLE.ID, TYPE, SOURCE, COMPOUND) %>% 
  summarize(PROPORTION = sum(PROPORTION)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(TYPE, SOURCE, COMPOUND), names_from = "SAMPLE.ID", values_from = "PROPORTION") %>% 
  mutate(COMPOUND.ID = paste("Compound", row_number(), sep = "_")) %>% 
  select(TYPE, SOURCE, COMPOUND, COMPOUND.ID, everything()) %>% 
  pivot_longer(-c(COMPOUND.ID, COMPOUND, TYPE, SOURCE), names_to = "SAMPLE.ID", values_to = "PROPORTION") %>% 
  # Add sample information
  inner_join(py.gcms.info.tbl, ., by = "SAMPLE.ID") %>% 
  # Update stand and plot columns
  mutate(STAND = paste("Stand", STAND, sep = "_"), 
         PLOT = ifelse(is.na(PLOT), PLOT, paste("Plot", PLOT, sep = "_")))

# Get compound information
py.gcms.compound.info.tbl <- py.gcms.som.tbl %>% 
  select(COMPOUND.ID, COMPOUND, TYPE, SOURCE) %>% 
  distinct()

#### 3. Format soil organic matter data ####

# Calculate SOM source level proportion
py.gcms.som.source.tbl <- py.gcms.som.tbl %>% 
  # Group
  group_by(STAND, PLOT, SOURCE) %>% 
  # Calculate source totals
  summarise(RELABUND = sum(PROPORTION)) %>% 
  # Ungroup
  ungroup() %>% 
  # Wide format
  pivot_wider(names_from = "SOURCE", values_from = "RELABUND") %>% 
  # Update column names
  rename_all(., toupper) %>% 
  # Update column names
  rename_all(~ str_replace_all(., '\\-', '\\.')) %>% 
  # Update column names
  rename_all(~ str_replace_all(., ' ', '\\.'))
