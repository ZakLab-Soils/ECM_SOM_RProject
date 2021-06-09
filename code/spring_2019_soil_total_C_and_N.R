#### 1. Set-up ####

# This script calculates soil total C and N for soil collected in spring 2019

# Load libraries
library(tidyverse)
library(data.table)

#### 2. Read in and format LECO data ####

# Read in LECO data, run 1
soil.leco.1.tbl <- as_tibble(fread("data/Manistee_env_data/leco/WA_LECO_run_from_19_07_25.csv", sep = ",", header = FALSE)) %>% 
  # Rename columns
  rename(RUN.ID = "V1", 
         SAMPLE.ID = "V2", 
         SAMPLE.MASS = "V3", 
         TOTAL.N = "V4", 
         TOTAL.C = "V5", 
         RUN.DATE.TIME = "V6") %>% 
  # Trim
  dplyr::select(SAMPLE.ID, TOTAL.N, TOTAL.C) %>% 
  # Filter
  filter(SAMPLE.ID != "Blank" & SAMPLE.ID != "Soil standard")

# Read in LECO data, run 2
soil.leco.2.tbl <- as_tibble(fread("data/Manistee_env_data/leco/WA_LECO_run_from_19_07_29.csv", sep = ",", header = FALSE)) %>% 
  # Rename columns
  rename(RUN.ID = "V1", 
         SAMPLE.ID = "V2", 
         SAMPLE.MASS = "V3", 
         TOTAL.N = "V4", 
         TOTAL.C = "V5", 
         RUN.DATE.TIME = "V6") %>% 
  # Trim
  dplyr::select(SAMPLE.ID, TOTAL.N, TOTAL.C) %>% 
  # Filter
  filter(SAMPLE.ID != "Blank" & SAMPLE.ID != "Soil standard")

#### 3. Combine LECO data and calculate masses ####

# Combined
soil.leco.tbl <- bind_rows(soil.leco.1.tbl, soil.leco.2.tbl) %>% 
  # Format sample ID
  separate(SAMPLE.ID, c("SUBSTRATE", "STAND", "PLOT", "REPLICATE"), sep = "_", convert = TRUE) %>% 
  # Update stand and plot names
  mutate(STAND = paste("Stand", STAND, sep = "_"), 
         PLOT = paste("Plot", PLOT, sep = "_")) %>% 
  # Group
  group_by(STAND, PLOT) %>% 
  # Calculate means
  summarise(TOTAL.C = mean(TOTAL.C), 
            TOTAL.N = mean(TOTAL.N)) %>% 
  # Ungroup
  ungroup() %>% 
  # Calculate C:N
  mutate(SOIL.CN = TOTAL.C * (1 / TOTAL.N)) %>% 
  # Calculate specific C and N
  mutate(SPEC.TOTAL.C = ((TOTAL.C * 1/100) * (1000)), 
         SPEC.TOTAL.N = ((TOTAL.N * 1/100) * (1000))) %>% 
  select(PLOT, STAND, SPEC.TOTAL.C, SPEC.TOTAL.N, SOIL.CN)
