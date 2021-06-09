#### 1. Set-up ####

# This script calculates net N mineralization for soil collected in spring 2019

# Load libraries
library(tidyverse)

# Source files
source("code/functions.R")
source("code/spring_2019_net_N_mineralization.R")
source("code/spring_2019_soil_total_C_and_N.R")
source("code/som_pre_root_py_gcms.R")
source("code/plot_coordinates.R")

# Functional group data
soil.fungal.functional.groups.tbl <- readRDS(file = "data/working/soil.fungal.functional.groups.tbl.RData")
roots.fungal.functional.groups.tbl <- readRDS(file = "data/working/roots.fungal.functional.groups.tbl.RData")

#### 2. Compile data for soil ####

# N min data
soil.compiled.data.tbl <- N.min.tbl %>% 
  select(-STAND) %>% 
  # Add soil C and N
  inner_join(., soil.leco.tbl, by = "PLOT") %>% 
  select(-STAND) %>% 
  # Add py-GC/MS data
  inner_join(., py.gcms.som.source.tbl, by = "PLOT") %>% 
  # Add plot coordinates
  inner_join(., plot.coordinates.tbl, by = "PLOT") %>% 
  rename(LAT = "lat", 
         LON = "lon") %>% 
  # Add functional group abundances
  inner_join(., soil.fungal.functional.groups.tbl, by = "PLOT") %>% 
  select(PLOT, STAND, everything())

#### 3. Compile data for roots ####

# N min data
roots.compiled.data.tbl <- N.min.tbl %>% 
  select(-STAND) %>% 
  # Add soil C and N
  inner_join(., soil.leco.tbl, by = "PLOT") %>% 
  select(-STAND) %>% 
  # Add py-GC/MS data
  inner_join(., py.gcms.som.source.tbl, by = "PLOT") %>% 
  # Add plot coordinates
  inner_join(., plot.coordinates.tbl, by = "PLOT") %>% 
  rename(LAT = "lat", 
         LON = "lon") %>% 
  # Add functional group abundances
  inner_join(., roots.fungal.functional.groups.tbl, by = "PLOT") %>% 
  select(PLOT, STAND, everything())

#### 4. Save compiled data ####

saveRDS(soil.compiled.data.tbl, 
        file = "data/working/soil.compiled.data.tbl.RData")
saveRDS(roots.compiled.data.tbl, 
        file = "data/working/roots.compiled.data.tbl.RData")

