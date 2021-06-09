#### 1. Set-up ####

# Load libraries
library(tidyverse)
library(TITAN2)

source("code/functions.R")

# Read in data
soil.compiled.data.tbl <- readRDS(file = "data/working/soil.compiled.data.tbl.RData")
source("code/calculate_genus_abundances.R")

#### 2. TITAN ####

# N mineralization gradient
env.N.min.df <- soil.compiled.data.tbl %>% 
  select(PLOT, N.MIN) %>% 
  arrange(PLOT) %>% 
  column_to_rownames(var = "PLOT") %>% 
  as.data.frame(.)

# Roots, genus
roots.genus.hlr.trim.titan.in.df <- roots.genus.trim.tbl %>% 
  select(PLOT, GENUS, HELLINGER) %>% 
  arrange(GENUS) %>% 
  pivot_wider(id_cols = PLOT, names_from = "GENUS", values_from = "HELLINGER") %>% 
  arrange(PLOT) %>% 
  column_to_rownames(var = "PLOT") %>% 
  as.data.frame(.)

# Soil, genus
soil.genus.hlr.trim.titan.in.df <- soil.genus.trim.tbl %>% 
  select(PLOT, GENUS, HELLINGER) %>% 
  arrange(GENUS) %>% 
  pivot_wider(id_cols = PLOT, names_from = "GENUS", values_from = "HELLINGER") %>% 
  arrange(PLOT) %>% 
  column_to_rownames(var = "PLOT") %>% 
  as.data.frame(.)

# Roots, TITAN
# roots.genus.hlr.trim.titan <- titan(env.N.min.df, roots.genus.hlr.trim.titan.in.df,
#                                     minSplt = 5,
#                                     numPerm = 1000,
#                                     boot = TRUE,
#                                     nBoot = 1000,
#                                     imax = FALSE,
#                                     ivTot = FALSE,
#                                     pur.cut = 0.95,
#                                     rel.cut = 0.95,
#                                     ncpus = 4,
#                                     memory = TRUE)

# Soil, TITAN
# soil.genus.hlr.trim.titan <- titan(env.N.min.df, soil.genus.hlr.trim.titan.in.df,
#                                    minSplt = 5,
#                                    numPerm = 1000,
#                                    boot = TRUE,
#                                    nBoot = 1000,
#                                    imax = FALSE,
#                                    ivTot = FALSE,
#                                    pur.cut = 0.95,
#                                    rel.cut = 0.95,
#                                    ncpus = 4,
#                                    memory = TRUE)

# Save data to reduce computational time later
# saveRDS(roots.genus.hlr.trim.titan, file = "data/working/roots.genus.hlr.trim.titan.rData")
# saveRDS(soil.genus.hlr.trim.titan, file = "data/working/soil.genus.hlr.trim.titan.rData")
