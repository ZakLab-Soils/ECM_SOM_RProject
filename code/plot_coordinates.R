#### 1. Set-up ####

# Load libraries
library(tidyverse)
library(data.table)
library(geosphere)

# Functions
source("code/functions.R")

#### 2. Calculate coordinates for each plot

# Read in sample information
plot.locations.tbl <- as_tibble(fread("data/Manistee_env_data/spatial/plot.locations.updated.2019.10.19.and.20.V2.txt", sep = "\t", header = TRUE)) %>% 
  dplyr::select(Site, Plot, `Bag angle`, `Bag distance (m)`, `Site latitude degrees (N)`, `Site latitude minutes (N)`, `Site longitude degrees (W)`, `Site longitude minutes (W)`) %>% 
  dplyr::rename(STAND = "Site", 
                PLOT = "Plot", 
                PLOT.BEARING = "Bag angle", 
                PLOT.DIST = "Bag distance (m)", 
                STAND.LAT.DEG = "Site latitude degrees (N)", 
                STAND.LAT.MIN = "Site latitude minutes (N)", 
                STAND.LON.DEG = "Site longitude degrees (W)", 
                STAND.LON.MIN = "Site longitude minutes (W)") %>% 
  dplyr::filter(PLOT <= 72) %>% 
  mutate(PLOT = paste("Plot", PLOT, sep = "_")) %>% 
  # Decimal degrees
  mutate(STAND.LAT = STAND.LAT.DEG + (STAND.LAT.MIN / 60), 
         STAND.LON = -1 * (STAND.LON.DEG + (STAND.LON.MIN / 60)))

# Format matrix of stand latitude and longitude
plot.locations.df <- as.data.frame(plot.locations.tbl)
row.names(plot.locations.df) <- plot.locations.df$PLOT
plot.locations.mat <- as.matrix(plot.locations.df[c("STAND.LON", "STAND.LAT")])

# Format plot bearings and distances
plot.bearing.mat <- as.matrix(plot.locations.df["PLOT.BEARING"])
plot.dist.mat <- as.matrix(plot.locations.df["PLOT.DIST"])

# Calculate coordinates for each plot
plot.coordinates.df <- as.data.frame(destPoint(plot.locations.mat, plot.bearing.mat, plot.dist.mat))
row.names(plot.coordinates.df) <- row.names(plot.locations.df)
plot.coordinates.tbl <- as_tibble(plot.coordinates.df, rownames = "PLOT")
