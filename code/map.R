#### 1. Set-up ####

# Load libraries
library(tidyverse)
library(data.table)
library(viridis)
library(ggpubr)
library(usmap)
library(rgdal)
library(ggrepel)

#### 4. Map of study sites ####

stand.locations.df <- as_tibble(fread("data/Manistee_env_data/spatial/plot.locations.updated.2019.10.19.and.20.V2.txt", sep = "\t", header = TRUE)) %>% 
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
  mutate(STAND = paste("Stand", STAND, sep = "_")) %>% 
  # Decimal degrees
  mutate(STAND.LAT = STAND.LAT.DEG + (STAND.LAT.MIN / 60), 
         STAND.LON = -1 * (STAND.LON.DEG + (STAND.LON.MIN / 60))) %>% 
  distinct(STAND, .keep_all = TRUE) %>% 
  select(STAND, STAND.LON, STAND.LAT) %>% 
  as.data.frame(.)

row.names(stand.locations.df) <- stand.locations.df$STAND
stand.locations.df <- stand.locations.df[2 : ncol(stand.locations.df)]

stand.locations.data <- usmap_transform(stand.locations.df)
stand.locations.data$STAND <- row.names(stand.locations.df)
stand.locations.data <- merge(stand.locations.data, soil.CN.stand.means.tbl, by = "STAND")

mi.map.plot <- plot_usmap(include = "MI") + 
  geom_point(data = stand.locations.data, 
             aes(x = STAND.LON.1, y = STAND.LAT.1), 
             size = 1, 
             alpha = 0.75) + 
  
  geom_text_repel(data = stand.locations.data, 
            aes(x = STAND.LON.1, y = STAND.LAT.1, label = STAND), 
            size = 2, 
            colour = "blue")

ggsave("figures/FigureS1.png", mi.map.plot, device = "png", width = 5, height = 5, units = "in", dpi = 600)
ggsave("manistee.map.fig.pdf", mi.map.plot, device = "pdf", width = 5, height = 5, units = "in")
