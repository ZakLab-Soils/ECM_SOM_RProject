#### 1. Set-up ####

# Load libraries
library(tidyverse)
library(ggpubr)
library(viridis)

source("code/functions.R")

# Format data
SOM.data.tbl <-readRDS("data/working/soil.compiled.data.tbl.RData") %>% 
  select(PLOT, AROMATIC, LIGNIN, LIPID, N.BEARING, PHENOL, POLYSACCHARIDE, PROTEIN, UNKNOWN.ORIGIN) %>% 
  pivot_longer(cols = AROMATIC : UNKNOWN.ORIGIN, 
               names_to = "SOURCE", 
               values_to = "PROP")

# Means
SOM.data.mean.tbl <- SOM.data.tbl %>% 
  group_by(SOURCE) %>% 
  summarise(PROP.MEAN = mean(PROP), 
            PROP.SE = se(PROP)) %>% 
  ungroup() %>% 
  mutate(LSE = PROP.MEAN - PROP.SE, 
         USE = PROP.MEAN + PROP.SE) %>% 
  arrange(desc(PROP.MEAN)) %>% 
  mutate(SOURCE.LABEL = NA, 
         SOURCE.LABEL = ifelse(SOURCE == "AROMATIC", "Aromatic", SOURCE.LABEL), 
         SOURCE.LABEL = ifelse(SOURCE == "LIGNIN", "Lignin", SOURCE.LABEL), 
         SOURCE.LABEL = ifelse(SOURCE == "LIPID", "Lipids", SOURCE.LABEL), 
         SOURCE.LABEL = ifelse(SOURCE == "N.BEARING", "N-bearing", SOURCE.LABEL), 
         SOURCE.LABEL = ifelse(SOURCE == "PHENOL", "Phenol", SOURCE.LABEL), 
         SOURCE.LABEL = ifelse(SOURCE == "POLYSACCHARIDE", "Polysaccharides", SOURCE.LABEL), 
         SOURCE.LABEL = ifelse(SOURCE == "PROTEIN", "Proteins", SOURCE.LABEL), 
         SOURCE.LABEL = ifelse(SOURCE == "UNKNOWN.ORIGIN", "Unknown origin", SOURCE.LABEL)) %>% 
  mutate(SOURCE.LABEL = factor(SOURCE.LABEL, levels = unique(SOURCE.LABEL)))

SOM.data.tbl <- SOM.data.tbl %>% 
  inner_join(., SOM.data.mean.tbl, by = "SOURCE") %>% 
  arrange(desc(PROP.MEAN)) %>% 
  mutate(SOURCE.LABEL = NA, 
         SOURCE.LABEL = ifelse(SOURCE == "AROMATIC", "Aromatic", SOURCE.LABEL), 
         SOURCE.LABEL = ifelse(SOURCE == "LIGNIN", "Lignin", SOURCE.LABEL), 
         SOURCE.LABEL = ifelse(SOURCE == "LIPID", "Lipids", SOURCE.LABEL), 
         SOURCE.LABEL = ifelse(SOURCE == "N.BEARING", "N-bearing", SOURCE.LABEL), 
         SOURCE.LABEL = ifelse(SOURCE == "PHENOL", "Phenol", SOURCE.LABEL), 
         SOURCE.LABEL = ifelse(SOURCE == "POLYSACCHARIDE", "Polysaccharides", SOURCE.LABEL), 
         SOURCE.LABEL = ifelse(SOURCE == "PROTEIN", "Proteins", SOURCE.LABEL), 
         SOURCE.LABEL = ifelse(SOURCE == "UNKNOWN.ORIGIN", "Unknown origin", SOURCE.LABEL)) %>% 
  mutate(SOURCE.LABEL = factor(SOURCE.LABEL, levels = unique(SOURCE.LABEL)))

# Figure
SOM.barplot <- ggplot() + 
  
  # Add plot level data points
  geom_point(data = SOM.data.tbl, 
             aes(x = SOURCE.LABEL, y = PROP, colour = SOURCE.LABEL), 
             position = position_jitter(0.2), 
             size = 1, 
             alpha = 0.5) + 
  
  # Bars
  geom_bar(data = SOM.data.mean.tbl, 
           aes(x = SOURCE.LABEL, y = PROP.MEAN, fill = SOURCE.LABEL, colour = SOURCE.LABEL), 
           stat = "identity", 
           alpha = 0.5) + 
  
  # Add error bars
  geom_errorbar(data = SOM.data.mean.tbl, 
                aes(x = SOURCE.LABEL, ymin = LSE, ymax = USE, colour = SOURCE.LABEL), 
                width = 0.2) + 
  
  # Format color for bar fills and colors for bar outlines, error bars, and plot-level data points
  scale_fill_viridis(discrete = TRUE, 
                     guide = FALSE, 
                     option = "plasma", 
                     begin = 0, 
                     end = 0.9) + 
  
  scale_color_viridis(discrete = TRUE, 
                      guide = FALSE, 
                      option = "plasma", 
                      begin = 0, 
                      end = 0.9) + 
  
  # Remove gap between x axis and bars
  scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
  
  # Titles
  labs(x = "Compound class", 
       y = "Relative abundance") + 
  
  # Format panel and strip
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(colour = "black", size = 10), 
        axis.title.y = element_text(colour = "black", size = 10), 
        # x axis text
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, size = 6), 
        axis.text.y = element_text(colour = "black", size = 8), 
        strip.background = element_blank())

# Save figure
ggsave(filename = "figures/FigureS8.png",
       plot = SOM.barplot, 
       device = "png",
       dpi = 400,
       width = 6.5,
       height = 4,
       units = "in")
