#### 1. Set-up ####

# Load libraries
library(tidyverse)
library(mgcv)
library(viridis)

#### 2. N min data ####

# Read in
Nmin.aug.2018.tbl <- read_csv("data/complete.min.august.manistee.csv") %>% 
  rename(TREE = "Tree ID", 
         INCUBATION = "Time", 
         REP = "Rep", 
         NITRATE = "Nitrate", 
         AMMONIUM = "Ammonia") %>% 
  select(INCUBATION, TREE, REP, NITRATE, AMMONIUM) %>% 
  group_by(INCUBATION, TREE) %>% 
  summarise(NITRATE = sum(NITRATE), 
            AMMONIUM = sum(AMMONIUM)) %>% 
  ungroup() %>% 
  mutate(TOTAL = NITRATE + AMMONIUM) %>% 
  select(INCUBATION, TREE, TOTAL) %>% 
  pivot_wider(id_cols = TREE, 
              names_from = "INCUBATION", 
              values_from = "TOTAL") %>% 
  mutate(N.MIN = (Post - Pre) / 14) %>% 
  select(TREE, N.MIN)
  
#### 3. Root biomass ####

# Read in
root.biomass.tbl <- read_tsv("data/root.picking.data.complete.txt") %>% 
  # Add site, tree, and core columns
  separate(SAMPLE, c("STAND", "TREE", "CORE"), sep = "_") %>% 
  # Remove unnecessary columns
  dplyr::select(STAND, TREE, CORE, DRY.ROOT.MASS) %>% 
  # Drop rows with no mass measurements
  dplyr::filter(!(is.na(DRY.ROOT.MASS))) %>% 
  # Group by site and tree
  dplyr::group_by(STAND, TREE) %>% 
  # Sum by tree
  dplyr::summarise(DRY.ROOT.MASS = sum(DRY.ROOT.MASS)) %>% 
  # Ungroup
  dplyr::ungroup() %>% 
  unite(TREE, STAND : TREE, sep = "_", remove = TRUE) %>% 
  inner_join(., Nmin.aug.2018.tbl, by = "TREE")

# GAM
root.lm <- lm(DRY.ROOT.MASS ~ N.MIN, data = root.biomass.tbl)
root.lm.summary<- summary(root.lm)
root.lm.P <- tibble(PVAL.LABEL = paste("atop(R[adj.]^2 == ", round(root.lm.summary$adj.r.squared, 3), ", ", "italic(P) < 0.001)", sep = ""), 
                     XPOS = 0.9, 
                     YPOS = 22)

# Plot
root.plot <- ggplot() + 
  
  # Points
  geom_point(data = root.biomass.tbl, 
             aes(x = N.MIN, y = DRY.ROOT.MASS), 
             colour = "#7B02A8FF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Pvalues
  geom_text(data = root.lm.P, 
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL), 
            parse = TRUE, 
            stat = "identity", 
            size = 3) + 
  
  # Line
  geom_smooth(data = root.biomass.tbl, 
              aes(x = N.MIN, y = DRY.ROOT.MASS), 
              formula = y ~ x,  
              method = "lm", 
              colour = "#7B02A8FF", 
              fill = "#7B02A8FF") + 
  
  # Titles
  labs(title = NULL, 
       x = bquote('Inorganic N ('*mu*g~N%.%g^-1%.%d^-1*')'), 
       y = "Fine root biomass (g)") + 
  
  # Format
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(colour = "black", size = 7), 
        axis.title.y = element_text(colour = "black", size = 7),  
        # x axis text
        axis.text.x = element_text(colour = "black", size = 6), 
        axis.text.y = element_text(colour = "black", size = 6))

# Save Figure S6
ggsave(filename = "figures/FigureS6.png",
       root.plot,
       device = "png",
       dpi = 600,
       width = 6.5,
       height = 4,
       units = "in")
