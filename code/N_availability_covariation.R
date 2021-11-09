#### 1. Set-up ####

# Load libraries
library(tidyverse)
library(mgcv)
library(viridis)
library(itsadug)
library(ggpubr)

# Functions and other data
source("code/functions.R")
source("code/compile_data.R")

#### 2. N min data, spring 2018 ####

# Read in
Nmin.aug.2018.mean.tbl <- read_csv("data/complete.min.august.manistee.csv") %>% 
  rename(TREE = "Tree ID", 
         INCUBATION = "Time", 
         REP = "Rep", 
         NITRATE = "Nitrate", 
         AMMONIUM = "Ammonia") %>% 
  select(INCUBATION, TREE, REP, NITRATE, AMMONIUM) %>% 
  mutate(NITRATE = ifelse(NITRATE < 0, 0, NITRATE), 
         AMMONIUM = ifelse(AMMONIUM < 0, 0, AMMONIUM)) %>% 
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
  select(TREE, N.MIN) %>% 
  separate(TREE, into = c("STAND", "TREE"), sep = "_") %>% 
  mutate(STAND = paste("Stand", STAND, sep = "_")) %>% 
  filter(!is.na(N.MIN)) %>% 
  group_by(STAND) %>% 
  summarise(NMIN.MEAN = mean(N.MIN), 
            NMIN.SE = se(N.MIN)) %>% 
  ungroup()

#### 3. Root C/N data ####

# Root CN
root.CN.tbl <- read_csv(file = "data/Manistee_env_data/roots/Zak_WA_pre_roots_3-26-19_1.csv") %>% 
  filter(!is.na(Sample) & Sample != "Blank" & Sample != "Standard") %>% 
  separate(Sample, into = c("STAND", "REP1", "REP2"), sep = "_") %>% 
  mutate(STAND = paste("Stand", STAND, sep = "_")) %>% 
  group_by(STAND, REP1) %>% 
  summarise(NITROGEN = mean(Nitrogen), 
            CARBON = mean(Carbon)) %>% 
  ungroup() %>% 
  mutate(ROOT.CN = CARBON / NITROGEN) %>% 
  group_by(STAND) %>% 
  summarise(ROOTCN.MEAN = mean(ROOT.CN), 
            ROOTCN.SE = se(ROOT.CN)) %>% 
  ungroup() %>% 
  inner_join(Nmin.aug.2018.mean.tbl, ., by = "STAND")

# LM
root.CN.lm <- lm(ROOTCN.MEAN ~ NMIN.MEAN, data = root.CN.tbl)
root.CN.lm.summary <- summary(root.CN.lm)
root.CN.lm.Pval = format(round(root.CN.lm.summary$coefficients[2, 4], 3), scientific = FALSE)
root.CN.lm.P <- tibble(PVAL.LABEL = paste("atop(R[adj.]^2 == ", 
                                          round(root.CN.lm.summary$adj.r.squared, 3), 
                                          ", ", 
                                          "italic(P) == ", 
                                          root.CN.lm.Pval, 
                                          ")", 
                                          sep = ""), 
                    XPOS = 0.65, 
                    YPOS = 55)

# Figure
root.CN.plot <- ggplot() + 
  
  # Vertical error bars
  geom_errorbar(data = root.CN.tbl, 
                aes(x = NMIN.MEAN, y = ROOTCN.MEAN, ymin = ROOTCN.MEAN - ROOTCN.SE, ymax = ROOTCN.MEAN + ROOTCN.SE), 
                colour = "#7B02A8FF") + 
  
  # Horizontal error bars
  geom_errorbar(data = root.CN.tbl, 
                aes(x = NMIN.MEAN, y = ROOTCN.MEAN, xmin = NMIN.MEAN - NMIN.SE, xmax = NMIN.MEAN + NMIN.SE), 
                colour = "#7B02A8FF") + 
  
  # Points
  geom_point(data = root.CN.tbl, 
             aes(x = NMIN.MEAN, y = ROOTCN.MEAN), 
             colour = "#7B02A8FF", 
             size = 1.5) + 
  
  # Pvalues
  geom_text(data = root.CN.lm.P, 
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL), 
            parse = TRUE, 
            stat = "identity", 
            size = 3) + 
  
  # Line
  geom_smooth(data = root.CN.tbl, 
              aes(x = NMIN.MEAN, y = ROOTCN.MEAN), 
              formula = y ~ x,  
              method = "lm", 
              colour = "#7B02A8FF", 
              fill = "#7B02A8FF") + 
  
  # Titles
  labs(title = "(a)", 
       x = bquote('Inorganic N ('*mu*g~N%.%g^-1%.%d^-1*')'), 
       y = "Fine root C/N") + 
  
  # Format
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 12, hjust = -0.1), 
        axis.title.x = element_text(colour = "black", size = 7), 
        axis.title.y = element_text(colour = "black", size = 7),  
        # x axis text
        axis.text.x = element_text(colour = "black", size = 6), 
        axis.text.y = element_text(colour = "black", size = 6))

#### 4. Soil CN and Nmin ####

# GAM
soilCN.Nmin.gam <- gam(SOIL.CN ~ s(N.MIN, k = -1), 
                         data = soil.compiled.data.tbl)
soilCN.Nmin.gam.summary <- summary(soilCN.Nmin.gam)
soilCN.Nmin.gam.P <- tibble(PVAL.LABEL = paste("atop(R[adj.]^2 == ", 
                                          round(soilCN.Nmin.gam.summary$r.sq, 3), 
                                          ", ", 
                                          "italic(P) < 0.001)", 
                                          sep = ""), 
                       XPOS = 0.65, 
                       YPOS = 30)

# Figure
soil.CN.plot <- ggplot() + 
  
  # Points
  geom_point(data = soil.compiled.data.tbl, 
             aes(x = N.MIN, y = SOIL.CN), 
             colour = "#FEBE2AFF", 
             size = 1.5) + 
  
  # Pvalues
  geom_text(data = soilCN.Nmin.gam.P, 
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL), 
            parse = TRUE, 
            stat = "identity", 
            size = 3) + 
  
  # Line
  geom_smooth(data = soil.compiled.data.tbl, 
              aes(x = N.MIN, y = SOIL.CN), 
              formula = y ~ s(x),  
              method = "gam", 
              colour = "#FEBE2AFF", 
              fill = "#FEBE2AFF") + 
  
  # Titles
  labs(title = "(b)", 
       x = bquote('Inorganic N ('*mu*g~N%.%g^-1%.%d^-1*')'), 
       y = "Soil C/N") + 
  
  # Format
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 12, hjust = -0.1), 
        axis.title.x = element_text(colour = "black", size = 7), 
        axis.title.y = element_text(colour = "black", size = 7),  
        # x axis text
        axis.text.x = element_text(colour = "black", size = 6), 
        axis.text.y = element_text(colour = "black", size = 6))


# Figure S2
FigureS2 <- ggarrange(root.CN.plot, 
                      soil.CN.plot, 
                      ncol = 2, 
                      nrow = 1, 
                      align = "hv")

# Save Figure S2
ggsave(filename = "figures/FigureS2.png",
       FigureS2,
       device = "png",
       dpi = 600,
       width = 6.5,
       height = 3.25,
       units = "in")

