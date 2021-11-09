#### 1. Set-up ####

# This script calculates fine root mass loss

# Load libraries
library(tidyverse)
library(mgcv)
library(itsadug)
library(ggpubr)
library(viridis)

source("code/functions.R")
source("code/environmental_summary.R")

# Read in data
soil.compiled.data.tbl <- readRDS(file = "data/working/soil.compiled.data.tbl.RData")

# Read in initial moisture data
decomp.root.initial.tbl <- readr::read_tsv(file = "data/root_mass_loss/manistee_roots_decomposed_masses_july_2020.txt") %>% 
  
  # Rename columns
  dplyr::rename(PLOT = "Plot", 
                TIN.MASS = "Tin (g)", 
                FRESH.TIN.MASS = "Fresh + tin (g)", 
                DRY.TIN.MASS = "Dry + tin (g)") %>% 
  
  # Calculate fresh and dry root mass, add measurement ID
  dplyr::mutate(FRESH.ROOT.MASS = FRESH.TIN.MASS - TIN.MASS, 
                DRY.ROOT.MASS = DRY.TIN.MASS - TIN.MASS, 
                MEASUREMENT = "Initial") %>% 
  
  # Trim table, remove test plots, and drop lost sample
  dplyr::select(PLOT, 
                MEASUREMENT, 
                FRESH.ROOT.MASS, 
                DRY.ROOT.MASS) %>% 
  dplyr::filter(PLOT != "73_A" & 
                  PLOT != "73_B" & 
                  PLOT != "74_A" & 
                  PLOT != "74_B") %>% 
  dplyr::mutate(PLOT = paste("Plot", PLOT, sep = "_"), 
                FRESH.ROOT.MASS = ifelse(PLOT == "Plot_57", 0, FRESH.ROOT.MASS), 
                DRY.ROOT.MASS = ifelse(PLOT == "Plot_57", 0, DRY.ROOT.MASS))

# Read in additional moisture data
decomp.root.moisture.tbl <- readr::read_tsv(file = "data/root_mass_loss/decomposed_roots_additional_moisture_data.txt") %>% 
  
  # Calculate fresh and dry root mass, add measurement ID
  dplyr::mutate(FRESH.ROOT.MASS = FRESH.ROOTS.TIN - TIN, 
                DRY.ROOT.MASS = DRY.ROOTS.TIN - TIN, 
                MEASUREMENT = "Additional") %>% 
  
  # Trim table and format plots
  dplyr::select(PLOT, 
                MEASUREMENT, 
                FRESH.ROOT.MASS, 
                DRY.ROOT.MASS) %>% 
  dplyr::mutate(PLOT = paste("Plot", PLOT, sep = "_")) %>% 
  
  # Add initial moisture and calculate combined masses
  dplyr::bind_rows(decomp.root.initial.tbl, .) %>% 
  dplyr::group_by(PLOT) %>% 
  dplyr::summarise(FRESH.ROOT.COMBINED.MASS = sum(FRESH.ROOT.MASS), 
                   DRY.ROOT.COMBINED.MASS = sum(DRY.ROOT.MASS)) %>% 
  dplyr::ungroup() %>% 
  
  # Calculate moisture content
  dplyr::mutate(MOISTURE.CONTENT = (FRESH.ROOT.COMBINED.MASS - DRY.ROOT.COMBINED.MASS) / FRESH.ROOT.COMBINED.MASS) %>% 
  dplyr::select(PLOT, 
                MOISTURE.CONTENT)

# Read in pre-incubated ash data
pre.incubated.ash.tbl <- readr::read_tsv(file = "data/root_mass_loss/pre_incubated_root_ash_data.txt") %>% 
  
  # Split sample ID into replicate and site, format stand
  tidyr::separate(col = SAMPLE, 
                  into = c("SITE", "REPL")) %>% 
  dplyr::mutate(SITE = paste("Site", SITE, sep = "_")) %>% 
  
  # Calculate ash-free proportion as site means
  dplyr::mutate(ASH.FREE.PROP = ((ROOT.CRUCIBLE - EMPTY.CRUCIBLE) - (ROOT.CRUCIBLE.ASHED - EMPTY.CRUCIBLE)) / 
                  (ROOT.CRUCIBLE - EMPTY.CRUCIBLE)) %>% 
  dplyr::select(SITE, 
                REPL, 
                ASH.FREE.PROP) %>% 
  
  # Format sample ID
  tidyr::unite(SITE : REPL, 
               col = SAMPLE.ID, 
               sep = "_", 
               remove = FALSE) %>% 
  dplyr::mutate(SAMPLE.ID = gsub("Site_", "", SAMPLE.ID))

# Calculate site means, pre-incubated roots
pre.incubated.ash.mean.tbl <- pre.incubated.ash.tbl %>% 
  dplyr::group_by(SITE) %>% 
  dplyr::summarise(MEAN.ASH.FREE.PROP = mean(ASH.FREE.PROP)) %>% 
  dplyr::ungroup()

# Read in decomposing root ash data
decomp.ash.tbl <- readr::read_tsv(file = "data/root_mass_loss/incubated_root_ash_data.txt") %>% 
  
  # Reformat plot column
  dplyr::rename(PLOT = "SAMPLE") %>% 
  dplyr::mutate(PLOT = paste("Plot", PLOT, sep = "_")) %>% 
  
  # Calculate ash-free proportion, trim table
  dplyr::mutate(ASH.FREE.PROP = ((ROOT.CRUCIBLE - EMPTY.CRUCIBLE) - (ROOT.CRUCIBLE.ASHED - EMPTY.CRUCIBLE)) / 
                  (ROOT.CRUCIBLE - EMPTY.CRUCIBLE)) %>% 
  dplyr::select(PLOT, ASH.FREE.PROP)


# Read in pre-incubated root masses
pre.incubated.ash.free.dry.mass.tbl <- readr::read_tsv(file = "data/root_mass_loss/manistee_undecomposed_root_masses.txt") %>% 
  
  # Format plot and site, drop test roots
  dplyr::rename(PLOT = "Plot", 
                SITE = "Stand", 
                ROOT.UNCORRECTED.MASS = "Root mass, mg") %>% 
  dplyr::filter(PLOT != "73" & 
                  PLOT != "73_A" & 
                  PLOT != "73_B" & 
                  PLOT != "74" & 
                  PLOT != "74_A" & 
                  PLOT != "74_B") %>% 
  dplyr::mutate(PLOT = paste("Plot", PLOT, sep = "_"), 
                SITE = paste("Site", SITE, sep = "_")) %>% 
  
  # Calculate ash-free dry mass, trim table, combine by plot
  dplyr::inner_join(., pre.incubated.ash.mean.tbl, by = "SITE") %>% 
  dplyr::mutate(PRE.INCUBATED.ASH.FREE.DRY.MASS = (ROOT.UNCORRECTED.MASS * MEAN.ASH.FREE.PROP) / 1000) %>% 
  dplyr::select(PLOT, 
                SITE, 
                PRE.INCUBATED.ASH.FREE.DRY.MASS) %>% 
  dplyr::group_by(SITE, PLOT) %>% 
  dplyr::summarise(PRE.INCUBATED.ROOT.MASS = sum(PRE.INCUBATED.ASH.FREE.DRY.MASS)) %>% 
  dplyr::ungroup()

# Read in decomposed root masses
decomp.ash.free.dry.mass.tbl <- readr::read_tsv(file = "data/root_mass_loss/manistee_roots_decomposed_masses_july_2020.txt") %>% 
  
  # Rename and format columns, and trim table
  dplyr::rename(PLOT = "Plot", 
                FRESH.MASS = "Total mass (g)") %>% 
  dplyr::select(PLOT, FRESH.MASS) %>% 
  dplyr::filter(PLOT != "73_A" & 
                  PLOT != "73_B" & 
                  PLOT != "74_A" & 
                  PLOT != "74_B") %>% 
  dplyr::mutate(PLOT = paste("Plot", PLOT, sep = "_")) %>% 
  
  # Calculate dry mass
  dplyr::inner_join(., decomp.root.moisture.tbl, by = "PLOT") %>% 
  dplyr::mutate(DECOMP.DRY.MASS = FRESH.MASS * (1 - MOISTURE.CONTENT)) %>% 
  
  # Calculate ash-free dry mass, trim table
  dplyr::inner_join(., decomp.ash.tbl, by = "PLOT") %>% 
  dplyr::mutate(DECOMP.ROOT.MASS = DECOMP.DRY.MASS * ASH.FREE.PROP) %>% 
  dplyr::select(PLOT, 
                DECOMP.ROOT.MASS)

# Read in decomposed ash-free dry mass
root.mass.loss.tbl <- decomp.ash.free.dry.mass.tbl %>% 
  
  # Combine masses, calculate mass loss, and trim table
  dplyr::inner_join(pre.incubated.ash.free.dry.mass.tbl, ., by = "PLOT") %>% 
  dplyr::mutate(PROP.MASS.LOSS = (PRE.INCUBATED.ROOT.MASS - DECOMP.ROOT.MASS) / PRE.INCUBATED.ROOT.MASS) %>% 
  dplyr::select(PLOT, 
                PROP.MASS.LOSS)

# Add to compiled data
root.mass.loss.full.tbl <- soil.compiled.data.tbl %>% 

  # Add root mass loss
  dplyr::inner_join(., root.mass.loss.tbl, by = "PLOT") %>% 
  mutate(LOG10.SOILC = log10(SPEC.TOTAL.C + 1))

# lignin~root mass loss
massloss.Nmin.gamm <- gamm(PROP.MASS.LOSS ~ s(N.MIN, k = -1), 
                           correlation = corSpatial(form = ~LON + LAT), 
                           data = root.mass.loss.full.tbl)

# lignin~root mass loss
lignin.massloss.gamm <- gamm(LIGNIN ~ s(PROP.MASS.LOSS, k = -1), 
                             correlation = corSpatial(form = ~LON + LAT), 
                             data = root.mass.loss.full.tbl)

# lignin~root mass loss
soilC.massloss.gamm <- gamm(LOG10.SOILC ~ s(PROP.MASS.LOSS, k = -1), 
                            correlation = corSpatial(form = ~LON + LAT), 
                            data = root.mass.loss.full.tbl)

# Pvalue tables
massloss.Nmin.pval.tbl <- tibble(XPOS = 0.25, 
                                 YPOS = 0.4, 
                                 PVAL.LABEL = "italic(n.s.)")

lignin.massloss.pval.tbl <- tibble(XPOS = 0.35, 
                                   YPOS = 0.2, 
                                   PVAL.LABEL = "italic(n.s.)")

soilC.massloss.pval.tbl <- tibble(XPOS = 0.35, 
                                  YPOS = 1.7, 
                                  PVAL.LABEL = "italic(n.s.)")

# Create plot
massloss.Nmin.plot <- ggplot() + 
  
  # Points
  geom_point(data = root.mass.loss.full.tbl, 
             aes(x = N.MIN, y = PROP.MASS.LOSS), 
             colour = "#240691FF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Pvalues
  geom_text(data = massloss.Nmin.pval.tbl,
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL),
            parse = TRUE,
            stat = "identity",
            size = 3) +
  
  # Titles
  labs(title = "(a)", 
       x = bquote('Inorganic N ('*mu*g~N%.%g^-1%.%d^-1*')'), 
       y = "Fine root mass loss\n(proportion of initial mass)") + 
  
  # Format
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 12, hjust = -0.3), 
        axis.title.x = element_text(colour = "black", size = 7), 
        axis.title.y = element_text(colour = "black", size = 7),  
        # x axis text
        axis.text.x = element_text(colour = "black", size = 6), 
        axis.text.y = element_text(colour = "black", size = 6))

# Create plot
lignin.massloss.plot <- ggplot() + 
  
  # Points
  geom_point(data = root.mass.loss.full.tbl, 
             aes(x = PROP.MASS.LOSS, y = LIGNIN), 
             colour = "#BE3885FF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Pvalues
  geom_text(data = lignin.massloss.pval.tbl,
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL),
            parse = TRUE,
            stat = "identity",
            size = 3) +

  # Titles
  labs(title = "(b)", 
       x = "Fine root mass loss\n(proportion of initial mass)", 
       y = "Lignin-derived SOM (prop. abund.)") + 
  
  # Format
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 12, hjust = -0.3), 
        axis.title.x = element_text(colour = "black", size = 7), 
        axis.title.y = element_text(colour = "black", size = 7), 
        # x axis text
        axis.text.x = element_text(colour = "black", size = 6), 
        axis.text.y = element_text(colour = "black", size = 6))

# Create plot
soilC.massloss.plot <- ggplot() + 
  
  # Points
  geom_point(data = root.mass.loss.full.tbl, 
             aes(x = PROP.MASS.LOSS, y = LOG10.SOILC), 
             colour = "#FEBE2AFF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Pvalues
  geom_text(data = soilC.massloss.pval.tbl,
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL),
            parse = TRUE,
            stat = "identity",
            size = 3) +

  # Titles
  labs(title = "(c)", 
       x = "Fine root mass loss\n(proportion of initial mass)", 
       y = bquote('Soil C ('*log[10]*'['*mg~C%.%g^-1*'])')) + 
  
  # Format
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 12, hjust = -0.3), 
        axis.title.x = element_text(colour = "black", size = 7), 
        axis.title.y = element_text(colour = "black", size = 7), 
        # x axis text
        axis.text.x = element_text(colour = "black", size = 6), 
        axis.text.y = element_text(colour = "black", size = 6))

# Create fine root mass loss
massloss.plot <- ggarrange(massloss.Nmin.plot, 
                           lignin.massloss.plot, 
                           soilC.massloss.plot,  
                           ncol = 3, 
                           nrow = 1, 
                           align = "hv")

# Save Figure 5
ggsave(filename = "figures/root.massloss.plots.png",
       massloss.plot,
       device = "png",
       dpi = 600,
       width = 6.5,
       height = 2.167,
       units = "in")
