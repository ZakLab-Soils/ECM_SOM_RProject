#### 1. Set-up ####

# This script calculates changes in fungal community size (qPCR of ITS1 region)

# Load libraries
library(tidyverse)
library(mgcv)
library(itsadug)
library(ggpubr)
library(viridis)

source("code/functions.R")
source("code/fine_root_mass_loss.R")
source("code/spring_2019_net_N_mineralization.R")

# Get soil DNA extraction masses
soil.DNA.extraction.masses.tbl <- read_tsv("data/qPCR/soil_extraction_mass.txt") %>% 
  
  # Format
  dplyr::rename(PLOT = "Sample number (plot)", 
                REPL = "Replicate (~0.25 g each)", 
                MASS = "Fresh mass of soil") %>% 
  dplyr::mutate(PLOT = gsub("\\ ", "_", PLOT)) %>% 
  
  # Calculate total mass
  dplyr::group_by(PLOT) %>% 
  dplyr::summarise(MASS = sum(MASS)) %>% 
  dplyr::ungroup() %>% 
  
  # Add percent moisture and calculate dry mass, trim
  dplyr::inner_join(., soil.moisture.tbl, by = "PLOT") %>% 
  dplyr::mutate(SOIL.DRY.MASS = MASS * (1 - GRAV)) %>% 
  dplyr::select(PLOT, 
                SOIL.DRY.MASS)

# Get root DNA extraction masses
roots.DNA.extraction.masses.tbl <- read_tsv("data/qPCR/root_extraction_mass.txt") %>% 
  
  # Format
  dplyr::rename(PLOT = "Sample number (plot)", 
                REPL = "Replicate", 
                MASS = "sample mass") %>% 
  dplyr::mutate(PLOT = gsub("\\ ", "_", PLOT)) %>% 
  
  # Calculate total mass
  dplyr::group_by(PLOT) %>% 
  dplyr::summarise(MASS = sum(MASS)) %>% 
  dplyr::ungroup() %>% 
  
  # Add root moisture and ash content, calculate dry mass
  dplyr::inner_join(., decomp.root.moisture.tbl, by = "PLOT") %>% 
  dplyr::inner_join(., decomp.ash.tbl, by = "PLOT") %>% 
  dplyr::mutate(DRY.MASS = (MASS * (1 - MOISTURE.CONTENT) * ASH.FREE.PROP))
  
# Read in data
soil.compiled.data.tbl <- readRDS(file = "data/working/soil.compiled.data.tbl.RData")
roots.compiled.data.tbl <- readRDS(file = "data/working/roots.compiled.data.tbl.RData")

# Read in qPCR data
soil.ITS.copies.tbl1 <- read_tsv(file = "data/qPCR/ITS1_MS_qPCR1_individual_full_curve_data.txt") %>% 
  
  # Trim
  dplyr::filter(`Well Type` == "Unknown" & Threshold != "Reference") %>% 
  dplyr::mutate(Quantity = as.numeric(Quantity), 
                PLOT = c(1 : nrow(.)), 
                PLOT = paste("Plot", PLOT, sep = "_"), 
                RUN = "Run_1") %>% 
  rename(INITIAL.COPIES = "Quantity") %>% 
  dplyr::select(PLOT, 
                RUN, 
                INITIAL.COPIES)

# Read in qPCR data
soil.ITS.copies.tbl2 <- read_tsv(file = "data/qPCR/ITS1_MS_qPCR2_individual_full_curve_data.txt") %>% 
  
  # Trim
  dplyr::filter(`Well Type` == "Unknown" & Threshold != "Reference") %>% 
  dplyr::mutate(Quantity = as.numeric(Quantity), 
                PLOT = c(1 : nrow(.)), 
                PLOT = paste("Plot", PLOT, sep = "_"), 
                RUN = "Run_2") %>% 
  rename(INITIAL.COPIES = "Quantity") %>% 
  dplyr::select(PLOT, 
                RUN, 
                INITIAL.COPIES)

# Read in qPCR data
roots.ITS.copies.tbl1 <- read_tsv(file = "data/qPCR/ITS1_R_qPCR1_individual_full_curve_data.txt") %>% 
  
  # Trim
  dplyr::filter(`Well Type` == "Unknown" & Threshold != "Reference") %>% 
  dplyr::mutate(Quantity = as.numeric(Quantity), 
                PLOT = c(1 : nrow(.)), 
                PLOT = paste("Plot", PLOT, sep = "_"), 
                RUN = "Run_1") %>% 
  rename(INITIAL.COPIES = "Quantity") %>% 
  dplyr::select(PLOT, 
                RUN, 
                INITIAL.COPIES)

# Read in qPCR data
roots.ITS.copies.tbl2 <- read_tsv(file = "data/qPCR/ITS1_R_qPCR2_individual_full_curve_data.txt") %>% 
  
  # Trim
  dplyr::filter(`Well Type` == "Unknown" & Threshold != "Reference") %>% 
  dplyr::mutate(Quantity = as.numeric(Quantity), 
                PLOT = c(1 : nrow(.)), 
                PLOT = paste("Plot", PLOT, sep = "_"), 
                RUN = "Run_2") %>% 
  rename(INITIAL.COPIES = "Quantity") %>% 
  dplyr::select(PLOT, 
                RUN, 
                INITIAL.COPIES)

# Combine and calculate copy number
soil.ITS.copies.tbl <- dplyr::bind_rows(soil.ITS.copies.tbl1, 
                                        soil.ITS.copies.tbl2) %>% 
  dplyr::group_by(PLOT) %>% 
  dplyr::summarise(SOIL.INITIAL.COPIES = mean(INITIAL.COPIES)) %>% 
  dplyr::ungroup() %>% 
  
  # Add soil mass
  dplyr::inner_join(., soil.DNA.extraction.masses.tbl, by = "PLOT") %>% 
  dplyr::mutate(VOLUME = 1, 
                DILUTION.FACTOR = 10, 
                EXTR.VOLUME = 400) %>% 
  dplyr::mutate(EXTR.VOLUME = ifelse(PLOT == "Plot_23" | 
                                       PLOT == "Plot_37" | 
                                       PLOT == "Plot_56" | 
                                       PLOT == "Plot_61", 200, EXTR.VOLUME)) %>% 
  
  # Calculate copy number
  dplyr::mutate(SOIL.COPY.NUMBER = (SOIL.INITIAL.COPIES * VOLUME * DILUTION.FACTOR * EXTR.VOLUME * (1 / SOIL.DRY.MASS))) %>% 
  
  # Add env data 
  dplyr::inner_join(., soil.compiled.data.tbl, by = "PLOT")

# Combine
roots.ITS.copies.tbl <- dplyr::bind_rows(roots.ITS.copies.tbl1, 
                                         roots.ITS.copies.tbl2) %>% 
  dplyr::group_by(PLOT) %>% 
  dplyr::summarise(ROOTS.INITIAL.COPIES = mean(INITIAL.COPIES)) %>% 
  dplyr::ungroup() %>% 
  
  # Add root mass
  dplyr::inner_join(., roots.DNA.extraction.masses.tbl, by = "PLOT") %>% 
  dplyr::mutate(VOLUME = 1, 
                DILUTION.FACTOR = 10, 
                CLEANUP.OUTPUT = 100, 
                CLEANUP.INPUT = 150, 
                EXTR.VOLUME = 600) %>% 
  
  # Calculate copy number
  dplyr::mutate(ROOTS.COPY.NUMBER = (ROOTS.INITIAL.COPIES * VOLUME * DILUTION.FACTOR * CLEANUP.OUTPUT * (1 / CLEANUP.INPUT) * EXTR.VOLUME * (1 / DRY.MASS))) %>% 
  
  # Add env data 
  dplyr::inner_join(., roots.compiled.data.tbl, by = "PLOT")

# soil ITS Nmin
soil.ITS.Nmin.gamm <- gamm(SOIL.COPY.NUMBER ~ s(N.MIN, k = -1), 
                           data = soil.ITS.copies.tbl)

roots.ITS.Nmin.gamm <- gamm(ROOTS.COPY.NUMBER ~ s(ECM.LIG, k = -1), 
                            data = roots.ITS.copies.tbl)

# Pvalue tables
soil.ITS.Nmin.pval.tbl <- tibble(XPOS = 0.8, 
                                 YPOS = 175000000, 
                                 PVAL.LABEL = "italic(P) < 0.001")

roots.ITS.Nmin.pval.tbl <- tibble(XPOS = 0.8, 
                                  YPOS = 5100000000, 
                                  PVAL.LABEL = "italic(P) == 0.009")

# Create plot
soil.ITS.Nmin.plot <- ggplot() + 
  
  # Points
  geom_point(data = soil.ITS.copies.tbl, 
             aes(x = N.MIN, y = SOIL.COPY.NUMBER), 
             colour = "#240691FF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Pvalues
  geom_text(data = soil.ITS.Nmin.pval.tbl,
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL),
            parse = TRUE,
            stat = "identity",
            size = 3) +
  
  # Scale
  scale_y_continuous(breaks = c(0, 50000000, 100000000, 150000000, 200000000), 
                     labels = expression(0, 5%*%10^7, 1.0%*%10^8, 1.5%*%10^8, 2.0%*%10^8),
                     limits = c(0, 200000000)) + 
  
  # Titles
  labs(title = "(a)", 
       x = bquote('Inorganic N ('*mu*g~N%.%g^-1%.%d^-1*')'), 
       y = bquote(atop('Fungal community size', '(ITS copies'%.%g^-1*')'))) + 
  
  # Line
  geom_smooth(data = soil.ITS.copies.tbl, 
              aes(x = N.MIN, y = SOIL.COPY.NUMBER), 
              formula = y ~ s(x), 
              method = "gam", 
              colour = "#240691FF", 
              fill = "#240691FF") + 
  
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
roots.ITS.Nmin.plot <- ggplot() + 
  
  # Points
  geom_point(data = roots.ITS.copies.tbl, 
             aes(x = N.MIN, y = ROOTS.COPY.NUMBER), 
             colour = "#FEBE2AFF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Pvalues
  geom_text(data = roots.ITS.Nmin.pval.tbl,
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL),
            parse = TRUE,
            stat = "identity",
            size = 3) +
  
  # Scale
  scale_y_continuous(breaks = c(0, 1500000000, 3000000000, 4500000000, 6000000000),
                     labels = expression(0, 1.5%*%10^9, 3.0%*%10^9, 4.5%*%10^9, 6.0%*%10^9),
                     limits = c(0, 6000000000)) +
  
  # Titles
  labs(title = "(b)", 
       x = bquote('Inorganic N ('*mu*g~N%.%g^-1%.%d^-1*')'), 
       y = bquote(atop('Fungal community size', '(ITS copies'%.%g^-1*')'))) + 
  
  # Line
  geom_smooth(data = roots.ITS.copies.tbl, 
              aes(x = N.MIN, y = ROOTS.COPY.NUMBER), 
              formula = y ~ s(x), 
              method = "gam", 
              colour = "#FEBE2AFF", 
              fill = "#FEBE2AFF") + 
  
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


# Create final fig
qPCR.plot <- ggarrange(soil.ITS.Nmin.plot, 
                       roots.ITS.Nmin.plot, 
                       ncol = 2, 
                       nrow = 1, 
                       align = "hv")

# Save Figure S8
ggsave(filename = "figures/FigS8.png",
       qPCR.plot,
       device = "png",
       dpi = 600,
       width = 6.5,
       height = 3.5,
       units = "in")

#### Prop of ECM seqs ####

soil.gamm.input.tbl2 <- soil.gamm.input.tbl %>% 
  
  dplyr::mutate(PROP.ECM.LIG = ECM.LIG * (1 / (ECM.LIG + ECM.NONLIG)))

soil.ECM.prop.gamm <- gamm(PROP.ECM.LIG ~ s(N.MIN, k = -1), 
                           correlation = corSpatial(form = ~LON + LAT), 
                           data = soil.gamm.input.tbl2)

soil.ECM.prop.P.tbl <- tibble(XPOS = 0.95, 
                              YPOS = 0.9, 
                              PVAL.LABEL = "italic(n.s.)")

# Figure
soil.ECM.prop.plot <- ggplot() + 
  
  # Points
  geom_point(data = soil.gamm.input.tbl2, 
             aes(x = N.MIN, y = PROP.ECM.LIG), 
             colour = "#E56B5DFF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Pvalues
  geom_text(data = soil.ECM.prop.P.tbl,
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL),
            parse = TRUE,
            stat = "identity",
            size = 3) +
  
  # Titles
  labs(y = "Proportion of ECM fungal sequences", 
       x = bquote('Inorganic N ('*mu*g~N%.%g^-1%.%d^-1*')')) + 
  
  # Format
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 12, hjust = -0.08), 
        axis.title.x = element_text(colour = "black", size = 8), 
        axis.title.y = element_text(colour = "black", size = 8), 
        # x axis text
        axis.text.x = element_text(colour = "black", size = 7), 
        axis.text.y = element_text(colour = "black", size = 7), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, colour = "black"))

ggsave(filename = "figures/FigS9.png",
       soil.ECM.prop.plot,
       device = "png",
       dpi = 600,
       width = 6.5,
       height = 3.5,
       units = "in")
