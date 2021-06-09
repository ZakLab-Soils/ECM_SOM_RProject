#### 1. Set-up ####

# Load libraries
library(tidyverse)
library(viridis)
library(ggpubr)
library(TITAN2)

source("code/functions.R")
source("code/calculate_genus_abundances.R")

# Read in data
roots.genus.hlr.trim.titan <- readRDS(file = "data/working/roots.genus.hlr.trim.titan.rData")
soil.genus.hlr.trim.titan <- readRDS(file = "data/working/soil.genus.hlr.trim.titan.rData")

#### 2. Format data ####

# Read in soil functional groups
soil.fungal.functional.groups.for.titan.tbl <- as_tibble(read_tsv("data/working/soil.fungal.functional.groups.classified.0.1.txt")) %>% 
  mutate(NUTRITIONAL.MODE = ifelse(is.na(NUTRITIONAL.MODE), "Uncertain", NUTRITIONAL.MODE), 
         DECAY = ifelse(is.na(DECAY), "Uncertain", DECAY)) %>% 
  mutate(FUNGAL.GROUP = ifelse(NUTRITIONAL.MODE == "Ectomycorrhizal" & DECAY == "Putatively ligninolytic", 
                               "ECM.LIG", NA)) %>% 
  mutate(FUNGAL.GROUP = ifelse(NUTRITIONAL.MODE == "Ectomycorrhizal" & DECAY != "Putatively ligninolytic", 
                               "ECM.NONLIG", FUNGAL.GROUP)) %>% 
  mutate(FUNGAL.GROUP = ifelse(NUTRITIONAL.MODE == "Mycorrhizal" | NUTRITIONAL.MODE == "Arbuscular mycorrhizal" | 
                                 NUTRITIONAL.MODE == "Putatively mycorrhizal" | 
                                 NUTRITIONAL.MODE == "Ericoid mycorrhizal" | NUTRITIONAL.MODE == "Orchid mycorrhizal", 
                               "OTHER.MYCORRHIZA", FUNGAL.GROUP)) %>% 
  mutate(FUNGAL.GROUP = ifelse(NUTRITIONAL.MODE == "Saprotrophic" & DECAY == "White-rot or ligninolytic litter decay", 
                               "SAP.LIG", FUNGAL.GROUP)) %>% 
  mutate(FUNGAL.GROUP = ifelse(NUTRITIONAL.MODE == "Saprotrophic" & DECAY != "White-rot or ligninolytic litter decay", 
                               "SAP.NONLIG", FUNGAL.GROUP)) %>% 
  mutate(FUNGAL.GROUP = ifelse(is.na(FUNGAL.GROUP), "OTHER", FUNGAL.GROUP)) %>% 
  select(GENUS, FUNGAL.GROUP) %>% 
  # Add genus abundances
  inner_join(., soil.genus.trim.tbl, by = "GENUS") %>% 
  select(GENUS, FUNGAL.GROUP) %>% 
  distinct()

# Read in roots functional groups
roots.fungal.functional.groups.for.titan.tbl <- as_tibble(read_tsv("data/working/roots.fungal.functional.groups.classified.0.1.txt")) %>% 
  mutate(NUTRITIONAL.MODE = ifelse(is.na(NUTRITIONAL.MODE), "Uncertain", NUTRITIONAL.MODE), 
         DECAY = ifelse(is.na(DECAY), "Uncertain", DECAY)) %>% 
  mutate(FUNGAL.GROUP = ifelse(NUTRITIONAL.MODE == "Ectomycorrhizal" & DECAY == "Putatively ligninolytic", 
                               "ECM.LIG", NA)) %>% 
  mutate(FUNGAL.GROUP = ifelse(NUTRITIONAL.MODE == "Ectomycorrhizal" & DECAY != "Putatively ligninolytic", 
                               "ECM.NONLIG", FUNGAL.GROUP)) %>% 
  mutate(FUNGAL.GROUP = ifelse(NUTRITIONAL.MODE == "Mycorrhizal" | NUTRITIONAL.MODE == "Arbuscular mycorrhizal" | 
                                 NUTRITIONAL.MODE == "Putatively mycorrhizal" | 
                                 NUTRITIONAL.MODE == "Ericoid mycorrhizal" | NUTRITIONAL.MODE == "Orchid mycorrhizal", 
                               "OTHER.MYCORRHIZA", FUNGAL.GROUP)) %>% 
  mutate(FUNGAL.GROUP = ifelse(NUTRITIONAL.MODE == "Saprotrophic" & DECAY == "White-rot or ligninolytic litter decay", 
                               "SAP.LIG", FUNGAL.GROUP)) %>% 
  mutate(FUNGAL.GROUP = ifelse(NUTRITIONAL.MODE == "Saprotrophic" & DECAY != "White-rot or ligninolytic litter decay", 
                               "SAP.NONLIG", FUNGAL.GROUP)) %>% 
  mutate(FUNGAL.GROUP = ifelse(is.na(FUNGAL.GROUP), "OTHER", FUNGAL.GROUP)) %>% 
  select(GENUS, FUNGAL.GROUP) %>% 
  # Add genus abundances
  inner_join(., roots.genus.trim.tbl, by = "GENUS") %>% 
  select(GENUS, FUNGAL.GROUP) %>% 
  distinct()

# Soil titan data
soil.titan.fig.input.tbl <- as_tibble(soil.genus.hlr.trim.titan$sppmax, 
                                      rownames = "GENUS") %>% 
  inner_join(., soil.fungal.functional.groups.for.titan.tbl, by = "GENUS") %>% 
  dplyr::rename(RESPONSE.GROUP = "filter") %>% 
  select(GENUS, FUNGAL.GROUP, RESPONSE.GROUP) %>% 
  distinct() %>% 
  group_by(FUNGAL.GROUP) %>% 
  summarise(POSITIVE.RESPONSE = length(FUNGAL.GROUP[RESPONSE.GROUP == 2]), 
            NEGATIVE.RESPONSE = length(FUNGAL.GROUP[RESPONSE.GROUP == 1]), 
            NO.RESPONSE = length(FUNGAL.GROUP[RESPONSE.GROUP == 0]), 
            TOTAL.GENERA = length(FUNGAL.GROUP)) %>% 
  ungroup() %>% 
  mutate(FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with\nperoxidases", NA), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without\nperoxidases", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other\nmycorrhizas", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic\nsaprotrophs", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic\nsaprotrophs", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Other or\nuncertain ecology", FUNGAL.GROUP.LABEL)) %>% 
  mutate(FUNGAL.GROUP.LABEL = factor(FUNGAL.GROUP.LABEL, levels = c("ECM with\nperoxidases", 
                                                                    "ECM without\nperoxidases", 
                                                                    "Other\nmycorrhizas", 
                                                                    "Ligninolytic\nsaprotrophs", 
                                                                    "Non-ligninolytic\nsaprotrophs", 
                                                                    "Other or\nuncertain ecology"))) %>% 
  select(FUNGAL.GROUP.LABEL, POSITIVE.RESPONSE, NEGATIVE.RESPONSE, TOTAL.GENERA) %>% 
  pivot_longer(cols = POSITIVE.RESPONSE : NEGATIVE.RESPONSE, names_to = "RESPONSE", values_to = "COUNT") %>% 
  mutate(PROPORTION = COUNT / TOTAL.GENERA) %>% 
  mutate(PROPORTION = ifelse(RESPONSE == "NEGATIVE.RESPONSE", -1 * PROPORTION, PROPORTION)) %>% 
  mutate(SIG.LABEL = paste("(", COUNT, "/", TOTAL.GENERA, ")", sep = "")) %>% 
  mutate(LABEL.POSITION = ifelse(RESPONSE == "POSITIVE.RESPONSE", 0.2 + PROPORTION, PROPORTION - 0.2))

# Roots titan data
roots.titan.fig.input.tbl <- as_tibble(roots.genus.hlr.trim.titan$sppmax, 
                                       rownames = "GENUS") %>% 
  inner_join(., roots.fungal.functional.groups.for.titan.tbl, by = "GENUS") %>% 
  dplyr::rename(RESPONSE.GROUP = "filter") %>% 
  select(GENUS, FUNGAL.GROUP, RESPONSE.GROUP) %>% 
  distinct() %>% 
  group_by(FUNGAL.GROUP) %>% 
  summarise(POSITIVE.RESPONSE = length(FUNGAL.GROUP[RESPONSE.GROUP == 2]), 
            NEGATIVE.RESPONSE = length(FUNGAL.GROUP[RESPONSE.GROUP == 1]), 
            NO.RESPONSE = length(FUNGAL.GROUP[RESPONSE.GROUP == 0]), 
            TOTAL.GENERA = length(FUNGAL.GROUP)) %>% 
  ungroup() %>% 
  mutate(FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with\nperoxidases", NA), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without\nperoxidases", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other\nmycorrhizas", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic\nsaprotrophs", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic\nsaprotrophs", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Other or\nuncertain ecology", FUNGAL.GROUP.LABEL)) %>% 
  mutate(FUNGAL.GROUP.LABEL = factor(FUNGAL.GROUP.LABEL, levels = c("ECM with\nperoxidases", 
                                                                    "ECM without\nperoxidases", 
                                                                    "Other\nmycorrhizas", 
                                                                    "Ligninolytic\nsaprotrophs", 
                                                                    "Non-ligninolytic\nsaprotrophs", 
                                                                    "Other or\nuncertain ecology"))) %>% 
  select(FUNGAL.GROUP.LABEL, POSITIVE.RESPONSE, NEGATIVE.RESPONSE, TOTAL.GENERA) %>% 
  pivot_longer(cols = POSITIVE.RESPONSE : NEGATIVE.RESPONSE, names_to = "RESPONSE", values_to = "COUNT") %>% 
  mutate(PROPORTION = COUNT / TOTAL.GENERA) %>% 
  mutate(PROPORTION = ifelse(RESPONSE == "NEGATIVE.RESPONSE", -1 * PROPORTION, PROPORTION)) %>% 
  mutate(SIG.LABEL = paste("(", COUNT, "/", TOTAL.GENERA, ")", sep = "")) %>% 
  mutate(LABEL.POSITION = ifelse(RESPONSE == "POSITIVE.RESPONSE", 0.2 + PROPORTION, PROPORTION - 0.2))

#### 3. Summary data for Figure 2 ####

# Soil
soil.titan.responses.tbl <- as_tibble(soil.genus.hlr.trim.titan$sppmax, rownames = NA) %>% 
  rownames_to_column(var = "GENUS") %>% 
  select(GENUS, zenv.cp, `5%`, `95%`, z.median, filter, maxgrp) %>% 
  rename(ZENV.CP = "zenv.cp", 
         LCI = "5%", 
         UCI = "95%", 
         Z.MEDIAN = "z.median", 
         GROUP = "filter", 
         MAXGRP = "maxgrp") %>% 
  mutate(Z.MEDIAN = ifelse(MAXGRP == 1, -1 * Z.MEDIAN, Z.MEDIAN)) %>% 
  inner_join(., soil.fungal.functional.groups.for.titan.tbl, by = "GENUS") %>% 
  # Add functional group
  mutate(FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with\nperoxidases", NA), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without\nperoxidases", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other\nmycorrhizas", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic\nsaprotrophs", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic\nsaprotrophs", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Other or\nuncertain ecology", FUNGAL.GROUP.LABEL)) %>% 
  mutate(FUNGAL.GROUP.LABEL = factor(FUNGAL.GROUP.LABEL, levels = c("ECM with\nperoxidases", 
                                                                    "ECM without\nperoxidases", 
                                                                    "Other\nmycorrhizas", 
                                                                    "Ligninolytic\nsaprotrophs", 
                                                                    "Non-ligninolytic\nsaprotrophs", 
                                                                    "Other or\nuncertain ecology"))) %>% 
  arrange(FUNGAL.GROUP.LABEL, Z.MEDIAN) %>% 
  mutate(GENUS = factor(GENUS, levels = GENUS)) %>% 
  mutate(SIG = ifelse(GROUP > 0, "Significant", "Not significant")) %>% 
  mutate(SIG = factor(SIG, levels = c("Significant", "Not significant")))

# Roots
roots.titan.responses.tbl <- as_tibble(roots.genus.hlr.trim.titan$sppmax, rownames = NA) %>% 
  rownames_to_column(var = "GENUS") %>% 
  select(GENUS, zenv.cp, `5%`, `95%`, z.median, filter, maxgrp) %>% 
  rename(ZENV.CP = "zenv.cp", 
         LCI = "5%", 
         UCI = "95%", 
         Z.MEDIAN = "z.median", 
         GROUP = "filter", 
         MAXGRP = "maxgrp") %>% 
  mutate(Z.MEDIAN = ifelse(MAXGRP == 1, -1 * Z.MEDIAN, Z.MEDIAN)) %>% 
  inner_join(., roots.fungal.functional.groups.for.titan.tbl, by = "GENUS") %>% 
  # Add functional group
  mutate(FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with\nperoxidases", NA), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without\nperoxidases", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other\nmycorrhizas", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic\nsaprotrophs", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic\nsaprotrophs", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Other or\nuncertain ecology", FUNGAL.GROUP.LABEL)) %>% 
  mutate(FUNGAL.GROUP.LABEL = factor(FUNGAL.GROUP.LABEL, levels = c("ECM with\nperoxidases", 
                                                                    "ECM without\nperoxidases", 
                                                                    "Other\nmycorrhizas", 
                                                                    "Ligninolytic\nsaprotrophs", 
                                                                    "Non-ligninolytic\nsaprotrophs", 
                                                                    "Other or\nuncertain ecology"))) %>% 
  arrange(FUNGAL.GROUP.LABEL, Z.MEDIAN) %>% 
  mutate(GENUS = factor(GENUS, levels = GENUS)) %>% 
  mutate(SIG = ifelse(GROUP > 0, "Significant", "Not significant")) %>% 
  mutate(SIG = factor(SIG, levels = c("Significant", "Not significant")))

#### 5. Figure 3 ####

# Soil
soil.titan.responses.plot <- ggplot() + 
  
  # Bars
  geom_bar(data = soil.titan.responses.tbl, 
           aes(x = GENUS, y = Z.MEDIAN, fill = FUNGAL.GROUP.LABEL, alpha = SIG), 
           stat = "identity", 
           width = 0.75) + 
  
  # Scale fill
  scale_fill_viridis_d(name = NULL, 
                       option = "plasma", 
                       begin = 0.1, 
                       end = 0.9) + 
  
  # Scale alpha
  scale_alpha_manual(name = NULL, 
                     values = c(1, 0.3)) + 
  
  # Line
  geom_hline(yintercept = 0, 
             linetype = 2, 
             size = 0.5, 
             colour = "black") + 
  
  # Titles
  labs(title = "(a)", 
       x = "Genus", 
       y = "Median Z-score\n(response to inorganic N availability)") + 
  
  coord_flip() + 
  
  scale_x_discrete(limits = rev(levels(soil.titan.responses.tbl$GENUS))) + 
  
  # Format panel
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 12, colour = "black", hjust = -0.5), 
        axis.title.x = element_text(size = 7, colour = "black"), 
        axis.title.y = element_text(size = 7, colour = "black"), 
        # Axis text
        axis.text.x = element_text(size = 6, colour = "black"), 
        axis.text.y = element_text(size = 6, colour = "black"), 
        legend.title = element_text(size = 7, colour = "black"), 
        legend.text = element_text(size = 6, colour = "black"), 
        legend.key = element_blank())

# Roots
roots.titan.responses.plot <- ggplot() + 
  
  # Bars
  geom_bar(data = roots.titan.responses.tbl, 
           aes(x = GENUS, y = Z.MEDIAN, fill = FUNGAL.GROUP.LABEL, alpha = SIG), 
           stat = "identity", 
           width = 0.75) + 
  
  # Scale fill
  scale_fill_viridis_d(name = NULL, 
                       option = "plasma", 
                       begin = 0.1, 
                       end = 0.9) + 
  
  # Scale alpha
  scale_alpha_manual(name = NULL, 
                     values = c(1, 0.3)) + 
  
  # Line
  geom_hline(yintercept = 0, 
             linetype = 2, 
             size = 0.5, 
             colour = "black") + 
  
  # Titles
  labs(title = "(b)", 
       x = "Genus", 
       y = "Median Z-score\n(response to inorganic N availability)") + 
  
  coord_flip() + 
  
  scale_x_discrete(limits = rev(levels(roots.titan.responses.tbl$GENUS))) + 
  
  # Format panel
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 12, colour = "black", hjust = -0.5), 
        axis.title.x = element_text(size = 7, colour = "black"), 
        axis.title.y = element_text(size = 7, colour = "black"), 
        # Axis text
        axis.text.x = element_text(size = 6, colour = "black"), 
        axis.text.y = element_text(size = 6, colour = "black"), 
        legend.title = element_text(size = 7, colour = "black"), 
        legend.text = element_text(size = 6, colour = "black"), 
        legend.key = element_blank())

# Combined figure
Fig3 <- ggarrange(soil.titan.responses.plot, 
                  roots.titan.responses.plot, 
                  ncol = 2, nrow = 1, 
                  align = "hv", 
                  common.legend = TRUE, 
                  legend = "bottom")

# Save figure
ggsave(filename = "figures/Fig_3.pdf",
       plot = Fig3, 
       device = "pdf",
       width = 6.5,
       height = 6,
       units = "in")
