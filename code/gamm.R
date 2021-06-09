#### 1. Set-up ####

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
roots.compiled.data.tbl <- readRDS(file = "data/working/roots.compiled.data.tbl.RData")

#### 2. Format GAMM input data ####

# Format GAMM input data, soil
soil.gamm.input.tbl <- env.input.tbl %>% 
  select(PLOT, SOIL.PH, VWC.MEAN, TEMP.MEAN) %>% 
  inner_join(soil.compiled.data.tbl, ., by = "PLOT") %>% 
  mutate(LOG10.SOILC = log10(SPEC.TOTAL.C + 1))

# Format GAMM input data, roots
roots.gamm.input.tbl <- env.input.tbl %>% 
  select(PLOT, SOIL.PH, VWC.MEAN, TEMP.MEAN) %>% 
  inner_join(roots.compiled.data.tbl, ., by = "PLOT") %>% 
  mutate(LOG10.SOILC = log10(SPEC.TOTAL.C + 1))

#### 3. Biogeochemistry GAMMs ####

# lignin~Nmin
lignin.Nmin.gamm <- gamm(LIGNIN ~ s(N.MIN, k = -1), 
                         correlation = corSpatial(form = ~LON + LAT), 
                         data = soil.gamm.input.tbl)

# soilC~lignin
soilC.lignin.gamm <- gamm(LOG10.SOILC ~ s(LIGNIN, k = -1), 
                          correlation = corSpatial(form = ~LON + LAT), 
                          data = soil.gamm.input.tbl)

# soilC~Nmin
soilC.Nmin.gamm <- gamm(LOG10.SOILC ~ s(N.MIN, k = -1), 
                        correlation = corSpatial(form = ~LON + LAT), 
                        data = soil.gamm.input.tbl)

#### 4. Fungi ~ Nmin GAMMs, soil ####

# Write function to collect summary data from GAMM
fungi_Nmin_GAMM_summary <- function(x) {
  
  # Get summary
  tmp1 <- summary(x$gam)
  
  # Get R2, F, and P
  tmp2 <- tmp1$r.sq
  tmp3 <- tmp1$s.table[3]
  tmp4 <- tmp1$s.table[4]
  
  # Create tibble
  tmp5 <- tibble(R2 = tmp2, 
                 FVAL = tmp3, 
                 PVAL = tmp4)
  
  # Return
  return(tmp5)
  
}

# ECM with peroxidases
soil.ECMLIG.Nmin.gamm <- gamm(ECM.LIG ~ s(N.MIN, k = -1), 
                              correlation = corSpatial(form = ~LON + LAT), 
                              data = soil.gamm.input.tbl)

# ECM without peroxidases
soil.ECMNONLIG.Nmin.gamm <- gamm(ECM.NONLIG ~ s(N.MIN, k = -1), 
                                 correlation = corSpatial(form = ~LON + LAT), 
                                 data = soil.gamm.input.tbl)

# Other mycorrhizas
soil.OTHERMYCO.Nmin.gamm <- gamm(OTHER.MYCORRHIZA ~ s(N.MIN, k = -1), 
                                 correlation = corSpatial(form = ~LON + LAT), 
                                 data = soil.gamm.input.tbl)
# Ligninolytic saprotrophs
soil.SAPLIG.Nmin.gamm <- gamm(SAP.LIG ~ s(N.MIN, k = -1), 
                              correlation = corSpatial(form = ~LON + LAT), 
                              data = soil.gamm.input.tbl)

# Non-ligninolytic saprotrophs
soil.SAPNONLIG.Nmin.gamm <- gamm(SAP.NONLIG ~ s(N.MIN, k = -1), 
                                 correlation = corSpatial(form = ~LON + LAT), 
                                 data = soil.gamm.input.tbl)

# Other ecology
soil.OTHER.Nmin.gamm <- gamm(OTHER ~ s(N.MIN, k = -1), 
                             correlation = corSpatial(form = ~LON + LAT), 
                             data = soil.gamm.input.tbl)

#### 5. Fungi ~ Nmin GAMMs, roots ####

# ECM with peroxidases
roots.ECMLIG.Nmin.gamm <- gamm(ECM.LIG ~ s(N.MIN, k = -1), 
                               correlation = corSpatial(form = ~LON + LAT), 
                               data = roots.gamm.input.tbl)

# ECM without peroxidases
roots.ECMNONLIG.Nmin.gamm <- gamm(ECM.NONLIG ~ s(N.MIN, k = -1), 
                                  correlation = corSpatial(form = ~LON + LAT), 
                                  data = roots.gamm.input.tbl)

# Other mycorrhizas
roots.OTHERMYCO.Nmin.gamm <- gamm(OTHER.MYCORRHIZA ~ s(N.MIN, k = -1), 
                                  correlation = corSpatial(form = ~LON + LAT), 
                                  data = roots.gamm.input.tbl)

# Ligninolytic saprotrophs
roots.SAPLIG.Nmin.gamm <- gamm(SAP.LIG ~ s(N.MIN, k = -1), 
                               correlation = corSpatial(form = ~LON + LAT), 
                               data = roots.gamm.input.tbl)

# Non-ligninolytic saprotrophs
roots.SAPNONLIG.Nmin.gamm <- gamm(SAP.NONLIG ~ s(N.MIN, k = -1), 
                                  correlation = corSpatial(form = ~LON + LAT), 
                                  data = roots.gamm.input.tbl)

# Other ecology
roots.OTHER.Nmin.gamm <- gamm(OTHER ~ s(N.MIN, k = -1), 
                              correlation = corSpatial(form = ~LON + LAT), 
                              data = roots.gamm.input.tbl)

#### 7. lignin ~ fungi, soil ####

# GAMM
soil.lignin.fungi.gamm <- gamm(LIGNIN ~ s(ECM.LIG, k = -1) + s(ECM.NONLIG, k = -1) + s(OTHER.MYCORRHIZA, k = -1) + 
                                 s(SAP.LIG, k = -1) + s(SAP.NONLIG, k = -1) + s(OTHER, k = -1), 
                               correlation = corSpatial(form = ~LON + LAT), 
                               data = soil.gamm.input.tbl)

#### 8. soilC ~ fungi, soil ####

# GAMM
soil.soilC.fungi.gamm <- gamm(LOG10.SOILC ~ s(ECM.LIG, k = -1) + s(ECM.NONLIG, k = -1) + s(OTHER.MYCORRHIZA, k = -1) + 
                                s(SAP.LIG, k = -1) + s(SAP.NONLIG, k = -1) + s(OTHER, k = -1), 
                              correlation = corSpatial(form = ~LON + LAT), 
                              data = soil.gamm.input.tbl)

#### 9. lignin ~ fungi, roots ####

# GAMM
roots.lignin.fungi.gamm <- gamm(LIGNIN ~ s(ECM.LIG, k = -1) + s(ECM.NONLIG, k = -1) + s(OTHER.MYCORRHIZA, k = -1) + 
                                  s(SAP.LIG, k = -1) + s(SAP.NONLIG, k = -1) + s(OTHER, k = -1), 
                                correlation = corSpatial(form = ~LON + LAT), 
                                data = roots.gamm.input.tbl)

#### 10. soilC ~ fungi, roots ####

# GAMM
roots.soilC.fungi.gamm <- gamm(LOG10.SOILC ~ s(ECM.LIG, k = -1) + s(ECM.NONLIG, k = -1) + s(OTHER.MYCORRHIZA, k = -1) + 
                                 s(SAP.LIG, k = -1) + s(SAP.NONLIG, k = -1) + s(OTHER, k = -1), 
                               correlation = corSpatial(form = ~LON + LAT), 
                               data = roots.gamm.input.tbl)

#### 11. lignin-Nmin figure ####

# Summary
lignin.Nmin.summary.input.tbl <- fungi_Nmin_GAMM_summary(lignin.Nmin.gamm) %>% 
  mutate(PVAL.LABEL = rep(NA, nrow(.))) %>% 
  mutate(XPOS = 0.99, 
         YPOS = 0.06, 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL < 0.001, paste("atop(R[adj.]^2 == ", round(R2, 3), ", ", "italic(P) < 0.001)", sep = ""), PVAL.LABEL), 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL <= 0.05 & PVAL >= 0.001, paste("atop(R[adj.]^2 == ", round(R2, 3), ", ", "italic(P) == ", round(PVAL, 3), ")", sep = ""), PVAL.LABEL), 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL > 0.05, "atop(italic(n.s.), )", PVAL.LABEL))

# Get predicted values
lignin.Nmin.gamm.pred <- get_predictions(lignin.Nmin.gamm$gam, 
                                         cond = list(N.MIN = seq(min(soil.gamm.input.tbl$N.MIN, na.rm = T),
                                                                 max(soil.gamm.input.tbl$N.MIN, na.rm = T), 
                                                                 length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1) %>% 
  as_tibble()

# Create plot
lignin.Nmin.plot <- ggplot() + 
  
  # Points
  geom_point(data = soil.gamm.input.tbl, 
             aes(x = N.MIN, y = LIGNIN), 
             colour = "#240691FF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Confidence intervals
  geom_ribbon(data = lignin.Nmin.gamm.pred, 
              aes(x = N.MIN, y = yval, ymin = fit - CI, ymax = fit + CI), 
              show.legend = FALSE, 
              fill = "#240691FF", 
              alpha = 0.3) + 
  
  # Line
  geom_line(data = lignin.Nmin.gamm.pred, 
            aes(x = N.MIN, y = fit), 
            show.legend = FALSE, 
            color = "#240691FF", 
            size = 1) + 
  
  # Pvalues
  geom_text(data = lignin.Nmin.summary.input.tbl, 
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL), 
            parse = TRUE, 
            stat = "identity", 
            size = 3) + 
  
  # Titles
  labs(title = "(a)", 
       x = bquote('Inorganic N ('*mu*g~N%.%g^-1%.%d^-1*')'), 
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

#### 14. soilC-lignin figure ####

# Summary
soilC.lignin.summary.input.tbl <- fungi_Nmin_GAMM_summary(soilC.lignin.gamm) %>% 
  mutate(PVAL.LABEL = rep(NA, nrow(.))) %>% 
  mutate(XPOS = 0.2, 
         YPOS = 1.135, 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL < 0.001, paste("atop(R[adj.]^2 == ", round(R2, 3), ", ", "italic(P) < 0.001)", sep = ""), PVAL.LABEL), 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL <= 0.05 & PVAL >= 0.001, paste("atop(R[adj.]^2 == ", round(R2, 3), ", ", "italic(P) == ", round(PVAL, 3), ")", sep = ""), PVAL.LABEL), 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL > 0.05, "atop(italic(n.s.), )", PVAL.LABEL))

# Get predicted values
soilC.lignin.gamm.pred <- get_predictions(soilC.lignin.gamm$gam, 
                                          cond = list(LIGNIN = seq(min(soil.gamm.input.tbl$LIGNIN, na.rm = T),
                                                                   max(soil.gamm.input.tbl$LIGNIN, na.rm = T), 
                                                                   length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1) %>% 
  as_tibble()

# Create plot
soilC.lignin.plot <- ggplot() + 
  
  # Points
  geom_point(data = soil.gamm.input.tbl, 
             aes(x = LIGNIN, y = LOG10.SOILC), 
             colour = "#BE3885FF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Confidence intervals
  geom_ribbon(data = soilC.lignin.gamm.pred, 
              aes(x = LIGNIN, y = yval, ymin = fit - CI, ymax = fit + CI), 
              show.legend = FALSE, 
              fill = "#BE3885FF", 
              alpha = 0.3) + 
  
  # Line
  geom_line(data = soilC.lignin.gamm.pred, 
            aes(x = LIGNIN, y = fit), 
            show.legend = FALSE, 
            color = "#BE3885FF", 
            size = 1) + 
  
  # Pvalues
  geom_text(data = soilC.lignin.summary.input.tbl, 
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL), 
            parse = TRUE, 
            stat = "identity", 
            size = 3) + 
  
  # Titles
  labs(title = "(b)", 
       x = "Lignin-derived SOM (prop. abund.)", 
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

#### 15. soilC-Nmin figure ####

# Summary
soilC.Nmin.summary.input.tbl <- fungi_Nmin_GAMM_summary(soilC.Nmin.gamm) %>% 
  mutate(PVAL.LABEL = rep(NA, nrow(.))) %>% 
  mutate(XPOS = 0.99, 
         YPOS = 1.135, 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL < 0.001, paste("atop(R[adj.]^2 == ", round(R2, 3), ", ", "italic(P) < 0.001)", sep = ""), PVAL.LABEL), 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL <= 0.05 & PVAL >= 0.001, paste("atop(R[adj.]^2 == ", round(R2, 3), ", ", "italic(P) == ", round(PVAL, 3), ")", sep = ""), PVAL.LABEL), 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL > 0.05, "atop(italic(n.s.), )", PVAL.LABEL))

# Get predicted values
soilC.Nmin.gamm.pred <- get_predictions(soilC.Nmin.gamm$gam, 
                                        cond = list(N.MIN = seq(min(soil.gamm.input.tbl$N.MIN, na.rm = T),
                                                                max(soil.gamm.input.tbl$N.MIN, na.rm = T), 
                                                                length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1) %>% 
  as_tibble()

# Create plot
soilC.Nmin.plot <- ggplot() + 
  
  # Points
  geom_point(data = soil.gamm.input.tbl, 
             aes(x = N.MIN, y = LOG10.SOILC), 
             colour = "#FEBE2AFF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Confidence intervals
  geom_ribbon(data = soilC.Nmin.gamm.pred, 
              aes(x = N.MIN, y = yval, ymin = fit - CI, ymax = fit + CI), 
              show.legend = FALSE, 
              fill = "#FEBE2AFF", 
              alpha = 0.3) + 
  
  # Line
  geom_line(data = soilC.Nmin.gamm.pred, 
            aes(x = N.MIN, y = fit), 
            show.legend = FALSE, 
            color = "#FEBE2AFF", 
            size = 1) + 
  
  # Pvalues
  geom_text(data = soilC.Nmin.summary.input.tbl, 
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL), 
            parse = TRUE, 
            stat = "identity", 
            size = 3) + 
  
  # Titles
  labs(title = "(c)", 
       x = bquote('Inorganic N ('*mu*g~N%.%g^-1%.%d^-1*')'), 
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

#### 16. Figure 5 - biogeochemistry ####

# Create figure 5
Figure5 <- ggarrange(lignin.Nmin.plot, 
                     soilC.lignin.plot, 
                     soilC.Nmin.plot, 
                     ncol = 3, 
                     nrow = 1, 
                     align = "hv")

# Save Figure 5
# ggsave(filename = "figures/Figure5.png",
#        Figure5,
#        device = "png",
#        dpi = 600,
#        width = 6.5,
#        height = 2.167,
#        units = "in")

# Save figure
ggsave(filename = "figures/Fig_5.pdf",
       plot = Figure5, 
       device = "pdf",
       width = 6.5,
       height = 2.167,
       units = "in")

#### Table S9 ####

write_tsv(lignin.Nmin.summary.input.tbl, 
          path = "data/working/lignin.N.min.gamm.summary.TableS9.txt")
write_tsv(soilC.lignin.summary.input.tbl, 
          path = "data/working/soilC.lignin.gamm.summary.TableS9.txt")
write_tsv(soilC.Nmin.summary.input.tbl, 
          path = "data/working/soilC.N.min.gamm.summary.TableS9.txt")

#### 17. Fungi ~ Nmin, soil ####

# Get summary data
fungi.N.min.summary.input.list <- list(soil.ECMLIG.Nmin.gamm, 
                                       soil.ECMNONLIG.Nmin.gamm, 
                                       soil.OTHERMYCO.Nmin.gamm, 
                                       soil.SAPLIG.Nmin.gamm, 
                                       soil.SAPNONLIG.Nmin.gamm, 
                                       soil.OTHER.Nmin.gamm, 
                                       roots.ECMLIG.Nmin.gamm, 
                                       roots.ECMNONLIG.Nmin.gamm, 
                                       roots.OTHERMYCO.Nmin.gamm, 
                                       roots.SAPLIG.Nmin.gamm, 
                                       roots.SAPNONLIG.Nmin.gamm, 
                                       roots.OTHER.Nmin.gamm)

fungi.N.min.summary.list <- lapply(fungi.N.min.summary.input.list, 
                                  fungi_Nmin_GAMM_summary) %>% 
  bind_rows() %>% 
  mutate(SUBSTRATE = c(rep("soil", 6), rep("decomp.roots", 6)), 
         FUNGAL.GROUP = rep(c("ECM.LIG", "ECM.NONLIG", "OTHER.MYCORRHIZA", 
                              "SAP.LIG", "SAP.NONLIG", "OTHER"), 2)) %>% 
  # B-H FDR correction
  mutate(PVAL.ADJ = p.adjust(PVAL, method = "BH")) %>% 
  group_by(SUBSTRATE) %>% 
  group_split()

soil.fungi.N.min.summary.tbl <- fungi.N.min.summary.list[[2]]
roots.fungi.N.min.summary.tbl <- fungi.N.min.summary.list[[1]]

#### Table S2 ####

write_tsv(soil.fungi.N.min.summary.tbl, 
          path = "data/working/soil.fungi.Nmin.gamm.summary.TableS2.txt")
write_tsv(roots.fungi.N.min.summary.tbl, 
          path = "data/working/roots.fungi.Nmin.gamm.summary.TableS2.txt")

# Get predicted values, ECM with peroxidases
soil.ECMLIG.Nmin.gamm.pred <- get_predictions(soil.ECMLIG.Nmin.gamm$gam, 
                                              cond = list(N.MIN = seq(min(soil.gamm.input.tbl$N.MIN, na.rm = T),
                                                                      max(soil.gamm.input.tbl$N.MIN, na.rm = T), 
                                                                      length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1, 
         FUNGAL.GROUP = "ECM.LIG") %>% 
  as_tibble()

# Get predicted values, ECM without peroxidases
soil.ECMNONLIG.Nmin.gamm.pred <- get_predictions(soil.ECMNONLIG.Nmin.gamm$gam, 
                                                 cond = list(N.MIN = seq(min(soil.gamm.input.tbl$N.MIN, na.rm = T),
                                                                         max(soil.gamm.input.tbl$N.MIN, na.rm = T), 
                                                                         length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1, 
         FUNGAL.GROUP = "ECM.NONLIG") %>% 
  as_tibble()

# Get predicted values, other mycorrhizas
soil.OTHERMYCO.Nmin.gamm.pred <- get_predictions(soil.OTHERMYCO.Nmin.gamm$gam, 
                                                 cond = list(N.MIN = seq(min(soil.gamm.input.tbl$N.MIN, na.rm = T),
                                                                         max(soil.gamm.input.tbl$N.MIN, na.rm = T), 
                                                                         length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1, 
         FUNGAL.GROUP = "OTHER.MYCORRHIZA") %>% 
  as_tibble()

# Get predicted values, ligninolytic saprotrophs
soil.SAPLIG.Nmin.gamm.pred <- get_predictions(soil.SAPLIG.Nmin.gamm$gam, 
                                              cond = list(N.MIN = seq(min(soil.gamm.input.tbl$N.MIN, na.rm = T),
                                                                      max(soil.gamm.input.tbl$N.MIN, na.rm = T), 
                                                                      length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1, 
         FUNGAL.GROUP = "SAP.LIG") %>% 
  as_tibble()

# Get predicted values, non-ligninolytic saprotrophs
soil.SAPNONLIG.Nmin.gamm.pred <- get_predictions(soil.SAPNONLIG.Nmin.gamm$gam, 
                                                 cond = list(N.MIN = seq(min(soil.gamm.input.tbl$N.MIN, na.rm = T),
                                                                         max(soil.gamm.input.tbl$N.MIN, na.rm = T), 
                                                                         length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1, 
         FUNGAL.GROUP = "SAP.NONLIG") %>% 
  as_tibble()

# Get predicted values, other
soil.OTHER.Nmin.gamm.pred <- get_predictions(soil.OTHER.Nmin.gamm$gam, 
                                             cond = list(N.MIN = seq(min(soil.gamm.input.tbl$N.MIN, na.rm = T),
                                                                     max(soil.gamm.input.tbl$N.MIN, na.rm = T), 
                                                                     length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1, 
         FUNGAL.GROUP = "OTHER") %>% 
  as_tibble()

# Combine original data for Figure 1
soil.fungi.N.min.input.tbl <- soil.gamm.input.tbl %>% 
  select(PLOT, N.MIN, ECM.LIG, ECM.NONLIG, OTHER.MYCORRHIZA, SAP.LIG, SAP.NONLIG, OTHER) %>% 
  pivot_longer(cols = ECM.LIG : OTHER, 
               names_to = "FUNGAL.GROUP", 
               values_to = "PROP") %>% 
  inner_join(., soil.fungi.N.min.summary.tbl, by = "FUNGAL.GROUP") %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Other mycorrhizas", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs", 
                                                      "Uncertain ecology")))

# Combine fit data for Figure 1
soil.fungi.Nmin.gamm.pred <- bind_rows(soil.ECMLIG.Nmin.gamm.pred, 
                                       soil.ECMNONLIG.Nmin.gamm.pred, 
                                       soil.OTHERMYCO.Nmin.gamm.pred, 
                                       soil.SAPLIG.Nmin.gamm.pred, 
                                       soil.SAPNONLIG.Nmin.gamm.pred, 
                                       soil.OTHER.Nmin.gamm.pred) %>% 
  inner_join(., soil.fungi.N.min.summary.tbl, by = "FUNGAL.GROUP") %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Other mycorrhizas", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs", 
                                                      "Uncertain ecology"))) %>% 
  # Remove non-significant fits
  mutate(fit = ifelse(PVAL.ADJ > 0.05, NA, fit))

# Pvals
soil.fungi.N.min.summary.input.tbl <- soil.fungi.N.min.summary.tbl %>% 
  #slice(rep(1:n(), each = 68)) %>% 
  mutate(PVAL.LABEL = rep(NA, nrow(.))) %>% 
  mutate(PVAL.LABEL = rep(NA, nrow(.))) %>% 
  mutate(XPOS = rep(0.9, 6), 
         YPOS = rep(0.65, 6), 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL.ADJ < 0.001, paste("atop(R[adj.]^2 == ", round(R2, 3), ", ", "italic(P[adj.]) < 0.001)", sep = ""), PVAL.LABEL), 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL.ADJ <= 0.05 & PVAL.ADJ >= 0.001, paste("atop(R[adj.]^2 == ", round(R2, 3), ", ", "italic(P[adj.]) == ", round(PVAL.ADJ, 3), ")", sep = ""), PVAL.LABEL), 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL.ADJ > 0.05, "atop(italic(n.s.), )", PVAL.LABEL)) %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Other mycorrhizas", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs", 
                                                      "Uncertain ecology")))

# Panel (a) for figure 1
soil.fungi.Nmin.plot <- ggplot() + 
  
  # Points
  geom_point(data = soil.fungi.N.min.input.tbl, 
             aes(x = N.MIN, y = PROP, colour = FACET.LABEL), 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Confidence intervals
  geom_ribbon(data = soil.fungi.Nmin.gamm.pred, 
              aes(x = N.MIN, y = yval, ymin = fit - CI, ymax = fit + CI, fill = FACET.LABEL), 
              show.legend = FALSE, 
              alpha = 0.3) + 
  
  # Line
  geom_line(data = soil.fungi.Nmin.gamm.pred, 
            aes(x = N.MIN, y = fit, colour = FACET.LABEL), 
            show.legend = FALSE, 
            size = 1) + 
  
  # Scale colors
  scale_colour_viridis_d(begin = 0.1, 
                         end = 0.44, 
                         option = "plasma", 
                         guide = FALSE) + 
  
  scale_fill_viridis_d(begin = 0.1, 
                       end = 0.44, 
                       option = "plasma", 
                       guide = FALSE) + 
  
  scale_y_continuous(limits = c(-0.125, 0.8)) + 
  
  # Pvalues
  geom_text(data = soil.fungi.N.min.summary.input.tbl, 
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL), 
            parse = TRUE, 
            stat = "identity", 
            size = 3) + 
  
  # Facet
  facet_wrap(~ FACET.LABEL, 
             nrow = 2, 
             ncol = 3, 
             scales = "free") + 
  
  # Titles
  labs(title = "(a)", 
       x = bquote('Inorganic N ('*mu*g~N%.%g^-1%.%d^-1*')'), 
       y = "Functional group relative abundance") + 
  
  # Format
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 12, hjust= -0.075), 
        axis.title.x = element_text(colour = "black", size = 8), 
        axis.title.y = element_text(colour = "black", size = 8), 
        # x axis text
        axis.text.x = element_text(colour = "black", size = 7), 
        axis.text.y = element_text(colour = "black", size = 7), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 8, colour = "black"))

#### 18. Fungi ~ Nmin, roots ####

# Get predicted values, ECM with peroxidases
roots.ECMLIG.Nmin.gamm.pred <- get_predictions(roots.ECMLIG.Nmin.gamm$gam, 
                                               cond = list(N.MIN = seq(min(roots.gamm.input.tbl$N.MIN, na.rm = T),
                                                                       max(roots.gamm.input.tbl$N.MIN, na.rm = T), 
                                                                       length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1, 
         FUNGAL.GROUP = "ECM.LIG") %>% 
  as_tibble()

# Get predicted values, ECM without peroxidases
roots.ECMNONLIG.Nmin.gamm.pred <- get_predictions(roots.ECMNONLIG.Nmin.gamm$gam, 
                                                  cond = list(N.MIN = seq(min(roots.gamm.input.tbl$N.MIN, na.rm = T),
                                                                          max(roots.gamm.input.tbl$N.MIN, na.rm = T), 
                                                                          length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1, 
         FUNGAL.GROUP = "ECM.NONLIG") %>% 
  as_tibble()

# Get predicted values, other mycorrhizas
roots.OTHERMYCO.Nmin.gamm.pred <- get_predictions(roots.OTHERMYCO.Nmin.gamm$gam, 
                                                  cond = list(N.MIN = seq(min(roots.gamm.input.tbl$N.MIN, na.rm = T),
                                                                          max(roots.gamm.input.tbl$N.MIN, na.rm = T), 
                                                                          length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1, 
         FUNGAL.GROUP = "OTHER.MYCORRHIZA") %>% 
  as_tibble()

# Get predicted values, ligninolytic saprotrophs
roots.SAPLIG.Nmin.gamm.pred <- get_predictions(roots.SAPLIG.Nmin.gamm$gam, 
                                               cond = list(N.MIN = seq(min(roots.gamm.input.tbl$N.MIN, na.rm = T),
                                                                       max(roots.gamm.input.tbl$N.MIN, na.rm = T), 
                                                                       length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1, 
         FUNGAL.GROUP = "SAP.LIG") %>% 
  as_tibble()

# Get predicted values, non-ligninolytic saprotrophs
roots.SAPNONLIG.Nmin.gamm.pred <- get_predictions(roots.SAPNONLIG.Nmin.gamm$gam, 
                                                  cond = list(N.MIN = seq(min(roots.gamm.input.tbl$N.MIN, na.rm = T),
                                                                          max(roots.gamm.input.tbl$N.MIN, na.rm = T), 
                                                                          length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1, 
         FUNGAL.GROUP = "SAP.NONLIG") %>% 
  as_tibble()

# Get predicted values, other
roots.OTHER.Nmin.gamm.pred <- get_predictions(roots.OTHER.Nmin.gamm$gam, 
                                              cond = list(N.MIN = seq(min(roots.gamm.input.tbl$N.MIN, na.rm = T),
                                                                      max(roots.gamm.input.tbl$N.MIN, na.rm = T), 
                                                                      length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  mutate(yval = 1, 
         FUNGAL.GROUP = "OTHER") %>% 
  as_tibble()

# Combine original data for Figure 1
roots.fungi.N.min.input.tbl <- roots.gamm.input.tbl %>% 
  select(PLOT, N.MIN, ECM.LIG, ECM.NONLIG, OTHER.MYCORRHIZA, SAP.LIG, SAP.NONLIG, OTHER) %>% 
  pivot_longer(cols = ECM.LIG : OTHER, 
               names_to = "FUNGAL.GROUP", 
               values_to = "PROP") %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Other mycorrhizas", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs", 
                                                      "Uncertain ecology")))

# Combine fit data for Figure 1
roots.fungi.Nmin.gamm.pred <- bind_rows(roots.ECMLIG.Nmin.gamm.pred, 
                                        roots.ECMNONLIG.Nmin.gamm.pred, 
                                        roots.OTHERMYCO.Nmin.gamm.pred, 
                                        roots.SAPLIG.Nmin.gamm.pred, 
                                        roots.SAPNONLIG.Nmin.gamm.pred, 
                                        roots.OTHER.Nmin.gamm.pred) %>% 
  inner_join(., roots.fungi.N.min.summary.tbl, by = "FUNGAL.GROUP") %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Other mycorrhizas", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs", 
                                                      "Uncertain ecology"))) %>% 
  # Remove non-significant fits
  mutate(fit = ifelse(PVAL.ADJ > 0.05, NA, fit))

# Pvals
roots.fungi.N.min.summary.input.tbl <- roots.fungi.N.min.summary.tbl %>% 
  mutate(PVAL.LABEL = rep(NA, nrow(.))) %>% 
  mutate(XPOS = rep(0.9, 6), 
         YPOS = rep(0.65, 6), 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL.ADJ < 0.001, paste("atop(R[adj.]^2 == ", round(R2, 3), ", ", "italic(P[adj.]) < 0.001)", sep = ""), PVAL.LABEL), 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL.ADJ <= 0.05 & PVAL.ADJ >= 0.001, paste("atop(R[adj.]^2 == ", round(R2, 3), ", ", "italic(P[adj.]) == ", round(PVAL.ADJ, 3), ")", sep = ""), PVAL.LABEL), 
         PVAL.LABEL = ifelse(!is.na(XPOS) & PVAL.ADJ > 0.05, "atop(italic(n.s.), )", PVAL.LABEL)) %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Other mycorrhizas", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs", 
                                                      "Uncertain ecology")))

# Panel (b) for figure 1
roots.fungi.Nmin.plot <- ggplot() + 
  
  # Points
  geom_point(data = roots.fungi.N.min.input.tbl, 
             aes(x = N.MIN, y = PROP, colour = FACET.LABEL), 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Confidence intervals
  geom_ribbon(data = roots.fungi.Nmin.gamm.pred, 
              aes(x = N.MIN, y = yval, ymin = fit - CI, ymax = fit + CI, fill = FACET.LABEL), 
              show.legend = FALSE, 
              alpha = 0.3) + 
  
  # Line
  geom_line(data = roots.fungi.Nmin.gamm.pred, 
            aes(x = N.MIN, y = fit, colour = FACET.LABEL), 
            show.legend = FALSE, 
            size = 1) + 
  
  # Scale colors
  scale_colour_viridis_d(begin = 0.5, 
                         end = 0.9, 
                         option = "plasma", 
                         guide = FALSE) + 
  
  scale_fill_viridis_d(begin = 0.5, 
                       end = 0.9, 
                       option = "plasma", 
                       guide = FALSE) + 
  
  scale_y_continuous(limits = c(-0.125, 0.8)) + 
  
  # Pvalues
  geom_text(data = roots.fungi.N.min.summary.input.tbl, 
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL), 
            parse = TRUE, 
            stat = "identity", 
            size = 3) + 
  
  # Facet
  facet_wrap(~ FACET.LABEL, 
             nrow = 2, 
             ncol = 3, 
             scales = "free") + 
  
  # Titles
  labs(title = "(b)", 
       x = bquote('Inorganic N ('*mu*g~N%.%g^-1%.%d^-1*')'), 
       y = "Functional group relative abundance") + 
  
  # Format
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 12, hjust= -0.075), 
        axis.title.x = element_text(colour = "black", size = 8), 
        axis.title.y = element_text(colour = "black", size = 8), 
        # x axis text
        axis.text.x = element_text(colour = "black", size = 6), 
        axis.text.y = element_text(colour = "black", size = 6), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 8, colour = "black"))

#### 19. Figure 2, fungi ~ Nmin ####

# Create figure 2
Figure2 <- ggarrange(soil.fungi.Nmin.plot, 
                     roots.fungi.Nmin.plot, 
                     ncol = 1, 
                     nrow = 2, 
                     align = "hv")

# Save figure
ggsave(filename = "figures/Fig_2.pdf",
       plot = Figure2, 
       device = "pdf",
       width = 6,
       height = 8,
       units = "in")

#### 20. lignin ~ fungi soil fig ####

# Summary
soil.lignin.fungi.summary <- summary(soil.lignin.fungi.gamm$gam)
soil.lignin.fungi.summary.tbl <- tibble(FUNGAL.GROUP = c("ECM.LIG", "ECM.NONLIG", "OTHER.MYCORRHIZA", "SAP.LIG", "SAP.NONLIG", "OTHER"), 
                                        FVAL = soil.lignin.fungi.summary$s.table[, 3], 
                                        PVAL = soil.lignin.fungi.summary$s.table[, 4]) %>% 
  mutate(PVAL.LABEL = rep(NA, nrow(.))) %>% 
  mutate(XPOS = c(0.45, 0.33, NA, 0.12, 0.55, NA), 
         YPOS = 0.0125, 
         PVAL.LABEL = ifelse(PVAL < 0.001, "italic(P) < 0.001)", PVAL.LABEL), 
         PVAL.LABEL = ifelse(PVAL <= 0.05 & PVAL >= 0.001, paste("italic(P) == ", round(PVAL, 3), sep = ""), PVAL.LABEL), 
         PVAL.LABEL = ifelse(PVAL > 0.05, "italic(n.s.)", PVAL.LABEL)) %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Other mycorrhizas", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs", 
                                                      "Uncertain ecology")))

soil.lignin.fungi.summary.input.tbl <- soil.lignin.fungi.summary.tbl %>% 
  filter(FUNGAL.GROUP != "OTHER.MYCORRHIZA" & FUNGAL.GROUP != "OTHER")

# Get predicted values, lignin ~ ECM with peroxidases
soil.lignin.ECMLIG.gamm.pred <- get_predictions(soil.lignin.fungi.gamm$gam, 
                                                cond = list(ECM.LIG = seq(min(soil.gamm.input.tbl$ECM.LIG, na.rm = T),
                                                                          max(soil.gamm.input.tbl$ECM.LIG, na.rm = T), 
                                                                          length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "ECM.LIG") %>% 
  select(-ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "ECM.LIG")

# Get predicted values, lignin ~ ECM without peroxidases
soil.lignin.ECMNONLIG.gamm.pred <- get_predictions(soil.lignin.fungi.gamm$gam, 
                                                   cond = list(ECM.NONLIG = seq(min(soil.gamm.input.tbl$ECM.NONLIG, na.rm = T),
                                                                                max(soil.gamm.input.tbl$ECM.NONLIG, na.rm = T), 
                                                                                length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "ECM.NONLIG") %>% 
  select(-ECM.LIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "ECM.NONLIG")

# Get predicted values, lignin ~ ligninolytic saprotrophs
soil.lignin.SAPLIG.gamm.pred <- get_predictions(soil.lignin.fungi.gamm$gam, 
                                                cond = list(SAP.LIG = seq(min(soil.gamm.input.tbl$SAP.LIG, na.rm = T),
                                                                          max(soil.gamm.input.tbl$SAP.LIG, na.rm = T), 
                                                                          length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "SAP.LIG") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "SAP.LIG")

# Get predicted values, lignin ~ non-ligninolytic saprotrophs
soil.lignin.SAPNONLIG.gamm.pred <- get_predictions(soil.lignin.fungi.gamm$gam, 
                                                   cond = list(SAP.NONLIG = seq(min(soil.gamm.input.tbl$SAP.NONLIG, na.rm = T),
                                                                                max(soil.gamm.input.tbl$SAP.NONLIG, na.rm = T), 
                                                                                length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "SAP.NONLIG") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "SAP.NONLIG")

# Combine initial data
soil.lignin.fungi.tbl <- soil.gamm.input.tbl %>% 
  select(PLOT, LIGNIN, ECM.LIG, ECM.NONLIG, SAP.LIG, SAP.NONLIG) %>% 
  pivot_longer(cols = ECM.LIG : SAP.NONLIG, 
               names_to = "INDEP.VARIABLE", 
               values_to = "INDEP.VALUE") %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs")))

# Combine fit data
soil.lignin.gamm.pred <- bind_rows(soil.lignin.ECMLIG.gamm.pred, 
                                   soil.lignin.ECMNONLIG.gamm.pred, 
                                   soil.lignin.SAPLIG.gamm.pred, 
                                   soil.lignin.SAPNONLIG.gamm.pred) %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs"))) %>% 
  # Remove non-significant fits
  mutate(fit = ifelse(INDEP.VARIABLE == "ECM.NONLIG", NA, fit), 
         fit = ifelse(INDEP.VARIABLE == "SAP.LIG", NA, fit), 
         fit = ifelse(INDEP.VARIABLE == "SAP.NONLIG", NA, fit))

# Figure 4a
soil.lignin.fungi.plot <- ggplot() + 
  
  # Points
  geom_point(data = soil.lignin.fungi.tbl, 
             aes(x = INDEP.VALUE, y = LIGNIN), 
             colour = "#270592FF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Confidence intervals
  geom_ribbon(data = soil.lignin.gamm.pred, 
              aes(x = INDEP.VALUE, y = yval, ymin = fit - CI, ymax = fit + CI), 
              show.legend = FALSE, 
              fill = "#270592FF", 
              alpha = 0.3) + 
  
  # Line
  geom_line(data = soil.lignin.gamm.pred, 
            aes(x = INDEP.VALUE, y = fit), 
            show.legend = FALSE, 
            color = "#270592FF", 
            size = 1) + 
  
  # Pvalues
  geom_text(data = soil.lignin.fungi.summary.input.tbl,
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL),
            parse = TRUE,
            stat = "identity",
            size = 3) +
  
  scale_y_continuous(limits = c(0, 0.25)) + 
  
  # Facet
  facet_wrap(~ FACET.LABEL, 
             nrow = 1, 
             ncol = 4, 
             scales = "free") + 
  
  # Titles
  labs(title = "(a)", 
       x = "Relative abundance of fungal functional group in soil", 
       y = "Lignin-derived SOM (prop. abund.)") + 
  
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

#### 21. soilC ~ fungi soil fig ####

# Summary
soil.soilC.fungi.summary <- summary(soil.soilC.fungi.gamm$gam)
soil.soilC.fungi.summary.tbl <- tibble(FUNGAL.GROUP = c("ECM.LIG", "ECM.NONLIG", "OTHER.MYCORRHIZA", "SAP.LIG", "SAP.NONLIG", "OTHER"), 
                                        FVAL = soil.soilC.fungi.summary$s.table[, 3], 
                                        PVAL = soil.soilC.fungi.summary$s.table[, 4]) %>% 
  mutate(PVAL.LABEL = rep(NA, nrow(.))) %>% 
  mutate(XPOS = c(0.45, 0.33, NA, 0.12, 0.55, NA), 
         YPOS = 1.65, 
         PVAL.LABEL = ifelse(PVAL < 0.001, "italic(P) < 0.001)", PVAL.LABEL), 
         PVAL.LABEL = ifelse(PVAL <= 0.05 & PVAL >= 0.001, paste("italic(P) == ", round(PVAL, 3), sep = ""), PVAL.LABEL), 
         PVAL.LABEL = ifelse(PVAL > 0.05, "italic(n.s.)", PVAL.LABEL)) %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Other mycorrhizas", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs", 
                                                      "Uncertain ecology")))

soil.soilC.fungi.summary.input.tbl <- soil.soilC.fungi.summary.tbl %>% 
  filter(FUNGAL.GROUP != "OTHER.MYCORRHIZA" & FUNGAL.GROUP != "OTHER")

# Get predicted values, lignin ~ ECM with peroxidases
soil.soilC.ECMLIG.gamm.pred <- get_predictions(soil.soilC.fungi.gamm$gam, 
                                                cond = list(ECM.LIG = seq(min(soil.gamm.input.tbl$ECM.LIG, na.rm = T),
                                                                          max(soil.gamm.input.tbl$ECM.LIG, na.rm = T), 
                                                                          length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "ECM.LIG") %>% 
  select(-ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "ECM.LIG")

# Get predicted values, soilC ~ ECM without peroxidases
soil.soilC.ECMNONLIG.gamm.pred <- get_predictions(soil.soilC.fungi.gamm$gam, 
                                                   cond = list(ECM.NONLIG = seq(min(soil.gamm.input.tbl$ECM.NONLIG, na.rm = T),
                                                                                max(soil.gamm.input.tbl$ECM.NONLIG, na.rm = T), 
                                                                                length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "ECM.NONLIG") %>% 
  select(-ECM.LIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "ECM.NONLIG")

# Get predicted values, lignin ~ ligninolytic saprotrophs
soil.soilC.SAPLIG.gamm.pred <- get_predictions(soil.soilC.fungi.gamm$gam, 
                                                cond = list(SAP.LIG = seq(min(soil.gamm.input.tbl$SAP.LIG, na.rm = T),
                                                                          max(soil.gamm.input.tbl$SAP.LIG, na.rm = T), 
                                                                          length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "SAP.LIG") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "SAP.LIG")

# Get predicted values, lignin ~ non-ligninolytic saprotrophs
soil.soilC.SAPNONLIG.gamm.pred <- get_predictions(soil.soilC.fungi.gamm$gam, 
                                                   cond = list(SAP.NONLIG = seq(min(soil.gamm.input.tbl$SAP.NONLIG, na.rm = T),
                                                                                max(soil.gamm.input.tbl$SAP.NONLIG, na.rm = T), 
                                                                                length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "SAP.NONLIG") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "SAP.NONLIG")

# Combine initial data
soil.soilC.fungi.tbl <- soil.gamm.input.tbl %>% 
  select(PLOT, LOG10.SOILC, ECM.LIG, ECM.NONLIG, SAP.LIG, SAP.NONLIG) %>% 
  pivot_longer(cols = ECM.LIG : SAP.NONLIG, 
               names_to = "INDEP.VARIABLE", 
               values_to = "INDEP.VALUE") %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs")))

# Combine fit data
soil.soilC.gamm.pred <- bind_rows(soil.soilC.ECMLIG.gamm.pred, 
                                  soil.soilC.ECMNONLIG.gamm.pred, 
                                  soil.soilC.SAPLIG.gamm.pred, 
                                  soil.soilC.SAPNONLIG.gamm.pred) %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs"))) %>% 
  # Remove non-significant fits
  mutate(fit = ifelse(INDEP.VARIABLE == "ECM.NONLIG", NA, fit), 
         fit = ifelse(INDEP.VARIABLE == "SAP.LIG", NA, fit), 
         fit = ifelse(INDEP.VARIABLE == "SAP.NONLIG", NA, fit))

# Figure 4b
soil.soilC.fungi.plot <- ggplot() + 
  
  # Points
  geom_point(data = soil.soilC.fungi.tbl, 
             aes(x = INDEP.VALUE, y = LOG10.SOILC), 
             colour = "#9C179EFF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Confidence intervals
  geom_ribbon(data = soil.soilC.gamm.pred, 
              aes(x = INDEP.VALUE, y = yval, ymin = fit - CI, ymax = fit + CI), 
              show.legend = FALSE, 
              fill = "#9C179EFF", 
              alpha = 0.3) + 
  
  # Line
  geom_line(data = soil.soilC.gamm.pred, 
            aes(x = INDEP.VALUE, y = fit), 
            show.legend = FALSE, 
            color = "#9C179EFF", 
            size = 1) + 
  
  # Pvalues
  geom_text(data = soil.soilC.fungi.summary.input.tbl,
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL),
            parse = TRUE,
            stat = "identity",
            size = 3) +
  
  scale_y_continuous(limits = c(1, 1.8)) + 
  
  # Facet
  facet_wrap(~ FACET.LABEL, 
             nrow = 1, 
             ncol = 4, 
             scales = "free") + 
  
  # Titles
  labs(title = "(b)", 
       x = "Relative abundance of fungal functional group in soil", 
       y = bquote('Soil C ('*log[10]*'['*mg~C%.%g^-1*'])')) + 
  
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

#### 22. lignin ~ fungi roots fig ####

# Summary
roots.lignin.fungi.summary <- summary(roots.lignin.fungi.gamm$gam)
roots.lignin.fungi.summary.tbl <- tibble(FUNGAL.GROUP = c("ECM.LIG", "ECM.NONLIG", "OTHER.MYCORRHIZA", "SAP.LIG", "SAP.NONLIG", "OTHER"), 
                                        FVAL = roots.lignin.fungi.summary$s.table[, 3], 
                                        PVAL = roots.lignin.fungi.summary$s.table[, 4]) %>% 
  mutate(PVAL.LABEL = rep(NA, nrow(.))) %>% 
  mutate(XPOS = c(0.55, 0.225, NA, 0.62, 0.32, NA), 
         YPOS = 0.025, 
         PVAL.LABEL = ifelse(PVAL < 0.001, "italic(P) < 0.001)", PVAL.LABEL), 
         PVAL.LABEL = ifelse(PVAL <= 0.05 & PVAL >= 0.001, paste("italic(P) == ", round(PVAL, 3), sep = ""), PVAL.LABEL), 
         PVAL.LABEL = ifelse(PVAL > 0.05, "italic(n.s.)", PVAL.LABEL)) %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Other mycorrhizas", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs", 
                                                      "Uncertain ecology")))

roots.lignin.fungi.summary.input.tbl <- roots.lignin.fungi.summary.tbl %>% 
  filter(FUNGAL.GROUP != "OTHER.MYCORRHIZA" & FUNGAL.GROUP != "OTHER")

# Get predicted values, lignin ~ ECM with peroxidases
roots.lignin.ECMLIG.gamm.pred <- get_predictions(roots.lignin.fungi.gamm$gam, 
                                                cond = list(ECM.LIG = seq(min(roots.gamm.input.tbl$ECM.LIG, na.rm = T),
                                                                          max(roots.gamm.input.tbl$ECM.LIG, na.rm = T), 
                                                                          length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "ECM.LIG") %>% 
  select(-ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "ECM.LIG")

# Get predicted values, lignin ~ ECM without peroxidases
roots.lignin.ECMNONLIG.gamm.pred <- get_predictions(roots.lignin.fungi.gamm$gam, 
                                                   cond = list(ECM.NONLIG = seq(min(roots.gamm.input.tbl$ECM.NONLIG, na.rm = T),
                                                                                max(roots.gamm.input.tbl$ECM.NONLIG, na.rm = T), 
                                                                                length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "ECM.NONLIG") %>% 
  select(-ECM.LIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "ECM.NONLIG")

# Get predicted values, lignin ~ ligninolytic saprotrophs
roots.lignin.SAPLIG.gamm.pred <- get_predictions(roots.lignin.fungi.gamm$gam, 
                                                cond = list(SAP.LIG = seq(min(roots.gamm.input.tbl$SAP.LIG, na.rm = T),
                                                                          max(roots.gamm.input.tbl$SAP.LIG, na.rm = T), 
                                                                          length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "SAP.LIG") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "SAP.LIG")

# Get predicted values, lignin ~ non-ligninolytic saprotrophs
roots.lignin.SAPNONLIG.gamm.pred <- get_predictions(roots.lignin.fungi.gamm$gam, 
                                                   cond = list(SAP.NONLIG = seq(min(roots.gamm.input.tbl$SAP.NONLIG, na.rm = T),
                                                                                max(roots.gamm.input.tbl$SAP.NONLIG, na.rm = T), 
                                                                                length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "SAP.NONLIG") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "SAP.NONLIG")

# Combine initial data
roots.lignin.fungi.tbl <- roots.gamm.input.tbl %>% 
  select(PLOT, LIGNIN, ECM.LIG, ECM.NONLIG, SAP.LIG, SAP.NONLIG) %>% 
  pivot_longer(cols = ECM.LIG : SAP.NONLIG, 
               names_to = "INDEP.VARIABLE", 
               values_to = "INDEP.VALUE") %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs")))

# Combine fit data
roots.lignin.gamm.pred <- bind_rows(roots.lignin.ECMLIG.gamm.pred, 
                                    roots.lignin.ECMNONLIG.gamm.pred, 
                                    roots.lignin.SAPLIG.gamm.pred, 
                                    roots.lignin.SAPNONLIG.gamm.pred) %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs"))) %>% 
  # Remove non-significant fits
  mutate(fit = ifelse(INDEP.VARIABLE == "ECM.NONLIG", NA, fit), 
         fit = ifelse(INDEP.VARIABLE == "SAP.LIG", NA, fit), 
         fit = ifelse(INDEP.VARIABLE == "ECM.LIG", NA, fit))

# Figure 4c
roots.lignin.fungi.plot <- ggplot() + 
  
  # Points
  geom_point(data = roots.lignin.fungi.tbl, 
             aes(x = INDEP.VALUE, y = LIGNIN), 
             colour = "#E56B5DFF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Confidence intervals
  geom_ribbon(data = roots.lignin.gamm.pred, 
              aes(x = INDEP.VALUE, y = yval, ymin = fit - CI, ymax = fit + CI), 
              show.legend = FALSE, 
              fill = "#E56B5DFF", 
              alpha = 0.3) + 
  
  # Line
  geom_line(data = roots.lignin.gamm.pred, 
            aes(x = INDEP.VALUE, y = fit), 
            show.legend = FALSE, 
            color = "#E56B5DFF", 
            size = 1) + 
  
  
  # Pvalues
  geom_text(data = roots.lignin.fungi.summary.input.tbl,
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL),
            parse = TRUE,
            stat = "identity",
            size = 3) +
  
  scale_y_continuous(limits = c(0, 0.25)) + 
  
  # Facet
  facet_wrap(~ FACET.LABEL, 
             nrow = 1, 
             ncol = 4, 
             scales = "free") + 
  
  # Titles
  labs(title = "(c)", 
       x = "Relative abundance of fungal functional group in roots", 
       y = "Lignin-derived SOM (prop. abund.)") + 
  
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

#### 23. soilC ~ fungi roots fig ####

# Summary
roots.soilC.fungi.summary <- summary(roots.soilC.fungi.gamm$gam)
roots.soilC.fungi.summary.tbl <- tibble(FUNGAL.GROUP = c("ECM.LIG", "ECM.NONLIG", "OTHER.MYCORRHIZA", "SAP.LIG", "SAP.NONLIG", "OTHER"), 
                                       FVAL = roots.soilC.fungi.summary$s.table[, 3], 
                                       PVAL = roots.soilC.fungi.summary$s.table[, 4]) %>% 
  mutate(PVAL.LABEL = rep(NA, nrow(.))) %>% 
  mutate(XPOS = c(0.55, 0.225, NA, 0.62, 0.4, NA), 
         YPOS = 1.8, 
         PVAL.LABEL = ifelse(PVAL < 0.001, "italic(P) < 0.001)", PVAL.LABEL), 
         PVAL.LABEL = ifelse(PVAL <= 0.05 & PVAL >= 0.001, paste("italic(P) == ", round(PVAL, 3), sep = ""), PVAL.LABEL), 
         PVAL.LABEL = ifelse(PVAL > 0.05, "italic(n.s.)", PVAL.LABEL)) %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Other mycorrhizas", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs", 
                                                      "Uncertain ecology")))

roots.soilC.fungi.summary.input.tbl <- roots.soilC.fungi.summary.tbl %>% 
  filter(FUNGAL.GROUP != "OTHER.MYCORRHIZA" & FUNGAL.GROUP != "OTHER")

# Get predicted values, soilC ~ ECM with peroxidases
roots.soilC.ECMLIG.gamm.pred <- get_predictions(roots.soilC.fungi.gamm$gam, 
                                                 cond = list(ECM.LIG = seq(min(roots.gamm.input.tbl$ECM.LIG, na.rm = T),
                                                                           max(roots.gamm.input.tbl$ECM.LIG, na.rm = T), 
                                                                           length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "ECM.LIG") %>% 
  select(-ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "ECM.LIG")

# Get predicted values, soilC ~ ECM without peroxidases
roots.soilC.ECMNONLIG.gamm.pred <- get_predictions(roots.soilC.fungi.gamm$gam, 
                                                    cond = list(ECM.NONLIG = seq(min(roots.gamm.input.tbl$ECM.NONLIG, na.rm = T),
                                                                                 max(roots.gamm.input.tbl$ECM.NONLIG, na.rm = T), 
                                                                                 length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "ECM.NONLIG") %>% 
  select(-ECM.LIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "ECM.NONLIG")

# Get predicted values, soilC ~ ligninolytic saprotrophs
roots.soilC.SAPLIG.gamm.pred <- get_predictions(roots.soilC.fungi.gamm$gam, 
                                                 cond = list(SAP.LIG = seq(min(roots.gamm.input.tbl$SAP.LIG, na.rm = T),
                                                                           max(roots.gamm.input.tbl$SAP.LIG, na.rm = T), 
                                                                           length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "SAP.LIG") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "SAP.LIG")

# Get predicted values, soilC ~ non-ligninolytic saprotrophs
roots.soilC.SAPNONLIG.gamm.pred <- get_predictions(roots.soilC.fungi.gamm$gam, 
                                                    cond = list(SAP.NONLIG = seq(min(roots.gamm.input.tbl$SAP.NONLIG, na.rm = T),
                                                                                 max(roots.gamm.input.tbl$SAP.NONLIG, na.rm = T), 
                                                                                 length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "SAP.NONLIG") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "SAP.NONLIG")

# Combine initial data
roots.soilC.fungi.tbl <- roots.gamm.input.tbl %>% 
  select(PLOT, LOG10.SOILC, ECM.LIG, ECM.NONLIG, SAP.LIG, SAP.NONLIG) %>% 
  pivot_longer(cols = ECM.LIG : SAP.NONLIG, 
               names_to = "INDEP.VARIABLE", 
               values_to = "INDEP.VALUE") %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs")))

# Combine fit data
roots.soilC.gamm.pred <- bind_rows(roots.soilC.ECMLIG.gamm.pred, 
                                   roots.soilC.ECMNONLIG.gamm.pred, 
                                   roots.soilC.SAPLIG.gamm.pred, 
                                   roots.soilC.SAPNONLIG.gamm.pred) %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.LIG", "ECM with peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "ECM.NONLIG", "ECM without peroxidases", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.LIG", "Ligninolytic saprotrophs", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "SAP.NONLIG", "Non-ligninolytic saprotrophs", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("ECM with peroxidases", 
                                                      "ECM without peroxidases", 
                                                      "Ligninolytic saprotrophs", 
                                                      "Non-ligninolytic saprotrophs")))

# Figure 4d
roots.soilC.fungi.plot <- ggplot() + 
  
  # Points
  geom_point(data = roots.soilC.fungi.tbl, 
             aes(x = INDEP.VALUE, y = LOG10.SOILC), 
             colour = "#FBD424FF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Pvalues
  geom_text(data = roots.soilC.fungi.summary.input.tbl,
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL),
            parse = TRUE,
            stat = "identity",
            size = 3) +
  
  scale_y_continuous(limits = c(1, 1.8)) + 
  
  # Facet
  facet_wrap(~ FACET.LABEL, 
             nrow = 1, 
             ncol = 4, 
             scales = "free") + 
  
  # Titles
  labs(title = "(d)", 
       x = "Relative abundance of fungal functional group in roots", 
       y = bquote('Soil C ('*log[10]*'['*mg~C%.%g^-1*'])')) + 
  
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

#### Create figure 4 ####
Figure4 <- ggarrange(soil.lignin.fungi.plot, 
                     soil.soilC.fungi.plot, 
                     roots.lignin.fungi.plot, 
                     roots.soilC.fungi.plot, 
                     ncol = 1, 
                     nrow = 4, 
                     align = "hv")

# Save Figure 4
ggsave(filename = "figures/Fig_4.pdf",
       Figure4,
       device = "pdf",
       width = 6.5,
       height = 8.5,
       units = "in")

#### Table S8 ####

write_tsv(soil.lignin.fungi.summary.tbl, 
          path = "data/working/soil.lignin.fungi.summary.TableS8.txt")
write_tsv(soil.soilC.fungi.summary.tbl, 
          path = "data/working/soil.soilC.fungi.summary.TableS8.txt")
write_tsv(roots.lignin.fungi.summary.tbl, 
          path = "data/working/roots.lignin.fungi.summary.TableS8.txt")
write_tsv(roots.soilC.fungi.summary.tbl, 
          path = "data/working/roots.soilC.fungi.summary.TableS8.txt")

#### 24. Supplemental GAMM figure ####

# Summary
soil.lignin.fungi.summary.tbl <- soil.lignin.fungi.summary.tbl %>% 
  mutate(XPOS = c(0.55, 0.225, 0.35, 0.62, 0.4, 0.15), 
         YPOS = c(0.24))

soil.lignin.fungi.summary.supp.tbl <- soil.lignin.fungi.summary.tbl %>% 
  filter(FUNGAL.GROUP == "OTHER.MYCORRHIZA" | FUNGAL.GROUP == "OTHER")

# Get predicted values, slignin ~ other myco
soil.lignin.OTHER.MYCO.gamm.pred <- get_predictions(soil.lignin.fungi.gamm$gam, 
                                                    cond = list(OTHER.MYCORRHIZA = seq(min(soil.gamm.input.tbl$OTHER.MYCORRHIZA, na.rm = T),
                                                                                       max(soil.gamm.input.tbl$OTHER.MYCORRHIZA, na.rm = T), 
                                                                                       length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "OTHER.MYCORRHIZA") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -SAP.LIG, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "OTHER.MYCORRHIZA")

# Get predicted values, slignin ~ other
soil.lignin.OTHER.gamm.pred <- get_predictions(soil.lignin.fungi.gamm$gam, 
                                               cond = list(OTHER = seq(min(soil.gamm.input.tbl$OTHER, na.rm = T),
                                                                       max(soil.gamm.input.tbl$OTHER, na.rm = T), 
                                                                       length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "OTHER") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -SAP.NONLIG) %>% 
  rename(INDEP.VALUE = "OTHER")

# Combine initial data
soil.lignin.fungi.supp.tbl <- soil.gamm.input.tbl %>% 
  select(PLOT, LIGNIN, OTHER.MYCORRHIZA, OTHER) %>% 
  pivot_longer(cols = OTHER.MYCORRHIZA : OTHER, 
               names_to = "INDEP.VARIABLE", 
               values_to = "INDEP.VALUE") %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("Other mycorrhizas", 
                                                      "Uncertain ecology")))

# Combine fit data
soil.lignin.gamm.supp.pred <- bind_rows(soil.lignin.OTHER.MYCO.gamm.pred, 
                                        soil.lignin.OTHER.gamm.pred) %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("Other mycorrhizas", 
                                                      "Uncertain ecology")))

# Figure a
soil.lignin.fungi.supp.plot <- ggplot() + 
  
  # Points
  geom_point(data = soil.lignin.fungi.supp.tbl, 
             aes(x = INDEP.VALUE, y = LIGNIN), 
             colour = "#270592FF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Pvalues
  geom_text(data = soil.lignin.fungi.summary.supp.tbl,
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL),
            parse = TRUE,
            stat = "identity",
            size = 3) +
  
  scale_y_continuous(limits = c(0, 0.25)) + 
  
  # Facet
  facet_wrap(~ FACET.LABEL, 
             nrow = 1, 
             ncol = 2, 
             scales = "free") + 
  
  # Titles
  labs(title = "(a)", 
       x = "Relative abundance of fungal functional group in soil", 
       y = "Lignin-derived SOM (prop. abund.)") + 
  
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

# Summary
soil.soilC.fungi.summary.supp.tbl <- soil.soilC.fungi.summary.tbl %>% 
  mutate(XPOS = c(0.45, 0.33, 0.35, 0.12, 0.55, 0.15), 
         YPOS = 1.65)

soil.soilC.fungi.summary.supp.tbl <- soil.soilC.fungi.summary.supp.tbl %>% 
  filter(FUNGAL.GROUP == "OTHER.MYCORRHIZA" | FUNGAL.GROUP == "OTHER")

# Get predicted values, lignin ~ other myco
soil.soilC.OTHER.MYCO.gamm.supp.pred <- get_predictions(soil.soilC.fungi.gamm$gam, 
                                                        cond = list(OTHER.MYCORRHIZA = seq(min(soil.gamm.input.tbl$OTHER.MYCORRHIZA, na.rm = T),
                                                                                           max(soil.gamm.input.tbl$OTHER.MYCORRHIZA, na.rm = T), 
                                                                                           length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "OTHER.MYCORRHIZA") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -SAP.LIG, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "OTHER.MYCORRHIZA")

# Get predicted values, soilC ~ other
soil.soilC.OTHER.gamm.supp.pred <- get_predictions(soil.soilC.fungi.gamm$gam, 
                                                   cond = list(OTHER = seq(min(soil.gamm.input.tbl$OTHER, na.rm = T),
                                                                           max(soil.gamm.input.tbl$OTHER, na.rm = T), 
                                                                           length.out = nrow(soil.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "OTHER") %>% 
  select(-ECM.LIG,
         -ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -SAP.NONLIG) %>% 
  rename(INDEP.VALUE = "OTHER")

# Combine initial data
soil.soilC.fungi.supp.tbl <- soil.gamm.input.tbl %>% 
  select(PLOT, LOG10.SOILC, OTHER.MYCORRHIZA, OTHER) %>% 
  pivot_longer(cols = OTHER.MYCORRHIZA : OTHER, 
               names_to = "INDEP.VARIABLE", 
               values_to = "INDEP.VALUE") %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("Other mycorrhizas", 
                                                      "Uncertain ecology")))

# Combine fit data
soil.soilC.gamm.supp.pred <- bind_rows(soil.soilC.OTHER.gamm.supp.pred, 
                                       soil.soilC.OTHER.MYCO.gamm.supp.pred) %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("Other mycorrhizas", 
                                                      "Uncertain ecology")))

# Figure b
soil.soilC.fungi.supp.plot <- ggplot() + 
  
  # Points
  geom_point(data = soil.soilC.fungi.supp.tbl, 
             aes(x = INDEP.VALUE, y = LOG10.SOILC), 
             colour = "#9C179EFF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Pvalues
  geom_text(data = soil.soilC.fungi.summary.supp.tbl,
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL),
            parse = TRUE,
            stat = "identity",
            size = 3) +
  
  scale_y_continuous(limits = c(1, 1.8)) + 
  
  # Facet
  facet_wrap(~ FACET.LABEL, 
             nrow = 1, 
             ncol = 4, 
             scales = "free") + 
  
  # Titles
  labs(title = "(b)", 
       x = "Relative abundance of fungal functional group in soil", 
       y = bquote('Soil C ('*log[10]*'['*mg~C%.%g^-1*'])')) + 
  
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

# Summary
roots.lignin.fungi.summary.supp.tbl <- roots.lignin.fungi.summary.tbl %>% 
  mutate(XPOS = c(0.55, 0.225, 0.1, 0.62, 0.32, 0.2), 
         YPOS = 0.025) %>% 
  filter(FUNGAL.GROUP == "OTHER.MYCORRHIZA" | FUNGAL.GROUP == "OTHER")

# Get predicted values, lignin ~ other myco
roots.lignin.OTHER.MYCO.gamm.pred <- get_predictions(roots.lignin.fungi.gamm$gam, 
                                                     cond = list(OTHER.MYCORRHIZA = seq(min(roots.gamm.input.tbl$OTHER.MYCORRHIZA, na.rm = T),
                                                                                        max(roots.gamm.input.tbl$OTHER.MYCORRHIZA, na.rm = T), 
                                                                                        length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "OTHER.MYCORRHIZA") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -SAP.LIG, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "OTHER.MYCORRHIZA")

# Get predicted values, lignin ~ other
roots.lignin.OTHER.gamm.pred <- get_predictions(roots.lignin.fungi.gamm$gam, 
                                                cond = list(OTHER = seq(min(roots.gamm.input.tbl$OTHER, na.rm = T),
                                                                        max(roots.gamm.input.tbl$OTHER, na.rm = T), 
                                                                        length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "OTHER") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -SAP.NONLIG) %>% 
  rename(INDEP.VALUE = "OTHER")

# Combine initial data
roots.lignin.fungi.supp.tbl <- roots.gamm.input.tbl %>% 
  select(PLOT, LIGNIN, OTHER.MYCORRHIZA, OTHER) %>% 
  pivot_longer(cols = OTHER.MYCORRHIZA : OTHER, 
               names_to = "INDEP.VARIABLE", 
               values_to = "INDEP.VALUE") %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("Other mycorrhizas", 
                                                      "Uncertain ecology")))

# Combine fit data
roots.lignin.gamm.supp.pred <- bind_rows(roots.lignin.OTHER.MYCO.gamm.pred, 
                                    roots.lignin.OTHER.gamm.pred) %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("Other mycorrhizas", 
                                                      "Uncertain ecology"))) %>% 
  # Remove non-significant fits
  mutate(fit = ifelse(INDEP.VARIABLE == "OTHER.MYCORRHIZA", NA, fit))

# Figure c
roots.lignin.fungi.supp.plot <- ggplot() + 
  
  # Points
  geom_point(data = roots.lignin.fungi.supp.tbl, 
             aes(x = INDEP.VALUE, y = LIGNIN), 
             colour = "#E56B5DFF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Confidence intervals
  geom_ribbon(data = roots.lignin.gamm.supp.pred, 
              aes(x = INDEP.VALUE, y = yval, ymin = fit - CI, ymax = fit + CI), 
              show.legend = FALSE, 
              fill = "#E56B5DFF", 
              alpha = 0.3) + 
  
  # Line
  geom_line(data = roots.lignin.gamm.supp.pred, 
            aes(x = INDEP.VALUE, y = fit), 
            show.legend = FALSE, 
            color = "#E56B5DFF", 
            size = 1) + 
  
  
  # Pvalues
  geom_text(data = roots.lignin.fungi.summary.supp.tbl,
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL),
            parse = TRUE,
            stat = "identity",
            size = 3) +
  
  scale_y_continuous(limits = c(0, 0.25)) + 
  
  # Facet
  facet_wrap(~ FACET.LABEL, 
             nrow = 1, 
             ncol = 2, 
             scales = "free") + 
  
  # Titles
  labs(title = "(c)", 
       x = "Relative abundance of fungal functional group in roots", 
       y = "Lignin-derived SOM (prop. abund.)") + 
  
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

#### 23. soilC ~ fungi roots fig ####

# Summary
roots.soilC.fungi.summary.supp.tbl <- roots.soilC.fungi.summary.tbl %>% 
  mutate(XPOS = c(0.55, 0.225, 0.1, 0.62, 0.4, 0.25), 
         YPOS = 1.8) %>% 
  filter(FUNGAL.GROUP == "OTHER.MYCORRHIZA" | FUNGAL.GROUP == "OTHER")

# Get predicted values, soilC ~ other mycorrhiza
roots.soilC.OTHER.MYCO.gamm.pred <- get_predictions(roots.soilC.fungi.gamm$gam, 
                                                    cond = list(OTHER.MYCORRHIZA = seq(min(roots.gamm.input.tbl$OTHER.MYCORRHIZA, na.rm = T),
                                                                                       max(roots.gamm.input.tbl$OTHER.MYCORRHIZA, na.rm = T), 
                                                                                       length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "OTHER.MYCORRHIZA") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -SAP.LIG, 
         -SAP.NONLIG, 
         -OTHER) %>% 
  rename(INDEP.VALUE = "OTHER.MYCORRHIZA")

# Get predicted values, soilC ~ other
roots.soilC.OTHER.gamm.pred <- get_predictions(roots.soilC.fungi.gamm$gam, 
                                               cond = list(OTHER = seq(min(roots.gamm.input.tbl$OTHER, na.rm = T),
                                                                       max(roots.gamm.input.tbl$OTHER, na.rm = T), 
                                                                       length.out = nrow(roots.gamm.input.tbl)), se = T)) %>% 
  as_tibble(.) %>% 
  mutate(yval = 1, 
         INDEP.VARIABLE = "OTHER") %>% 
  select(-ECM.LIG, 
         -ECM.NONLIG, 
         -OTHER.MYCORRHIZA, 
         -SAP.LIG, 
         -SAP.NONLIG) %>% 
  rename(INDEP.VALUE = "OTHER")

# Combine initial data
roots.soilC.fungi.supp.tbl <- roots.gamm.input.tbl %>% 
  select(PLOT, LOG10.SOILC, OTHER.MYCORRHIZA, OTHER) %>% 
  pivot_longer(cols = OTHER.MYCORRHIZA : OTHER, 
               names_to = "INDEP.VARIABLE", 
               values_to = "INDEP.VALUE") %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("Other mycorrhizas", 
                                                      "Uncertain ecology")))

# Combine fit data
roots.soilC.gamm.supp.pred <- bind_rows(roots.soilC.OTHER.MYCO.gamm.pred, 
                                   roots.soilC.OTHER.gamm.pred) %>% 
  # Add label
  mutate(FACET.LABEL = NA) %>% 
  mutate(FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER.MYCORRHIZA", "Other mycorrhizas", FACET.LABEL), 
         FACET.LABEL = ifelse(INDEP.VARIABLE == "OTHER", "Uncertain ecology", FACET.LABEL)) %>% 
  # Factor
  mutate(FACET.LABEL = factor(FACET.LABEL, levels = c("Other mycorrhizas", 
                                                      "Uncertain ecology")))

# Figure d
roots.soilC.fungi.supp.plot <- ggplot() + 
  
  # Points
  geom_point(data = roots.soilC.fungi.supp.tbl, 
             aes(x = INDEP.VALUE, y = LOG10.SOILC), 
             colour = "#FBD424FF", 
             alpha = 0.5, 
             size = 1.5) + 
  
  # Pvalues
  geom_text(data = roots.soilC.fungi.summary.supp.tbl,
            aes(x = XPOS, y = YPOS, label = PVAL.LABEL),
            parse = TRUE,
            stat = "identity",
            size = 3) +
  
  scale_y_continuous(limits = c(1, 1.8)) + 
  
  # Facet
  facet_wrap(~ FACET.LABEL, 
             nrow = 1, 
             ncol = 2, 
             scales = "free") + 
  
  # Titles
  labs(title = "(d)", 
       x = "Relative abundance of fungal functional group in roots", 
       y = bquote('Soil C ('*log[10]*'['*mg~C%.%g^-1*'])')) + 
  
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

# Create figure S7
FigureS7 <- ggarrange(soil.lignin.fungi.supp.plot, 
                      soil.soilC.fungi.supp.plot, 
                      roots.lignin.fungi.supp.plot, 
                      roots.soilC.fungi.supp.plot, 
                      ncol = 2, 
                      nrow = 2, 
                      align = "hv")

# Save Figure S7
ggsave(filename = "figures/FigureS7.png",
       FigureS4,
       device = "png",
       dpi = 600,
       width = 10,
       height = 6.5,
       units = "in")

#### 24. Potential covariates with N min (Table S3) ####

lignin.cov.gamm <- gamm(LIGNIN ~ s(N.MIN, k = -1) + s(SOIL.PH, k = -1) + s(VWC.MEAN, k = -1) + s(TEMP.MEAN, k = -1), 
                        correlation = corSpatial(form = ~LON + LAT), 
                        data = soil.gamm.input.tbl)
lignin.cov.gamm.summary <- summary(lignin.cov.gamm$gam)

soilC.cov.gamm <- gamm(LOG10.SOILC ~ s(N.MIN, k = -1) + s(SOIL.PH, k = -1) + s(VWC.MEAN, k = -1) + s(TEMP.MEAN, k = -1), 
                       correlation = corSpatial(form = ~LON + LAT), 
                       data = soil.gamm.input.tbl)
soilC.cov.gamm.summary <- summary(soilC.cov.gamm$gam)

soil.ECMLIG.cov.gamm <- gamm(ECM.LIG ~ s(N.MIN, k = -1) + s(SOIL.PH, k = -1) + s(VWC.MEAN, k = -1) + s(TEMP.MEAN, k = -1), 
                             correlation = corSpatial(form = ~LON + LAT), 
                             data = soil.gamm.input.tbl)
soil.ECMLIG.cov.gamm.summary <- summary(soil.ECMLIG.cov.gamm$gam)

soil.SAPLIG.cov.gamm <- gamm(SAP.LIG ~ s(N.MIN, k = -1) + s(SOIL.PH, k = -1) + s(VWC.MEAN, k = -1) + s(TEMP.MEAN, k = -1), 
                             correlation = corSpatial(form = ~LON + LAT), 
                             data = soil.gamm.input.tbl)
soil.SAPLIG.cov.gamm.summary <- summary(soil.SAPLIG.cov.gamm$gam)

roots.ECMLIG.cov.gamm <- gamm(ECM.LIG ~ s(N.MIN, k = -1) + s(SOIL.PH, k = -1) + s(VWC.MEAN, k = -1) + s(TEMP.MEAN, k = -1), 
                              correlation = corSpatial(form = ~LON + LAT), 
                              data = roots.gamm.input.tbl)
roots.ECMLIG.cov.gamm.summary <- summary(roots.ECMLIG.cov.gamm$gam)

roots.SAPLIG.cov.gamm <- gamm(SAP.LIG ~ s(N.MIN, k = -1) + s(SOIL.PH, k = -1) + s(VWC.MEAN, k = -1) + s(TEMP.MEAN, k = -1), 
                              correlation = corSpatial(form = ~LON + LAT), 
                              data = roots.gamm.input.tbl)
roots.SAPLIG.cov.gamm.summary <- summary(roots.SAPLIG.cov.gamm$gam)

