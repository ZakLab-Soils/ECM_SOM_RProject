#### 1. Set up ####

# Load required packages
library(tidyverse)
library(viridis)

# Load data
soil.pH.tbl <- read_tsv(file = "data/Manistee_env_data/pH/Manistee_spring_2019_soil_pH.txt") %>% 
  as_tibble(.) %>% 
  mutate(PLOT = paste("Plot", PLOT, sep = "_"), 
         SOIL.PH = as.numeric(SOIL.PH))
source("code/functions.R")
source("code/spring_2019_net_N_mineralization.R")
source("code/spring_2019_soil_total_C_and_N.R")
source("code/som_pre_root_py_gcms.R")

# N min data
soil.env.data.tbl <- N.min.tbl %>% 
  select(-STAND) %>% 
  # Add soil C and N
  inner_join(., soil.leco.tbl, by = "PLOT") %>% 
  select(-STAND) %>% 
  # Add py-GC/MS data
  inner_join(., py.gcms.som.source.tbl, by = "PLOT")

#### 2. Soil water and temp data ####

# Water and temp
vwc.temp.tbl <- read_tsv("data/vwc.temp.txt") %>% 
  mutate(VWC.PLOT.MEAN = VWC.PLOT.MEAN / 100)

# Correction factor list
vwc.temp.corr.list <- vwc.temp.tbl %>% 
  filter(PLOT == "Plot_75" | PLOT == "Plot_76" | PLOT == "Plot_77" | PLOT == "Plot_78" | 
           PLOT == "Plot_79" | PLOT == "Plot_80" | PLOT == "Plot_81" | PLOT == "Plot_82" | 
           PLOT == "Plot_83" | PLOT == "Plot_84" | PLOT == "Plot_85" | PLOT == "Plot_86") %>% 
  group_by(STAND) %>% 
  group_split()
names(vwc.temp.corr.list) <- unique(vwc.temp.tbl$STAND)

# Write function to get coefficients from linear regressions to correct plot level data to match logger readings
vwc_temp_corr_coeff <- function(x) {
  
  # Get site
  tmp1 <- unique(x$STAND)
  
  # VWC regression
  tmp2 <- lm(VWC.LOGGER ~ VWC.PLOT.MEAN, data = x)
  # Summary
  tmp3 <- summary(tmp2)
  # Coefficients
  tmp4 <- tmp3$coefficients
  
  # Temp regression
  tmp5 <- lm(TEMP.LOGGER ~ TEMP.PLOT.MEAN, data = x)
  # Summary
  tmp6 <- summary(tmp5)
  # Coefficients
  tmp7 <- tmp6$coefficients
  
  # Compile data
  tmp8 <- tibble(STAND = tmp1, 
                 VWC.CORR.INTERCEPT = tmp4[1, 1], 
                 VWC.CORR.SLOPE = tmp4[2, 1], 
                 VWC.CORR.INTERCEPT.P = tmp4[1, 4], 
                 VWC.CORR.SLOPE.P = tmp4[2, 4], 
                 VWC.CORR.R2 = tmp3$r.squared, 
                 TEMP.CORR.INTERCEPT = tmp7[1, 1], 
                 TEMP.CORR.SLOPE = tmp7[2, 1], 
                 TEMP.CORR.INTERCEPT.P = tmp7[1, 4], 
                 TEMP.CORR.SLOPE.P = tmp7[2, 4], 
                 TEMP.CORR.R2 = tmp6$r.squared)
  
  # Return
  return(tmp8)
  
}

# Calculate regression coefficients for correction factors, and correct plot point measurements
vwc.temp.corr.tbl <- lapply(vwc.temp.corr.list, vwc_temp_corr_coeff) %>% 
  bind_rows() %>% 
  # Add plot level data
  inner_join(vwc.temp.tbl, ., by = "STAND") %>% 
  filter(PLOT != "Plot_75" & PLOT != "Plot_76" & PLOT != "Plot_77" & PLOT != "Plot_78" & 
           PLOT != "Plot_79" & PLOT != "Plot_80" & PLOT != "Plot_81" & PLOT != "Plot_82" & 
           PLOT != "Plot_83" & PLOT != "Plot_84" & PLOT != "Plot_85" & PLOT != "Plot_86") %>% 
  mutate(VWC.PLOT.MEAN.CORR = (VWC.CORR.SLOPE * VWC.PLOT.MEAN) + VWC.CORR.INTERCEPT, 
         TEMP.PLOT.MEAN.CORR = (TEMP.CORR.SLOPE * TEMP.PLOT.MEAN) + TEMP.CORR.INTERCEPT)

# Write function to get coefficients from linear regressions of plot vs. logger readings
vwc_temp_coeff <- function(x) {
  
  # Get plot
  tmp1 <- unique(x$PLOT)
  # Get site
  tmp1b <- unique(x$STAND)
  
  # VWC regression
  tmp2 <- lm(VWC.PLOT.MEAN.CORR ~ VWC.LOGGER, data = x)
  # Summary
  tmp3 <- summary(tmp2)
  # Coefficients
  tmp4 <- tmp3$coefficients
  
  # Temp regression
  tmp5 <- lm(TEMP.PLOT.MEAN.CORR ~ TEMP.LOGGER, data = x)
  # Summary
  tmp6 <- summary(tmp5)
  # Coefficients
  tmp7 <- tmp6$coefficients
  
  # Compile data
  tmp8 <- tibble(PLOT = tmp1, 
                 STAND = tmp1b, 
                 VWC.INTERCEPT = tmp4[1, 1], 
                 VWC.SLOPE = tmp4[2, 1], 
                 VWC.INTERCEPT.P = tmp4[1, 4], 
                 VWC.SLOPE.P = tmp4[2, 4], 
                 VWC.R2 = tmp3$r.squared, 
                 TEMP.INTERCEPT = tmp7[1, 1], 
                 TEMP.SLOPE = tmp7[2, 1], 
                 TEMP.INTERCEPT.P = tmp7[1, 4], 
                 TEMP.SLOPE.P = tmp7[2, 4], 
                 TEMP.R2 = tmp6$r.squared)
  
  # Return
  return(tmp8)
  
}

# Input list
vwc.temp.list <- vwc.temp.corr.tbl %>% 
  group_by(PLOT) %>% 
  group_split()
names(vwc.temp.list) <- unique(vwc.temp.corr.tbl$PLOT)

# Get coefficients for regressions predicting plot values from logger values
vwc.temp.coeff.tbl <- lapply(vwc.temp.list, vwc_temp_coeff) %>% 
  bind_rows()

# Read in hourly logger data
vwc.temp.plot.tbl <- read_tsv("data/vwc.temp.logger.txt") %>% 
  select(STAND, VWC, TEMP, DATE, TIME, JULIAN) %>% 
  rename(VWC.LOGGER.HOURLY = "VWC", 
         TEMP.LOGGER.HOURLY = "TEMP") %>% 
  # Add coefficients to calculate plot level hourly data
  inner_join(., vwc.temp.coeff.tbl, by = "STAND") %>% 
  # Calculate hourly values for plots
  mutate(VWC.PLOT.HOURLY = (VWC.SLOPE * VWC.LOGGER.HOURLY) + VWC.INTERCEPT, 
         TEMP.PLOT.HOURLY = (TEMP.SLOPE * TEMP.LOGGER.HOURLY) + TEMP.INTERCEPT) %>% 
  # Calculate means by plot over growing season
  group_by(PLOT) %>% 
  summarise(VWC.MEAN = mean(VWC.PLOT.HOURLY), 
            TEMP.MEAN = mean(TEMP.PLOT.HOURLY)) %>% 
  ungroup()

#### 3. Combine data ####

# Get data
env.input.tbl <- soil.env.data.tbl %>% 
  # Add soil pH
  inner_join(., soil.pH.tbl, by = "PLOT") %>% 
  inner_join(., vwc.temp.plot.tbl, by = "PLOT") %>% 
  select(PLOT, N.MIN, SPEC.TOTAL.C, SOIL.CN, SPEC.TOTAL.N, SOIL.PH, LIGNIN, N.BEARING, PROTEIN, VWC.MEAN, TEMP.MEAN) %>% 
  arrange(PLOT)

#### 4. Run PCA ####

# Convert to data frame
env.PCA.input.df <- env.input.tbl %>% 
  filter(PLOT != "Plot_6" & PLOT != "Plot_19" & PLOT != "Plot_43" & PLOT != "Plot_60") %>% 
  select(-PROTEIN, -N.BEARING, -LIGNIN) %>% 
  column_to_rownames(var = "PLOT") %>% 
  as.data.frame(.)

# Add row names
row.names(env.PCA.input.df) <- env.PCA.input.df$PLOT

# Run PCA
env.pca <- prcomp(env.PCA.input.df, scale. = TRUE)

# Summarize PCA
env.pca.summary <- summary(env.pca)

# Get variance explained
env.pca.pc1.var <- as.character(round((100 * env.pca.summary$importance[2, 1]), 1))
env.pca.pc2.var <- as.character(round((100 * env.pca.summary$importance[2, 2]), 1))

# Get biplot loadings
env.pca.biplot.tbl <- as_tibble(as.data.frame(env.pca.summary$rotation), rownames = "VARIABLE") %>% 
  mutate(VECTOR.LABEL = NA, 
         VECTOR.LABEL = ifelse(VARIABLE == "N.MIN", "Inorg. N", VECTOR.LABEL), 
         VECTOR.LABEL = ifelse(VARIABLE == "SPEC.TOTAL.C", "Soil C", VECTOR.LABEL), 
         VECTOR.LABEL = ifelse(VARIABLE == "SPEC.TOTAL.N", "Soil total N", VECTOR.LABEL), 
         VECTOR.LABEL = ifelse(VARIABLE == "SOIL.CN", "Soil C/N", VECTOR.LABEL), 
         VECTOR.LABEL = ifelse(VARIABLE == "SOIL.PH", "Soil pH", VECTOR.LABEL), 
         VECTOR.LABEL = ifelse(VARIABLE == "VWC.MEAN", "Soil water", VECTOR.LABEL), 
         VECTOR.LABEL = ifelse(VARIABLE == "TEMP.MEAN", "Soil temp.", VECTOR.LABEL))

# Get plot loadings
env.pca.points.tbl <- as_tibble(as.data.frame(env.pca.summary$x))

# Create plot
env.pca.plot <- ggplot() + 
  
  # Points
  geom_point(data = env.pca.points.tbl, 
             aes(x = PC1, y = PC2), 
             size = 3, 
             colour = "grey30", 
             alpha = 0.6) + 
  
  # Vectors
  geom_segment(data = env.pca.biplot.tbl, 
               aes(x = 0, y = 0, 
                   xend = PC1 * 9, yend = PC2 * 9), 
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = "black") + 
  
  # Vector label
  geom_text(data = env.pca.biplot.tbl, 
            aes(x = PC1 * 9.4, y = PC2 * 9.4, label = VECTOR.LABEL), 
            size = 2.5, 
            hjust = -0.25) + 
  
  # Axis titles
  labs(y = paste("PC 2 (", env.pca.pc2.var, "% of variance)", sep = ""), 
       x = paste("PC 1 (", env.pca.pc1.var, "% of variance)", sep = "")) + 
  
  # Format
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(colour = "black", size = 8), 
        axis.title.y = element_text(colour = "black", size = 8),  
        # x axis text
        axis.text.x = element_text(colour = "black", size = 6), 
        axis.text.y = element_text(colour = "black", size = 6))

# Save Figure S6
ggsave(filename = "figures/FigureS6.png",
       env.pca.plot,
       device = "png",
       dpi = 600,
       width = 4.5,
       height = 4.5,
       units = "in")

#### 5. Bar plot (Fig. S4) ####

# Stands
stand.tbl <- soil.env.data.tbl %>% 
  select(PLOT, STAND)

# Site means and SE
FigS4.barplot.mean.tbl <- env.input.tbl %>% 
  inner_join(., stand.tbl, by = "PLOT") %>% 
  select(PLOT, STAND, N.MIN, SOIL.CN, SPEC.TOTAL.C, SPEC.TOTAL.N, PROTEIN, N.BEARING, LIGNIN, SOIL.PH, VWC.MEAN, TEMP.MEAN) %>% 
  group_by(STAND) %>% 
  mutate(N.MIN.MEAN = mean(N.MIN)) %>% 
  ungroup() %>% 
  pivot_longer(cols = N.MIN : TEMP.MEAN, names_to = "VARIABLE", values_to = "VALUE") %>% 
  group_by(STAND, VARIABLE, N.MIN.MEAN) %>% 
  summarise(VARIABLE.MEAN = mean(VALUE), 
            VARIABLE.SE = se(VALUE)) %>% 
  ungroup() %>% 
  arrange(N.MIN.MEAN) %>% 
  mutate(STAND.LABEL = gsub("Stand_", "Stand ", STAND)) %>% 
  # Add a factor for variable label
  mutate(VARIABLE.LABEL = rep(NA, nrow(.))) %>% 
  mutate(VARIABLE.LABEL = ifelse(VARIABLE == "N.MIN", "Net~N~mineralization", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "SOIL.CN", "Soil~C/N", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "SOIL.PH", "Soil~pH", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "SPEC.TOTAL.C", "Soil~C", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "SPEC.TOTAL.N", "Soil~N", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "LIGNIN", "`Lignin-derived`~SOM", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "PROTEIN", "Proteinaceous~SOM", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "N.BEARING", "`Non-proteinaceous`~`N-SOM`", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "VWC.MEAN", "Soil~water~content", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "TEMP.MEAN", "Soil~temperature", VARIABLE.LABEL)) %>% 
  mutate(VARIABLE.LABEL = factor(VARIABLE.LABEL, levels = c("Net~N~mineralization", "Soil~C/N", "Soil~C", "Soil~N", "`Lignin-derived`~SOM", "Proteinaceous~SOM", "`Non-proteinaceous`~`N-SOM`", "Soil~pH", "Soil~water~content", "Soil~temperature"))) %>% 
  # Add a second line for variable label with units
  mutate(VARIABLE.LABEL2 = rep(NA, nrow(.))) %>% 
  mutate(VARIABLE.LABEL2 = ifelse(VARIABLE == "N.MIN", "(mu*g~N%.%g^-1%.%d^-1)", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "SOIL.CN", " ", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "SOIL.PH", " ", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "SPEC.TOTAL.C", "(mg~C%.%g^-1)", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "SPEC.TOTAL.N", "(mg~N%.%g^-1)", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "LIGNIN", "(prop.~abund.)", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "PROTEIN", "(prop.~abund.)", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "N.BEARING", "(prop.~abund.)", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "VWC.MEAN", "(theta)", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "TEMP.MEAN", "(degree*C)", VARIABLE.LABEL2)) %>% 
  mutate(STAND.LABEL = factor(STAND.LABEL, levels = unique(STAND.LABEL)))

# Plot level data
FigS4.barplot.tbl <- env.input.tbl %>% 
  inner_join(., stand.tbl, by = "PLOT") %>% 
  select(PLOT, STAND, N.MIN, SOIL.CN, SPEC.TOTAL.C, SPEC.TOTAL.N, PROTEIN, N.BEARING, LIGNIN, SOIL.PH, VWC.MEAN, TEMP.MEAN) %>% 
  group_by(STAND) %>% 
  mutate(N.MIN.MEAN = mean(N.MIN)) %>% 
  ungroup() %>% 
  pivot_longer(cols = N.MIN : TEMP.MEAN, names_to = "VARIABLE", values_to = "VALUE") %>% 
  mutate(STAND.LABEL = gsub("Stand_", "Stand ", STAND)) %>% 
  # Add a factor for variable label
  mutate(VARIABLE.LABEL = rep(NA, nrow(.))) %>% 
  mutate(VARIABLE.LABEL = ifelse(VARIABLE == "N.MIN", "Net~N~mineralization", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "SOIL.CN", "Soil~C/N", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "SOIL.PH", "Soil~pH", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "SPEC.TOTAL.C", "Soil~C", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "SPEC.TOTAL.N", "Soil~N", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "LIGNIN", "`Lignin-derived`~SOM", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "PROTEIN", "Proteinaceous~SOM", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "N.BEARING", "`Non-proteinaceous`~`N-SOM`", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "VWC.MEAN", "Soil~water~content", VARIABLE.LABEL), 
         VARIABLE.LABEL = ifelse(VARIABLE == "TEMP.MEAN", "Soil~temperature", VARIABLE.LABEL)) %>% 
  mutate(VARIABLE.LABEL = factor(VARIABLE.LABEL, levels = c("Net~N~mineralization", "Soil~C/N", "Soil~C", "Soil~N", "`Lignin-derived`~SOM", "Proteinaceous~SOM", "`Non-proteinaceous`~`N-SOM`", "Soil~pH", "Soil~water~content", "Soil~temperature"))) %>% 
  # Add a second line for variable label with units
  mutate(VARIABLE.LABEL2 = rep(NA, nrow(.))) %>% 
  mutate(VARIABLE.LABEL2 = ifelse(VARIABLE == "N.MIN", "(mu*g~N%.%g^-1%.%d^-1)", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "SOIL.CN", " ", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "SOIL.PH", " ", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "SPEC.TOTAL.C", "(mg~C%.%g^-1)", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "SPEC.TOTAL.N", "(mg~N%.%g^-1)", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "LIGNIN", "(prop.~abund.)", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "PROTEIN", "(prop.~abund.)", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "N.BEARING", "(prop.~abund.)", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "VWC.MEAN", "(theta)", VARIABLE.LABEL2), 
         VARIABLE.LABEL2 = ifelse(VARIABLE == "TEMP.MEAN", "(degree*C)", VARIABLE.LABEL2)) %>% 
  # Add a column of y axis plot positions for outliers
  mutate(OUTLIER.POSITION = rep(NA, nrow(.))) %>% 
  mutate(OUTLIER.POSITION = ifelse(PLOT == "Plot_6" & VARIABLE == "LIGNIN", VALUE, OUTLIER.POSITION), 
         OUTLIER.POSITION = ifelse(PLOT == "Plot_43" & VARIABLE == "LIGNIN", VALUE, OUTLIER.POSITION),
         OUTLIER.POSITION = ifelse(PLOT == "Plot_60" & VARIABLE == "PROTEIN", VALUE, OUTLIER.POSITION), 
         OUTLIER.POSITION = ifelse(PLOT == "Plot_19" & VARIABLE == "N.MIN", VALUE, OUTLIER.POSITION)) %>% 
  # Calculate the maximum value for each variable; used for scaling outlier label locations
  group_by(VARIABLE) %>% 
  mutate(MAX.VALUE = max(VALUE)) %>% 
  ungroup() %>% 
  # Column of outlier label positions
  mutate(LABEL.POSITION = rep(NA, nrow(.))) %>% 
  mutate(LABEL.POSITION = ifelse(PLOT == "Plot_6" & VARIABLE == "LIGNIN", VALUE + 0.09 * MAX.VALUE, LABEL.POSITION), 
         LABEL.POSITION = ifelse(PLOT == "Plot_43" & VARIABLE == "LIGNIN", VALUE + 0.09 * MAX.VALUE, LABEL.POSITION),
         LABEL.POSITION = ifelse(PLOT == "Plot_60" & VARIABLE == "PROTEIN", VALUE + 0.09 * MAX.VALUE, LABEL.POSITION), 
         LABEL.POSITION = ifelse(PLOT == "Plot_19" & VARIABLE == "N.MIN", VALUE + 0.09 * MAX.VALUE, LABEL.POSITION)) %>% 
  # Create outlier labels
  mutate(OUTLIER.LABEL = rep(NA, nrow(.))) %>% 
  mutate(OUTLIER.LABEL = ifelse(!is.na(OUTLIER.POSITION), gsub("\\_", " ", PLOT), OUTLIER.LABEL)) %>% 
  arrange(N.MIN.MEAN) %>% 
  mutate(STAND.LABEL = factor(STAND.LABEL, levels = unique(STAND.LABEL)))

# Plot
FigS4.barplot <- ggplot() + 
  
  # Add plot level data points
  geom_point(data = FigS4.barplot.tbl, 
             aes(x = STAND.LABEL, y = VALUE, colour = VARIABLE.LABEL), 
             size = 1, 
             alpha = 0.5) + 
  
  # Bars
  geom_bar(data = FigS4.barplot.mean.tbl, 
           aes(x = STAND.LABEL, y = VARIABLE.MEAN, fill = VARIABLE.LABEL, colour = VARIABLE.LABEL), 
           stat = "identity", 
           alpha = 0.5) + 
  
  # Add error bars
  geom_errorbar(data = FigS4.barplot.mean.tbl, 
                aes(x = STAND.LABEL, ymin = VARIABLE.MEAN - VARIABLE.SE, ymax = VARIABLE.MEAN + VARIABLE.SE, colour = VARIABLE.LABEL), 
                width = 0.2) + 
  
  # Add red circles for outliers
  geom_point(data = FigS4.barplot.tbl, 
             aes(x = STAND.LABEL, y = OUTLIER.POSITION), 
             pch = 21, 
             fill = NA, 
             colour = "red", 
             size = 3) + 
  
  # Add outlier labels
  geom_text(data = FigS4.barplot.tbl, 
            aes(x = STAND.LABEL, y = LABEL.POSITION, label = OUTLIER.LABEL), 
            nudge_x = -0.5, 
            size = 2) + 
  
  # Split into panels
  facet_wrap(~ VARIABLE.LABEL + VARIABLE.LABEL2, 
             ncol = 4, 
             nrow = 3, 
             scales = "free",  
             labeller = label_parsed) + 
  
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
  
  # X axis
  #scale_x_continuous(expand = expansion(mult = c(0.025, 0.025))) + 
  
  # Remove gap between x axis and bars
  scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
  
  # Format panel and strip
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        # x axis text
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, size = 6), 
        axis.text.y = element_text(colour = "black", size = 6), 
        strip.text = element_text(colour = "black", size = 7), 
        strip.background = element_blank())

# Save figure
ggsave(filename = "figures/FigureS4.png",
       plot = FigS4.barplot,
       device = "png",
       dpi = 600,
       width = 6.5,
       height = 6.5,
       units = "in")

