#### 1. Set-up ####

# Load libraries
library(tidyverse)
library(phyloseq)
library(viridis)

# Functions
source("code/functions.R")

# Read in data
fnFs.qual.plots <- readRDS("data/Manistee_MiSeq_data/quality_output/fnFs.qual.rds")
track.tbl <- as_tibble(readRDS("data/Manistee_MiSeq_data/seq_summary/dada2.seq.summary.rds"), rownames = "SAMPLE.ID")
ps.tax.tbl <- readRDS("data/working/ps.tax.tbl.RData")
ASV.all.tbl <- readRDS("data/working/ASV.all.tbl.Rdata")
ps <- readRDS("data/Manistee_MiSeq_data/asv_data/phyloseq.rds")

# Genus-level and functional group data
source("code/calculate_genus_abundances.R")
source("code/calculate_functional_group_abundances.R")

#### 2. Seq summary (Table S4) ####

# Calculate raw reads
raw.reads.tbl <- tibble(RAW = unlist(lapply(fnFs.qual.plots, raw_seq_counts))) %>% 
  mutate(SAMPLE.ID = names(fnFs.qual.plots)) %>% 
  select(SAMPLE.ID, RAW)

# Get sample data from phyloseq object
manistee.track.tbl <- as.data.frame(sample_data(ps)) %>% 
  as_tibble(.) %>% 
  inner_join(., track.tbl, by = "SAMPLE.ID") %>% 
  inner_join(raw.reads.tbl, ., by = "SAMPLE.ID") %>% 
  filter(TYPE == "environmental") %>% 
  select(SAMPLE.ID, TYPE, PLOT, STAND, SUBSTRATE, RUN.TYPE, everything()) %>% 
  group_by(SUBSTRATE, STAND, PLOT) %>% 
  summarise(RAW = sum(RAW), 
            INPUT = sum(input), 
            FILTERED = sum(filtered), 
            DENOISEDF = sum(denoisedF), 
            NONCHIM = sum(nonchim)) %>% 
  ungroup() %>% 
  select(-STAND) %>% 
  # Drop unused plots (outliers from SOM data)
  filter(PLOT != "Plot_6" & PLOT != "Plot_19" & PLOT != "Plot_43" & PLOT != "Plot_60") %>% 
  # Calculate totals by substrate
  group_by(SUBSTRATE) %>% 
  summarise(RAW = sum(RAW), 
            HQ = sum(NONCHIM)) %>% 
  ungroup() %>% 
  mutate(PERC.REMAINING = 100 * (HQ / RAW))

#### 3. ASV summary (Table S5) ####

ASV.summary.tbl <- ASV.all.tbl %>% 
  filter(PLOT != "Plot_6" & PLOT != "Plot_19" & PLOT != "Plot_43" & PLOT != "Plot_60") %>% 
  group_by(SUBSTRATE, ASV) %>% 
  mutate(TOTAL.COUNT = sum(COUNT)) %>% 
  ungroup() %>% 
  filter(TOTAL.COUNT > 0) %>% 
  select(ASV) %>% 
  distinct()
length(ASV.summary.tbl$ASV)

ASV.substrate.summary.tbl <- ASV.all.tbl %>% 
  filter(PLOT != "Plot_6" & PLOT != "Plot_19" & PLOT != "Plot_43" & PLOT != "Plot_60") %>% 
  select(SUBSTRATE, ASV, COUNT) %>% 
  group_by(SUBSTRATE, ASV) %>% 
  mutate(TOTAL.COUNT = sum(COUNT)) %>% 
  ungroup() %>% 
  filter(TOTAL.COUNT > 0) %>% 
  select(SUBSTRATE, ASV) %>% 
  distinct() %>% 
  group_by(SUBSTRATE) %>% 
  summarise(ASVs = length(ASV))

#### 4. Phylum summary (Table S6) ####

# Phylum
phylum.tbl <- ASV.all.tbl %>% 
  filter(PLOT != "Plot_6" & PLOT != "Plot_19" & PLOT != "Plot_43" & PLOT != "Plot_60") %>% 
  inner_join(ps.tax.modified.tbl, ., by = "ASV") %>% 
  # Calculate by plot
  group_by(PHYLUM, SUBSTRATE) %>% 
  summarize(COUNT = sum(COUNT)) %>% 
  ungroup() %>% 
  select(SUBSTRATE, PHYLUM, COUNT) %>% 
  group_by(SUBSTRATE) %>% 
  mutate(TOTAL = sum(COUNT)) %>% 
  ungroup() %>% 
  mutate(PERC = 100 * (COUNT / TOTAL)) %>% 
  arrange(SUBSTRATE, desc(PERC))

write_tsv(phylum.tbl, 
          path = "data/working/phylum.summary.tableS6.txt")

#### 5. Genus summary (Table S7) ####

tax.genus.tbl <- ps.tax.modified.tbl %>% 
  select(PHYLUM, CLASS, GENUS) %>% 
  distinct()

# Genus
genus.tbl <- bind_rows(soil.genus.trim.tbl, roots.genus.trim.tbl) %>% 
  select(GENUS) %>% 
  inner_join(tax.genus.tbl, ., by = "GENUS") %>% 
  select(PHYLUM, CLASS, GENUS) %>% 
  distinct() %>% 
  arrange(PHYLUM, CLASS, GENUS)

write_tsv(genus.tbl, 
          path = "data/working/genus.tbl.tableS7.txt")

#### 6. Functional group summary (Fig. S5) ####

# Soil
soil.fungal.functional.groups.summary.tbl <- as_tibble(read_tsv("data/working/soil.fungal.functional.groups.classified.0.1.txt")) %>% 
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
  select(GENUS, FUNGAL.GROUP, GENUS.PERC) %>% 
  distinct() %>% 
  # Calculate functional group percentages
  group_by(FUNGAL.GROUP) %>% 
  summarise(PERC = sum(GENUS.PERC)) %>% 
  ungroup() %>% 
  arrange(desc(PERC)) %>% 
  mutate(FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with\nperoxidases", NA), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without\nperoxidases", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other\nmycorrhizas", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic\nsaprotrophs", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic\nsaprotrophs", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Other or\nuncertain ecology", FUNGAL.GROUP.LABEL)) %>% 
  mutate(FUNGAL.GROUP.LABEL = factor(FUNGAL.GROUP.LABEL, levels = FUNGAL.GROUP.LABEL)) %>% 
  mutate(SUBSTRATE = rep("Soil", nrow(.)))

# Roots
roots.fungal.functional.groups.summary.tbl <- as_tibble(read_tsv("data/working/roots.fungal.functional.groups.classified.0.1.txt")) %>% 
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
  select(GENUS, FUNGAL.GROUP, GENUS.PERC) %>% 
  distinct() %>% 
  # Calculate functional group percentages
  group_by(FUNGAL.GROUP) %>% 
  summarise(PERC = sum(GENUS.PERC)) %>% 
  ungroup() %>% 
  arrange(desc(PERC)) %>% 
  mutate(FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "ECM.LIG", "ECM with\nperoxidases", NA), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "ECM.NONLIG", "ECM without\nperoxidases", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "OTHER.MYCORRHIZA", "Other\nmycorrhizas", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "SAP.LIG", "Ligninolytic\nsaprotrophs", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "SAP.NONLIG", "Non-ligninolytic\nsaprotrophs", FUNGAL.GROUP.LABEL), 
         FUNGAL.GROUP.LABEL = ifelse(FUNGAL.GROUP == "OTHER", "Other or\nuncertain ecology", FUNGAL.GROUP.LABEL)) %>% 
  mutate(FUNGAL.GROUP.LABEL = factor(FUNGAL.GROUP.LABEL, levels = FUNGAL.GROUP.LABEL)) %>% 
  mutate(SUBSTRATE = rep("Decaying fine roots", nrow(.)))

# Combine
FigS5.data.tbl <- bind_rows(soil.fungal.functional.groups.summary.tbl, roots.fungal.functional.groups.summary.tbl) %>% 
  arrange(desc(PERC)) %>% 
  mutate(SUBSTRATE = factor(SUBSTRATE, levels = c("Soil", "Decaying fine roots"))) %>% 
  mutate(PERC.LABEL = format(round(PERC, 1)))

# Figure
FigS5 <- ggplot() + 
  
  # Bars
  geom_bar(data = FigS5.data.tbl, 
           aes(x = FUNGAL.GROUP.LABEL, y = PERC, colour = SUBSTRATE, fill = SUBSTRATE), 
           stat = "identity", 
           alpha = 0.5) + 
  
  # Split
  facet_wrap(~ SUBSTRATE, 
             scales = "free") + 
  
  # Scale color
  scale_colour_viridis(name = "Fungal group", 
                       discrete = TRUE, 
                       begin = 0.2, 
                       end = 0.8, 
                       option = "plasma", 
                       guide = FALSE) + 
  
  # Scale fill
  scale_fill_viridis(name = "Fungal group", 
                     discrete = TRUE, 
                     begin = 0.2, 
                     end = 0.8, 
                     option = "plasma", 
                     guide = FALSE) + 
  
  # Bar labels
  geom_text(data = FigS5.data.tbl, 
            aes(x = FUNGAL.GROUP.LABEL, y = PERC + 3, label = PERC.LABEL)) + 
  
  # Axes
  scale_y_continuous(limits = c(0, 36), 
                     expand = expansion(c(0, 0.05))) + 
  
  # Labels
  labs(y = "% of sequences") + 
  
  # Format panel
  theme(axis.line = element_line(colour = "black", size = 0.5), 
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 7, colour = "black"), 
        # Axis text
        axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 6, colour = "black"), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 12, colour = "black"))

# Save figure
ggsave(filename = "figures/FigureS5.png",
       plot = FigS5, 
       device = "png",
       dpi = 300,
       width = 173,
       height = 87,
       units = "mm")

