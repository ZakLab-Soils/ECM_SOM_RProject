#### 1. Set-up ####

# Load libraries
library(tidyverse)

# Source files
source("code/functions.R")
source("code/calculate_genus_abundances.R")

#### 2. Add funguild data to genera ####

# Format funguild data
funguild.db.tbl <- as_tibble(read_tsv("data/working/funguild_db_test.txt",col_names = FALSE)) %>% 
  pivot_longer(everything(), names_to = "ENTRY.NUMBER", values_to = "ENTRY") %>% 
  separate(col = ENTRY, into = c("ID", "TAX", "TAX.LEVEL", "TROPHIC.MODE", "GUILD", "CONF", "FORM", "TRAIT", "NOTES", "CITATION"), sep = " , ") %>% 
  mutate(ENTRY.NUMBER = str_replace(ENTRY.NUMBER, "\\X", "")) %>% 
  mutate(ID = str_replace_all(ID, "[^[:alnum:]]", "")) %>% 
  mutate(ID = str_replace(ID, "idoid", "")) %>% 
  mutate(TAX = str_replace(TAX, '\\ : ', "")) %>% 
  mutate(TAX = str_replace_all(TAX, '\\"', "")) %>% 
  mutate(TAX = str_replace(TAX, "taxon", "")) %>% 
  mutate(TAX.LEVEL = str_replace(TAX.LEVEL, "taxonomicLevel", "")) %>% 
  mutate(TAX.LEVEL = str_replace_all(TAX.LEVEL, "[^[:alnum:]]", "")) %>% 
  mutate(TAX.LEVEL = as.numeric(TAX.LEVEL)) %>% 
  mutate(TROPHIC.MODE = str_replace(TROPHIC.MODE, "trophicMode", "")) %>% 
  mutate(TROPHIC.MODE = str_replace_all(TROPHIC.MODE, '\\"', "")) %>% 
  mutate(TROPHIC.MODE = str_replace_all(TROPHIC.MODE, "\\ : ", "")) %>% 
  mutate(GUILD = str_replace(GUILD, "guild", "")) %>% 
  mutate(GUILD = str_replace_all(GUILD, '\\"', "")) %>% 
  mutate(GUILD = str_replace_all(GUILD, "\\ : ", "")) %>% 
  mutate(CONF = str_replace(CONF, "confidenceRanking", "")) %>% 
  mutate(CONF = str_replace_all(CONF, '\\"', "")) %>% 
  mutate(CONF = str_replace_all(CONF, "\\ : ", "")) %>% 
  mutate(FORM = str_replace(FORM, "growthForm", "")) %>% 
  mutate(FORM = str_replace_all(FORM, '\\"', "")) %>% 
  mutate(FORM = str_replace_all(FORM, "\\ : ", "")) %>% 
  mutate(TRAIT = str_replace(TRAIT, "trait", "")) %>% 
  mutate(TRAIT = str_replace_all(TRAIT, '\\"', "")) %>% 
  mutate(TRAIT = str_replace_all(TRAIT, "\\ : ", "")) %>% 
  mutate(NOTES = str_replace(NOTES, "notes", "")) %>% 
  mutate(NOTES = str_replace_all(NOTES, '\\"', "")) %>% 
  mutate(NOTES = str_replace_all(NOTES, "\\ : ", "")) %>% 
  mutate(CITATION = str_replace(CITATION, "citationSource", "")) %>% 
  mutate(CITATION = str_replace_all(CITATION, '\\"', "")) %>% 
  mutate(CITATION = str_replace_all(CITATION, "\\ : ", "")) %>% 
  filter(TAX.LEVEL == 13) %>% 
  select(-ENTRY.NUMBER, -ID, -TAX.LEVEL) %>% 
  dplyr::rename_with(.fn = ~ paste(., "funguild", sep = "_")) %>% 
  dplyr::rename(GENUS = "TAX_funguild")

# Summary table for ecological classification, soil
soil.genera.to.classify.tbl <- soil.genus.trim.tbl %>% 
  select(GENUS) %>% 
  distinct() %>% 
  inner_join(., soil.genus.perc.tbl, by = "GENUS") %>% 
  left_join(., funguild.db.tbl, by = "GENUS") %>% 
  inner_join(ps.tax.modified.tbl, ., by = "GENUS") %>% 
  select(-ASV, -SPECIES) %>% 
  distinct() %>% 
  arrange(desc(GENUS.PERC))

# Summary table for ecological classification, roots
roots.genera.to.classify.tbl <- roots.genus.trim.tbl %>% 
  select(GENUS) %>% 
  distinct() %>% 
  inner_join(., roots.genus.perc.tbl, by = "GENUS") %>% 
  left_join(., funguild.db.tbl, by = "GENUS") %>% 
  inner_join(ps.tax.modified.tbl, ., by = "GENUS") %>% 
  select(-ASV, -SPECIES) %>% 
  distinct() %>% 
  arrange(desc(GENUS.PERC))

# Write out tables for classification
# write_tsv(soil.genera.to.classify.tbl,
#           path = "data/working/soil.genera.to.classify.tbl.txt")
# 
# write_tsv(roots.genera.to.classify.tbl,
#           path = "data/working/roots.genera.to.classify.tbl.txt")

#### 3. Format function group data ####

# Read in soil functional groups
soil.fungal.functional.groups.tbl <- as_tibble(read_tsv("data/working/soil.fungal.functional.groups.classified.0.1.txt")) %>% 
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
  select(GENUS, FUNGAL.GROUP, PLOT, PROP, GENUS.PERC) %>% 
  # Calculate abundance of functional groups
  group_by(PLOT, FUNGAL.GROUP) %>% 
  summarise(PROP = sum(PROP)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = PLOT, 
              names_from = "FUNGAL.GROUP", 
              values_from = "PROP")

# Read in roots functional groups
roots.fungal.functional.groups.tbl <- as_tibble(read_tsv("data/working/roots.fungal.functional.groups.classified.0.1.txt")) %>% 
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
  select(GENUS, FUNGAL.GROUP, PLOT, PROP, GENUS.PERC) %>% 
  # Calculate abundance of functional groups
  group_by(PLOT, FUNGAL.GROUP) %>% 
  summarise(PROP = sum(PROP)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = PLOT, 
              names_from = "FUNGAL.GROUP", 
              values_from = "PROP")

# Save data
saveRDS(soil.fungal.functional.groups.tbl, 
        file = "data/working/soil.fungal.functional.groups.tbl.RData")
saveRDS(roots.fungal.functional.groups.tbl, 
        file = "data/working/roots.fungal.functional.groups.tbl.RData")

