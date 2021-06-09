#### 1. Set-up ####

# This script calculates net N mineralization for soil collected in spring 2019

# Load libraries
library(tidyverse)
library(data.table)

#### 2. Format soil moisture data ####

# Read in soil moisture data
soil.moisture.tbl <- as_tibble(fread("data/Manistee_env_data/aq2/spring_2019_soil_moisture.txt", sep = "\t", header = TRUE)) %>% 
  # Calculate percent moisture
  mutate(GRAV = (FRESH - (DRY.TOT - TIN )) * (1 / FRESH)) %>% 
  # Calculate mean by plot
  group_by(STAND, PLOT) %>% 
  summarize(GRAV = mean(GRAV)) %>% 
  ungroup() %>% 
  # Rename stands and plots
  mutate(STAND = paste("Stand", STAND, sep = "_"), 
         PLOT = paste("Plot", PLOT, sep = "_")) %>% 
  # Drop stand
  dplyr::select(-STAND)

#### 3. Format extraction and incubation masses ####

# Read in extraction and incubation fresh masses
extrac.incub.mass.tbl <- as_tibble(fread("data/Manistee_env_data/aq2/spring_2019_N_mineralization_incubation_and_pre_extraction_mass.txt", sep = "\t", header = TRUE)) %>% 
  # Select only fresh masses
  filter(TYPE == "F") %>% 
  # Rename stands and plots
  mutate(STAND = paste("Stand", STAND, sep = "_"), 
         PLOT = paste("Plot", PLOT, sep = "_")) %>% 
  # Add gravimetric water contents
  inner_join(., soil.moisture.tbl, by = "PLOT") %>% 
  # Calculate pre-incubation extraction dry mass of soil
  mutate(PRE.INCUB.EXTRACT.DRY.MASS = (1 - GRAV) * PRE.INCUB.EXTRACT.MASS) %>% 
  # Calculate post-incubation extraction dry mass of soil
  mutate(POST.INCUB.EXTRACT.DRY.MASS = (1 - GRAV) * INCUB.MASS) %>% 
  # Add column
  mutate(POST.INCUB.EXTRACT.MASS = rep(NA, length(INCUB.MASS))) %>% 
  # Update masses
  mutate(POST.INCUB.EXTRACT.MASS = ifelse(ADD.H2O == "Yes", INCUB.MASS.H2O.ADDED, INCUB.MASS)) %>% 
  # Calculate adjusted gravimetric water content
  mutate(GRAV.CORR = (POST.INCUB.EXTRACT.MASS - POST.INCUB.EXTRACT.DRY.MASS) * (1 / POST.INCUB.EXTRACT.MASS)) %>% 
  # Trim
  dplyr::select(PLOT, GRAV, GRAV.CORR, PRE.INCUB.EXTRACT.DRY.MASS, POST.INCUB.EXTRACT.DRY.MASS)

# Read in incubation times
incub.schedule.tbl <- as_tibble(fread("data/Manistee_env_data/aq2/spring_2019_N_mineralization_incubation_schedule.txt", sep = "\t", header = TRUE)) %>% 
  # Rename columns
  rename(STAND = "Stand", 
         PLOT = "Plot", 
         START.DATE = "Date started", 
         START.TIME = "Time started", 
         END.DATE = "Date finished", 
         END.TIME = "Rounded time finished") %>% 
  # Fresh samples only
  filter(Type == "F") %>% 
  # Combine dates and times
  mutate(START = paste(START.DATE, START.TIME, sep = " "), 
         END = paste(END.DATE, END.TIME, sep = " ")) %>% 
  # Convert dates
  mutate(START = as.POSIXct(strptime(START, format = "%Y-%m-%d %H:%M", tz = "EST5EDT")), 
         END = as.POSIXct(strptime(END, format = "%Y-%m-%d %H:%M", tz = "EST5EDT"))) %>% 
  # Time interval
  mutate(TIME.INTERVAL = difftime(END, START, tz = "EST5EDT", units = "days")) %>% 
  # Convert to numeric
  mutate(TIME.INTERVAL = as.numeric(TIME.INTERVAL)) %>% 
  # Rename stands and plots
  mutate(PLOT = paste("Plot", PLOT, sep = "_")) %>% 
  # Trim
  dplyr::select(PLOT, TIME.INTERVAL)

#### 4. Format AQ2 data ####

# Calculate blank 1 correction factor
aq2.blank.tbl <- as_tibble(fread("data/Manistee_env_data/aq2/WA_aq2_run_from_19_07_23.csv", sep = ",", header = TRUE)) %>% 
  # Redo test column
  mutate(TEST = ifelse(Test == "Ammonia 10 KCl", "Ammonium", "Nitrate")) %>% 
  # Remove standards and analytical blanks
  filter(`Sample ID` != "Standard 1" & 
           `Sample ID` != "Standard 90" & 
           `Sample ID` != "Standard 91" & 
           `Sample ID` != "Standard 92" & 
           `Sample ID` != "Standard 93" & 
           `Sample ID` != "Standard 94" & 
           `Sample ID` != "Standard 0" & 
           `Sample ID` != "CCV" & 
           `Sample ID` != "CCB" & 
           `Sample ID` != "NO2 1 mg L" & 
           `Sample ID` != "NO3 1 mg L" & 
           `Sample ID` != "KCl analytical") %>% 
  # Remove unnecessary columns
  dplyr::select(-`Sample Details`, -Test, -Units, -Absorbance, -`QC Pro result`, -Operator, 
         -`Man Dil Factor`, -`Auto Dil Factor`, -`Date and Time`) %>% 
  # Rename columns
  rename(SAMPLE.ID = `Sample ID`, 
         RESULTS = "Results") %>% 
  # Select blanks
  filter(SAMPLE.ID == "KCl blank 1 m" | 
           SAMPLE.ID == "KCl blank 1 n" | 
           SAMPLE.ID == "KCl blank 2 m" | 
           SAMPLE.ID == "KCl blank 1 o" | 
           SAMPLE.ID == "KCl blank 1test a" | 
           SAMPLE.ID == "KCl blank 1test b" | 
           SAMPLE.ID == "KCl blank 1test c" | 
           SAMPLE.ID == "KCl blank 1test d") %>% 
  # Separate samples
  separate(SAMPLE.ID, into = c("STAND", "PLOT", "EXTRAC", "REPLICATE"), sep = " ", convert = TRUE) %>% 
  # Group
  group_by(STAND, PLOT, EXTRAC, TEST) %>% 
  # Calculate mean values by sample
  summarise(RESULTS = mean(RESULTS)) %>% 
  # Ungroup
  ungroup()

# Calculate ammonium blank ratio
amm.blank.ratio <- unlist(aq2.blank.tbl[aq2.blank.tbl$EXTRAC == "1test" & aq2.blank.tbl$TEST == "Ammonium", "RESULTS"], use.names = FALSE) / unlist(aq2.blank.tbl[aq2.blank.tbl$EXTRAC == "1" & aq2.blank.tbl$TEST == "Ammonium", "RESULTS"], use.names = FALSE)

# Calculate nitrate blank ratio
nit.blank.ratio <- unlist(aq2.blank.tbl[aq2.blank.tbl$EXTRAC == "1test" & aq2.blank.tbl$TEST == "Nitrate", "RESULTS"], use.names = FALSE) / unlist(aq2.blank.tbl[aq2.blank.tbl$EXTRAC == "1" & aq2.blank.tbl$TEST == "Nitrate", "RESULTS"], use.names = FALSE)

# Get blanks, re-run
aq2.7.blanks.tbl <- as_tibble(fread("data/Manistee_env_data/aq2/WA_aq2_run_from_19_07_23.csv", sep = ",", header = TRUE)) %>% 
  # Redo test column
  mutate(TEST = ifelse(Test == "Ammonia 10 KCl", "Ammonium", "Nitrate")) %>% 
  # Remove standards and analytical blanks
  filter(`Sample ID` != "Standard 1" & 
           `Sample ID` != "Standard 90" & 
           `Sample ID` != "Standard 91" & 
           `Sample ID` != "Standard 92" & 
           `Sample ID` != "Standard 93" & 
           `Sample ID` != "Standard 94" & 
           `Sample ID` != "Standard 0" & 
           `Sample ID` != "CCV" & 
           `Sample ID` != "CCB" & 
           `Sample ID` != "NO2 1 mgL" & 
           `Sample ID` != "NO3 1 mgL" & 
           `Sample ID` != "KCl analytical") %>% 
  # Remove unnecessary columns
  dplyr::select(-`Sample Details`, -Test, -Units, -Absorbance, -`QC Pro result`, -Operator, 
         -`Man Dil Factor`, -`Auto Dil Factor`, -`Date and Time`) %>% 
  # Rename columns
  rename(SAMPLE.ID = `Sample ID`, 
         RESULTS = "Results") %>% 
  # Get blanks only
  filter(SAMPLE.ID == "KCl blank 1 m" | 
           SAMPLE.ID == "KCl blank 1 n" | 
           SAMPLE.ID == "KCl blank 1 o" | 
           SAMPLE.ID == "KCl blank 2 m") %>% 
  # Update blank names
  mutate(SAMPLE.ID = gsub("\\ ", ".", SAMPLE.ID)) %>% 
  # Replace all periods
  mutate(SAMPLE.ID = gsub("\\.", "_", SAMPLE.ID)) %>% 
  # Separate samples
  separate(SAMPLE.ID, into = c("STAND", "PLOT", "EXTRAC", "REPLICATE"), sep = "_", convert = TRUE) %>% 
  # Group
  group_by(STAND, PLOT, EXTRAC, TEST) %>% 
  # Calculate mean values by sample
  summarise(RESULTS = mean(RESULTS)) %>% 
  # Ungroup
  ungroup()

# Read in AQ2 data, re-run
aq2.7.tbl <- as_tibble(fread("data/Manistee_env_data/aq2/WA_aq2_run_from_19_07_23.csv", sep = ",", header = TRUE)) %>% 
  # Redo test column
  mutate(TEST = ifelse(Test == "Ammonia 10 KCl", "Ammonium", "Nitrate")) %>% 
  # Remove standards and analytical blanks
  filter(`Sample ID` != "Standard 1" & 
           `Sample ID` != "Standard 90" & 
           `Sample ID` != "Standard 91" & 
           `Sample ID` != "Standard 92" & 
           `Sample ID` != "Standard 93" & 
           `Sample ID` != "Standard 94" & 
           `Sample ID` != "Standard 0" & 
           `Sample ID` != "CCV" & 
           `Sample ID` != "CCB" & 
           `Sample ID` != "NO2 1 mgL" & 
           `Sample ID` != "NO3 1 mgL" & 
           `Sample ID` != "KCl analytical") %>% 
  # Remove unnecessary columns
  dplyr::select(-`Sample Details`, -Test, -Units, -Absorbance, -`QC Pro result`, -Operator, 
         -`Man Dil Factor`, -`Auto Dil Factor`, -`Date and Time`) %>% 
  # Rename columns
  rename(SAMPLE.ID = `Sample ID`, 
         RESULTS = "Results") %>% 
  # Keep only samples that were re-run
  filter((SAMPLE.ID %in% c("22.31.post.c", "22.31.post.d", "24.37.post.c", "24.37.post.d", 
                           "41.49.post.c", "41.49.post.d", "100.67.pre.c", "100.67.pre.d",
                           "100.67.post.c", "100.67.post.d", "50.55.pre.c", "50.55.pre.d",
                           "50.55.post.c", "50.55.post.d", "58.61.pre.c", "58.61.pre.d", 
                           "58.61.post.c", "58.61.post.d"))) %>% 
  # Replace all periods
  mutate(SAMPLE.ID = gsub("\\.", "_", SAMPLE.ID)) %>% 
  # Separate samples
  separate(SAMPLE.ID, into = c("STAND", "PLOT", "EXTRAC", "REPLICATE"), sep = "_") %>% 
  # Create a new column for KCl values
  mutate(KCl.BLANK = rep(NA, length(RESULTS))) %>% 
  # Add KCl blank pre ammonium
  mutate(KCl.BLANK = ifelse((EXTRAC == "pre" & TEST == "Ammonium"), 
                            unlist(aq2.7.blanks.tbl[aq2.7.blanks.tbl$EXTRAC == "1" & aq2.7.blanks.tbl$TEST == "Ammonium", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank post ammonium
  mutate(KCl.BLANK = ifelse((EXTRAC == "post" & TEST == "Ammonium"), 
                            unlist(aq2.7.blanks.tbl[aq2.7.blanks.tbl$EXTRAC == "2" & aq2.7.blanks.tbl$TEST == "Ammonium", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank pre nitrate
  mutate(KCl.BLANK = ifelse((EXTRAC == "pre" & TEST == "Nitrate"), 
                            unlist(aq2.7.blanks.tbl[aq2.7.blanks.tbl$EXTRAC == "1" & aq2.7.blanks.tbl$TEST == "Nitrate", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank post nitrate
  mutate(KCl.BLANK = ifelse((EXTRAC == "post" & TEST == "Nitrate"), 
                            unlist(aq2.7.blanks.tbl[aq2.7.blanks.tbl$EXTRAC == "2" & aq2.7.blanks.tbl$TEST == "Nitrate", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK))

# Get blanks, run 1
aq2.1.blanks.tbl <- as_tibble(fread("data/Manistee_env_data/aq2/WA_aq2_run_from_19_07_15.csv", sep = ",", header = TRUE)) %>% 
  # Redo test column
  mutate(TEST = ifelse(Test == "Ammonia 10 KCl", "Ammonium", "Nitrate")) %>% 
  # Remove standards and analytical blanks
  filter(`Sample ID` != "Standard 1" & 
           `Sample ID` != "Standard 90" & 
           `Sample ID` != "Standard 91" & 
           `Sample ID` != "Standard 92" & 
           `Sample ID` != "Standard 93" & 
           `Sample ID` != "Standard 94" & 
           `Sample ID` != "Standard 0" & 
           `Sample ID` != "CCV" & 
           `Sample ID` != "CCB" & 
           `Sample ID` != "NO2 1 mgL" & 
           `Sample ID` != "NO3 1 mgL" & 
           `Sample ID` != "KCl analytical") %>% 
  # Remove unnecessary columns
  dplyr::select(-`Sample Details`, -Test, -Units, -Absorbance, -`QC Pro result`, -Operator, 
         -`Man Dil Factor`, -`Auto Dil Factor`, -`Date and Time`) %>% 
  # Rename columns
  rename(SAMPLE.ID = `Sample ID`, 
         RESULTS = "Results") %>% 
  # Get blanks only
  filter(SAMPLE.ID == "KCl blank 1 a" | 
           SAMPLE.ID == "KCl blank 2 a" | 
           SAMPLE.ID == "KCl blank 2 b") %>% 
  # Update blank names
  mutate(SAMPLE.ID = gsub("\\ ", ".", SAMPLE.ID)) %>% 
  # Replace all periods
  mutate(SAMPLE.ID = gsub("\\.", "_", SAMPLE.ID)) %>% 
  # Separate samples
  separate(SAMPLE.ID, into = c("STAND", "PLOT", "EXTRAC", "REPLICATE"), sep = "_", convert = TRUE) %>% 
  # Group
  group_by(STAND, PLOT, EXTRAC, TEST) %>% 
  # Calculate mean values by sample
  summarise(RESULTS = mean(RESULTS)) %>% 
  # Ungroup
  ungroup()

# Read in AQ2 data, run 1
aq2.1.tbl <- as_tibble(fread("data/Manistee_env_data/aq2/WA_aq2_run_from_19_07_15.csv", sep = ",", header = TRUE)) %>% 
  # Redo test column
  mutate(TEST = ifelse(Test == "Ammonia 10 KCl", "Ammonium", "Nitrate")) %>% 
  # Remove standards and analytical blanks
  filter(`Sample ID` != "Standard 1" & 
           `Sample ID` != "Standard 90" & 
           `Sample ID` != "Standard 91" & 
           `Sample ID` != "Standard 92" & 
           `Sample ID` != "Standard 93" & 
           `Sample ID` != "Standard 94" & 
           `Sample ID` != "Standard 0" & 
           `Sample ID` != "CCV" & 
           `Sample ID` != "CCB" & 
           `Sample ID` != "NO2 1 mgL" & 
           `Sample ID` != "NO3 1 mgL" & 
           `Sample ID` != "KCl analytical") %>% 
  # Remove unnecessary columns
  dplyr::select(-`Sample Details`, -Test, -Units, -Absorbance, -`QC Pro result`, -Operator, 
         -`Man Dil Factor`, -`Auto Dil Factor`, -`Date and Time`) %>% 
  # Rename columns
  rename(SAMPLE.ID = `Sample ID`, 
         RESULTS = "Results") %>% 
  # Remove samples that need to be re-run
  filter(!((SAMPLE.ID %in% c("22.31.post.a", "22.31.post.b", "24.37.post.a", "24.37.post.b", 
                             "41.49.post.a", "41.49.post.b", "KCl blank 1 b", "100.67.pre.b", 
                             "100.67.post.b", "100.67.post.a", "50.55.pre.b", "50.55.post.b", 
                             "58.61.pre.b", "58.61.post.b")) & (TEST == "Nitrate"))) %>% 
  # Update blank names
  mutate(SAMPLE.ID = gsub("\\ ", ".", SAMPLE.ID)) %>% 
  # Replace all periods
  mutate(SAMPLE.ID = gsub("\\.", "_", SAMPLE.ID)) %>% 
  # Separate samples
  separate(SAMPLE.ID, into = c("STAND", "PLOT", "EXTRAC", "REPLICATE"), sep = "_", convert = TRUE) %>% 
  # Create a new column for KCl values
  mutate(KCl.BLANK = rep(NA, length(RESULTS))) %>% 
  # Add KCl blank pre ammonium
  mutate(KCl.BLANK = ifelse((EXTRAC == "pre" & TEST == "Ammonium"), 
                            unlist(aq2.1.blanks.tbl[aq2.1.blanks.tbl$EXTRAC == "1" & aq2.1.blanks.tbl$TEST == "Ammonium", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank post ammonium
  mutate(KCl.BLANK = ifelse((EXTRAC == "post" & TEST == "Ammonium"), 
                            unlist(aq2.1.blanks.tbl[aq2.1.blanks.tbl$EXTRAC == "2" & aq2.1.blanks.tbl$TEST == "Ammonium", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank pre nitrate
  mutate(KCl.BLANK = ifelse((EXTRAC == "pre" & TEST == "Nitrate"), 
                            unlist(aq2.1.blanks.tbl[aq2.1.blanks.tbl$EXTRAC == "1" & aq2.1.blanks.tbl$TEST == "Nitrate", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank post nitrate
  mutate(KCl.BLANK = ifelse((EXTRAC == "post" & TEST == "Nitrate"), 
                            unlist(aq2.1.blanks.tbl[aq2.1.blanks.tbl$EXTRAC == "2" & aq2.1.blanks.tbl$TEST == "Nitrate", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Trim off blanks
  filter(STAND != "KCl") %>% 
  # Add re-run samples
  bind_rows(., aq2.7.tbl) %>% 
  # Correct blank readings for ammonium
  mutate(KCl.BLANK.CORR = ifelse((EXTRAC == "pre" & TEST == "Ammonium"), KCl.BLANK * amm.blank.ratio, KCl.BLANK)) %>% 
  # Correct blank readings for nitrate
  mutate(KCl.BLANK.CORR = ifelse((EXTRAC == "pre" & TEST == "Nitrate"), KCl.BLANK * nit.blank.ratio, KCl.BLANK)) %>% 
  # Correct the blank readings
  mutate(KCl.BLANK.FINAL = ifelse((EXTRAC == "pre"), KCl.BLANK.CORR, KCl.BLANK)) %>% 
  # Correct samples for blank readings
  mutate(CORR.RESULTS = RESULTS - KCl.BLANK.FINAL) %>% 
  # Group
  group_by(STAND, PLOT, EXTRAC, TEST) %>% 
  # Calculate mean values by sample
  summarise(CORR.RESULTS = mean(RESULTS)) %>% 
  # Ungroup
  ungroup() %>% 
  # Convert negative values to 0
  mutate(CORR.RESULTS = ifelse(CORR.RESULTS < 0, 0, CORR.RESULTS))

# Read in AQ2 data, run 2
aq2.2.tbl <- as_tibble(fread("data/Manistee_env_data/aq2/WA_aq2_run_from_19_07_16.csv", sep = ",", header = TRUE)) %>% 
  # Redo test column
  mutate(TEST = ifelse(Test == "Ammonia 10 KCl", "Ammonium", "Nitrate")) %>% 
  # Remove standards and analytical blanks
  filter(`Sample ID` != "Standard 1" & 
           `Sample ID` != "Standard 90" & 
           `Sample ID` != "Standard 91" & 
           `Sample ID` != "Standard 92" & 
           `Sample ID` != "Standard 93" & 
           `Sample ID` != "Standard 94" & 
           `Sample ID` != "Standard 0" & 
           `Sample ID` != "CCV" & 
           `Sample ID` != "CCB" & 
           `Sample ID` != "NO2 1 mg L" & 
           `Sample ID` != "NO3 1 mg L" & 
           `Sample ID` != "KCl analytical") %>% 
  # Remove unnecessary columns
  dplyr::select(-`Sample Details`, -Test, -Units, -Absorbance, -`QC Pro result`, -Operator, 
         -`Man Dil Factor`, -`Auto Dil Factor`, -`Date and Time`) %>% 
  # Rename columns
  rename(SAMPLE.ID = `Sample ID`, 
         RESULTS = "Results") %>% 
  # Update blank names
  mutate(SAMPLE.ID = gsub("\\ ", ".", SAMPLE.ID)) %>% 
  # Replace all periods
  mutate(SAMPLE.ID = gsub("\\.", "_", SAMPLE.ID)) %>% 
  # Separate samples
  separate(SAMPLE.ID, into = c("STAND", "PLOT", "EXTRAC", "REPLICATE"), sep = "_", convert = TRUE) %>% 
  # Group
  group_by(STAND, PLOT, EXTRAC, TEST) %>% 
  # Calculate mean values by sample
  summarise(RESULTS = mean(RESULTS)) %>% 
  # Ungroup
  ungroup() %>% 
  # Create a new column for KCl values
  mutate(KCl.BLANK = rep(NA, length(RESULTS))) %>% 
  # Add KCl blank pre ammonium
  mutate(KCl.BLANK = ifelse((EXTRAC == "pre" & TEST == "Ammonium"), 
                            unlist(.[.$EXTRAC == "1" & .$TEST == "Ammonium", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank post ammonium
  mutate(KCl.BLANK = ifelse((EXTRAC == "post" & TEST == "Ammonium"), 
                            unlist(.[.$EXTRAC == "2" & .$TEST == "Ammonium", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank pre nitrate
  mutate(KCl.BLANK = ifelse((EXTRAC == "pre" & TEST == "Nitrate"), 
                            unlist(.[.$EXTRAC == "1" & .$TEST == "Nitrate", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank post nitrate
  mutate(KCl.BLANK = ifelse((EXTRAC == "post" & TEST == "Nitrate"), 
                            unlist(.[.$EXTRAC == "2" & .$TEST == "Nitrate", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Trim off blanks
  filter(STAND != "KCl") %>% 
  # Correct blank readings for ammonium
  mutate(KCl.BLANK.CORR = ifelse((EXTRAC == "pre" & TEST == "Ammonium"), KCl.BLANK * amm.blank.ratio, KCl.BLANK)) %>% 
  # Correct blank readings for nitrate
  mutate(KCl.BLANK.CORR = ifelse((EXTRAC == "pre" & TEST == "Nitrate"), KCl.BLANK * nit.blank.ratio, KCl.BLANK)) %>% 
  # Correct the blank readings
  mutate(KCl.BLANK.FINAL = ifelse((EXTRAC == "pre"), KCl.BLANK.CORR, KCl.BLANK)) %>% 
  # Correct samples for blank readings
  mutate(CORR.RESULTS = RESULTS - KCl.BLANK.FINAL) %>% 
  # Convert negative values to 0
  mutate(CORR.RESULTS = ifelse(CORR.RESULTS < 0, 0, CORR.RESULTS)) %>% 
  # Trim
  dplyr::select(STAND, PLOT, EXTRAC, TEST, CORR.RESULTS)

# Read in AQ2 data, run 3
aq2.3.tbl <- as_tibble(fread("data/Manistee_env_data/aq2/WA_aq2_run_from_19_07_17.csv", sep = ",", header = TRUE)) %>% 
  # Redo test column
  mutate(TEST = ifelse(Test == "Ammonia 10 KCl", "Ammonium", "Nitrate")) %>% 
  # Remove standards and analytical blanks
  filter(`Sample ID` != "Standard 1" & 
           `Sample ID` != "Standard 90" & 
           `Sample ID` != "Standard 91" & 
           `Sample ID` != "Standard 92" & 
           `Sample ID` != "Standard 93" & 
           `Sample ID` != "Standard 94" & 
           `Sample ID` != "Standard 0" & 
           `Sample ID` != "CCV" & 
           `Sample ID` != "CCB" & 
           `Sample ID` != "NO2 1 mg L" & 
           `Sample ID` != "NO3 1 mg L" & 
           `Sample ID` != "KCl analytical") %>% 
  # Remove unnecessary columns
  dplyr::select(-`Sample Details`, -Test, -Units, -Absorbance, -`QC Pro result`, -Operator, 
         -`Man Dil Factor`, -`Auto Dil Factor`, -`Date and Time`) %>% 
  # Rename columns
  rename(SAMPLE.ID = `Sample ID`, 
         RESULTS = "Results") %>% 
  # Update blank names
  mutate(SAMPLE.ID = gsub("\\ ", ".", SAMPLE.ID)) %>% 
  # Replace all periods
  mutate(SAMPLE.ID = gsub("\\.", "_", SAMPLE.ID)) %>% 
  # Separate samples
  separate(SAMPLE.ID, into = c("STAND", "PLOT", "EXTRAC", "REPLICATE"), sep = "_", convert = TRUE) %>% 
  # Group
  group_by(STAND, PLOT, EXTRAC, TEST) %>% 
  # Calculate mean values by sample
  summarise(RESULTS = mean(RESULTS)) %>% 
  # Ungroup
  ungroup() %>% 
  # Create a new column for KCl values
  mutate(KCl.BLANK = rep(NA, length(RESULTS))) %>% 
  # Add KCl blank pre ammonium
  mutate(KCl.BLANK = ifelse((EXTRAC == "pre" & TEST == "Ammonium"), 
                            unlist(.[.$EXTRAC == "1" & .$TEST == "Ammonium", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank post ammonium
  mutate(KCl.BLANK = ifelse((EXTRAC == "post" & TEST == "Ammonium"), 
                            unlist(.[.$EXTRAC == "2" & .$TEST == "Ammonium", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank pre nitrate
  mutate(KCl.BLANK = ifelse((EXTRAC == "pre" & TEST == "Nitrate"), 
                            unlist(.[.$EXTRAC == "1" & .$TEST == "Nitrate", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank post nitrate
  mutate(KCl.BLANK = ifelse((EXTRAC == "post" & TEST == "Nitrate"), 
                            unlist(.[.$EXTRAC == "2" & .$TEST == "Nitrate", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Trim off blanks
  filter(STAND != "KCl") %>% 
  # Correct blank readings for ammonium
  mutate(KCl.BLANK.CORR = ifelse((EXTRAC == "pre" & TEST == "Ammonium"), KCl.BLANK * amm.blank.ratio, KCl.BLANK)) %>% 
  # Correct blank readings for nitrate
  mutate(KCl.BLANK.CORR = ifelse((EXTRAC == "pre" & TEST == "Nitrate"), KCl.BLANK * nit.blank.ratio, KCl.BLANK)) %>% 
  # Correct the blank readings
  mutate(KCl.BLANK.FINAL = ifelse((EXTRAC == "pre"), KCl.BLANK.CORR, KCl.BLANK)) %>% 
  # Correct samples for blank readings
  mutate(CORR.RESULTS = RESULTS - KCl.BLANK.FINAL) %>% 
  # Convert negative values to 0
  mutate(CORR.RESULTS = ifelse(CORR.RESULTS < 0, 0, CORR.RESULTS)) %>% 
  # Trim
  dplyr::select(STAND, PLOT, EXTRAC, TEST, CORR.RESULTS)

# Read in AQ2 data, run 4
aq2.4.tbl <- as_tibble(fread("data/Manistee_env_data/aq2/WA_aq2_run_from_19_07_18.csv", sep = ",", header = TRUE)) %>% 
  # Redo test column
  mutate(TEST = ifelse(Test == "Ammonia 10 KCl", "Ammonium", "Nitrate")) %>% 
  # Remove standards and analytical blanks
  filter(`Sample ID` != "Standard 1" & 
           `Sample ID` != "Standard 90" & 
           `Sample ID` != "Standard 91" & 
           `Sample ID` != "Standard 92" & 
           `Sample ID` != "Standard 93" & 
           `Sample ID` != "Standard 94" & 
           `Sample ID` != "Standard 0" & 
           `Sample ID` != "CCV" & 
           `Sample ID` != "CCB" & 
           `Sample ID` != "NO2 1 mg L" & 
           `Sample ID` != "NO3 1 mg L" & 
           `Sample ID` != "KCl analytical") %>% 
  # Remove unnecessary columns
  dplyr::select(-`Sample Details`, -Test, -Units, -Absorbance, -`QC Pro result`, -Operator, 
         -`Man Dil Factor`, -`Auto Dil Factor`, -`Date and Time`) %>% 
  # Rename columns
  rename(SAMPLE.ID = `Sample ID`, 
         RESULTS = "Results") %>% 
  # Update blank names
  mutate(SAMPLE.ID = gsub("\\ ", ".", SAMPLE.ID)) %>% 
  # Replace all periods
  mutate(SAMPLE.ID = gsub("\\.", "_", SAMPLE.ID)) %>% 
  # Separate samples
  separate(SAMPLE.ID, into = c("STAND", "PLOT", "EXTRAC", "REPLICATE"), sep = "_", convert = TRUE) %>% 
  # Group
  group_by(STAND, PLOT, EXTRAC, TEST) %>% 
  # Calculate mean values by sample
  summarise(RESULTS = mean(RESULTS)) %>% 
  # Ungroup
  ungroup() %>% 
  # Create a new column for KCl values
  mutate(KCl.BLANK = rep(NA, length(RESULTS))) %>% 
  # Add KCl blank pre ammonium
  mutate(KCl.BLANK = ifelse((EXTRAC == "pre" & TEST == "Ammonium"), 
                            unlist(.[.$EXTRAC == "1" & .$TEST == "Ammonium", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank post ammonium
  mutate(KCl.BLANK = ifelse((EXTRAC == "post" & TEST == "Ammonium"), 
                            unlist(.[.$EXTRAC == "2" & .$TEST == "Ammonium", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank pre nitrate
  mutate(KCl.BLANK = ifelse((EXTRAC == "pre" & TEST == "Nitrate"), 
                            unlist(.[.$EXTRAC == "1" & .$TEST == "Nitrate", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank post nitrate
  mutate(KCl.BLANK = ifelse((EXTRAC == "post" & TEST == "Nitrate"), 
                            unlist(.[.$EXTRAC == "2" & .$TEST == "Nitrate", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Trim off blanks
  filter(STAND != "KCl") %>% 
  # Correct blank readings for ammonium
  mutate(KCl.BLANK.CORR = ifelse((EXTRAC == "pre" & TEST == "Ammonium"), KCl.BLANK * amm.blank.ratio, KCl.BLANK)) %>% 
  # Correct blank readings for nitrate
  mutate(KCl.BLANK.CORR = ifelse((EXTRAC == "pre" & TEST == "Nitrate"), KCl.BLANK * nit.blank.ratio, KCl.BLANK)) %>% 
  # Correct the blank readings
  mutate(KCl.BLANK.FINAL = ifelse((EXTRAC == "pre"), KCl.BLANK.CORR, KCl.BLANK)) %>% 
  # Correct samples for blank readings
  mutate(CORR.RESULTS = RESULTS - KCl.BLANK.FINAL) %>% 
  # Convert negative values to 0
  mutate(CORR.RESULTS = ifelse(CORR.RESULTS < 0, 0, CORR.RESULTS))  %>% 
  # Trim
  dplyr::select(STAND, PLOT, EXTRAC, TEST, CORR.RESULTS)

# Read in AQ2 data, run 5
aq2.5.tbl <- as_tibble(fread("data/Manistee_env_data/aq2/WA_aq2_run_from_19_07_19.csv", sep = ",", header = TRUE)) %>% 
  # Redo test column
  mutate(TEST = ifelse(Test == "Ammonia 10 KCl", "Ammonium", "Nitrate")) %>% 
  # Remove standards and analytical blanks
  filter(`Sample ID` != "Standard 1" & 
           `Sample ID` != "Standard 90" & 
           `Sample ID` != "Standard 91" & 
           `Sample ID` != "Standard 92" & 
           `Sample ID` != "Standard 93" & 
           `Sample ID` != "Standard 94" & 
           `Sample ID` != "Standard 0" & 
           `Sample ID` != "CCV" & 
           `Sample ID` != "CCB" & 
           `Sample ID` != "NO2 1 mg L" & 
           `Sample ID` != "NO3 1 mg L" & 
           `Sample ID` != "KCl analytical") %>% 
  # Remove unnecessary columns
  dplyr::select(-`Sample Details`, -Test, -Units, -Absorbance, -`QC Pro result`, -Operator, 
         -`Man Dil Factor`, -`Auto Dil Factor`, -`Date and Time`) %>% 
  # Rename columns
  rename(SAMPLE.ID = `Sample ID`, 
         RESULTS = "Results") %>% 
  # Update blank names
  mutate(SAMPLE.ID = gsub("\\ ", ".", SAMPLE.ID)) %>% 
  # Replace all periods
  mutate(SAMPLE.ID = gsub("\\.", "_", SAMPLE.ID)) %>% 
  # Separate samples
  separate(SAMPLE.ID, into = c("STAND", "PLOT", "EXTRAC", "REPLICATE"), sep = "_", convert = TRUE) %>% 
  # Group
  group_by(STAND, PLOT, EXTRAC, TEST) %>% 
  # Calculate mean values by sample
  summarise(RESULTS = mean(RESULTS)) %>% 
  # Ungroup
  ungroup() %>% 
  # Create a new column for KCl values
  mutate(KCl.BLANK = rep(NA, length(RESULTS))) %>% 
  # Add KCl blank pre ammonium
  mutate(KCl.BLANK = ifelse((EXTRAC == "pre" & TEST == "Ammonium"), 
                            unlist(.[.$EXTRAC == "1" & .$TEST == "Ammonium", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank post ammonium
  mutate(KCl.BLANK = ifelse((EXTRAC == "post" & TEST == "Ammonium"), 
                            unlist(.[.$EXTRAC == "2" & .$TEST == "Ammonium", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank pre nitrate
  mutate(KCl.BLANK = ifelse((EXTRAC == "pre" & TEST == "Nitrate"), 
                            unlist(.[.$EXTRAC == "1" & .$TEST == "Nitrate", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank post nitrate
  mutate(KCl.BLANK = ifelse((EXTRAC == "post" & TEST == "Nitrate"), 
                            unlist(.[.$EXTRAC == "2" & .$TEST == "Nitrate", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Trim off blanks
  filter(STAND != "KCl") %>% 
  # Correct blank readings for ammonium
  mutate(KCl.BLANK.CORR = ifelse((EXTRAC == "pre" & TEST == "Ammonium"), KCl.BLANK * amm.blank.ratio, KCl.BLANK)) %>% 
  # Correct blank readings for nitrate
  mutate(KCl.BLANK.CORR = ifelse((EXTRAC == "pre" & TEST == "Nitrate"), KCl.BLANK * nit.blank.ratio, KCl.BLANK)) %>% 
  # Correct the blank readings
  mutate(KCl.BLANK.FINAL = ifelse((EXTRAC == "pre"), KCl.BLANK.CORR, KCl.BLANK)) %>% 
  # Correct samples for blank readings
  mutate(CORR.RESULTS = RESULTS - KCl.BLANK.FINAL) %>% 
  # Convert negative values to 0
  mutate(CORR.RESULTS = ifelse(CORR.RESULTS < 0, 0, CORR.RESULTS))  %>% 
  # Trim
  dplyr::select(STAND, PLOT, EXTRAC, TEST, CORR.RESULTS)

# Read in AQ2 data, run 5
aq2.6.tbl <- as_tibble(fread("data/Manistee_env_data/aq2/WA_aq2_run_from_19_07_22.csv", sep = ",", header = TRUE)) %>% 
  # Redo test column
  mutate(TEST = ifelse(Test == "Ammonia 10 KCl", "Ammonium", "Nitrate")) %>% 
  # Remove standards and analytical blanks
  filter(`Sample ID` != "Standard 1" & 
           `Sample ID` != "Standard 90" & 
           `Sample ID` != "Standard 91" & 
           `Sample ID` != "Standard 92" & 
           `Sample ID` != "Standard 93" & 
           `Sample ID` != "Standard 94" & 
           `Sample ID` != "Standard 0" & 
           `Sample ID` != "CCV" & 
           `Sample ID` != "CCB" & 
           `Sample ID` != "NO2 1 mg L" & 
           `Sample ID` != "NO3 1 mg L" & 
           `Sample ID` != "KCl analytical") %>% 
  # Remove unnecessary columns
  dplyr::select(-`Sample Details`, -Test, -Units, -Absorbance, -`QC Pro result`, -Operator, 
         -`Man Dil Factor`, -`Auto Dil Factor`, -`Date and Time`) %>% 
  # Rename columns
  rename(SAMPLE.ID = `Sample ID`, 
         RESULTS = "Results") %>% 
  # Update blank names
  mutate(SAMPLE.ID = gsub("\\ ", ".", SAMPLE.ID)) %>% 
  # Replace all periods
  mutate(SAMPLE.ID = gsub("\\.", "_", SAMPLE.ID)) %>% 
  # Separate samples
  separate(SAMPLE.ID, into = c("STAND", "PLOT", "EXTRAC", "REPLICATE"), sep = "_", convert = TRUE) %>% 
  # Group
  group_by(STAND, PLOT, EXTRAC, TEST) %>% 
  # Calculate mean values by sample
  summarise(RESULTS = mean(RESULTS)) %>% 
  # Ungroup
  ungroup() %>% 
  # Create a new column for KCl values
  mutate(KCl.BLANK = rep(NA, length(RESULTS))) %>% 
  # Add KCl blank pre ammonium
  mutate(KCl.BLANK = ifelse((EXTRAC == "pre" & TEST == "Ammonium"), 
                            unlist(.[.$EXTRAC == "1" & .$TEST == "Ammonium", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank post ammonium
  mutate(KCl.BLANK = ifelse((EXTRAC == "post" & TEST == "Ammonium"), 
                            unlist(.[.$EXTRAC == "2" & .$TEST == "Ammonium", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank pre nitrate
  mutate(KCl.BLANK = ifelse((EXTRAC == "pre" & TEST == "Nitrate"), 
                            unlist(.[.$EXTRAC == "1" & .$TEST == "Nitrate", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Add KCl blank post nitrate
  mutate(KCl.BLANK = ifelse((EXTRAC == "post" & TEST == "Nitrate"), 
                            unlist(.[.$EXTRAC == "2" & .$TEST == "Nitrate", "RESULTS"], use.names = FALSE), 
                            KCl.BLANK)) %>% 
  # Trim off blanks
  filter(STAND != "KCl") %>% 
  # Correct blank readings for ammonium
  mutate(KCl.BLANK.CORR = ifelse((EXTRAC == "pre" & TEST == "Ammonium"), KCl.BLANK * amm.blank.ratio, KCl.BLANK)) %>% 
  # Correct blank readings for nitrate
  mutate(KCl.BLANK.CORR = ifelse((EXTRAC == "pre" & TEST == "Nitrate"), KCl.BLANK * nit.blank.ratio, KCl.BLANK)) %>% 
  # Correct the blank readings
  mutate(KCl.BLANK.FINAL = ifelse((EXTRAC == "pre"), KCl.BLANK.CORR, KCl.BLANK)) %>% 
  # Correct samples for blank readings
  mutate(CORR.RESULTS = RESULTS - KCl.BLANK.FINAL) %>% 
  # Convert negative values to 0
  mutate(CORR.RESULTS = ifelse(CORR.RESULTS < 0, 0, CORR.RESULTS))  %>% 
  # Trim
  dplyr::select(STAND, PLOT, EXTRAC, TEST, CORR.RESULTS)

#### 5. Calculate net N mineralization ####

# Combine data from all runs
N.min.tbl <- bind_rows(aq2.1.tbl, aq2.2.tbl, aq2.3.tbl, aq2.4.tbl, aq2.5.tbl, aq2.6.tbl) %>% 
  # Group
  group_by(STAND, PLOT, EXTRAC) %>% 
  # Calculate total inorganic N by plot and extraction
  summarise(INORG.N = sum(CORR.RESULTS)) %>% 
  # Ungroup
  ungroup() %>% 
  # Wide format
  pivot_wider(names_from = "EXTRAC", values_from = "INORG.N") %>% 
  # Rename columns
  rename(PRE.INCUB.N = "pre", 
         POST.INCUB.N = "post") %>% 
  # Reorder columns
  dplyr::select(STAND, PLOT, PRE.INCUB.N, POST.INCUB.N) %>% 
  # Convert plot to numeric
  mutate(PLOT = as.numeric(PLOT)) %>% 
  # Rename stands and plots
  mutate(STAND = paste("Stand", STAND, sep = "_"), 
         PLOT = paste("Plot", PLOT, sep = "_")) %>% 
  # Add extraction masses
  inner_join(., extrac.incub.mass.tbl, by = "PLOT") %>% 
  # Add time intervals
  inner_join(., incub.schedule.tbl, by = "PLOT") %>% 
  # Calculate N concentration per g dry soil
  mutate(PRE.N.PER.MASS.SOIL = PRE.INCUB.N * (1 / 1000) * 60 * (1 / PRE.INCUB.EXTRACT.DRY.MASS) * 1000, 
         POST.N.PER.MASS.SOIL = POST.INCUB.N * (1 / 1000) * 60 * (1 / POST.INCUB.EXTRACT.DRY.MASS) * 1000) %>% 
  # Calculate mineralization rate
  mutate(N.MIN = (POST.N.PER.MASS.SOIL - PRE.N.PER.MASS.SOIL) * (1 / TIME.INTERVAL)) %>% 
  rename(N.INORG = "PRE.N.PER.MASS.SOIL") %>% 
  dplyr::select(PLOT, STAND, N.MIN, N.INORG)

