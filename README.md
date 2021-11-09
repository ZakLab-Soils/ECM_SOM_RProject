README for ECM\_SOM\_RProject code accompanying the manuscript, *Decay
by ectomycorrhizal fungi couples soil organic matter to nitrogen
availability*
================
William A. Argiroff, Donald R. Zak, Peter T.
Pellitier, Rima A. Upchurch, and Julia P. Belke
6/7/2021

# I. Overview

`R` scripts for all analyses are contained in the `code` folder. All
environmental data are contained in the `Manistee_env_data` folder
within the `data` folder. Water content data are contained in the
`9_Water_content_data` folder in the `data` folder. After downloading
sequence files, these should be moved to the `Manistee_MiSeq_data`
folder in the `data` folder. Data produced at intermediate steps are
written to and read from the `working` folder in the `data` folder. All
figures are written to the `figures` folder in the main project
directory.

# II. Data processing

## A. Soil characteristics and SOM biochemistry

### 1. Net N mineralization (soil inorganic N availability)

We calculated net N mineralization using the code in
`spring_2019_net_N_mineralization.R`.

### 2. Soil total C and N

Code for formatting soil C and N LECO data is contained in
`spring_2019_soil_total_C_and_N.R`.

### 3. Soil water content, soil temperature, and soil pH

`environmental_summary.R` contains code to calculate soil water content
and temperature, and format soil pH data.

### 4. SOM biochemistry

`som_pre_root_pygcms.R` contains code to format SOM biochemistry data.

### 5. Plot locations

Geographic coordinates for each plot were formatted using the code in
`plot_coordinates.R`.

### 6. Identify outlier plots

Outlier plots were identified using the code in
`environmental_summary.R`.

## B. Fungal sequences

### 1. Quality filter sequences and quantify ASVs

We used the code in `manistee_dada2.R` to process ITS2 Illumina reads.
This process takes a while (overnight due to taxonomic classification of
reads). To shorten time needed for downstream analyses, this script
saves the filtered and classified sequence data as a `phyloseq` object
in a temporary folder (`asv_data`) in the `Manistee_Miseq_data` folder.

The script `calculate_ASV.R` contains code to produce an ASV table (plot
by ASV matrix with sequence counts for each ASV). These tables are saved
as `.Rdata` objects (`ps.tax.tbl.RData` and `ASV.all.tbl.Rdata`) to save
future computation time.

### 2. Genus-level abundances

We used the code in the `calculate_genus_abundances.R` script to
calculate the abundances of fungal genera from `ps.tax.tbl.RData` and
`ASV.all.tbl.Rdata`, drop outlier plots, and trim genera that were not
present in at least 5 plots and at least 0.1% of fungal sequences. These
calculations were made separately for soil and decaying fine root
litter.

### 3. Functional group abundances

We calculated functional groups for abundant genera using the code in
`calculate_genus_abundances.R`. This script reads in genus abundances
source from `calculate_genus_abundances.R`, a text file of the FUNGuild
database (`funguild_db_test.txt`), and manually-assigned functional
groups for each genus
(`roots.fungal.functional.groups.classified.0.1.txt` and
`soil.fungal.functional.groups.classified.0.1.txt`) as described in the
main text and Table **S1**.

## C. Compile data

We used the code in the script `compile_data.R` to combine soil data and
functional group data for GAMM analysis. Output is saved as
`soil.compiled.data.tbl.RData` and `roots.compiled.data.tbl.RData`.

# III. Statistical analyses

## A. GAMM

Code to run all GAMM is contained in `gamm.R`, which reads in
`soil.compiled.data.tbl.RData` and `roots.compiled.data.tbl.RData`.

## B. TITAN2

We ran TITAN2 analyses using the code contained in `titan.R`, which
reads in `soil.compiled.data.tbl.RData` and sources
`calculate_genus_abundances.R`. This script saves the output of TITAN
analyses as `roots.genus.hlr.trim.titan.rData` and
`soil.genus.hlr.trim.titan.rData` to save time in downstream analyses.

# IV. Figures and tables

## A. Main text figures

### 1. Fig. **2**

Code to produce Fig. **2** is contained in `gamm.R`.

### 2. Fig. **3**

We produced Fig. **3** using code contained in `titan_summary.R`.

### 3. Fig. **4**

Code to produce Fig. **4** is contained in `gamm.R`.

### 4. Fig. **5**

We produced Fig. **5** using code contained in `gamm.R`.

## B. Supplemental figures

### 1. Fig. **S1**

We produced Fig. **S1** using the code contained in `map.R`.

### 2. Fig. **S2**

Code to produce Fig. **S2** is contained in `environmental_summary.R`.

### 3. Fig. **S3**

We produced Fig. **S3** using the code contained in
`sequence_summary.R`.

### 4. Fig. **S4**

We produced Fig. **S4** using the code contained in
`environmental_summary.R`.

### 5. Fig. **S5**

Code to produce Fig. **S5** is contained in `sequence_summary.R`.

### 6. Fig. **S6**

We produced Fig. **S7** using code contained in `root_biomass.R`.

### 7. Fig. **S7**

We produced Fig. **S7** using code contained in `gamm.R`.

### 8. Fig **S8**

Code to produce Fig. **S8** is contained in `SOM_summary.R`.

## C. Supplemental tables

### 1. Table **S2**

Code to produce Table **S2** is contained in `gamm.R`.

### 2. Table **S3**

We produced Table **S3** using code contained in `gamm.R`.

### 3. Table **S4**

Code to produce Table **S4** is contained in `sequence_summary.R`.

### 4. Table **S5**

We produced Table **S5** using code contained in `sequence_summary.R`.

### 5. Table **S6**

Code to produce Table **S6** is contained in `sequence_summary.R`.

### 6. Table **S7**

We produced Table **S7** using code contained in `sequence_summary.R`.

### 7. Table **S8**

We produced Table **S8** using code contained in `gamm.R`.

### 8. Table **S9**

Code to produce Table **S9** is contained in `gamm.R`.
