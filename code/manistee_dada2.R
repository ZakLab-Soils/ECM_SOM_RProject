#### 1. Set up environment ####

# Load packages
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

# Load source scripts
source("r/functions.R")

#### 2. Check initial read quality ####

# Directory containing fastq files
path <- "data/Manistee_MiSeq_data"

# Get a list of files in the directory
list.files(path)

# Create vectors of paths to each file
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Check initial quality of reads
fnFs.qual <- vector(mode = "list", length = length(fnFs))
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
names(fnFs.qual) <- sample.names

#for(i in seq_along(fnFs)) {
#  
#  print(paste("Processing:", sample.names[i], sep = " "))
#  fnFs.qual[[i]] <- plotQualityProfile(fnFs[i])
#  
#}

# Save plots
#saveRDS(fnFs.qual, file = "data/Manistee_MiSeq_data_v2/quality_output/fnFs.qual.rds")

# Check initial quality of reads
fnRs.qual <- vector(mode = "list", length = length(fnRs))
names(fnRs.qual) <- sample.names

#for(i in seq_along(fnRs)) {
#  
#  print(paste("Processing:", sample.names[i], sep = " "))
#  fnRs.qual[[i]] <- plotQualityProfile(fnRs[i])
#  
#}

# Save plots
#saveRDS(fnRs.qual, file = "data/Manistee_MiSeq_data/quality_output/fnRs.qual.rds")

# Drop reverse reads from the rest of the analysis due to low quality

#### 3. Filter Ns and check for primer sequences ####

# Forward and reverse primer sequences
FWD <- "AGCCTCCGCTTATTGATATGCTTAART"
REV <- "AACTTTYRRCAAYGGATCWCT"

# Create function for all possible orientations of each primer
allOrients <- function(primer) {
  
  # Create all orientations of the input sequence
  require(Biostrings)
  # The Biostrings works w/ DNAString objects rather than character vectors
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  # Convert back to character vector
  return(sapply(orients, toString))
}

# Create all possible orientations of each primer
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# Pre-filter sequences to remove sequences with any ambiguous base calls; put N-filtered files in filtN/ subdirectory
fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
# Filter
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE)

# Tabulate primer hits
fnFs.filtN.primerHits <- vector(mode = "list", length = length(fnFs.filtN))
names(fnFs.filtN.primerHits) <- sample.names

#for(i in seq_along(fnFs.filtN.primerHits)) {
#  
#  print(paste("Processing:", sample.names[i], sep = " "))
#  fnFs.filtN.primerHits[[i]] <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[i]]), 
#                                      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[i]]))
#  
#}

# Save data
saveRDS(fnFs.filtN.primerHits, file = "data/Manistee_MiSeq_data/primer_hits_output/fnFs.primer.hits.rds")

#### 4. Remove primers using cutadapt ####

# Specify path to cutadapt
cutadapt <- "/Users/williamargiroff/.local/bin/cutadapt"
# Check that cutadapt is reached from R
system2(cutadapt, args = "--version")

# Create a directory to hold the cutadapt-trimmed sequences
path.cut <- file.path(path, "cutadapt_processed")
if(!dir.exists(path.cut)) dir.create(path.cut)

# Create names for the files in the directory
fnFs.cut <- file.path(path.cut, basename(fnFs))

# Reverse complement the reverse primer (as it would show up in the forward reads if the read sequenced into the reverse primer)
REV.RC <- dada2:::rc(REV)

# Create flags to trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 

# Run Cutadapt
#for(i in seq_along(fnFs.cut)) {
#  
  #-n 2 required to remove FWD and REV primers from reads
#  system2(cutadapt, args = c(R1.flags, "-n", 2,
                             # output files
#                             "-o", fnFs.cut[i],
                             # input files
#                             fnFs.filtN[i]))
#  
#}

# Tabulate primer hits
fnFs.cut.primerHits <- vector(mode = "list", length = length(fnFs.cut))
names(fnFs.cut.primerHits) <- sample.names

for(i in seq_along(fnFs.cut.primerHits)) {
  
  print(paste("Processing:", sample.names[i], sep = " "))
  fnFs.cut.primerHits[[i]] <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[i]]), 
                                    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[i]]))
  
}

# Save data
saveRDS(fnFs.cut.primerHits, file = "data/Manistee_MiSeq_data/primer_hits_output/fnFs.cut.primer.hits.rds")

# Forward read fastq filenames have the following format
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format from above
sample.names <- unname(sapply(cutFs, get.sample.name))

# Quality
cutFs.qual <- vector(mode = "list", length = length(cutFs))
names(cutFs.qual) <- sample.names

for(i in seq_along(cutFs)) {
  
  print(paste("Processing:", sample.names[i], sep = " "))
  cutFs.qual[[i]] <- plotQualityProfile(cutFs[i])
  
}

# Save plots
saveRDS(cutFs.qual, file = "data/Manistee_MiSeq_data/quality_output/cutFs.qual.rds")

#### 5. Filter reads based on errors and length; quantify error rate ####

# Create a directory for filtered reads
filtFs <- file.path(path.cut, "filtered", basename(cutFs))

# filter reads
filter.summary <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 2, 
                                truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

# Quantify error rate
errF <- learnErrors(filtFs, multithread = TRUE)

# Plot error rates
errF.plot <- plotErrors(errF, nominalQ = TRUE)

# Quality
filtFs.qual <- vector(mode = "list", length = length(filtFs))
names(filtFs.qual) <- sample.names

for(i in seq_along(filtFs)) {
  
  print(paste("Processing:", sample.names[i], sep = " "))
  filtFs.qual[[i]] <- plotQualityProfile(filtFs[i])
  
}

# Save plots
saveRDS(filtFs.qual, file = "data/Manistee_MiSeq_data/quality_output/filtFs.qual.rds")

#### 6. De-replicate seqs, run ASV inference, remove chimeras, ASV table ####

# De-replicate the sequences
derepFs <- derepFastq(filtFs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

# Run sample inference on dereplicated sequences
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)

# Create ASV table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

# Inspect distribution of sequence lengths
seqtab.table <- table(nchar(getSequences(seqtab)))
saveRDS(seqtab.table, file = "data/Manistee_MiSeq_data/seq_summary/seqtab.table.rds")

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
seqtab.nochim.table <- table(nchar(getSequences(seqtab.nochim)))
saveRDS(seqtab.nochim.table, file = "data/Manistee_MiSeq_data_v2/seq_summary/seqtab.nochim.table.rds")

#### 7. Summarize filtering output ####

# Write function to get uniques to start tracking progress through the pipeline
getN <- function(x) sum(getUniques(x))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- as.data.frame(cbind(filter.summary, sapply(dadaFs, getN), rowSums(seqtab.nochim)))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
# Percent remaining
track$PERC.REMAINING <- 100 * (track$nonchim / track$input)

# Track sequences remaining after different steps; save data
saveRDS(track, file = "data/Manistee_MiSeq_data/seq_summary/dada2.seq.summary.rds")

#### 8. Assign taxonomy using UNITE ####

# Assign taxonomy
unite.ref <- "data/Manistee_MiSeq_data/unite/sh_general_release_dynamic_s_04.02.2020.fasta"
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, 
                       multithread = TRUE, 
                       tryRC = TRUE, 
                       minBoot = 80, 
                       outputBootstraps = TRUE, 
                       verbose = TRUE)

saveRDS(taxa, file = "data/Manistee_MiSeq_data/unite/taxa.rds")

#### 9. Save data as a phyloseq object ####

# Construct a sample data frame
samples.out <- rownames(seqtab.nochim)
samdf <- data.frame(SAMPLE.ID = samples.out, stringsAsFactors = FALSE)
row.names(samdf) <- samples.out

# Get sample type
samdf$TYPE <- sapply(samdf$SAMPLE.ID, sample_type, simplify = TRUE)
# Get sample plot
samdf$PLOT <- mapply(sample_plot, samdf$SAMPLE.ID, samdf$TYPE, SIMPLIFY = TRUE)
# Get sample stand
samdf$STAND <- mapply(sample_stand, samdf$PLOT, samdf$TYPE, SIMPLIFY = TRUE)
# Update plot
samdf$PLOT <- mapply(sample_plot_update, samdf$PLOT, samdf$TYPE, SIMPLIFY = TRUE)
# Update stand
samdf$STAND <- mapply(sample_stand_update, samdf$STAND, samdf$TYPE, SIMPLIFY = TRUE)
# Add substrate
samdf$SUBSTRATE <- mapply(sample_substrate, samdf$SAMPLE.ID, samdf$TYPE, SIMPLIFY = TRUE)
# Add run
samdf$RUN.TYPE <- mapply(sample_seqRun, samdf$SAMPLE.ID, samdf$TYPE, SIMPLIFY = TRUE)

# Load phyloseq
library(phyloseq)

# Create phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
               sample_data(samdf), 
               tax_table(taxa$tax))

# Get refseqs
dna <- Biostrings::DNAStringSet(taxa_names(ps))
# Add seqs to phyloseq object
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
# Update ASV names
taxa_names(ps) <- paste0("ASV_", seq(ntaxa(ps)))

# Save data
saveRDS(ps, file = "data/Manistee_MiSeq_data/asv_data/phyloseq.rds")
