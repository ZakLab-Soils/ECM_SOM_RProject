#### Functions for SOM feedback analysis ####

#### Standard error function ####
se <- function(x) {
  # Calculate se
  tmp1 <- sd(x) * (1 / sqrt(length(x)))
  # Return the result
  return(tmp1)
}

#### Degrees to radians ####
deg_rad <- function(x) {
  # Calculate radians
  tmp1 <- x * (pi / 180)
  # Return result
  return(tmp1)
}

#### Better rounding function than R's base round (Marian Schmidt) ####
matround <- function(x) {
  tmp1 <- trunc(x + 0.5)
  return(tmp1)
}

#### Scale 0 to 1 ####
scale_01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#### Create and order two-element vectors ####
sort_vec <- function(x, y) {
  tmp1 <- x
  tmp2 <- y
  
  tmp3 <- c(x, y)
  tmp4 <- sort(tmp3)
  tmp5 <- toString(tmp4)
  tmp6 <- gsub("\\, ", ".", tmp5)
  
  return(tmp6)
}

#### Get P value from overall regression model (lm) ####
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Write a function to tabulate primer hits
primerHits <- function(primer, fn) {
  
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
  
}

#### Detach a package ####

# https://stackoverflow.com/questions/6979917/how-to-unload-a-package-without-restarting-r

detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

#### Write a function to get sample type ####
sample_type <- function(x) {
  
  if(x == "MOCK1" | x == "MOCK2" | x == "MOCK1nano" | x == "MOCK2nano") {
    tmp1 <- "mock"
  } else if(x == "NC1" | x == "NC2" | x == "NC1nano" | x == "NC2nano") {
    tmp1 <- "negative.control"
  } else {
    tmp1 <- "environmental"
  }
  
  return(tmp1)
  
}

#### Write a function to get plot number ####
sample_plot <- function(sample.id, sample.type) {
  
  if(sample.type == "environmental") {
    tmp1 <- gsub("[a-zA-Z]", "", sample.id)
  } else {
    tmp1 <- sample.id
  }
  
  return(tmp1)
  
}

#### Write a function to get stand number ####
sample_stand <- function(plot.number, sample.type) {
  
  if(sample.type == "environmental" & as.numeric(plot.number) >= 1 & as.numeric(plot.number) <= 6) {
    tmp1 <- 3
  } else if(sample.type == "environmental" & as.numeric(plot.number) >= 7 & as.numeric(plot.number) <= 12) {
    tmp1 <- 6
  } else if(sample.type == "environmental" & as.numeric(plot.number) >= 13 & as.numeric(plot.number) <= 18) {
    tmp1 <- 7
  } else if(sample.type == "environmental" & as.numeric(plot.number) >= 19 & as.numeric(plot.number) <= 24) {
    tmp1 <- 9
  } else if(sample.type == "environmental" & as.numeric(plot.number) >= 25 & as.numeric(plot.number) <= 30) {
    tmp1 <- 20
  } else if(sample.type == "environmental" & as.numeric(plot.number) >= 31 & as.numeric(plot.number) <= 36) {
    tmp1 <- 22
  } else if(sample.type == "environmental" & as.numeric(plot.number) >= 37 & as.numeric(plot.number) <= 42) {
    tmp1 <- 24
  } else if(sample.type == "environmental" & as.numeric(plot.number) >= 43 & as.numeric(plot.number) <= 48) {
    tmp1 <- 31
  } else if(sample.type == "environmental" & as.numeric(plot.number) >= 49 & as.numeric(plot.number) <= 54) {
    tmp1 <- 41
  } else if(sample.type == "environmental" & as.numeric(plot.number) >= 55 & as.numeric(plot.number) <= 60) {
    tmp1 <- 50
  } else if(sample.type == "environmental" & as.numeric(plot.number) >= 61 & as.numeric(plot.number) <= 66) {
    tmp1 <- 58
  } else if(sample.type == "environmental" & as.numeric(plot.number) >= 67 & as.numeric(plot.number) <= 72) {
    tmp1 <- 100
  } else if(sample.type == "mock") {
    tmp1 <- plot.number
  } else if(sample.type == "negative.control") {
    tmp1 <- plot.number
  } else {
    tmp1 <- NA
  }
  
  return(tmp1)
}

#### Write function to update plot ####
sample_plot_update <- function(plot.number, sample.type) {
  
  if(sample.type == "environmental") {
    tmp1 <- paste("Plot", plot.number, sep = "_")
  } else {
    tmp1 <- plot.number
  }
  
  return(tmp1)
  
}

#### Write function to update stand ####
sample_stand_update <- function(stand.number, sample.type) {
  
  if(sample.type == "environmental") {
    tmp1 <- paste("Stand", stand.number, sep = "_")
  } else {
    tmp1 <- stand.number
  }
  
  return(tmp1)
  
}

#### Write function to get substrate ####
sample_substrate <- function(sample.id, sample.type) {
  
  if(sample.type == "environmental") {
    tmp1 <- gsub("[[:digit:]]", "", sample.id)
    tmp2 <- ifelse(tmp1 == "S" | tmp1 == "Snano", "soil", "decomp.roots")
  } else {
    tmp2 <- sample.id
  }
  
  return(tmp2)
  
}

#### Write a function to get sequence run ####
sample_seqRun <- function(sample.id, sample.type) {
  
  if(sample.type == "environmental") {
    tmp1 <- gsub("[[:digit:]]", "", sample.id)
    tmp2 <- ifelse(tmp1 == "S" | tmp1 == "R" | tmp1 == "RA" | tmp1 == "RB", "full", "nano")
  } else if(sample.id == "MOCK1" | sample.id == "MOCK2" | sample.id == "NC1" | sample.id == "NC2") {
    tmp2 <- "full"
  } else {
    tmp2 <- "nano"
  }
  
  return(tmp2)
  
}

#### Write function to clean up taxa ####
tax_clean <- function(x) {
  
  if(!is.na(x)) {
    tmp1 <- gsub(".*_", "", x)
  } else {
    tmp1 <- x
  }
  
  return(tmp1)
  
}

#### Write a function to handle unclassified species ####
s_unclass <- function(kingdom, phylum, class, order, family, genus, species) {
  
  if(!is.na(species)) {
    tmp1 <- species
  } else if(is.na(species) & !is.na(genus) & genus != "sedis") {
    tmp1 <- paste("unclassified", genus, sep = " ")
  } else if(is.na(species) & is.na(genus) & !is.na(family) & family != "sedis") {
    tmp1 <- paste("unclassified", family, sep = " ")
  } else if(is.na(species) & is.na(genus) & is.na(family) & !is.na(order) & order != "sedis") {
    tmp1 <- paste("unclassified", order, sep = " ")
  } else if(is.na(species) & is.na(genus) & is.na(family) & is.na(order) & !is.na(class) & class != "sedis") {
    tmp1 <- paste("unclassified", class, sep = " ")
  } else if(is.na(species) & is.na(genus) & is.na(family) & is.na(order) & is.na(class) & !is.na(phylum) & phylum != "sedis") {
    tmp1 <- paste("unclassified", phylum, sep = " ")
  } else if(is.na(species) & is.na(genus) & is.na(family) & is.na(order) & is.na(class) & is.na(phylum) & !is.na(kingdom) & kingdom != "sedis") {
    tmp1 <- paste("unclassified", kingdom, sep = " ")
  } else {
    tmp1 <- NA
  }

  return(tmp1)
  
}

#### Write a function to handle unclassified genera ####
g_unclass <- function(kingdom, phylum, class, order, family, genus) {
  
  if(!is.na(genus)) {
    tmp1 <- genus
  } else if(is.na(genus) & !is.na(family) & family != "sedis") {
    tmp1 <- paste("unclassified", family, sep = " ")
  } else if(is.na(genus) & is.na(family) & !is.na(order) & order != "sedis") {
    tmp1 <- paste("unclassified", order, sep = " ")
  } else if(is.na(genus) & is.na(family) & is.na(order) & !is.na(class) & class != "sedis") {
    tmp1 <- paste("unclassified", class, sep = " ")
  } else if(is.na(genus) & is.na(family) & is.na(order) & is.na(class) & !is.na(phylum) & phylum != "sedis") {
    tmp1 <- paste("unclassified", phylum, sep = " ")
  } else if(is.na(genus) & is.na(family) & is.na(order) & is.na(class) & is.na(phylum) & !is.na(kingdom) & kingdom != "sedis") {
    tmp1 <- paste("unclassified", kingdom, sep = " ")
  } else {
    tmp1 <- NA
  }
  
  return(tmp1)
  
}

#### Write a function to handle unclassified families ####
f_unclass <- function(kingdom, phylum, class, order, family) {
  
  if(!is.na(family)) {
    tmp1 <- family
  } else if(is.na(family) & !is.na(order) & order != "sedis") {
    tmp1 <- paste("unclassified", order, sep = " ")
  } else if(is.na(family) & is.na(order) & !is.na(class) & class != "sedis") {
    tmp1 <- paste("unclassified", class, sep = " ")
  } else if(is.na(family) & is.na(order) & is.na(class) & !is.na(phylum) & phylum != "sedis") {
    tmp1 <- paste("unclassified", phylum, sep = " ")
  } else if(is.na(family) & is.na(order) & is.na(class) & is.na(phylum) & !is.na(kingdom) & kingdom != "sedis") {
    tmp1 <- paste("unclassified", kingdom, sep = " ")
  } else {
    tmp1 <- NA
  }
  
  return(tmp1)
  
}

#### Write a function to handle unclassified orders ####
o_unclass <- function(kingdom, phylum, class, order) {
  
  if(!is.na(order)) {
    tmp1 <- order
  } else if(is.na(order) & !is.na(class) & class != "sedis") {
    tmp1 <- paste("unclassified", class, sep = " ")
  } else if(is.na(order) & is.na(class) & !is.na(phylum) & phylum != "sedis") {
    tmp1 <- paste("unclassified", phylum, sep = " ")
  } else if(is.na(order) & is.na(class) & is.na(phylum) & !is.na(kingdom) & kingdom != "sedis") {
    tmp1 <- paste("unclassified", kingdom, sep = " ")
  } else {
    tmp1 <- NA
  }
  
  return(tmp1)
  
}

#### Write a function to handle unclassified classes ####
c_unclass <- function(kingdom, phylum, class) {
  
  if(!is.na(class)) {
    tmp1 <- class
  } else if(is.na(class) & !is.na(phylum) & phylum != "sedis") {
    tmp1 <- paste("unclassified", phylum, sep = " ")
  } else if(is.na(class) & is.na(phylum) & !is.na(kingdom) & kingdom != "sedis") {
    tmp1 <- paste("unclassified", kingdom, sep = " ")
  } else {
    tmp1 <- NA
  }
  
  return(tmp1)
  
}

#### Write a function to handle unclassified phyla ####
p_unclass <- function(kingdom, phylum) {
  
  if(!is.na(phylum)) {
    tmp1 <- phylum
  } else if(is.na(phylum) & !is.na(kingdom) & kingdom != "sedis") {
    tmp1 <- paste("unclassified", kingdom, sep = " ")
  } else {
    tmp1 <- NA
  }
  
  return(tmp1)
  
}

#### Get raw sequence counts ####
raw_seq_counts <- function(x) {
  
  tmp1 <- sum(x$data$Count) / 251
  
  return(tmp1)
  
}
