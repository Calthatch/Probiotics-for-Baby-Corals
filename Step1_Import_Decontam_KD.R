# install required packages with the command:

#install.packages("Package name")

# If that doesn't work, do the following:

# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("version = "3.15""Package name)


# Load required librairies
library(phyloseq)
library(decontam)
library(ggplot2)
library(ape)
library(dplyr)
library(vegan)

# Command to remove all variables in R:
rm(list=ls())

# Set working directory and check we are in the right place
setwd("/Users/katarinadamjanovic/Work AIMS/AmpliconAnalysis")
getwd()
list.files(path = ".")


###################################################################################
########################## IMPORTING and CHECKING DATA ############################
###################################################################################

# Import ASV Table
ASV <- read.table(
  "ImportR/feature-table.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

# Import taxonomic table
TAXA <- read.table(
  "ImportR/taxonomy-edited.txt",
  header = TRUE,
  sep = "\t",
  # fill = TRUE, needed if we had left blank cells in the table
  row.names = 1
)

# Import metadata file
META <- read.table(
  "ImportR/ProbioticsMetadata.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

# Convert ASV and taxonomic tables to matrices for phyloseq
ASV <- as.matrix(ASV)
TAXA <- as.matrix(TAXA)


# Combine ASV, taxa and metadata files into a phyloseq object
phy <- phyloseq(
  otu_table(ASV,taxa_are_rows = T),
  tax_table(TAXA),
  sample_data(META)
)



# Explore phyloseq object to make sure everything is coherent:
phy # There are 4714 ASVs and 265 samples

head(sample_data(phy)) # Displays the first few lines of the metadata file
str(get_variable(phy)) # Summarizes the type of variables in the metadata
head(otu_table(phy)) # Displays the first few lines of the ASV table
head(tax_table(phy)) # Displays the first few lines of the taxonomic table
rank_names(phy) # Taxonomic ranks are called Domain, Phylum, Class, Family, Genus.

# Testing if there remains any Mitochondria or Chloroplasts in the dataset (we removed them in qiime2)
any(tax_table(phy)[, "Family"] == "Endozoicomonadaceae") # TRUE
any(tax_table(phy)[, "Family"] == "Mitochondria") # FALSE
any(tax_table(phy)[, "Order"] == "Chloroplast") # FALSE

# If needed, here is the code to remove particular taxa once already in R:
#phy <- phy %>%
#  subset_taxa(
#    Domain == "Bacteria" &
#      Family  != "Mitochondria" &
#      Order   != "Chloroplast"
#  )


# See number of reads (in ascending order) for each sample:
sort(sample_sums(phy))

# Statistics for the dataset:
summary(sample_sums(phy)) # Min: 33 reads; Median: 39326 reads; Max: 64032 reads
sum(sample_sums(phy)) # There is a total of 10'214'118 reads



###################################################################################
########################## DECONTAMINATION ########################################
###################################################################################

# The code to decontaminate the dataset using package "decontam" is written below
# Before running it, let's first identify the main ASVs present in the bacterial inoculation
# samples. If some sequences cross-contaminated other samples, including blanks, we need
# to make sure they are not labelled as contaminants and entirely removed from the dataset 

#### Check ASV identities of our inoculation solutions ####

inocAten <- subset_samples(phy, Sp_Age_Tr == "At_inoc_bact")
inocPdae <- subset_samples(phy, Sp_Age_Tr == "Pd_inoc_bact")

# Remove zero sum ASVs
inocAten <- prune_taxa((taxa_sums(inocAten) > 0), inocAten)
inocPdae <- prune_taxa((taxa_sums(inocPdae) > 0), inocPdae)

# Interrogate the objects 
inocAten # 38 ASVs, 2 samples
inocPdae # 39 ASVs, 2 samples

# Write the ASV table and taxonomic information into two .csv files, which we can combine
# to have the taxonomy and read counts in one place. The files are combined in a new
# file called "incoAten_Info.xlsx" and "incoPdae_Info.xlsx" and merged correctly in Excel
# by sorting both tables (the otu and taxonomy ones) by alphabetic order based on the ASV identifier:

write.csv(x = otu_table(inocAten), row.names = TRUE, file = "Contaminants/inocAtenCounts.csv")
write.csv(x = tax_table(inocAten), row.names = TRUE, file = "Contaminants/inocAtenTaxas.csv")

write.csv(x = otu_table(inocPdae), row.names = TRUE, file = "Contaminants/inocPdaeCounts.csv")
write.csv(x = tax_table(inocPdae), row.names = TRUE, file = "Contaminants/inocPdaeTaxas.csv")

# Reminder, the following taxa were used as inocula for A.tenuis and P. daedalea:

             ###### A.tenuis #######        ###### P.daedalea #######
                
                  # E. coli                    # E. coli 
                  # Pseudovibrio               # Neptunomonas
                  # Ruegeria                   # Shimia
                  # Neptunomonas               # Ruegeria
                  # Halomonas                  # Pseudovibrio
                  # Endozoicomonas             # Pseudoalteromonas 
                  # Roseivivax                 # Amphitrea
                  # Pseudoalteromonas          # Bacterioplanes
                  # Thallasobius               # Colwellia
                  # V. coralliilyticus         # V. coralliilyticus 


# We observe that somne of these taxa can appear several times in the inoculation data - these "different"
# ASVs could represent the same bacteria but come from distinct copies of their 16S rRNA genes. This could
# be misleading as we can treat the ASVs either as different bacteria, or as the same bacteria with internal
# sequence variations within their 16S rRNA genes... For example, there are several Escherichia ASVs, or several
# "Colwelliaceae" in the P. daedalea dataset

# We also observe that the taxonomy in our table can be incomplete (this is due to the SILVA database
# being incomplete, as it was based mostly on human/earth research and not on marine organisms). So
# for example we might only have information up to the family level

# To be certain that we don't remove essential ASVs during the decontamination process, we can create
# a list of the most abundant bacteria present in our inoculation samples, and ensure they are kept
# in the dataset. Note: the ASVs that have zero counts in one of the incoulation replicate samples
# are probably spurious and we won't inlcude them in our inoculation list

# 1) Removing from our inoculation list the ASVs that have zero count in at least one replicate 

# Convert the asv table of these objects into a dataframe
Aten.asv.table <- data.frame(otu_table(inocAten))
Pdae.asv.table <- data.frame(otu_table(inocPdae))

# Check what it looks like
head(Aten.asv.table)
dim(Aten.asv.table)

# Extract all rows for which there is a zero in any column:
Aten.zeros <- Aten.asv.table %>% filter_all(any_vars(. %in% c(0)))
Pdae.zeros <- Pdae.asv.table %>% filter_all(any_vars(. %in% c(0)))

# Get a list of these ASVs that you want to eliminate
Aten.bad.ASVs <- rownames(Aten.zeros)
Pdae.bad.ASVs <- rownames(Pdae.zeros)

# This is a list of all ASVs in the original phyloseq object:
Aten.allTaxa = taxa_names(inocAten)
Pdae.allTaxa = taxa_names(inocPdae)

# List of ASVs that you want to keep (the opposite of what you want to elimninate:
Aten.goodTaxa <- Aten.allTaxa[!(Aten.allTaxa %in% Aten.bad.ASVs)]
Pdae.goodTaxa <- Pdae.allTaxa[!(Pdae.allTaxa %in% Pdae.bad.ASVs)]

# Eliminate the bad ASVs from original phyloseq object (keep the good ones)
inocAten2 <- prune_taxa(Aten.goodTaxa, inocAten)
inocPdae2 <- prune_taxa(Pdae.goodTaxa, inocPdae)

# Check what these inocualtion samples look like now:
otu_table(inocAten2)
otu_table(inocPdae2)

# Create list of "inoculation ASVs", selecting only the 20 most abundant ones:
listAten <- names(sort(taxa_sums(inocAten2), TRUE)[1:20])
listPdae <- names(sort(taxa_sums(inocPdae2), TRUE)[1:20])
listInoc <- c(listAten, listPdae) # Joins the two lists
listInoc2 <- unique(listInoc) # removes duplicates


########################### Code using the decontam package ######################
##################################################################################

# Identify contaminants in two steps, first from PCR negatives, then extraction blanks.

#### PCR negatives ####
consList <- isContaminant(seqtab = phy, neg = "PCR_control", method = "prevalence")

# Pull out the names of contaminants
cons <- rownames(consList)[consList$contaminant=="TRUE"]
cons <- as.character(cons) #3 contaminants identified from PCR negatives

# To get info on the contaminants (taxas and relative abundance), run
# the following code ON THE FILE WITH THE CONTAMINANT ASVs IN-SITU,
# then combine the consPCR.csv and taxonomy.csv file data
# This code was provided by Ashley Dungan

# - - - - - - - - - - - - - - - - - - - - - - - - - #
# PCR controls

# subset the non neg-control samples
vvv1 <- subset_samples(phy, PCR_control == "FALSE")
# merge the samples into one so that all ASV abundances are summed
yyy1 <- merge_samples(vvv1, "PCR_control", fun = sum)
# transform counts to percentages
yyy1 <- transform_sample_counts(yyy1, function(x) 100 * x/sum(x))
# extract the contaminants percentage data
zzz1 <- prune_taxa(x = yyy1, taxa = cons)
# write otu table to dataframe
xxx1 <- data.frame(t(zzz1@otu_table))
# write xxx to csv
write.csv(x = xxx1, row.names = TRUE, file = "ContaminantsR/consPCR_Control.csv")
# subset the contaminant ASVs
PhyCons1 <- prune_taxa(phy, taxa = cons)
# write the contaminants to a file for reference
contaminants1 <- PhyCons1@tax_table
contaxa1 <- contaminants1@.Data
write.csv(contaxa1, "ContaminantsR/contaxaPCR_Control.csv")

# In total, the 3 PCR controls in the dataset comprise ~2% of the reads in actual samples
# The most abundant contaminant is a Halomonas (2.17%). However we note that this particular Halomonas
# ASV (08a3e391fb401071943e7585a9679dd7) belongs to the Atenuis inoculation group.
# This might have been a sample cross contamination rather than a lab contamination.
# So we will make sure not to remove it from the dataset


#### Blank extractions ####
consList2 <- isContaminant(seqtab = phy, neg = "Extract_control", method = "prevalence")
# Not sure whether it is better to use with batch...
# consList2 <- isContaminant(seqtab = phy, neg = "Extract_control", method = "prevalence", batch = "ExtractionBatch") # Not sure

# Pull out the names of contaminants
cons2 <- rownames(consList2)[consList2$contaminant=="TRUE"]
cons2 <- as.character(cons2) # 117 contaminants identified from blank extractions

# To get info on the contaminants (taxas and rel. ab.), run
# the following code ON THE FILE WITH THE CONTAMINANT ASVs IN-SITU,
# then combine the consPCR.csv and taxonomy.csv file data
# Again, this code was provided by Ashley Dungan

# - - - - - - - - - - - - - - - - - - - - - - - - - #
# Blank extractions

# subset the non extraction-control samples
vvv2 <- subset_samples(phy, Extract_control == "FALSE")
# merge the samples
yyy2 <- merge_samples(vvv2, "Extract_control", fun = sum)
# transform counts to percentages
yyy2 <- transform_sample_counts(yyy2, function(x) 100 * x/sum(x))
# extract the cons percentage data
zzz2 <- prune_taxa(x = yyy2, taxa = cons2)
# write otu table to dataframe
xxx2 <- data.frame(t(zzz2@otu_table))
# write xxx to csv
write.csv(x = xxx2, row.names = TRUE, file = "ContaminantsR/consExtract_Control.csv")
# subset the contaminant ASVs
PhyCons2 <- prune_taxa(phy, taxa = cons2)
# write the contaminants to a file for reference
contaminants2 <- PhyCons2@tax_table
contaxa2 <- contaminants2@.Data
write.csv(contaxa2, "ContaminantsR/contaxaExtract_Control.csv")

# In total, the 117 blank extraction controls in thedataset comprise 3.9% of the reads


############# Removing negative control samples and contaminant taxas ############
##################################################################################

# List of all contaminants
consAll <- c(cons, cons2)

# Check which elements from our inoculation list are in the contaminant list:
falseContam <- table(consAll[match(listInoc2, consAll)]) # there are 2 ASVs 
falseContam

# 08a3e391fb401071943e7585a9679dd7: Halomonas (in Atenuis inoc)
# 876826f3480f487da8c7c723faf6c719: Ruegeria (in Pdaedalea inoc)

# Remove those 2 ASVs from contamination list:
consAll2 <- consAll[consAll %in% names(falseContam) == FALSE]

# Remove contaminant ASVs from phyloseq object
allTaxa = taxa_names(phy)
goodTaxa <- allTaxa[!(allTaxa %in% consAll2)]
phy.clean0 <- prune_taxa(goodTaxa, phy)

# Check that numbers match (we should have 4714 - 118 = 4596 taxas)
phy.clean0

# Remove all negative controls
phy.clean1 <- subset_samples(phy.clean0, PCR_control == "FALSE")
phy.clean2 <- subset_samples(phy.clean1, Extract_control == "FALSE")

# Remove zero-sum ASVs
phy.clean <- prune_taxa((taxa_sums(phy.clean2) > 0), phy.clean2)

# Final object:
phy.clean # 4377 taxa in 243 samples 

# Some ASVs were not identified as contaminants but occurred only in negative controls,
# and not in any other sample. Hence they were removed from the dataset when we removed
# the negative control samples. 

# Remove variables we won't use anymore
rm(ASV, META, TAXA, consList, consList2, contaxa1, allTaxa,listAten, listPdae,
   contaxa2, cons, cons2, consAll,consAll2, PhyCons1, PhyCons2, inocAten, inocPdae,
   phy.clean0, phy.clean1, phy.clean2, xxx1, xxx2, yyy1, yyy2, zzz1, zzz2,
   vvv1, vvv2, contaminants1, contaminants2, falseContam, goodTaxa, listInoc, listInoc2,
   Aten.allTaxa, Aten.asv.table, Aten.bad.ASVs, Aten.goodTaxa, Pdae.allTaxa, Pdae.asv.table,
   Pdae.bad.ASVs, Pdae.goodTaxa, Pdae.zeros, Aten.zeros, inocAten2, inocPdae2)


##################################################################################
##################################################################################


