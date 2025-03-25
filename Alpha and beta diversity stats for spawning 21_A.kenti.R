######## Final Script for modeling alpha and beta diversity stats on the spawning 21 amplicon data - A. kenti ONLY #######

#NOTES for github: At the start of this project the taxonomic assignment of the coral was Acropora tenuis. It has since 
# been changes to Acropora kenti, therefore At and Ak are interchangeable. 

library(phyloseq)
library(vegan)
library(dplyr)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(pairwiseAdonis)

# Set working directory and check we are in the right place
setwd("C:/Users/calla/Desktop/Amplicon analysis spawning 2021/Cal R files")

#Upload your probiotics_workspace which is an RDATA FIle made by Kat in this case, it should have the
#clean and rareified phy objects

# Create a table with alpha diversity indices
rich <- estimate_richness(phy.cl.r)

# Attach the sample metadata to this table of diversity indices
# First add a common column (sampleID) at the start of both matrices, then attach the metadata
sampleID <- rownames(rich)
rich <- cbind(sampleID, rich) %>% data.frame()
sampleID <- rownames(sample_data(phy.cl.r))
sample_data(phy.cl.r) <- cbind(sampleID, sample_data(phy.cl.r)) %>% data.frame()
rich_fac <- rich %>% left_join(data.frame(sample_data(phy.cl.r)))

# We can split the dataset into the two species:

rich_Aten <- subset(rich_fac, Experiment == "Atenuis")
rich_Pdae <- subset(rich_fac, Experiment == "Pdaedalea")

#Change variables to factors 

rich_Aten$Treatment <- as.factor(rich_Aten$Treatment)
rich_Aten$Age <- as.factor(rich_Aten$Age)

#Subset the rich_Aten data to eliminate factors that will confuse our model

rich_Aten_subset <- subset(rich_Aten, Treatment!="Zoox" & Treatment!="Culture" & Treatment!="Nothing" & Treatment!="T2_Rugeria")

# Build your model, we can use the rareified data in alpha diversity metrics 

model_Shannon_Aten <- glmmTMB(Shannon ~ Age*Treatment, data = rich_Aten_subset)

# Now we need to check the assumptions of normality and homogenity
# The following code does both at once

S.resid <- model_Shannon_Aten %>% simulateResiduals(plot=TRUE, integerResponse = TRUE)

#Want there to be no significance, in our case our P value is 0.6 so it's fine 
#The model works so now we look for significance in the placebo control compared to each treatment

model_Shannon_Aten %>% summary ()

#Run a post hoc

model_Shannon_Aten %>%
  emmeans(~Treatment | Age) %>%
  pairs(adjust = "BH")

# To see P values and adjusted P values

# Unadjusted
unadj <- model_Shannon_Aten %>%
  emmeans(~Treatment | Age) %>%
  pairs(adjust = "none") %>%
  as.data.frame() %>%
  rename(p_unadj = p.value)

# Adjusted (e.g., BH)
adj <- model_Shannon_Aten %>%
  emmeans(~Treatment | Age) %>%
  pairs(adjust = "BH") %>%
  as.data.frame() %>%
  select(contrast, Age, p.value) %>%
  rename(p_adj = p.value)

compare_p <- left_join(unadj, adj, by = c("contrast", "Age"))
compare_p

# We also want to see if treatments change over time i.e. Roseivivax Day 14 to 18

#Remove T2 data as there is none for day 18

rich_Aten_subset <- subset(rich_Aten, Treatment!="Zoox" & Treatment!="Culture" & Treatment!="Nothing" & Treatment!="T2_Rugeria")

#Build your model
model_Shannon_Aten_Age <- glmmTMB(Shannon ~ Treatment*Age, data = rich_Aten_subset)

# Check the residuals
S.resid <- model_Shannon_Aten_Age %>% simulateResiduals(plot=TRUE, integerResponse = TRUE)

#Run the model

model_Shannon_Aten_Age %>% 
  emmeans (~Age|Treatment) %>%
  pairs(adjust='BH')

# Unadjusted
unadj.age <- model_Shannon_Aten_Age %>%
  emmeans(~Age | Treatment) %>%
  pairs(adjust = "none") %>%
  as.data.frame() %>%
  rename(p_unadj = p.value)


###################################################################
############## Beta Diversity Metrics #############################
###################################################################

# Subset data by species (done above) then day 
# Building the model using adonis2 functions

# Split the data set into different species

Aten.clean <- subset_samples(phy.clean, Experiment == "Atenuis")
Pdae.clean <- subset_samples(phy.clean, Experiment == "Pdaedalea")

####### A. kenti Day 14 Beta diversity stats #####

Aten.clean.subset.D14 <- subset_samples(Aten.clean, Age %in% c('14 Days') & Treatment %in% c('C1_FSW', 'C2_Ecoli', 'T1_Pseudovibrio', 'T2_Rugeria', 'T3_Neptunomonas', 'T4_Halomonas', 'T5_Endozoicomonas', 'T6_Roseivivax', 'T7_Pseudoalteromonas', 'T8_Thallasobius', 'T9_Vcoralliilyticus')) 

DistBC_Aten_subset.D14 = phyloseq::distance(Aten.clean.subset.D14, method = "bray") 

dat_Aten_subset.D14 = as(sample_data(Aten.clean.subset.D14), "data.frame")

adonis2(DistBC_Aten_subset.D14 ~ Treatment, data = dat_Aten_subset.D14, permutations = 9999)

pairwise.adonis2(DistBC_Aten_subset.D14 ~ Treatment, data = dat_Aten_subset.D14, nperm= 9999)

#Check variances using beta disper

beta_Aten.subset.D14 <- betadisper(DistBC_Aten_subset.D14, sample_data(Aten.clean.subset.D14)$Treatment)

disper.test.Aten.subset.D14 = permutest(beta_Aten.subset.D14, permutations =9999)

disper.test.Aten.subset.D14

#Correcting for multiple comparisons as requested by Kat

# Use "BH" the Benjamin-Hochberg correction as it is commonly usesd in beta analyses

# Run pairwise Adonis2 test
pairwise_results_AkD13 <- pairwise.adonis2(DistBC_Aten_subset.D14 ~ Treatment, data = dat_Aten_subset.D14, nperm= 9999)

# Extract p-values from the pairwise comparisons
p_values <- pairwise_results_AkD13$p.value  # Ensure this corresponds to the column containing p-values

# Apply multiple testing correction (Benjamini-Hochberg FDR)
pairwise_results_AkD13$p.adjusted <- p.adjust(p_values, method = "BH")

# Print the updated results with corrected p-values
print(pairwise_results_AkD13)


# Safely extract p-values from pairwise_results_AkD13
p_values <- sapply(pairwise_results_AkD13, function(x) {
  if (is.data.frame(x) && "Pr(>F)" %in% colnames(x)) {  # Ensure x is a data frame and has the p-value column
    return(x$`Pr(>F)`[1])  # Extract the first row's p-value
  } else {
    return(NA)  # Return NA if not applicable
  }
})

# Remove NA values (if any)
p_values <- p_values[!is.na(p_values)]

# Print extracted p-values
print(p_values)


# Apply Benjamini-Hochberg FDR correction
p_adjusted <- p.adjust(p_values, method = "BH")

# Print adjusted p-values
print(p_adjusted)

##### A. kenti Day 18 Beta Stats #####

#Repeat the above steps for Atenuis Day 18 however excluding T2 as it doesn't exist

Aten.clean.subset.D18 <- subset_samples(Aten.clean, Age %in% c('18 Days') & Treatment %in% c('C1_FSW', 'C2_Ecoli', 'T1_Pseudovibrio', 'T3_Neptunomonas', 'T4_Halomonas', 'T5_Endozoicomonas', 'T6_Roseivivax', 'T7_Pseudoalteromonas', 'T8_Thallasobius', 'T9_Vcoralliilyticus')) 

DistBC_Aten_subset.D18 = phyloseq::distance(Aten.clean.subset.D18, method = "bray") 

dat_Aten_subset.D18 = as(sample_data(Aten.clean.subset.D18), "data.frame")

adonis2(DistBC_Aten_subset.D18 ~ Treatment, data = dat_Aten_subset.D18, permutations = 9999)

pairwise.adonis2(DistBC_Aten_subset.D18 ~ Treatment, data = dat_Aten_subset.D18, nperm= 9999)

#Check variances on this new subset of data

beta_Aten.subset.D18 <- betadisper(DistBC_Aten_subset.D18, sample_data(Aten.clean.subset.D18)$Treatment)

disper.test.Aten.subset.D18 = permutest(beta_Aten.subset.D18, permutations =9999)

disper.test.Aten.subset.D18

#

pairwise_results_AKD17 <- pairwise.adonis2(DistBC_Aten_subset.D18 ~ Treatment, data = dat_Aten_subset.D18, nperm= 9999)


#### Check for multiple testing using "BH" correction ####

p_values_AKD17 <- sapply(pairwise_results_AKD17, function(x) {
  if (is.data.frame(x) && "Pr(>F)" %in% colnames(x)) {  # Ensure x is a data frame and has the p-value column
    return(x$`Pr(>F)`[1])  # Extract the first row's p-value
  } else {
    return(NA)  # Return NA if not applicable
  }
})

# Remove NA values (if any)
p_values_AkD17 <- p_values_AKD17[!is.na(p_values_AKD17)]

# Print extracted p-values
print(p_values_AKD17)

# Apply Benjamini-Hochberg FDR correction
p_adjusted_AKD17 <- p.adjust(p_values_AKD17, method = "BH")

# Print adjusted p-values
print(p_adjusted_AKD17)

