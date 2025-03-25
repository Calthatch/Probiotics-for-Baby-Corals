# Load required librairies
library(phyloseq)
library(ggplot2)
library(ape)
library(dplyr)
library(vegan)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(RVAideMemoire)

# Set working directory and check we are in the right place
setwd("/Users/katarinadamjanovic/Work AIMS/AmpliconAnalysis")
getwd()

# Get a set of colours to represent the different treatments:
getPalette <- colorRampPalette(brewer.pal(9, "Set1")) # See https://rdrr.io/cran/RColorBrewer/man/ColorBrewer.html
TreatmentPalette <- getPalette(13)

# The following code will represent PCoA (based on actual distances between samples),
# but we can also represent the data using NMDS plots and analyses, which are based
# on ranking

###################################################################################
########################## PCoA ###################################################
###################################################################################

phy.cl.r # Decontaminated and rarefied dataset: 3786 ASVs and 242 samples

# Create different variables to avoid putting all samples in a single plot:

# Atenuis recruits
Aten.cl.r <- subset_samples(phy.cl.r, Experiment == "Atenuis" & LifeStage == "Recruit")
Aten.cl.r <- prune_taxa((taxa_sums(Aten.cl.r)>0), Aten.cl.r)

# Pdaedalea recruits
Pdae.cl.r <- subset_samples(phy.cl.r, Experiment == "Pdaedalea" & LifeStage == "Recruit")
Pdae.cl.r <- prune_taxa((taxa_sums(Pdae.cl.r)>0), Pdae.cl.r)

# Both species, before inoculation recruits
Naive.cl.r <- subset_samples(phy.cl.r, Age %in% c("3 Days", "5 Days") & LifeStage == "Recruit")
Naive.cl.r <- prune_taxa((taxa_sums(Naive.cl.r)>0), Naive.cl.r)

# Both species, with C1 and C2
Controls.cl.r <- subset_samples(phy.cl.r, Treatment %in% c("C1_FSW", "C2_Ecoli") & LifeStage == "Recruit")
Controls.cl.r <- prune_taxa((taxa_sums(Controls.cl.r)>0), Controls.cl.r)

#### Bray-Curtis ####

# Note: Bray-curtis is based on relative abundances. Using Jaccard will provide results 
# based on presence/absence, while UniFrac is based on phylogeny

## Atenuis ## 

Aten.cl.r.pcoa <- ordinate(
  physeq = Aten.cl.r, 
  method = "PCoA", 
  distance = "bray")

plot_ordination(
  physeq = Aten.cl.r,
  ordination = Aten.cl.r.pcoa,
  color = "Treatment",
  shape = "Age",
  title = "PCoA of bacterial communities in A. tenuis") + 
  geom_point(colour = "black", size = 6) +
  geom_point(aes(fill = Treatment), alpha = 1, size = 5) +
  scale_colour_manual(values = TreatmentPalette) +
  theme_bw(base_size = 20)

## Pdaedalea ##

Pdae.cl.r.pcoa <- ordinate(
  physeq = Pdae.cl.r, 
  method = "PCoA", 
  distance = "bray")

plot_ordination(
  physeq = Pdae.cl.r,
  ordination = Pdae.cl.r.pcoa,
  color = "Treatment",
  shape = "Age",
  title = "PCoA of bacterial communities in P. daedalea") + 
  geom_point(colour = "black", size = 6) +
  geom_point(aes(fill = Treatment), alpha = 1, size = 5) +
  scale_colour_manual(values = TreatmentPalette) +
  theme_bw(base_size = 20)


# Colour palette:
my_colors1 <- RColorBrewer::brewer.pal(10, "BrBG")[2:4]
my_colors2 <- RColorBrewer::brewer.pal(10, "BrBG")[8:10]
my_colors <- c(my_colors1, my_colors2)

## Naive and controls ##

Naive.cl.r.pcoa <- ordinate(
  physeq = Naive.cl.r, 
  method = "PCoA", 
  distance = "bray")

pNaive <- plot_ordination(
  physeq = Naive.cl.r,
  ordination = Naive.cl.r.pcoa,
  color = "Species",
  shape = "Age",
  title = "A. tenuis and P. daedalea before inoculation") + 
  geom_point(colour = "black", size = 6) +
  geom_point(aes(fill = Species), alpha = 1, size = 5) +
  scale_colour_manual(values = c(my_colors[2], my_colors[5])) +
  theme_bw(base_size = 20)

Controls.cl.r.pcoa <- ordinate(
  physeq = Controls.cl.r, 
  method = "PCoA", 
  distance = "bray")

pControl <- plot_ordination(
  physeq = Controls.cl.r,
  ordination = Controls.cl.r.pcoa,
  color = "Sp_Age",
  shape = "Treatment",
  title = "Control treatments") + 
  geom_point(colour = "black", size = 6) +
  geom_point(aes(fill = Sp_Age), alpha = 1, size = 5) +
  scale_colour_manual(values = c(my_colors[1], my_colors[3], my_colors[4], my_colors[6])) +
  theme_bw(base_size = 20)



grid.arrange(pNaive, pControl, ncol = 1, nrow = 2 )


rm(Aten.cl.r, Aten.cl.r.pcoa, Pdae.cl.r, Pdae.cl.r.pcoa, TreatmentPalette, Naive.cl.r, Control.cl.r,
   Controls.cl.r, Controls.cl.r.pcoa, Naive.cl.r, Naive.cl.r.pcoa, pControl, pNaive, my_colors, my_colors1, my_colors2)





# Then comparisons between sample groups can be made using different models and pairwise comparisons:
# The following functions could be used, but has to be checked with a statistician...
# Also not sure if this has to be applied on rartefied or non rarefied data

# Calculate dissimilarity matrix
# DistBC = phyloseq::distance(Aten.cl.r, method = "bray") # Relative abundance
# dat = as(sample_data(Aten.cl.r), "data.frame")
# adonis2(DistBC ~ Treatment, strata = dat$Age, data = dat, permutations = 9999) # Treatment seems to play a role

# Check for homogeneity of variance
# beta <- betadisper(DistBC, sample_data(Aten.cl.r)$Treatment)
# disper.test = permutest(beta, permutations =9999)
# disper.test # Ok if p > 0.05, here it is significant -> variances are not homogenous (dispersion is different btwn treatments)

# We shouldn't then be allowed to do that:
# testingBC = pairwise.perm.manova(DistBC, sample_data(Aten.cl.r)$Treatment,
#                                       nperm=9999, p.method = "BH")
# testingBC$p.value















