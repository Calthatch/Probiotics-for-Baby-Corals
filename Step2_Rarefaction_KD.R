# Load required librairies
library(phyloseq)
library(ggplot2)
library(ape)
library(dplyr)
library(vegan)
library(RColorBrewer)


# Set working directory and check we are in the right place
setwd("/Users/katarinadamjanovic/Work AIMS/AmpliconAnalysis")
getwd()

# Colours for plots
coloursAten <- c("burlywood1", "burlywood2", "burlywood3", "burlywood4", "lightseagreen", "turquoise4")
coloursPdae <- c("wheat2", "wheat3", "lightgoldenrod3", "lightgoldenrod4", "lightseagreen", "turquoise4")

# Colours for rarefaction curves
my_colors1 <- RColorBrewer::brewer.pal(10, "BrBG")[2:4]
my_colors2 <- RColorBrewer::brewer.pal(10, "BrBG")[8:10]
my_colors <- c(my_colors1, my_colors2)


###################################################################################
########################## CHECK SEQUENCING DEPTH #################################
###################################################################################

# Stats about the "cleaned" phyloseq object

phy.clean # 4377 taxa in 243 samples 

# Statistics for the dataset:
summary(sample_sums(phy.clean)) # Min: 108 reads; Median: 39393 reads; Max: 63532 reads
sum(sample_sums(phy.clean)) # There is a total of 9'614'958 reads

# Make bar plots with the number of reads per sample
# Hard to see so many samples in one plot, so split the phyloseq object in two:
Aten.clean <- subset_samples(phy.clean, Experiment == "Atenuis")
Pdae.clean <- subset_samples(phy.clean, Experiment == "Pdaedalea")

Aten.clean.melt <- psmelt(Aten.clean)
Pdae.clean.melt <- psmelt(Pdae.clean)

Aten.clean.plot <- ggplot(Aten.clean.melt, aes(x = SampleName, y = Abundance, fill = Sp_Age )) +
  geom_bar(stat = "identity", colour = "black") + 
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) +
  scale_fill_manual(values = coloursAten) +
  labs(title = "Read counts for A. tenuis", x = "Sample ID", y = "Number of reads", fill = "Sample group")

Pdae.clean.plot <- ggplot(Pdae.clean.melt, aes(x = SampleName, y = Abundance, fill = Sp_Age )) +
  geom_bar(stat = "identity", colour = "black") + 
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) +
  scale_fill_manual(values = coloursPdae) +
  labs(title = "Read counts for P. daedalea", x = "Sample ID", y = "Number of reads", fill = "Sample group")

# Display the plots and save them as .PDF using the Export tool and selecting "A4, Landscape":
Aten.clean.plot
Pdae.clean.plot

# Let's have a look at the read counts for all samples:
sort(sample_sums(phy.clean))

# We will be losing one Pdaedalea sample (Pd_D14_T2_1), which failed to sequence properly as
# it has only 108 reads. The next samples with lowest counts are Bact_At_inoc_1 (18600) and
# Bact_Pd_inoc_1 (22844). If we accept to lose this one bacterial inoculation replicate sample
# for each experiment, we can push up the limit of reads to the next sample (At_D18_T3_1) and
# rarefy at 23160 reads. Otherwise, if we prefer keeping all samples and have the two inoculation
# replicates, we can rarefy at 18600 reads

set.seed(1)
rare_depth <- 18600
phy.cl.r = rarefy_even_depth(phy.clean, sample.size = rare_depth)
# 591 ASVs were removed because they are no longer present in any sample after random subsampling

phy.cl.r # 3786 taxa and 242 samples

rm(Aten.clean, Pdae.clean, Aten.clean.plot, Pdae.clean.plot, rare_depth, coloursAten, coloursPdae,
   Aten.clean.melt, Pdae.clean.melt)




###################################################################################
########################## PLOT RAREFACTION CURVES ################################
###################################################################################

# Load function for OPTION 1:
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/load-extra-functions.R")

# Load function for OPTION 2: https://github.com/joey711/phyloseq/issues/143
set.seed(1)
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}


###### Subsetting data before applying rarefaction plotting functions

# Let's draw rarefaction curves to see where they plateau and decide on a sequencing depth.
# We will separate the datasets according to species so that the curves are less packed:
# For each one, elimitae ASVs which are completely absent from the dataset

Aten.clean <- subset_samples(phy.clean, Experiment == "Atenuis")
Aten.clean <- prune_taxa((taxa_sums(Aten.clean)>0), Aten.clean)

Pdae.clean <- subset_samples(phy.clean, Experiment == "Pdaedalea")
Pdae.clean <- prune_taxa((taxa_sums(Pdae.clean)>0), Pdae.clean)


# OPTION 1

Aten.r.plot <- ggrare(Aten.clean, step = 10, plot = TRUE, color = "Sp_Age", label = NULL, se = FALSE)

Aten.r.plot +
  scale_colour_manual(values = my_colors) +
  labs(title = "Rarefaction curves for A. tenuis", x = "Read counts", y = "Richness", color = "Sample group") +
  theme_bw(base_size = 16)

Pdae.r.plot <- ggrare(Pdae.clean, step = 10, plot = TRUE, color = "Sp_Age", label = NULL, se = FALSE)

Pdae.r.plot +
  scale_colour_manual(values = my_colors) +
  labs(title = "Rarefaction curves for P. daedalea", x = "Read counts", y = "Richness", color = "Sample group") +
  theme_bw(base_size = 16)
  



# OPTION 2

# Rarefaction curves for Atenuis:
rarefaction_curve_data.Aten <- calculate_rarefaction_curves(Aten.clean, c('Observed', 'Shannon'),
                                                          rep(c(1:100 * 10, 1:100 * 1000), each = 10))
summary(rarefaction_curve_data.Aten)

# Sumarise alpha diversity:
rarefaction_curve_data_summary.Aten<- ddply(rarefaction_curve_data.Aten, c('Depth', 'Sample', 'Measure'),
                                           summarise, Alpha_diversity_mean = mean(Alpha_diversity), 
                                           Alpha_diversity_sd = sd(Alpha_diversity))
# Add sample data:
rarefaction_curve_data_summary.Aten_verbose <- merge(rarefaction_curve_data_summary.Aten,
                                                   data.frame(sample_data(Aten.clean)),
                                                   by.x = 'Sample', by.y = 'row.names')

# Plot:(we see a huge diversity from the water samples)
ggplot(
  data = rarefaction_curve_data_summary.Aten_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = Sp_Age,
    group = Sample)) +
  geom_line() +
  facet_wrap(facets = ~ Measure, scales = 'free_y') +
  scale_colour_manual(values = my_colors) +
  ggtitle("Rarefaction curves for A. tenuis")


# Rarefaction curves for Pdaedalea:
rarefaction_curve_data.Pdae <- calculate_rarefaction_curves(Pdae.clean, c('Observed', 'Shannon'),
                                                            rep(c(1:100 * 10, 1:100 * 1000), each = 10))
summary(rarefaction_curve_data.Pdae)

# Sumarise alpha diversity:
rarefaction_curve_data_summary.Pdae<- ddply(rarefaction_curve_data.Pdae, c('Depth', 'Sample', 'Measure'),
                                            summarise, Alpha_diversity_mean = mean(Alpha_diversity), 
                                            Alpha_diversity_sd = sd(Alpha_diversity))
# Add sample data:
rarefaction_curve_data_summary.Pdae_verbose <- merge(rarefaction_curve_data_summary.Pdae,
                                                     data.frame(sample_data(Pdae.clean)),
                                                     by.x = 'Sample', by.y = 'row.names')

# Plot:(we see a huge diversity from the water samples)
ggplot(
  data = rarefaction_curve_data_summary.Pdae_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = Sp_Age,
    group = Sample)) +
  geom_line() +
  facet_wrap(facets = ~ Measure, scales = 'free_y') +
  scale_colour_manual(values = my_colors) +
  ggtitle("Rarefaction curves for P. daedalea")


rm(my_colors, my_colors1, my_colors2, Aten.clean, Pdae.clean, rarefaction_curve_data_summary.Pdae_verbose, rarefaction_curve_data.Pdae,
   rarefaction_curve_data_summary.Pdae, rarefaction_curve_data_summary.Aten_verbose, rarefaction_curve_data_summary.Aten, rarefaction_curve_data.Aten,
   coloursAten, coloursPdae, calculate_rarefaction_curves, scripts, url, urls)

rm(list=lsf.str()) # removes functions
