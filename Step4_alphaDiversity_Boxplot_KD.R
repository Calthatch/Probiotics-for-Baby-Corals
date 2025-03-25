# Load required librairies
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(ape)
library(dplyr)
library(vegan)
library(RColorBrewer)

# Set working directory and check we are in the right place
setwd("/Users/katarinadamjanovic/Work AIMS/AmpliconAnalysis")
getwd()


###################################################################################
########################## DIVERSITY METRICS ######################################
###################################################################################

phy.cl.r # Decontaminated and rarefied dataset: 3786 ASVs and 242 samples

# Create a table with alpha diversity indices
rich <- estimate_richness(phy.cl.r)

# Attach the sample metadata to this table of diversity indices
# First add a common column (sampleID) at the start of both matrices, then attach the metadata
sampleID <- rownames(rich)
rich <- cbind(sampleID, rich) %>% data.frame()
sampleID <- rownames(sample_data(phy.cl.r))
sample_data(phy.cl.r) <- cbind(sampleID, sample_data(phy.cl.r)) %>% data.frame()
rich_fac <- rich %>% left_join(data.frame(sample_data(phy.cl.r)))

head(rich_fac) # Check what the variable looks like


###################################################################################
########################## BOX PLOTS ##############################################
###################################################################################


# We can split the dataset into the two species:

rich_Aten <- subset(rich_fac, Experiment == "Atenuis")
rich_Pdae <- subset(rich_fac, Experiment == "Pdaedalea")

# Reorder levels to plot in a different order:

Sp_Age_Tr2 <- c("At_inoc_bact", "At_inoc_zoox", "At_D03", "At_D05_Zx", "At_D14_C1", "At_D18_C1", "At_D14_C2",
                "At_D18_C2", "At_D14_T1", "At_D18_T1", "At_D14_T2", "At_D18_T2", "At_D14_T3", "At_D18_T3","At_D14_T4", 
               "At_D18_T4", "At_D14_T5", "At_D18_T5",  "At_D14_T6", "At_D18_T6", "At_D14_T7", "At_D18_T7",
               "At_D14_T8", "At_D18_T8", "At_D14_T9",  "At_D18_T9")

rich_Aten$Sp_Age_Tr <- factor(rich_Aten$Sp_Age_Tr, 
                                      levels(rich_Aten$Sp_Age_Tr) <- Sp_Age_Tr2)


# Atenuis
plot_Aten_observed <- ggplot(rich_Aten, aes(x = Sp_Age_Tr , y = Observed, fill = Treatment)) +
                        geom_boxplot(position =position_dodge(width  =.8)) +
                        geom_point(colour = "black") +
                        theme_bw(base_size = 16) +
                        theme(legend.position="right", axis.title.x = element_blank(),
                        axis.text.x = element_text(angle = 90, hjust = 1)) 

plot_Aten_shannon <- ggplot(rich_Aten, aes(x = Sp_Age_Tr , y = Shannon, fill = Treatment)) +
                        geom_boxplot(position =position_dodge(width  =.8)) +
                        geom_point(colour = "black") +
                        theme_bw(base_size = 16) +
                        theme(legend.position="right", axis.title.x = element_blank(),
                              axis.text.x = element_text(angle = 90, hjust = 1)) 

boxplot_Aten <- ggarrange(plot_Aten_observed, plot_Aten_shannon,
                    common.legend = TRUE, legend = "right",
                    ncol = 1, nrow = 2)


boxplot_Aten
annotate_figure(boxplot_Aten, top = text_grob("Diversity indices in A. tenuis", 
                                      color = "black", face = "bold", size = 16), fig.lab.pos = "top.left")


# Pdaedalea
plot_Pdae_observed <- ggplot(rich_Pdae, aes(x = Sp_Age_Tr , y = Observed, fill = Treatment)) +
  geom_boxplot(position =position_dodge(width  =.8)) +
  geom_point(colour = "black") +
  theme_bw(base_size = 16) +
  theme(legend.position="right", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) 

plot_Pdae_shannon <- ggplot(rich_Pdae, aes(x = Sp_Age_Tr , y = Shannon, fill = Treatment)) +
  geom_boxplot(position =position_dodge(width  =.8)) +
  geom_point(colour = "black") +
  theme_bw(base_size = 16) +
  theme(legend.position="right", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) 

boxplot_Pdae <- ggarrange(plot_Pdae_observed, plot_Pdae_shannon,
                          common.legend = TRUE, legend = "right",
                          ncol = 1, nrow = 2)


boxplot_Pdae
annotate_figure(boxplot_Pdae, top = text_grob("Diversity indices in P. daedalea", 
                                              color = "black", face = "bold", size = 16))

rm(sampleID, rich_Aten, rich_Pdae, plot_Aten_observed, plot_Aten_shannon, boxplot_Aten,
   plot_Pdae_observed, plot_Pdae_shannon, boxplot_Pdae, rich, rich_fac)


