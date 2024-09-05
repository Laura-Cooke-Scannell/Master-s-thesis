# clear all
rm(list = ls())

# set working directory to bio609 on your desktop
setwd("/Users/Laura/Desktop/Thesis/Another redo of samples/Joost nopoop")
#libraries
library(devtools)
library(phyloseq)
library(qiime2R)
library(ggplot2)
library(vegan)

#create the phyloseq object
physeq <- qza_to_phyloseq(features = "samplesclean.qza", tree = "rooted-tree.qza",
                          taxonomy = "taxonomy.qza", metadata = "Joostmetadata.tsv")

# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = phyloseq::distance(physeq, method ="unifrac", weighted=F)
ordination = ordinate(physeq, method="PCoA", distance = wunifrac_dist)
plot_ordination(physeq, ordination, color="description", title="Sulfur Cave UniFrac results P = 0.33") + theme(aspect.ratio = 1)

#permanova significance test of unweighted unifrac
adonis2(wunifrac_dist ~ sample_data(physeq)$description)