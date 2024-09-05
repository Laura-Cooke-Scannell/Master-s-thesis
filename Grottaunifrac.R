# clear all
rm(list = ls())

# set working directory to another redo of samples on my desktop
setwd("/Users/Laura/Desktop/another redo of samples/nopoop")

#libraries
library(devtools)
library(phyloseq)
library(qiime2R)
library(ggplot2)
library(vegan)

#create the phyloseq object
physeq <- qza_to_phyloseq(features = "nopooptable.qza", tree = "rooted-tree.qza",
                          taxonomy = "taxonomynopoop.qza", metadata = "metadata.csv")

# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = phyloseq::distance(physeq, method ="unifrac", weighted=F)
ordination = ordinate(physeq, method="PCoA", distance = wunifrac_dist)
plot_ordination(physeq, ordination, color="Location", title = "Grotta Del Cane bacterial UniFrac results P = 0.563") + theme(aspect.ratio = 1)

#permanova significance test of unweighted unifrac
adonis2(wunifrac_dist ~ sample_data(physeq)$Location)

#PCoA plot using the weighted UniFrac as distance
uunifrac_dist = phyloseq::distance(physeq, method = "unifrac", weighted=T)
ordination2 = ordinate(physeq, method = "PCoA", distance = uunifrac_dist)
plot_ordination(physeq, ordination2, color="Location") + theme(aspect.ratio = 1)

#permanova significance test of weighted unifrac
adonis2(uunifrac_dist ~ sample_data(physeq)$Location)


