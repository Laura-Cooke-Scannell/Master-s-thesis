# clear all
rm(list = ls())

# set working directory to nopoop on my desktop
setwd("/Users/Laura/Desktop/Another redo of samples/Joost nopoop")

#load libraries
library(phyloseq)
library(qiime2R)
library(vegan)
library(stats)
library(dendextend)

#create phyloseq object
physeq <- qza_to_phyloseq(features = "nopoop.qza", tree = "rooted-tree.qza", 
                          taxonomy = "taxonomy.qza", metadata = "Joostmetadata.tsv")

#compute relative abundance
ps_rel_abund = phyloseq::transform_sample_counts(physeq, function(x){x /sum(x)})
phyloseq::otu_table(physeq)

#Extract OTU table and compute bray-curtis dissimilarity 
ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abund))
ps_rel_otu <- t(ps_rel_otu)
bc_dist <- vegan::vegdist(ps_rel_otu, method ="bray")
as.matrix(bc_dist)

#save as dendrogram
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))

#color codes
meta <- data.frame(phyloseq::sample_data(ps_rel_abund))
colorCode <- c(Interface = "hotpink3", Above.Interface= "goldenrod3", Below.Interface = "darkslategray4",
               Negative = "red", Positive = "blue")
labels_colors(ward) <- colorCode[meta$description][order.dendrogram(ward)]

#plot
plot(ward, main = "Sulfur Cave Bray-Curtis clustering results P = 0.07", 
     ylab = "Bray-Curtis distance", xlab = "SampleID")

#probability
adonis2(bc_dist ~ sample_data(physeq)$description)
