#clear console
rm(list = ls())

#Check if Rtools is working
Sys.which("make")

# set working directory to bio609 on your desktop
setwd("C:\\Users\\Laura\\Desktop\\Thesis\\Joost\\Nopoop")

# Install bioconductor packages
BiocManager::install("ANCOMBC")
BiocManager::install("ComplexHeatmap")
devtools::install_github("jbisanz/qiime2R")

#libraries
library(devtools)
library(phyloseq)
library(qiime2R)
library(vegan)
library(BiocManager)
library(ComplexHeatmap)
library(ANCOMBC)

#create the phyloseq object
physeq <- qza_to_phyloseq(features = "table-nopoop.qza", tree = "rooted-tree.qza",
                          taxonomy = "taxonomy.qza", metadata = "Joostmetadata.tsv")

#Factorize data
sample_data(physeq)$description <- as.factor(sample_data(physeq)$description)
physeq.taxa <- tax_glom(physeq, taxrank = 'Species', NArm = FALSE)

#pairwise comparison between interface and above
physeq.taxa.sub1 <- subset_samples(physeq.taxa, description %in% c("Interface", "Above.Interface"))

#pairwise comparison between interface and below
physeq.taxa.sub2 <- subset_samples(physeq.taxa, description %in% c("Interface", "Above.Interface"))

#ancombc 1
out1 <- ancombc(phyloseq = physeq.taxa.sub1, formula = "descripption",
                p_adj_method = "holm", zero_cut = 0, group = "description",
                struc_zero = TRUE, conserve = TRUE)
res1 <-out1$res1

#create realtive form of taxa comparisons 1
physeq.taxa.rel <- transform_sample_counts(physeq, function(x) x/sum(x)*100)

# Select the bottom 20 with lowest p values 1
res1.or_p <- rownames(res1$q_val[,"desciptionInterface"])[base:order(res$q_val[,"descriptionInterface"])]
taxa_sig1 <- res1.or_p[1:20]
physeq.taxa.rel.sig1 <- prune_taxa(taxa_sig1, physeq.taxa.rel)

#Only keep interface and samples from above
ps.taxa.rel.sig1 <- prune_samples(colnames(otu_table(physeq.taxa.sub1)))

#Create heatmap
matrix1 <- as.matrix(data.frame(otu_table(physeq.taxa.rel.sig1)))
rownames(matrix1) <- as.character(tax_table(physeq.taxa.rel.sig1)[, "Species"])
metadata_sub1 <- data.frame(sample_data(physeq.taxa.rel.sig1))

#color columns and rows
annotation_col1 = data.frame(
  Location = as.factor(metadata_sub1$description),
  check.names = FALSE
)

rownames(annotation_col1) = rownames(metadata_sub1)

annotation_row1 = data.frame(
  Phylum = as.factor(tax_table(physeq.taxa.rel.sig1)[, "Phylum"])
)

rownames(annotation_row1) = rownames(matrix1)

# ann_color should be named vectors
phylum_col1 = RColorBrewer::brewer.pal(length(levels(annotation_row1$Phylum)), "Paired")
names(phylum_col1) = levels(annotation_row1$Phylum)
ann_colors1 = list(
  Location = c('Interface' = "blue", 'Above.Interface' = "red"),
  Phylum = phylum_col1
)

ComplexHeatmap1::pheatmap(matrix1, scale = "row",
                          annotation_col = annotation_col1,
                          annotation_row = annotation_row1,
                          annotation_colors = ann_colors1)