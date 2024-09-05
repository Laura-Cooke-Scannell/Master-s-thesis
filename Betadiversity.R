# clear all
rm(list = ls())

# set working directory to fasta-qual-mapping-files-NOTRIM on my desktop
setwd("/home/laura/Desktop/Thesis/samplesredo")

# load packages 
library(vegan)
library(biomformat)
library(ecodist)

# working directory set to fasta-qual-mapping-files-NOTRIM on desktop, loading
#the biom files of my data 
ITS_data_bact <- read_biom("./exported-table-bact/feature-table.biom")
ITS_data_euk <- read_biom("./exported-table-euk/feature-table.biom")

#Turn the files into matrices for Vegan
ITS_data1 <- as.matrix(ITS_data_bact)
ITS_data2 <- as.matrix(ITS_data_euk)

#Turn labels into one dataframe and data into another for the bacterial data
ITS_data1_labels <- ITS_data1[1:8,1]
ITS_data1_data <- ITS_data1[1:8,2:40]

#Turn labels into one dataframe and data into another for the eukaryotes data
ITS_data2_labels <- ITS_data_euk[1:5,1]
ITS_data2_data <-ITS_data_euk[1:5,2:12]

#Calculate bray curtis dissimilarity index (standardized)
ITS1_diss_sta <- vegdist(ITS_data1_data, method = 'bray', binary = TRUE, diag = FALSE,
                         upper = FALSE)
ITS2_diss_sta <- vegdist(ITS_data2_data, method = 'bray', binary = TRUE, diag = FALSE,
                         upper = FALSE)


#Calculate bray curtis dissimalarity index (not standardized)
ITS1_diss <- vegdist(ITS_data1_data, method = 'bray', binary = FALSE, diag = FALSE,
                     upper = FALSE)
ITS2_diss <- vegdist(ITS_data2_data, method = 'bray', binary = FALSE, diag = FALSE,
                     upper = FALSE)

#Calculate Chi Squared (standardized)
ITS1_Chi_sta <- vegdist(ITS_data1_data, method = 'chisq', binary = TRUE, diag = FALSE,
                        upper = FALSE)
ITS2_Chi_sta <- vegdist(ITS_data2_data, method = 'chisq', binary = TRUE, diag = FALSE,
                        upper = FALSE)

#Calculate Chi Squared (not standardized)
ITS1_Chi <- vegdist(ITS_data1_data, method = 'chisq', binary = FALSE, diag = FALSE,
                    upper = FALSE)
ITS2_Chi <- vegdist(ITS_data2_data, method = 'chisq', binary = FALSE, diag = FALSE,
                    upper = FALSE)

#creating a PCoA plot  for the bray-curtis bacteria (standardized)
pcoa_ITS1_dis_sta <- pco(ITS1_diss_sta, negvals = "zero", dround = 0)
plot(pcoa_ITS1_dis_sta$vectors[,1], pcoa_ITS1_dis_sta$vectors[,2], type = "n", xlab = 
       "PCoA1, 0.629", ylab = "PCoA2, .5", axes = TRUE, main = "PCoA (ecodist on 16s data)")

text(pcoa_ITS1_dis_sta$vectors[,1], pcoa_ITS1_dis_sta$vectors[,2], labels= (ITS_data1_labels),
     cex = 0.9, xpd = TRUE)

#creating a PCoA plot  for the bray-curtis eukarotes (standardized)
pcoa_ITS2_dis_sta <- pco(ITS2_diss_sta, negvals = "zero", dround = 0)
plot(pcoa_ITS2_dis_sta$vectors[,1], pcoa_ITS2_dis_sta$vectors[,2], type = "n", xlab = 
       "PCoA1, 0.629", ylab = "PCoA2, .5", axes = TRUE, main = "PCoA (ecodist on 16s data)")

text(pcoa_ITS2_dis_sta$vectors[,1], pcoa_ITS2_dis_sta$vectors[,2], labels= (ITS_data2_labels),
     cex = 0.9, xpd = TRUE)

#creating a PCoA plot for the bray-curtis of the bacteria (not standardized)
pcoa_ITS1_dis <- pco(ITS1_diss, negvals = "zero", dround = 0)
plot(pcoa_ITS1_dis$vectors[,1], pcoa_ITS1_dis$vectors[,2], type = "n", xlab = 
       "PCoA1, 0.629", ylab = "PCoA2, .5", axes = TRUE, main = "PCoA (ecodist on 16s data)")

text(pcoa_ITS1_dis$vectors[,1], pcoa_ITS1_dis$vectors[,2], labels= (ITS_data1_labels),
     cex = 0.9, xpd = TRUE)

#creating a PCoA plot for the bray-curtis of the eukaryotes (not standardized)
pcoa_ITS2_dis <- pco(ITS2_diss, negvals = "zero", dround = 0)
plot(pcoa_ITS2_dis$vectors[,1], pcoa_ITS2_dis$vectors[,2], type = "n", xlab = 
       "PCoA1, 0.629", ylab = "PCoA2, .5", axes = TRUE, main = "PCoA (ecodist on 16s data)")

text(pcoa_ITS2_dis$vectors[,1], pcoa_ITS2_dis$vectors[,2], labels= (ITS_data2_labels),
     cex = 0.9, xpd = TRUE)

#creating a PCoA plot for the Chi squared of the bacteria (standardized)
pcoa_ITS1_chi_sta <- pco(ITS1_Chi_sta, negvals = "zero", dround = 0)
plot(pcoa_ITS1_chi_sta$vectors[,1], pcoa_ITS1_chi_sta$vectors[,2], type = "n", xlab = 
       "PCoA1, 0.629", ylab = "PCoA2, .5", axes = TRUE, main = "PCoA (ecodist on 16s data)")

text(pcoa_ITS1_chi_sta$vectors[,1], pcoa_ITS1_chi_sta$vectors[,2], labels= (ITS_data1_labels),
     cex = 0.9, xpd = TRUE)

#creating a PCoA plot for the Chi squared of the eukaryotes (standardized)
pcoa_ITS2_chi_sta <- pco(ITS2_Chi_sta, negvals = "zero", dround = 0)
plot(pcoa_ITS2_chi_sta$vectors[,1], pcoa_ITS2_chi_sta$vectors[,2], type = "n", xlab = 
       "PCoA1, 0.629", ylab = "PCoA2, .5", axes = TRUE, main = "PCoA (ecodist on 16s data)")

text(pcoa_ITS2_chi_sta$vectors[,1], pcoa_ITS2_chi_sta$vectors[,2], labels= (ITS_data2_labels),
     cex = 0.9, xpd = TRUE)

#Creating a PCoA plot for the Chi squared of the bacteria (not standardized)
pcoa_ITS1_chi <- pco(ITS1_Chi, negvals = "zero", dround = 0)
plot(pcoa_ITS1_chi$vectors[,1], pcoa_ITS1_chi$vectors[,2], type = "n", xlab = 
       "PCoA1, 0.629", ylab = "PCoA2, .5", axes = TRUE, main = "PCoA (ecodist on 16s data)")

text(pcoa_ITS1_chi$vectors[,1], pcoa_ITS1_chi$vectors[,2], labels= (ITS_data1_labels),
     cex = 0.9, xpd = TRUE)

#Creating a PCoA plot for the Chi squared of the eukaryotes (not standardized)
pcoa_ITS2_chi <- pco(ITS2_Chi, negvals = "zero", dround = 0)
plot(pcoa_ITS2_chi$vectors[,1], pcoa_ITS2_chi$vectors[,2], type = "n", xlab = 
       "PCoA1, 0.629", ylab = "PCoA2, .5", axes = TRUE, main = "PCoA (ecodist on 16s data)")

text(pcoa_ITS2_chi$vectors[,1], pcoa_ITS2_chi$vectors[,2], labels= (ITS_data2_labels),
     cex = 0.9, xpd = TRUE)