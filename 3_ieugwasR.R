######################################################################################################################
######## PROJECT: MR analysis of the gut microbiome on pancreatic cancer using PanScan + PanC4 data (IARC) and FinnGen
######## Script: Look up the SNPs associated with the gut microbial traits that are being followed up
######## Date: 05/06/23
######################################################################################################################
## Clear space
rm(list=ls())

## Set working directory
setwd("/path/to/directory/prostate_cancer/")

## INSTALL PACKAGES BEFOREHAND
#remotes::install_github("mrcieu/ieugwasr")
library(ieugwasr)
#install.packages("readxl")   # run once
library(readxl)
library(dplyr)


#######################################################################################
## READING IN THE EXPOSURE DATA AND RESTRICTING TO THE FOUR TRAITS BEING FOLLOWED UP ##
#######################################################################################
## Hughes
Hughes_exposure_dat <- read_exposure_data("./Data/Hughes_ovarian_exposure_data_main.txt", sep="\t", phenotype_col = "mt", snp_col = "rsid", beta_col = "beta", se_col = "se", pval_col = "P_value", eaf_col = "eaf", effect_allele_col = "allele_B", other_allele_col = "allele_A", samplesize_col = "n")
head(Hughes_exposure_dat)

## MiBioGen
MiBioGen_exposure_dat <- read_exposure_data("./Data/MiBioGen_ovarian_exposure_data_main.txt", sep="\t", phenotype_col = "exposure", snp_col = "SNP", beta_col = "beta.exposure", se_col = "se.exposure", pval_col = "pval.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure", samplesize_col = "samplesize.exposure")
head(MiBioGen_exposure_dat)

## Bind these together
all_exposures <- rbind(Hughes_exposure_dat, MiBioGen_exposure_dat)
dim(all_exposures)

## Restrict to only the four exposures upon which sensitivity analyses will focus
ofinterest <- as.vector(c("G_unclassified_P_Firmicutes_HB"))
SNPsofinterest <- as.vector(unique(all_exposures$SNP[all_exposures$exposure %in% ofinterest]))



############################################################################
## USE IEUGWASR FUNCTION TO LOOK UP THE TRAITS ASSOCIATED WITH THESE SNPS ##
############################################################################
## Query using the "phewas" function

phewas <- as.data.frame(phewas(SNPsofinterest, pval = 0.0001))
colnames(phewas)
head(phewas)

phewas_formatted <- phewas[,c(8,5,3,9,10,11,12,7,2,1,6,4)]
head(phewas_formatted)
write.table(phewas, "./IEUGWAS_phewas_results.txt", col.names = T, row.names = F, quote = F, sep = "\t")

