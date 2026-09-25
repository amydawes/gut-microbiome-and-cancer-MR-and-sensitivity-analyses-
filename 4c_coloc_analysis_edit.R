##############################################################################################################################
## PROJECT: MR analysis of the gut microbiome on prostate cancer using PRACTICAL, UK Biobank and FinnGen
## Script: Undertaking colocalisation analysis 
## Date: 25/09/26
##############################################################################################################################
## Clear space
rm(list=ls())

## Set working directory
setwd("/path/to/directory/prostate_cancer/")

## Load in the appropriate packages
#install.packages("dplyr")
#install.packages("openxlsx")
#install.packages("readxl")
#install.packages("TwoSampleMR")
#install.packages("coloc")
#install.packages("data.table")
library(dplyr)
library(openxlsx)
library(readxl)
library(TwoSampleMR)
library(coloc)
library(data.table)
#install.packages("MendelianRandomization")
library(MendelianRandomization)
library(dplyr)
library(ieugwasr)
#################################################### HUGHES ####################################################
### microbial traits: G_Parabacteroides_RNT, G_unclassified_F_Porphyromonadaceae_RNT, G_unclassified_P_Firmicutes_RNT

gwasinfo("ebi-a-GCST006085")
## SET ARGUMENTS 
N1 <- 2223
N2 <- 79148+61106 
type1 <- "quant"
type2 <- "cc"
exposure_name <- "Gut microbiome"
outcome_name <- "prostate cancer"

## Read in exposure data and format for colocalisation
exposure <- read.table("./Data/G_Parabacteroides_RNT_rs13207588_plusminus_1e6.txt", header = T, sep = "\t", fill = T)
exposure$eaf <- ifelse((exposure$all_BB < exposure$all_AA), exposure$all_maf, (1-exposure$all_maf))
head(exposure)
exposure_coloc <- format_data(exposure, type = "exposure", snps = NULL, header = T, snp_col = "rsid", beta_col = "frequentist_add_beta_1", se_col = "frequentist_add_se_1", eaf_col = "eaf", effect_allele_col = "alleleB", other_allele_col = "alleleA", pval_col = "frequentist_add_pvalue", chr_col = "chromosome", pos_col = "position")
exposure_coloc$exposure <- "Gut microbiome"
head(exposure_coloc)
exposure_coloc_snps <- as.vector(exposure_coloc$SNP)

## Read in outcome data and make sure all missing values are removed (big file so going to make sure that the SNPs that we need first are selected) then format for colocalisation analyses
prostate_data_for_coloc <- extract_outcome_data(
  snps = exposure_coloc_snps,
  outcomes = 'ebi-a-GCST006085') # 

prostate_data_for_coloc <- unique(prostate_data_for_coloc)
dim(prostate_data_for_coloc)
head(prostate_data_for_coloc)
prostate_data_for_coloc$rsid <- prostate_data_for_coloc$SNP

outcome_coloc <- format_data(prostate_data_for_coloc, type = "outcome", snps = NULL, header = T, snp_col = "rsid", beta_col = "beta.outcome", se_col = "se.outcome", eaf_col = "eaf.outcome", effect_allele_col = "effect_allele.outcome", other_allele_col = "other_allele.outcome", pval_col = "pval.outcome", chr_col = "chr", pos_col = "pos")# Duplicated SNPs removed in this step
outcome_coloc$outcome <- "Prostate cancer"
head(outcome_coloc)
dim(outcome_coloc)
## Remove any NAs
outcome_coloc_complete <- outcome_coloc[outcome_coloc$mr_keep.outcome == TRUE,] ## of those 5899, there are 5899 complete SNPs
dim(outcome_coloc_complete)
## Harmonise data (exposure and outcome data needs to be harmonised and then EAF converted to MAF)
df <- harmonise_data(exposure_dat = exposure_coloc, outcome_dat = outcome_coloc_complete, action = 2)
df <- df[!is.na(df$eaf.outcome), ]
## Get the MAF
for (j in 1:nrow(df)){
  if(df$eaf.exposure[j]>0.5){df$MAF1[j]=1-df$eaf.exposure[j]} else {df$MAF1[j]=df$eaf.exposure[j]}  
  if(df$eaf.outcome[j]>0.5){df$MAF2[j]=1-df$eaf.outcome[j]} else {df$MAF2[j]=df$eaf.outcome[j]}  
}
head(df)

## Perform analyses (note: that code assumes you have provided standard error and then squares this to get the variance)
coloc <- coloc.abf(dataset1=list(snp=df$SNP, beta=df$beta.exposure, 
                                 varbeta=df$se.exposure^2, MAF=df$MAF1, N=N1, type=type1),
                   dataset2=list(snp=df$SNP, beta=df$beta.outcome, 
                                 varbeta=df$se.outcome^2, MAF=df$MAF2, N=N2, type=type2),
                   p1=1e-04, p2=1e-04, p12=1e-05)

## Make table of results 
table <- data.frame(Exposure = exposure_name, Outcome = outcome_name,
                    NSNPs = coloc["summary"][[1]][1],
                    H0 = coloc["summary"][[1]][2],
                    H1 = coloc["summary"][[1]][3],
                    H2 = coloc["summary"][[1]][4],
                    H3 = coloc["summary"][[1]][5],
                    H4 = coloc["summary"][[1]][6])
write.table(table, "./Output/G_Parabacteroides_RNT_coloc_results.txt", col.names = T, row.names = F, sep = "\t", quote = F)

########## G_unclassified_F_Porphyromonadaceae_RNT


N1 <- 1420
N2 <- 79148+61106 
type1 <- "quant"
type2 <- "cc"
exposure_name <- "Gut microbiome"
outcome_name <- "prostate cancer"

## Read in exposure data and format for colocalisation
exposure <- read.table("./Data/G_unclassified_F_Porphyromonadaceae_RNT_rs35980751_plusminus_1e6.txt", header = T, sep = "\t", fill = T)
exposure$eaf <- ifelse((exposure$all_BB < exposure$all_AA), exposure$all_maf, (1-exposure$all_maf))
head(exposure)
exposure_coloc <- format_data(exposure, type = "exposure", snps = NULL, header = T, snp_col = "rsid", beta_col = "frequentist_add_beta_1", se_col = "frequentist_add_se_1", eaf_col = "eaf", effect_allele_col = "alleleB", other_allele_col = "alleleA", pval_col = "frequentist_add_pvalue", chr_col = "chromosome", pos_col = "position")
exposure_coloc$exposure <- "Gut microbiome"
head(exposure_coloc)
exposure_coloc_snps <- as.vector(exposure_coloc$SNP)

## Read in outcome data and make sure all missing values are removed (big file so going to make sure that the SNPs that we need first are selected) then format for colocalisation analyses
prostate_data_for_coloc <- extract_outcome_data(
  snps = exposure_coloc_snps,
  outcomes = 'ebi-a-GCST006085') # 

prostate_data_for_coloc <- unique(prostate_data_for_coloc)
dim(prostate_data_for_coloc)
head(prostate_data_for_coloc)
prostate_data_for_coloc$rsid <- prostate_data_for_coloc$SNP

outcome_coloc <- format_data(prostate_data_for_coloc, type = "outcome", snps = NULL, header = T, snp_col = "rsid", beta_col = "beta.outcome", se_col = "se.outcome", eaf_col = "eaf.outcome", effect_allele_col = "effect_allele.outcome", other_allele_col = "other_allele.outcome", pval_col = "pval.outcome", chr_col = "chr", pos_col = "pos")# Duplicated SNPs removed in this step
outcome_coloc$outcome <- "Prostate cancer"
head(outcome_coloc)
dim(outcome_coloc)
## Remove any NAs
outcome_coloc_complete <- outcome_coloc[outcome_coloc$mr_keep.outcome == TRUE,] ## 
dim(outcome_coloc_complete) # 6815 SNPs
## Harmonise data (exposure and outcome data needs to be harmonised and then EAF converted to MAF)
df <- harmonise_data(exposure_dat = exposure_coloc, outcome_dat = outcome_coloc_complete, action = 2)
df <- df[!is.na(df$eaf.outcome), ]
## Get the MAF
for (j in 1:nrow(df)){
  if(df$eaf.exposure[j]>0.5){df$MAF1[j]=1-df$eaf.exposure[j]} else {df$MAF1[j]=df$eaf.exposure[j]}  
  if(df$eaf.outcome[j]>0.5){df$MAF2[j]=1-df$eaf.outcome[j]} else {df$MAF2[j]=df$eaf.outcome[j]}  
}
head(df)

## Perform analyses (note: that code assumes you have provided standard error and then squares this to get the variance)
coloc <- coloc.abf(dataset1=list(snp=df$SNP, beta=df$beta.exposure, 
                                 varbeta=df$se.exposure^2, MAF=df$MAF1, N=N1, type=type1),
                   dataset2=list(snp=df$SNP, beta=df$beta.outcome, 
                                 varbeta=df$se.outcome^2, MAF=df$MAF2, N=N2, type=type2),
                   p1=1e-04, p2=1e-04, p12=1e-05)

## Make table of results 
table <- data.frame(Exposure = exposure_name, Outcome = outcome_name,
                    NSNPs = coloc["summary"][[1]][1],
                    H0 = coloc["summary"][[1]][2],
                    H1 = coloc["summary"][[1]][3],
                    H2 = coloc["summary"][[1]][4],
                    H3 = coloc["summary"][[1]][5],
                    H4 = coloc["summary"][[1]][6])
write.table(table, "./Output/G_unclassified_F_Porphyromonadaceae_RNT_rs35980751_coloc_results.txt", col.names = T, row.names = F, sep = "\t", quote = F)




########## G_unclassified_P_Firmicutes_RNT


N1 <- 1904
N2 <- 79148+61106 
type1 <- "quant"
type2 <- "cc"
exposure_name <- "Gut microbiome"
outcome_name <- "prostate cancer"

## Read in exposure data and format for colocalisation
exposure <- read.table("./Data/G_unclassified_P_Firmicutes_RNT_rs11788336_plusminus_1e6.txt", header = T, sep = "\t", fill = T)
exposure$eaf <- ifelse((exposure$all_BB < exposure$all_AA), exposure$all_maf, (1-exposure$all_maf))
head(exposure)
exposure_coloc <- format_data(exposure, type = "exposure", snps = NULL, header = T, snp_col = "rsid", beta_col = "frequentist_add_beta_1", se_col = "frequentist_add_se_1", eaf_col = "eaf", effect_allele_col = "alleleB", other_allele_col = "alleleA", pval_col = "frequentist_add_pvalue", chr_col = "chromosome", pos_col = "position")
exposure_coloc$exposure <- "Gut microbiome"
head(exposure_coloc)
exposure_coloc_snps <- as.vector(exposure_coloc$SNP)

## Read in outcome data and make sure all missing values are removed (big file so going to make sure that the SNPs that we need first are selected) then format for colocalisation analyses
prostate_data_for_coloc <- extract_outcome_data(
  snps = exposure_coloc_snps,
  outcomes = 'ebi-a-GCST006085') # 

prostate_data_for_coloc <- unique(prostate_data_for_coloc)
dim(prostate_data_for_coloc)
head(prostate_data_for_coloc)
prostate_data_for_coloc$rsid <- prostate_data_for_coloc$SNP

outcome_coloc <- format_data(prostate_data_for_coloc, type = "outcome", snps = NULL, header = T, snp_col = "rsid", beta_col = "beta.outcome", se_col = "se.outcome", eaf_col = "eaf.outcome", effect_allele_col = "effect_allele.outcome", other_allele_col = "other_allele.outcome", pval_col = "pval.outcome", chr_col = "chr", pos_col = "pos")# Duplicated SNPs removed in this step
outcome_coloc$outcome <- "Prostate cancer"
head(outcome_coloc)
dim(outcome_coloc)
## Remove any NAs
outcome_coloc_complete <- outcome_coloc[outcome_coloc$mr_keep.outcome == TRUE,] ## 
dim(outcome_coloc_complete) # 6815 SNPs
## Harmonise data (exposure and outcome data needs to be harmonised and then EAF converted to MAF)
df <- harmonise_data(exposure_dat = exposure_coloc, outcome_dat = outcome_coloc_complete, action = 2)
df <- df[!is.na(df$eaf.outcome), ]
## Get the MAF
for (j in 1:nrow(df)){
  if(df$eaf.exposure[j]>0.5){df$MAF1[j]=1-df$eaf.exposure[j]} else {df$MAF1[j]=df$eaf.exposure[j]}  
  if(df$eaf.outcome[j]>0.5){df$MAF2[j]=1-df$eaf.outcome[j]} else {df$MAF2[j]=df$eaf.outcome[j]}  
}
head(df)

## Perform analyses (note: that code assumes you have provided standard error and then squares this to get the variance)
coloc <- coloc.abf(dataset1=list(snp=df$SNP, beta=df$beta.exposure, 
                                 varbeta=df$se.exposure^2, MAF=df$MAF1, N=N1, type=type1),
                   dataset2=list(snp=df$SNP, beta=df$beta.outcome, 
                                 varbeta=df$se.outcome^2, MAF=df$MAF2, N=N2, type=type2),
                   p1=1e-04, p2=1e-04, p12=1e-05)

## Make table of results 
table <- data.frame(Exposure = exposure_name, Outcome = outcome_name,
                    NSNPs = coloc["summary"][[1]][1],
                    H0 = coloc["summary"][[1]][2],
                    H1 = coloc["summary"][[1]][3],
                    H2 = coloc["summary"][[1]][4],
                    H3 = coloc["summary"][[1]][5],
                    H4 = coloc["summary"][[1]][6])
write.table(table, "./Output/G_unclassified_P_Firmicutes_RNT_rs11788336_coloc_results.txt", col.names = T, row.names = F, sep = "\t", quote = F)


#################################################### MIBIOGEN ####################################################
### microbial traits: 

### genus.Allisonella.id.2174

## SET ARGUMENTS 
N1 <-  15316
N2 <- 79148+61106 
type1 <- "quant"
type2 <- "cc"
exposure_name <- "Gut microbiome"
outcome_name <- "Prostate cancer"

## Read in exposure data and format for colocalisation
exposure <- read.table("/Users/xh18454/Library/CloudStorage/OneDrive-UniversityofBristol/MR/snps_in_locus_gwshits_all_ancestry/genus.Allisonella.id.2174_rs602075_plusminus_1e6.txt", header = T, sep = "\t", fill = T)
head(exposure)
## create the alt.allele column
exposure <- exposure %>%
  mutate(alt.allele = ifelse(eff.allele == substr(alleles, 1, 1),
                             substr(alleles, 3, 3),
                             substr(alleles, 1, 1)))

head(exposure)


exposure_coloc <- format_data(exposure, type = "exposure", snps = NULL, header = T, snp_col = "SNP", beta_col = "beta", se_col = "SE", effect_allele_col = "eff.allele", other_allele_col = "alt.allele", pval_col = "P.wz", chr_col = "chromosome", pos_col = "position")
exposure_coloc$exposure <- "Gut microbiome"
head(exposure_coloc)
exposure_coloc <- exposure_coloc %>% arrange(SNP)
exposure_coloc_snps <- as.vector(exposure_coloc$SNP)


## Read in outcome data and make sure all missing values are removed (big file so going to make sure that the SNPs that we need first are selected) then format for colocalisation analyses
prostate_data_for_coloc <- extract_outcome_data(
  snps = exposure_coloc_snps,
  outcomes = 'ebi-a-GCST006085') # 

prostate_data_for_coloc <- unique(prostate_data_for_coloc)
dim(prostate_data_for_coloc)
head(prostate_data_for_coloc)
prostate_data_for_coloc$rsid <- prostate_data_for_coloc$SNP
prostate_data_for_coloc$MarkerName <- paste(prostate_data_for_coloc$chr, ":", prostate_data_for_coloc$pos, sep ="")

outcome_coloc <- format_data(prostate_data_for_coloc, type = "outcome", snps = NULL, header = T, snp_col = "MarkerName", beta_col = "beta.outcome", se_col = "se.outcome", eaf_col = "eaf.outcome", effect_allele_col = "effect_allele.outcome", other_allele_col = "other_allele.outcome", pval_col = "pval.outcome", chr_col = "chr", pos_col = "pos")# Duplicated SNPs removed in this step
outcome_coloc$outcome <- "Prostate cancer"
head(outcome_coloc)
dim(outcome_coloc)
## Remove any NAs
outcome_coloc_complete <- outcome_coloc[outcome_coloc$mr_keep.outcome == TRUE,] ## 
dim(outcome_coloc_complete)
outcome_coloc_complete <- outcome_coloc_complete %>% arrange(SNP)

## Harmonise data (exposure and outcome data needs to be harmonised and then EAF converted to MAF)
df <- harmonise_data(exposure_dat = exposure_coloc, outcome_dat = outcome_coloc_complete, action = 1)
head(df)
dim(df)

# extract the mafs and merge with df
all_mafs <- data.table::fread("/Users/xh18454/Library/CloudStorage/OneDrive-UniversityofBristol/MR/MiBioGen/all_mafs.txt.gz")
head(all_mafs)
colnames(all_mafs)[1] <- "SNP"
# Calculate the MAF1 (excluding the first column and ignoring NA values)
all_mafs$MAF1 <- rowMeans(all_mafs[,-1], na.rm = TRUE)
head(all_mafs)

# Select only the SNP and eaf columns from all_mafs
all_mafs_subset <- all_mafs[, c("SNP", "MAF1")]

# Merge the exposure dataframe with the subset of all_mafs on the SNP column
df <- merge(df, all_mafs_subset, by = "SNP")

# View the merged dataframe

head(df)
df <- df[!is.na(df$eaf.outcome), ]
dim(df) ## 7746 SNPs

df <- df[!is.na(df$MAF1), ]
dim(df) # 1 missing value in MAF1 column removed

## Get the MAF
for (j in 1:nrow(df)){
  #if(df$eaf.exposure[j]>0.5){df$MAF1[j]=1-df$eaf.exposure[j]} else {df$MAF1[j]=df$eaf.exposure[j]}  # remove this line
  if(df$eaf.outcome[j]>0.5){df$MAF2[j]=1-df$eaf.outcome[j]} else {df$MAF2[j]=df$eaf.outcome[j]}  
}
head(df)

## Perform analyses (note: that code assumes you have provided standard error and then squares this to get the variance)
coloc <- coloc.abf(dataset1=list(snp=df$SNP, beta=df$beta.exposure, 
                                 varbeta=df$se.exposure^2, MAF=df$MAF1, N=N1, type=type1),
                   dataset2=list(snp=df$SNP, beta=df$beta.outcome, 
                                 varbeta=df$se.outcome^2, MAF=df$MAF2, N=N2, type=type2),

                                      p1=1e-04, p2=1e-04, p12=1e-05)


## Make table of results 
table <- data.frame(Exposure = exposure_name, Outcome = outcome_name,
                    NSNPs = coloc["summary"][[1]][1],
                    H0 = coloc["summary"][[1]][2],
                    H1 = coloc["summary"][[1]][3],
                    H2 = coloc["summary"][[1]][4],
                    H3 = coloc["summary"][[1]][5],
                    H4 = coloc["summary"][[1]][6])
write.table(table, "./Output/genus.Allisonella.id.2174_rs602075_coloc_results.txt", col.names = T, row.names = F, sep = "\t", quote = F)



################### END #######################################
