####################################################################################################################################
######## PROJECT: MR analysis of the gut microbiome on prostate  cancer using PRACTICAL prostate cancer (any/overall prostate cancer) data
######## Script: Extract the microbiome-related SNPs for main analysis from the chromosomal-level data 
######## Date: 07/09/2026
#####################################################################################################################################
## Clear space
rm(list=ls())

## Set working directory
setwd("/path/to/directory/prostate_cancer/")

## INSTALL PACKAGES BEFOREHAND
#install.packages("openxlsx")
#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR") #to update the package
#devtools::install_github("MRCIEU/MRInstruments")
#install.packages("LDlinkR")
#install.packages("dplyr")
#install.packages("ieugwasr")
library(openxlsx)
library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(LDlinkR)
library(dplyr)
library(ieugwasr)

##edit r.environ file to ensure the ieu open GWAS API token for authentication is up to date
Sys.getenv("R_ENVIRON_USER")
usethis::edit_r_environ() ###You can check the location of your .Renviron file by running Sys.getenv("R_ENVIRON_USER") in R.
ieugwasr::get_opengwas_jwt()  ### To check that your token is being recognised, run ieugwasr::get_opengwas_jwt(). If it returns a long random string then you are authenticated.- see https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication

user() # check user
api_status()




#######################################################
## CREATING A VECTOR FOR ALL MICROBIOME-RELATED SNPS ##
#######################################################
### HUGHES
## Read in the file that has all the instrument information (from David Hughes paper)
Hughes_exposures <- read.table("Bug2DiseaseMRformat.txt", header = T, sep = "\t")
Hughes_exposures$MarkerName <- substr(Hughes_exposures$snpid, 1, nchar(Hughes_exposures$snpid)-4)
head(Hughes_exposures)
Hughes_instruments <- as.vector(Hughes_exposures$MarkerName)

### KURILSHIKOV
## Read in the file that has all the instrument information (downloaded from MiBioGen consortium - https://mibiogen.gcc.rug.nl/)
MiBioGen_exposures <- read.table("./Data/all_hits_from_website.txt", header = T, sep = "\t")

MiBioGen_hits <- MiBioGen_exposures[MiBioGen_exposures$P.weightedSumZ<=5e-8,]
MiBioGen_hits_formated <- format_data(MiBioGen_hits, type = "exposure", snps = NULL, header = TRUE, phenotype_col = "bac", snp_col = "rsID", beta_col = "beta", se_col = "SE", effect_allele_col = "eff.allele", other_allele_col = "ref.allele",  pval_col = "P.weightedSumZ", samplesize_col = "N", min_pval = 1e-200, z_col = "Z.weightedSumZ", chr_col = "chr", pos_col = "bp", log_pval = FALSE)
write.table(MiBioGen_hits_formated,"MiBioGen_GW_hits.txt", col.names = T, row.names = F, quote = F, sep = "\t")
MiBioGen_hits_clumped <- clump_data(MiBioGen_hits_formated, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, pop = "EUR") # there should be 30
MiBioGen_hits_clumped$MarkerName <- paste(MiBioGen_hits_clumped$chr.exposure, ":", MiBioGen_hits_clumped$pos.exposure, sep ="")
MiBioGen_hits_clumped <- MiBioGen_hits_clumped[order(MiBioGen_hits_clumped$SNP),]
head(MiBioGen_hits_clumped)
rm(MiBioGen_exposures,MiBioGen_hits,MiBioGen_hits_formated)

## Read in the supplementary table of results from the paper (downloaded - found in paper supplementary table)
MiBioGen_results <- read.xlsx("./Data/41588_2020_763_MOESM3_ESM.xlsx", sheet = "7_mbQTL_p<5x10-8", startRow = 3, colNames = T, rowNames = F)

colnames(MiBioGen_results)
MiBioGen_results <- MiBioGen_results[-c(43),-c(7, 13, 14, 15)]
MiBioGen_results <- MiBioGen_results[MiBioGen_results$Type=="Quant",]
dim(MiBioGen_results)
MiBioGen_results$MarkerName <- paste(MiBioGen_results$SNPChr, ":", MiBioGen_results$SNPChrPos, sep ="")
MiBioGen_results <- MiBioGen_results[order(MiBioGen_results$rsID),]
head(MiBioGen_results)

## Check that the SNPs are the same across the clumped data and those that are presented in the table # why did we need to compare supp table and website hits? 
colnames(MiBioGen_results)
colnames(MiBioGen_hits_clumped)
MiBioGen_results[,c(3,12,7,10)]
MiBioGen_hits_clumped[,c(3,16,11,8)]

## All 30 of them are in there but one of the SNPs is missing from the clumped data (genus.CandidatusSoleaferrea.id.11350)
MiBioGen_hits_clumped$SNP[MiBioGen_hits_clumped$exposure=="genus.CandidatusSoleaferrea.id.11350"] <- "rs830151"
head(MiBioGen_hits_clumped)

## Merge together based on Markername
MiBioGen_merged <- merge(MiBioGen_results, MiBioGen_hits_clumped, by.x = c("MarkerName", "TaxonName"), by.y = c("MarkerName", "exposure"))
colnames(MiBioGen_merged)
MiBioGen_merged <- MiBioGen_merged[order(MiBioGen_merged$SNP),]
MiBioGen_merged[,c(1,5,15,2,11,19)] ## Why selecting these columns?
MiBioGen_results[,c(3,12,7,10)]

## Two of the traits have different names so matching
MiBioGen_hits_clumped$exposure[MiBioGen_hits_clumped$exposure=="family.unknownfamily.id.1000001214"] <- "family.Gastranaerophilales.unknown.id.1000001214"
MiBioGen_hits_clumped$exposure[MiBioGen_hits_clumped$exposure=="genus.unknowngenus.id.1000001215"] <- "genus.Gastranaerophilales.unknown.id.1000001215"
tail(MiBioGen_hits_clumped)
write.table(MiBioGen_hits_clumped,"./MiBioGen_indep_GW_hits.txt", col.names = T, row.names = F, quote = F, sep = "\t")

MiBioGen_merged <- merge(MiBioGen_results, MiBioGen_hits_clumped, by.x = c("MarkerName", "TaxonName"), by.y = c("MarkerName", "exposure"))
colnames(MiBioGen_merged)
dim(MiBioGen_merged)
MiBioGen_merged <- MiBioGen_merged[order(MiBioGen_merged$SNP),]
MiBioGen_instruments <- as.vector(MiBioGen_merged$MarkerName)
rm(MiBioGen_merged) 

### ALL INSTRUMENTS TOGETHER
all_instruments <- as.vector(c(Hughes_instruments,MiBioGen_instruments))


################################################
## SEARCHING FOR THE MICROBIOME-SPECIFIC SNPS ##
################################################
## Search for the SNPs of interest in the lung cancer data

prostate_data <- extract_outcome_data(
  snps = all_instruments,
  outcomes = 'ebi-a-GCST006085'
)

View(prostate_data)
dim(prostate_data)
prostate_data$MarkerName <- paste(prostate_data$chr, ":", prostate_data$pos, sep = "")
head(prostate_data)
prostate_snp_data <- subset(prostate_data, MarkerName %in% all_instruments)
dim(prostate_snp_data) # should have 35 because there are 9 duplicates in the MiBioGen consortium data but only have 34
prostate_snp_data = unique(prostate_snp_data)
unique(prostate_snp_data$SNP)
## Check which ones are missing and whether this is due to using the MarkerName rather than rsid
all_snps <- as.vector(c(Hughes_exposures$rsid, MiBioGen_hits_clumped$SNP))
all_snps # there are 44 (i.e., 14 from Hughes and 30 from MiBioGen)
unique(all_snps) # though 35 of these are unique
prostate_test_snp_data <- subset(prostate_snp_data, SNP %in% all_snps) 
dim(prostate_test_snp_data) # still 34  so just not present - now to check which ones are not available
rm(prostate_test_snp_data)

## Now checking which SNPs aren't present
colnames(prostate_snp_data)[1] <- "rsid"
prostate_hughes_merge <- merge(Hughes_exposures, prostate_snp_data, by.x = c("rsid"), by.y = c("rsid"))
head(prostate_hughes_merge)
dim(prostate_hughes_merge) # only 11/14 SNPs available so may have to use proxies for 3 of these
unique(prostate_hughes_merge$mt) # representing 11 taxa - G_unclassified_P_Firmicutes_RNT and G_unclassified_P_Firmicutes_HB both included
write.table(prostate_hughes_merge, "./Data/Hughes_prostate_merge.txt", row.names = F, col.names = T, quote = F, sep = "\t")

prostate_MiBioGen_merge <- merge(MiBioGen_hits_clumped, prostate_snp_data, by.x = c("SNP"), by.y = c("rsid"))
View(prostate_MiBioGen_merge)
dim(prostate_MiBioGen_merge) # only 34 SNPs available so may have to use proxies - 34 SNPs 26 unique. 
unique(prostate_MiBioGen_merge$exposure) # representing 26 individual taxa
write.table(prostate_MiBioGen_merge, "./Data/MiBioGen_prostate_merge.txt", row.names = F, col.names = T, quote = F, sep = "\t")


###########################################
## SEARCHING FOR PROXIES SNPS FOR HUGHES ##
###########################################

## Can use LDLink (specifically, LDproxy_batch function)
dim(prostate_MiBioGen_merge)
head(prostate_MiBioGen_merge) # only 11/14 SNPs available so may have to use proxies for these
unique(prostate_MiBioGen_merge$exposure) # representing 27 unique taxa
Hughes_missing <- Hughes_exposures$rsid[-which(Hughes_exposures$rsid %in% prostate_MiBioGen_merge$rsid)]
Hughes_unique_missing <- unique(Hughes_missing) # no missing SNPs - don't need to find proxies 


New_Hughes_exposures = Hughes_exposures
New_Hughes_exposures$rsid <- as.character(New_Hughes_exposures$rsid)
New_Hughes_exposures$exposure_marker <- NULL
New_Hughes_exposures 

write.table(New_Hughes_exposures, "./Data/Hughes_prostate_exposure_data_main.txt", quote = F, sep = "\t", col.names = T, row.names = F)
## 14 SNPs, 11 of which were found originally, 3 of which were not; of those 3, 3 proxies were found; 1 of these were in the exposure data; 0 of these was in the outcome data 
## meaning only 11 results will be used (i.e., original 11)


#############################################
## SEARCHING FOR PROXIES SNPS FOR MIBIOGEN ##
#############################################
## Can use LDLink (specifically, LDproxy_batch function)
dim(prostate_MiBioGen_merge) # 30/30 SNPs available so may have to use proxies for these
prostate_MiBioGen_merge = unique(prostate_MiBioGen_merge)
unique(prostate_MiBioGen_merge$exposure)# representing 27 individual taxa
unique(prostate_MiBioGen_merge$SNP)
MiBioGen_missing <- MiBioGen_hits_clumped$SNP[-which(MiBioGen_hits_clumped$SNP %in% prostate_MiBioGen_merge$SNP)]
MiBioGen_unique_missing <- unique(MiBioGen_missing) # No missing snps - don't need to find proxies 

## MiBioGen data is the same as before
New_MiBioGen_exposures <- MiBioGen_hits_clumped
write.table(New_MiBioGen_exposures, "./Data/MiBioGen_prostate_exposure_data_main.txt", quote = F, sep = "\t", col.names = T, row.names = F)


###########################################################
## FINALISING NEW OUTCOME DATA FOR ANALYSIS WITH PROXIES ##
###########################################################
#rm(lung_snp_data, lung_hughes_merge, Hughes_instruments, MiBioGen_instruments)
Hughes_instruments <- as.vector(New_Hughes_exposures$MarkerName)
MiBioGen_instruments <- as.vector(New_MiBioGen_exposures$MarkerName)
all_instruments <- as.vector(c(Hughes_instruments,MiBioGen_instruments))
prostate_snp_data <- subset(prostate_data, MarkerName %in% all_instruments)
dim(prostate_snp_data) 
prostate_snp_data <- unique(prostate_snp_data)
dim(prostate_snp_data)
prostate_snp_data$rsid <- prostate_snp_data$SNP
## Merge back into the exposure datasets to get the rsid for panscan data and make the panscan outcome data for the two exposure files


prostate_hughes_merge <- merge(New_Hughes_exposures, prostate_snp_data, by.x = "rsid", by.y = "rsid", sort = FALSE)
head(prostate_hughes_merge)
dim(prostate_hughes_merge) # 11/14 SNPs available now
unique(prostate_hughes_merge$mt) # representing 11 individual taxa
Hughes_outcome_data <- prostate_hughes_merge[, c("rsid","chr", "pos", "effect_allele.outcome", "other_allele.outcome","beta.outcome","se.outcome","pval.outcome","eaf.outcome", "outcome", "samplesize.outcome")]
Hughes_outcome_data$outcome <- "prostate cancer"
head(Hughes_outcome_data)
dim(Hughes_outcome_data)
write.table(Hughes_outcome_data, "./Data/Hughes_prostate_cancer_outcome_data_main.txt", sep="\t", col.names = T, row.names = F, quote = F)

head(prostate_MiBioGen_merge)
dim(prostate_MiBioGen_merge) # 
prostate_MiBioGen_merge <- unique(prostate_MiBioGen_merge)
unique(prostate_MiBioGen_merge$exposure) # 30 SNPs representing 27 individual taxa
MiBioGen_outcome_data <- prostate_MiBioGen_merge[, c("SNP","chr", "pos", "effect_allele.outcome", "other_allele.outcome","beta.outcome","se.outcome","pval.outcome","eaf.outcome", "samplesize.outcome")]

MiBioGen_outcome_data$outcome <- "prostate cancer"
head(MiBioGen_outcome_data)
dim(MiBioGen_outcome_data)
write.table(MiBioGen_outcome_data, "./Data/MiBioGen_prostate_outcome_data_main.txt", sep="\t", col.names = T, row.names = F, quote = F)

##################################### END #####################################