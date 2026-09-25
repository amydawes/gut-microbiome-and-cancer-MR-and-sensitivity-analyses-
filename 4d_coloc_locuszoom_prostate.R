##############################################################################################################################
## PROJECT: MR analysis of the gut microbiome on prostate cancer in PRACTICAL, UK Biobank and FinnGen
## Script: Plotting locus zoom graphs from colocalisation analysis 
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
#devtools::install_github("boxiangliu/locuscomparer")
library(dplyr)
library(openxlsx)
library(readxl)
library(TwoSampleMR)
library(coloc)
library(data.table)
library(locuscomparer)


#################################################### HUGHES ####################################################
### microbial traits: G_Parabacteroides_RNT
#################################################
## READING IN AND FORMATTING DATA FOR PLOTTING ##
#################################################
### microbial traits: G_Parabacteroides_RNT, G_unclassified_F_Porphyromonadaceae_RNT, G_unclassified_P_Firmicutes_RNT

## Load in the data that were generated from colocalisation analysis (i.e., the "4c_coloc_analysis.R" script)
## This uses the same data as that previous script and note that this is assuming hg19

## Read in exposure data and format for colocalisation
exposure <- read.table("./Data/G_Parabacteroides_RNT_rs13207588_plusminus_1e6.txt", header = T, sep = "\t", fill = T)
head(exposure)
exposure_coloc <- format_data(exposure, type = "exposure", snps = NULL, header = T, snp_col = "rsid", beta_col = "frequentist_add_beta_1", se_col = "frequentist_add_se_1", eaf_col = "eaf", effect_allele_col = "alleleB", other_allele_col = "alleleA", pval_col = "frequentist_add_pvalue", chr_col = "chromosome", pos_col = "position")
exposure_coloc$exposure <- "Gut microbiome"
colnames(exposure_coloc)[1] <- c("rsid")
colnames(exposure_coloc)[7] <- c("pval")
head(exposure_coloc)
exposure_coloc_snps <- as.vector(exposure_coloc$rsid)

## Read in outcome data and make sure all missing values are removed (big file so going to make sure that the SNPs that we need first are selected) then format for colocalisation analyses
prostate_data_for_coloc <- extract_outcome_data(
  snps = exposure_coloc_snps,
  outcomes = 'ebi-a-GCST006085') # 

prostate_data_for_coloc <- unique(prostate_data_for_coloc)
dim(prostate_data_for_coloc)
head(prostate_data_for_coloc)
prostate_data_for_coloc$rsid <- prostate_data_for_coloc$SNP

outcome_coloc <- format_data(prostate_data_for_coloc, type = "outcome", snps = NULL, header = T, snp_col = "rsid", beta_col = "beta.outcome", se_col = "se.outcome", eaf_col = "eaf.outcome", effect_allele_col = "effect_allele.outcome", other_allele_col = "other_allele.outcome", pval_col = "pval.outcome", chr_col = "chr", pos_col = "pos")# Duplicated SNPs removed in this step
outcome_coloc$outcome <- "prostate cancer"
head(outcome_coloc)
dim(outcome_coloc)
## Remove any NAs
outcome_coloc_complete <- outcome_coloc[outcome_coloc$mr_keep.outcome == TRUE,] ## 
dim(outcome_coloc_complete)
head(outcome_coloc_complete)
colnames(outcome_coloc_complete)[9] <- c("rsid")
colnames(outcome_coloc_complete)[5] <- c("pval")
head(exposure_coloc)
head(outcome_coloc_complete)

## Formatting data
exposure_metal <- read_metal(exposure_coloc, marker_col = "rsid", pval_col = "pval")
exposure_pos <- get_position(exposure_metal, genome = c("hg19"))
outcome_metal <- read_metal(outcome_coloc_complete, marker_col = "rsid", pval_col = "pval")
outcome_pos <- get_position(outcome_metal, genome = c("hg19"))


##############
## PLOTTING ##
##############
## Need to make sure the rsid is the one associated with the gut microbiome trait (rs13207588 with G_Parabacteroides_RNT)
pdf("./Figures/G_Parabacteroides_RNT_rs13207588_locus_zoom.pdf", width=8, height=4)
  locus_compared <- locuscompare(exposure_pos, outcome_pos, marker_col1 = "rsid", pval_col1 = "pval", title1 = "gm", marker_col2 = "rsid", pval_col2 = "pval", title2 = "pc", snp = "rs13207588", population = "EUR", combine = TRUE, legend = TRUE, legend_position = "topright", lz_ylab_linebreak = FALSE, genome = c("hg19"))
  locus_compared
dev.off()

###############   G_unclassified_F_Porphyromonadaceae_RNT_rs35980751
## Read in exposure data and format for colocalisation
exposure <- read.table("./Data/G_unclassified_F_Porphyromonadaceae_RNT_rs35980751_plusminus_1e6.txt", header = T, sep = "\t", fill = T)
head(exposure)
exposure_coloc <- format_data(exposure, type = "exposure", snps = NULL, header = T, snp_col = "rsid", beta_col = "frequentist_add_beta_1", se_col = "frequentist_add_se_1", eaf_col = "eaf", effect_allele_col = "alleleB", other_allele_col = "alleleA", pval_col = "frequentist_add_pvalue", chr_col = "chromosome", pos_col = "position")
exposure_coloc$exposure <- "Gut microbiome"
colnames(exposure_coloc)[1] <- c("rsid")
colnames(exposure_coloc)[7] <- c("pval")
head(exposure_coloc)
exposure_coloc_snps <- as.vector(exposure_coloc$rsid)

## Read in outcome data and make sure all missing values are removed (big file so going to make sure that the SNPs that we need first are selected) then format for colocalisation analyses
prostate_data_for_coloc <- extract_outcome_data(
  snps = exposure_coloc_snps,
  outcomes = 'ebi-a-GCST006085') # 


prostate_data_for_coloc <- unique(prostate_data_for_coloc)
dim(prostate_data_for_coloc)
head(prostate_data_for_coloc)
prostate_data_for_coloc$rsid <- prostate_data_for_coloc$SNP

outcome_coloc <- format_data(prostate_data_for_coloc, type = "outcome", snps = NULL, header = T, snp_col = "rsid", beta_col = "beta.outcome", se_col = "se.outcome", eaf_col = "eaf.outcome", effect_allele_col = "effect_allele.outcome", other_allele_col = "other_allele.outcome", pval_col = "pval.outcome", chr_col = "chr", pos_col = "pos")# Duplicated SNPs removed in this step
outcome_coloc$outcome <- "prostate cancer"
head(outcome_coloc)
dim(outcome_coloc)
## Remove any NAs
outcome_coloc_complete <- outcome_coloc[outcome_coloc$mr_keep.outcome == TRUE,] ##
dim(outcome_coloc_complete)
head(outcome_coloc_complete)
colnames(outcome_coloc_complete)[9] <- c("rsid")
colnames(outcome_coloc_complete)[5] <- c("pval")
head(exposure_coloc)
head(outcome_coloc_complete)

## Formatting data
exposure_metal <- read_metal(exposure_coloc, marker_col = "rsid", pval_col = "pval")
exposure_pos <- get_position(exposure_metal, genome = c("hg19"))
outcome_metal <- read_metal(outcome_coloc_complete, marker_col = "rsid", pval_col = "pval")
outcome_pos <- get_position(outcome_metal, genome = c("hg19"))


##############
## PLOTTING ##
##############
## Need to make sure the rsid is the one associated with the gut microbiome trait (rs13207588 with G_Parabacteroides_RNT)
pdf("./Figures/G_unclassified_F_Porphyromonadaceae_RNT_rs35980751_locus_zoom.pdf", width=8, height=4)
locus_compared <- locuscompare(exposure_pos, outcome_pos, marker_col1 = "rsid", pval_col1 = "pval", title1 = "gm", marker_col2 = "rsid", pval_col2 = "pval", title2 = "pc", snp = "rs35980751", population = "EUR", combine = TRUE, legend = TRUE, legend_position = "topright", lz_ylab_linebreak = FALSE, genome = c("hg19"))
locus_compared
dev.off()



###############   G_unclassified_P_Firmicutes_RNT
## Read in exposure data and format for colocalisation
exposure <- read.table("./Data/G_unclassified_P_Firmicutes_RNT_rs11788336_plusminus_1e6.txt", header = T, sep = "\t", fill = T)
head(exposure)
exposure_coloc <- format_data(exposure, type = "exposure", snps = NULL, header = T, snp_col = "rsid", beta_col = "frequentist_add_beta_1", se_col = "frequentist_add_se_1", eaf_col = "eaf", effect_allele_col = "alleleB", other_allele_col = "alleleA", pval_col = "frequentist_add_pvalue", chr_col = "chromosome", pos_col = "position")
exposure_coloc$exposure <- "Gut microbiome"
colnames(exposure_coloc)[1] <- c("rsid")
colnames(exposure_coloc)[7] <- c("pval")
head(exposure_coloc)
exposure_coloc_snps <- as.vector(exposure_coloc$rsid)

## Read in outcome data and make sure all missing values are removed (big file so going to make sure that the SNPs that we need first are selected) then format for colocalisation analyses
prostate_data_for_coloc <- extract_outcome_data(
  snps = exposure_coloc_snps,
  outcomes = 'ebi-a-GCST006085') # 


prostate_data_for_coloc <- unique(prostate_data_for_coloc)
dim(prostate_data_for_coloc)
head(prostate_data_for_coloc)
prostate_data_for_coloc$rsid <- prostate_data_for_coloc$SNP

outcome_coloc <- format_data(prostate_data_for_coloc, type = "outcome", snps = NULL, header = T, snp_col = "rsid", beta_col = "beta.outcome", se_col = "se.outcome", eaf_col = "eaf.outcome", effect_allele_col = "effect_allele.outcome", other_allele_col = "other_allele.outcome", pval_col = "pval.outcome", chr_col = "chr", pos_col = "pos")# Duplicated SNPs removed in this step
outcome_coloc$outcome <- "prostate cancer"
head(outcome_coloc)
dim(outcome_coloc)
## Remove any NAs
outcome_coloc_complete <- outcome_coloc[outcome_coloc$mr_keep.outcome == TRUE,] ## 
dim(outcome_coloc_complete)
head(outcome_coloc_complete)
colnames(outcome_coloc_complete)[9] <- c("rsid")
colnames(outcome_coloc_complete)[5] <- c("pval")
head(exposure_coloc)
head(outcome_coloc_complete)

## Formatting data
exposure_metal <- read_metal(exposure_coloc, marker_col = "rsid", pval_col = "pval")
exposure_pos <- get_position(exposure_metal, genome = c("hg19"))
outcome_metal <- read_metal(outcome_coloc_complete, marker_col = "rsid", pval_col = "pval")
outcome_pos <- get_position(outcome_metal, genome = c("hg19"))


##############
## PLOTTING ##
##############
## Need to make sure the rsid is the one associated with the gut microbiome trait (rs13207588 with G_Parabacteroides_RNT)
pdf("./Figures/G_unclassified_P_Firmicutes_RNT_rs11788336_locus_zoom.pdf", width=8, height=4)
locus_compared <- locuscompare(exposure_pos, outcome_pos, marker_col1 = "rsid", pval_col1 = "pval", title1 = "gm", marker_col2 = "rsid", pval_col2 = "pval", title2 = "pc", snp = "rs11788336", population = "EUR", combine = TRUE, legend = TRUE, legend_position = "topright", lz_ylab_linebreak = FALSE, genome = c("hg19"))
locus_compared
dev.off()



#################################################### MiBioGen ####################################################

### microbial traits:  genus.Allisonella.id.2174

## Read in exposure data and format for colocalisation
exposure <- read.table("/Users/xh18454/Library/CloudStorage/OneDrive-UniversityofBristol/MR/snps_in_locus_gws_hits_all_ancestry/genus.Allisonella.id.2174_rs602075_plusminus_1e6.txt", header = T, sep = "\t", fill = T)
head(exposure)
## create the alt.allele column
exposure <- exposure %>%
  mutate(alt.allele = ifelse(eff.allele == substr(alleles, 1, 1),
                             substr(alleles, 3, 3),
                             substr(alleles, 1, 1)))

head(exposure)
exposure$exposure <- exposure$bac
dim(exposure) #  SNPs
#exposure_coloc_snps <- as.vector(exposure_coloc$SNP)
exposure_snps <- as.vector(exposure$SNP)

length(exposure_snps)
length(unique(exposure_snps))

sum(is.na(exposure_snps))
head(exposure_snps)
tail(exposure_snps)

########### received a 401 error in authentication of OpenGWAS a/c. testing on a small subset of the variants
check_reset()  #The purpose of check_reset() is specifically to determine whether OpenGWAS thinks your allowance has been exhausted and, if so, report when it resets.

test_snp <- exposure_snps[1:10]

test <- ieugwasr::associations(
  variants = test_snp,
  id = "ebi-a-GCST006085"
)

head(test). ## test worked


############ now going to complete the full 8992 snp extraction in batches of 1000 and combine them together

batch_size <- 1000

batches <- split(
  exposure_snps,
  ceiling(seq_along(exposure_snps) / batch_size)
)

prostate_batches <- vector("list", length(batches))

for (i in seq_along(batches)) {
  
  cat("Starting batch", i, "of", length(batches), 
      "(", length(batches[[i]]), "SNPs )\n")
  
  prostate_batches[[i]] <- extract_outcome_data(
    snps = batches[[i]],
    outcomes = "ebi-a-GCST006085"
  )
  
  cat("Completed batch", i, "\n")
}

############## then combine each batch ##################
prostate_data_for_coloc <- do.call(
  rbind,
  prostate_batches
)

###### The API failed after completing the first 1000 SNP batch, so I am trying to run them one at a time next...

dir.create(
  "./prostate_coloc_chunks",
  showWarnings = FALSE
)

chunk_size <- 500

snp_chunks <- split(
  exposure_snps,
  ceiling(seq_along(exposure_snps) / chunk_size)
)

length(snp_chunks)

###### Now run one at a time


for (i in seq_along(snp_chunks)) {
  
  outfile <- paste0(
    "./prostate_coloc_chunks/prostate_chunk_",
    i,
    ".rds"
  )
  
  # Skip chunks that have already been successfully saved
  if (file.exists(outfile)) {
    cat("Chunk", i, "already exists - skipping\n")
    next
  }
  
  cat(
    "\nStarting chunk", i,
    "of", length(snp_chunks),
    "(",
    length(snp_chunks[[i]]),
    "SNPs)\n"
  )
  
  result <- tryCatch(
    
    {
      associations(
        variants = snp_chunks[[i]],
        id = "ebi-a-GCST006085",
        proxies = 0
      )
    },
    
    error = function(e) {
      cat(
        "ERROR in chunk", i, ":\n",
        conditionMessage(e), "\n"
      )
      return(NULL)
    }
  )
  
  if (!is.null(result)) {
    
    saveRDS(
      result,
      outfile
    )
    
    cat(
      "Saved chunk", i,
      "with", nrow(result), "rows\n"
    )
  }
  
  # Small pause between requests
  Sys.sleep(2)
}

######### chunks 1,6,12 failed but the rest (other 15) chunks have snps extracted successfully and these files saved. 
##### Repeating chunks 1,6 and 12 only, then combining ...

failed_chunks <- c(1, 6, 12)

for (i in failed_chunks) {
  
  outfile <- paste0(
    "./prostate_coloc_chunks/prostate_chunk_",
    i,
    ".rds"
  )
  
  cat(
    "\nRetrying chunk", i,
    "of", length(snp_chunks[[i]]), "SNPs\n"
  )
  
  result <- tryCatch(
    
    {
      ieugwasr::associations(
        variants = snp_chunks[[i]],
        id = "ebi-a-GCST006085",
        proxies = 0
      )
    },
    
    error = function(e) {
      cat(
        "ERROR in chunk", i, ":\n",
        conditionMessage(e), "\n"
      )
      return(NULL)
    }
  )
  
  if (!is.null(result)) {
    
    saveRDS(
      result,
      outfile
    )
    
    cat(
      "Successfully saved chunk", i,
      "with", nrow(result), "rows\n"
    )
  }
  
  Sys.sleep(2)
}


############################## Now checking that all 18 chunks are present:
list.files(
  "./prostate_coloc_chunks",
  pattern = "^prostate_chunk_\\d+\\.rds$"
)
#### checking specifically that chunks 1,6,12 are present:
file.exists(
  paste0(
    "./prostate_coloc_chunks/prostate_chunk_",
    c(1, 6, 12),
    ".rds"
  )
)


############## combining chunks 1-18:
chunk_files <- list.files(
  "./prostate_coloc_chunks",
  pattern = "^prostate_chunk_\\d+\\.rds$",
  full.names = TRUE
)

prostate_data_for_coloc <- do.call(
  rbind,
  lapply(chunk_files, readRDS)
)
############ 
dim(prostate_data_for_coloc)
length(unique(prostate_data_for_coloc$rsid))
prostate_data_for_coloc$MarkerName <- paste(prostate_data_for_coloc$chr, ":", prostate_data_for_coloc$pos, sep ="")
head(prostate_data_for_coloc)
head(exposure)
MiBioGen_prostate_merge = merge(exposure, prostate_data_for_coloc, by.x = "SNP", by.y = "MarkerName")
head(MiBioGen_prostate_merge)


exposure_coloc <- format_data(MiBioGen_prostate_merge, type = "exposure", header = T, snp_col = "rsid", beta_col = "beta.x", se_col = "SE", effect_allele_col = "eff.allele", other_allele_col = "alt.allele", pval_col = "P.wz", chr_col = "chromosome", pos_col = "position.x", samplesize_col = "N")
exposure_coloc$exposure <- "Gut microbiome"
head(exposure_coloc)
dim(exposure_coloc) # 7752 SNPs available


outcome_coloc <- format_data(MiBioGen_prostate_merge, type = "outcome", snps = NULL, header = T, snp_col = "rsid", beta_col = "beta.y", se_col = "se", eaf_col = "eaf", effect_allele_col = "ea", other_allele_col = "nea", pval_col = "p", chr_col = "chr", pos_col = "position.y", samplesize_col = "n") # Duplicated SNPs removed in this step
outcome_coloc$outcome <- "prostate cancer"
head(outcome_coloc)
dim(outcome_coloc) # 7752 SNPs 
## Remove any NAs
outcome_coloc_complete <- outcome_coloc[outcome_coloc$mr_keep.outcome == TRUE,] ## of those 5899, there are 5899 complete SNPs
dim(outcome_coloc_complete)
head(outcome_coloc_complete)
colnames(outcome_coloc_complete)[10] <- c("pval")
colnames(exposure_coloc)[6] <- c("pval")
head(outcome_coloc_complete)
head(exposure_coloc)
exposure_coloc$rsid <- exposure_coloc$SNP
outcome_coloc_complete$rsid <- outcome_coloc_complete$SNP

# Extract positions from SNP column in exposure_coloc
exposure_pos <- data.frame(
  rsid = exposure_coloc$rsid,
  chr = exposure_coloc$chr.exposure,
  pos = exposure_coloc$pos.exposure,
  pval = exposure_coloc$pval,
  logp = -log10(exposure_coloc$pval)
)

# Extract positions from SNP column in outcome_coloc_complete
outcome_pos <- data.frame(
  rsid = outcome_coloc_complete$rsid,
  chr = outcome_coloc_complete$chr.outcome,
  pos = outcome_coloc_complete$pos.outcome,
  pval = outcome_coloc_complete$pval,
  logp = -log10(outcome_coloc_complete$pval)
)



# Check the structure of your data frames to confirm column names
str(exposure_pos)
str(outcome_pos)
unique(exposure_pos$chr)
unique(outcome_pos$chr)


##############
## PLOTTING ##
##############
Sys.setenv(MARIADB_TLS_DISABLE_PEER_VERIFICATION = "1")

## Need to make sure the rsid is the one associated with the gut microbiome trait 
pdf("./Figures/genus.Allisonella.id.2174_rs602075_locus_zoom.pdf", width=8, height=6)
locus_compared <- locuscompare(exposure_pos, outcome_pos, marker_col1 = "rsid", pval_col1 = "pval", title1 = "gut micriobiome", marker_col2 = "rsid", pval_col2 = "pval", title2 = "lung cancer", snp = "rs602075", population = "EUR", combine = TRUE, legend = TRUE, legend_position = "topright", lz_ylab_linebreak = FALSE, genome = c("hg19"))
print(locus_compared, res = 300)
dev.off()


########################################################################################

locus_compared <- locuscompare(
  exposure_pos,
  outcome_pos,
  marker_col1 = "rsid",
  pval_col1 = "pval",
  title1 = "gut microbiome",
  marker_col2 = "rsid",
  pval_col2 = "pval",
  title2 = "prostate cancer",
  snp = "rs602075",
  population = "EUR",
  combine = TRUE,
  legend = TRUE,
  legend_position = "bottomright",
  lz_ylab_linebreak = FALSE,
  genome = "hg19"
)

pdf("./Figures/genus.Allisonella.id.2174_rs602075_locus_zoom_v2.pdf",
    width = 8,
    height = 6)

print(locus_compared, res = 300)

dev.off()


