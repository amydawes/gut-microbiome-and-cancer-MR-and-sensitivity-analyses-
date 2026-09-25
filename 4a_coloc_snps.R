##############################################################################################################################
## PROJECT: MR analysis of the gut microbiome on lung cancer using PRACTICAL, UK Biobank and Finngen 
## Script: Creating the files that will be used to extract the summary-level data for colocalisation analyses
## Date: 25/09/26
##############################################################################################################################
## Clear space
rm(list=ls())

## Set working directory
setwd("/path/to/directory/prostate_cancer/")

### traits for colocalisation 

#"G_Parabacteroides_RNT", "G_unclassified_P_Firmicutes_RNT", "G_unclassified_F_Porphyromonadaceae_RNT", "genus.Allisonella.id.2174"
############################################################################
## SAVING FILES WITH THE SNPS WE WANT TO USE FOR COLOCALISATION ANALYSES  ##
############################################################################
## Read in exposure datasets that were prepared by previous step (i.e., "1_extraction_main.R" script) and identified as being related to pancreatic cancer (i.e., with the "2_mainMR_analysis_PanScan3+C4+FinnGen.R" script)
## Hughes
Hughes_exposure <- read.table("./Data/Hughes_prostate_exposure_data_main.txt", header = T, sep="\t")
head(Hughes_exposure)
Hughes_exposure$chr <- substr(Hughes_exposure$snpid, 1, regexpr(":", Hughes_exposure$snpid)-1)
Hughes_exposure$pos <- substr(Hughes_exposure$snpid, regexpr(":", Hughes_exposure$snpid)+1, nchar(Hughes_exposure$snpid)-4)
colnames(Hughes_exposure)
Hughes_snps <- Hughes_exposure[Hughes_exposure$mt %in% c("G_Parabacteroides_RNT", "G_unclassified_P_Firmicutes_RNT", "G_unclassified_F_Porphyromonadaceae_RNT"), c(1,2,12,13)]
write.table(Hughes_snps, "./Data/Hughes_coloc_snps.txt", sep="\t", quote = F, col.names = T, row.names = F)

# MiBioGen
MiBioGen_exposure <- read.table("./Data/MiBioGen_prostate_exposure_data_main.txt", header = T, sep="\t")
head(MiBioGen_exposure)
colnames(MiBioGen_exposure)
MiBioGen_snps <- MiBioGen_exposure[MiBioGen_exposure$exposure %in% c("genus.Allisonella.id.2174"), c(11,3,1,2)]
colnames(MiBioGen_snps)[1] = "chr"
colnames(MiBioGen_snps)[2] = "pos"
write.table(MiBioGen_snps, "./Data/MiBioGen_coloc_snps.txt", sep="\t", quote = F, col.names = T, row.names = F)

## Copy these text files across to HPC for extracting the summary-level data around these SNPs
## Then use Blue Pebble and "4b_coloc_extraction.R" script (submitted via a job using the "FGFP_mGWAS_coloc_extraction" job submission script)

#################################################### END ####################################################