######################################################################################################################
######## PROJECT: MR analysis of the gut microbiome on prostate cancer using PRACTICAL, UK Biobank and FinnGen
######## Script: Run the reverse MR analysis of prostate cancer and the gut microbiome 
######## Date: 25/09/26
######################################################################################################################
## Clear space
rm(list=ls())

## Set working directory
setwd("/path/to/directory/prostate_cancer/")

## INSTALL PACKAGES BEFOREHAND
#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR") #to update the package
#devtools::install_github("MRCIEU/MRInstruments")
#install.packages("plyr")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("xlsx")
#install.packages("png")
#install.packages("openxlsx")
#install.packages("ggrepel")
#install.packages("ggthemes")
#install.packages("ggforestplot")
library(ggrepel)
library(ggthemes)
library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(openxlsx)
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(png)
library(ggforestplot)



############################################################
## PREPARING EXPOSURE DATA FROM THE PRACTICAL DATA ##
############################################################
## 
prostate_data <- extract_instruments(outcomes = 'ebi-a-GCST006085', clump = TRUE, r2 = 0.001, kb= 10000, p1=5e-8, ) # This returns a set of LD clumped SNPs that are GWAS signficant for lung cancer
dim(prostate_data)
View(prostate_data)
prostate_data$MarkerName <- paste(prostate_data$chr.exposure, ":", prostate_data$pos.exposure, sep ="")
prostate_data <- prostate_data[order(prostate_data$SNP),]
head(prostate_data)
write.table(prostate_data,"./Data/PRACTICAL_GW_hits.txt", col.names = T, row.names = F, quote = F, sep = "\t")

## using the 16 genome wide significant SNPs 
prostate_cancer_snps <- as.data.frame(prostate_data$SNP)
colnames(prostate_cancer_snps) <- c("SNP")
write.table(prostate_cancer_snps,"./Data/prostate_cancer_snps.txt", col.names = F, row.names = F, quote = F, sep = "\t")
prostate_cancer_snps
dim(prostate_cancer_snps) # 142 snps

prostate_cancer_MarkerName <- as.data.frame(prostate_data$MarkerName)
colnames(prostate_cancer_MarkerName) <- c("MarkerName")
write.table(prostate_cancer_MarkerName,"./Data/prostate_cancer_MarkerName.txt", col.names = F, row.names = F, quote = F, sep = "\t")
prostate_cancer_MarkerName


######################################################################  HUGHES ET AL ######################################################################


###############################################################
## EXTRACTING LUNG CANCER SNPS FROM FGFP MGWAS RESULTS ##
###############################################################
## Use the "FGFP_mGWAS_SumStats_Extraction.sh" which is on HPC to extract this list (i.e., "prostate_cancer_snps.txt") from the FGFP full meta-analysis results.
## This outputs .txt files for each bug, each of which contains the number of SNPs available in the GWAS summary statistics
## These were copied across from HPC to the working directory

#########################################################################################################################
## READING IN AND HARMONIZING EXPOSURE AND OUTCOME DATA, RUNNING ALL ANALYSES, CREATING FIGURES AND OUTPUTTING RESULTS ##
#########################################################################################################################
## Exposure (prostate cancer) data is already formatted and has everything we need for MR analyses
head(prostate_data)
dim(prostate_data)
## Read in outcome (microbiome) datasets prepared by previous step (i.e., "FGFP_mGWAS_SumStats_Extraction_PanScan.sh" script with an interactive node on Blue Pebble)
## Assign names for loop to read in all the files 
bugs <- c("G_Parabacteroides_RNT", "G_unclassified_P_Firmicutes_RNT", "G_unclassified_F_Porphyromonadaceae_RNT")

## Make a list and combine all outcomes together, to read out MR results and results from heterogeneity statistics and pleiotropy tests
outcome_list <- list()
reverseMR_results_list <- list()
reverseMR_het_results_list <- list()
reverseMR_pleio_results_list <- list()

## Run the loop for analyses, sensitivity analyses and figures
for(i in bugs){
  
  # Read in outcome data
  location <- paste("./Data/reverse_MR_Hughes_results/", i, ".txt", sep = "")
  outcome <- read.table(location, sep="\t", header = T)  
  outcome$outcome <- i
  outcome$eaf = -9
  outcome$eaf <- ifelse((outcome$all_BB < outcome$all_AA), outcome$all_maf, (1-outcome$all_maf))
  outcome_restrict <- outcome[,c("rsid","snpid","chromosome","position","alleleA", "alleleB","frequentist_add_pvalue","frequentist_add_beta_1","frequentist_add_se_1","outcome", "eaf")]
  assign(paste("out_",i,sep=""), outcome_restrict)
  outcome_list[[i]] <- assign(i,outcome_restrict)

  # Format outcome data for analyses
  outcome_dat <- format_data(outcome_restrict, type = "outcome", snps = NULL, header = TRUE, snp_col = "rsid", phenotype_col = "outcome", beta_col = "frequentist_add_beta_1", se_col = "frequentist_add_se_1", effect_allele_col = "alleleB", other_allele_col = "alleleA",  pval_col = "frequentist_add_pvalue", eaf_col = "eaf", min_pval = 1e-200, chr_col = "chromosome", pos_col = "position", log_pval = FALSE)
  outcome_dat <- outcome_dat[outcome_dat$mr_keep.outcome == TRUE,] 
  prostate_data <- prostate_data[prostate_data$mr_keep.exposure == TRUE,] 
  prostate_data <- format_data(prostate_data, type = "exposure", snps = NULL, header = TRUE, snp_col = "SNP", phenotype_col = "exposure", beta_col = "beta.exposure", se_col = "se.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure",  pval_col = "pval.exposure", eaf_col = "eaf.exposure", min_pval = 1e-200, chr_col = "chr.exposure", pos_col = "pos.exposure", log_pval = FALSE)

  # Harmonize between exposure and outcome data
  prostate_microbiome_dat <- harmonise_data(prostate_data, outcome_dat, action =1)
  
  gwasinfo("ebi-a-GCST006085") # 29863 cases 55586 controls
  
  ## Add the prevalence to the binary phenotypes from FGFP supplementary data (Supplementary Table ) , where the continuous ones will have a prevalence of 1
  prostate_microbiome_dat$ncase.outcome <- NA
  prostate_microbiome_dat$ncontrol.outcome <- NA
  prostate_microbiome_dat$prevalence.outcome <- NA
  
  #G_Parabacteroides_RNT
  prostate_microbiome_dat$ncase.outcome[prostate_microbiome_dat$outcome == "G_Parabacteroides_RNT"] <- 2223
  prostate_microbiome_dat$ncontrol.outcome[prostate_microbiome_dat$outcome == "G_Parabacteroides_RNT"] <- 2223
  prostate_microbiome_dat$samplesize.outcome[prostate_microbiome_dat$outcome == "G_Parabacteroides_RNT"] <- 2223
  prostate_microbiome_dat$prevalence.outcome[prostate_microbiome_dat$outcome == "G_Parabacteroides_RNT"] <- 1
  
  #G_unclassified_P_Firmicutes_RNT
  prostate_microbiome_dat$ncase.outcome[prostate_microbiome_dat$outcome == "G_unclassified_P_Firmicutes_RNT"] <- 1904
  prostate_microbiome_dat$ncontrol.outcome[prostate_microbiome_dat$outcome == "G_unclassified_P_Firmicutes_RNT"] <- 1904
  prostate_microbiome_dat$samplesize.outcome[prostate_microbiome_dat$outcome == "G_unclassified_P_Firmicutes_RNT"] <- 1904
  prostate_microbiome_dat$prevalence.outcome[prostate_microbiome_dat$outcome == "G_unclassified_P_Firmicutes_RNT"] <- 1
  
  #G_unclassified_F_Porphyromonadaceae_RNT
  prostate_microbiome_dat$ncase.outcome[prostate_microbiome_dat$outcome == "G_unclassified_F_Porphyromonadaceae_RNT"] <- 1420
  prostate_microbiome_dat$ncontrol.outcome[prostate_microbiome_dat$outcome == "G_unclassified_F_Porphyromonadaceae_RNT"] <- 1420
  prostate_microbiome_dat$samplesize.outcome[prostate_microbiome_dat$outcome == "G_unclassified_F_Porphyromonadaceae_RNT"] <- 1420
  prostate_microbiome_dat$prevalence.outcome[prostate_microbiome_dat$outcome == "G_unclassified_F_Porphyromonadaceae_RNT"] <- 1
  
  prostate_microbiome_dat$ncase.exposure <- 79148
  prostate_microbiome_dat$ncontrol.exposure <- 61106
  prostate_microbiome_dat$prevalence.exposure <- 0.1
  prostate_microbiome_dat$units.exposure <- "log odds"
  
  
  Hughes_steiger <- steiger_filtering(prostate_microbiome_dat)

  Hughes_steiger$samplesize.exposure <- Hughes_steiger$ncase.exposure + Hughes_steiger$ncontrol.exposure
  
  # Step 2: Recalculate F for binary exposures
  Hughes_steiger$F_exposure <- (Hughes_steiger$rsq.exposure * (Hughes_steiger$samplesize.exposure - 1 - 1)) / 
    ((1 - Hughes_steiger$rsq.exposure) * 1)
  
  
  # Save the updated table
  write.table(Hughes_steiger, paste("./Output/reverseMR_",i,"_steiger_filtering.txt", sep=""), sep = "\t", 
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  #keep_snps <- as.vector(Hughes_steiger$steiger_dir == "TRUE")
  keep_snps <- Hughes_steiger$SNP[Hughes_steiger$steiger_dir == "TRUE"]
  exposure_steiger <- Hughes_steiger[Hughes_steiger$SNP %in% keep_snps,]
  lung_microbiome_mr_results_steiger <- mr(exposure_steiger, method_list=c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
  
  # Outputting each MR result into results_list
  a <- cbind.data.frame(lung_microbiome_mr_results_steiger$exposure, lung_microbiome_mr_results_steiger$outcome,lung_microbiome_mr_results_steiger$nsnp,lung_microbiome_mr_results_steiger$method,lung_microbiome_mr_results_steiger$b,lung_microbiome_mr_results_steiger$se,lung_microbiome_mr_results_steiger$pval)
  name <- paste("reverseMR_steiger_filter_",i,sep="")
  reverseMR_results_list[[paste("reverseMR_steiger_filter_",i,sep="")]] <- assign(name,a)
  rm(a,name)
  
  # Merge all results files together and input into text files
  all_reverseMR_results_steiger_filter <- do.call(rbind, reverseMR_results_list)
  rownames(all_reverseMR_results_steiger_filter) <- NULL
  write.table(all_reverseMR_results_steiger_filter, "./Output/ReverseMR_Hughes_results_steiger_filter.txt", row.names = F, col.names = T, quote = F, sep="\t")
  

  # Run the MR analyses, heterogeneity statistics and pleiotropy tests
  prostate_microbiome_mr_results <- mr(prostate_microbiome_dat, method_list=c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
  write.table(prostate_microbiome_mr_results, "./Output/Reverse_MR_Hughes_microbiome_PRACTICAL_MR_results.txt", sep = "\t", col.names = T, row.names = F, quote = F)
  
  # Outputting each MR result into results_list
  a <- cbind.data.frame(prostate_microbiome_mr_results$exposure, prostate_microbiome_mr_results$outcome,prostate_microbiome_mr_results$nsnp,prostate_microbiome_mr_results$method,prostate_microbiome_mr_results$b,prostate_microbiome_mr_results$se,prostate_microbiome_mr_results$pval)
  name <- paste("reverseMR_",i,sep="")
  reverseMR_results_list[[paste("reverseMR_",i,sep="")]] <- assign(name,a)
  rm(a,name)

  # Heterogeneity statistics
  prostate_microbiome_het <- mr_heterogeneity(prostate_microbiome_dat)
  b <- cbind.data.frame(prostate_microbiome_het$exposure, prostate_microbiome_het$outcome, prostate_microbiome_het$method, prostate_microbiome_het$Q_pval)
  name <- paste("reverseMR_het_",i,sep="")
  reverseMR_het_results_list[[paste("reverseMR_",i,sep="")]] <- assign(name, b)
  rm(b,name)

  # Pleiotropy tests
  prostate_microbiome_pleio <- mr_pleiotropy_test(prostate_microbiome_dat)
  c <- cbind.data.frame(prostate_microbiome_pleio$exposure, prostate_microbiome_pleio$outcome, prostate_microbiome_pleio$egger_intercept, prostate_microbiome_pleio$se, prostate_microbiome_pleio$pval)
  name <- paste("reverseMR_",i,sep="")
  reverseMR_pleio_results_list[[paste("reverseMR_",i,sep="")]] <- assign(name, c)
  rm(c,name) 

  # Scatter plot
  png(paste("./Figures/reverseMR_",i,"_scatter.png", sep=""))
  print(mr_scatter_plot(prostate_microbiome_mr_results, prostate_microbiome_dat))
  dev.off()

  # Forest plot
  reverseMR_single <- mr_singlesnp(prostate_microbiome_dat)
  png(paste("./Figures/reverseMR_",i,"_forest.png", sep=""), width = 2400, height = 5000, res = 300)
  print(mr_forest_plot(reverseMR_single))
  dev.off()
  
  # Funnel plot
  png(paste("./Figures/reverseMR_",i,"_funnel.png", sep=""))
  print(mr_funnel_plot(reverseMR_single))
  dev.off()
  
  # Leave one out analysis
  reverseMR_loo <- mr_leaveoneout(prostate_microbiome_dat)
  png(paste("./Figures/reverseMR_",i,"_loo.png", sep=""),width = 2400, height = 5000, res = 300)
  print(mr_leaveoneout_plot(reverseMR_loo))
  dev.off()
  
  # Remove objects before next iteration of loop
  rm(location, outcome, outcome_restrict, outcome_dat, prostate_microbiome_dat, prostate_microbiome_mr_results)
  
}

## Make a data frame for all outcomes
all_outcomes <- do.call(rbind, outcome_list)
head(all_outcomes)

# Merge all results files together and input into text files
all_reverseMR_results <- do.call(rbind, reverseMR_results_list)
rownames(all_reverseMR_results) <- NULL
write.table(all_reverseMR_results, "./Output/ReverseMR_Hughes_results.txt", row.names = F, col.names = T, quote = F, sep="\t")
all_reverseMR_het_results <- do.call(rbind, reverseMR_het_results_list)
rownames(all_reverseMR_het_results) <- NULL
write.table(all_reverseMR_het_results, "./Output/ReverseMR_het_Hughes_results.txt", row.names = F, col.names = T, quote = F, sep="\t")
all_reverseMR_pleio_results <- do.call(rbind, reverseMR_pleio_results_list)
rownames(all_reverseMR_pleio_results) <- NULL
write.table(all_reverseMR_pleio_results, "./Output/ReverseMR_pleio_Hughes_results.txt", row.names = F, col.names = T, quote = F, sep="\t")


#################################################
## ADDING FOREST PLOTS OF ALL RESULTS TOGETHER ##
#################################################
colnames(all_reverseMR_results) <- c("exposure","outcome","nsnp","method","beta","se","pval")
head(all_reverseMR_results)
fres <- all_reverseMR_results[(all_reverseMR_results$method=="Inverse variance weighted"),]
fres$lci <- fres$beta - (1.96*fres$se)
fres$uci <- fres$beta + (1.96*fres$se)
fres$sig <- fres$pval < 0.05
fres$category <-"<0.05"
fres$label <- fres$outcome
fres$observation <- 1:nrow(fres) 
head(fres)

ggplot(data=fres, aes(y=observation, x=beta, xmin=lci, xmax=uci, colour=method)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(fres), labels=fres$label) +
  labs(title="Effect of prostate cancer on the gut microbiome", x="Beta", y = "Microbial trait (Hughes)") +
  geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) +
  theme_classic()
ggsave("./Figures/forestplot_reverseMR_Hughes.png", width=7, height=8)

######################################################################  MIBIOGEN ######################################################################

####### Using MiBioGen data for the Europeans

#### copy 6_MiBioGen_Eur_SumStats_Extraction.sh and prostate_cancer_snps.txt to HPC and run this script. 
## Assign names for loop to read in all the files 
bugs <- c("genus.Allisonella.id.2174")



## Make a list and combine all outcomes together, to read out MR results and results from heterogeneity statistics and pleiotropy tests
outcome_list <- list()
reverseMR_results_list <- list()
reverseMR_het_results_list <- list()
reverseMR_pleio_results_list <- list()


### extract the lung data and format data 
prostate_data <- extract_instruments(outcomes = 'ebi-a-GCST006085', clump = TRUE, r2 = 0.001, kb= 10000, p1=5e-8, ) # This returns a set of LD clumped SNPs that are GWAS signficant for lung cancer
dim(prostate_data)
View(prostate_data)
prostate_data$MarkerName <- paste(prostate_data$chr.exposure, ":", prostate_data$pos.exposure, sep ="")
prostate_data <- prostate_data[prostate_data$mr_keep.exposure == TRUE,] 
prostate_data <- format_data(prostate_data, type = "exposure", snps = NULL, header = TRUE, snp_col = "MarkerName", phenotype_col = "exposure", beta_col = "beta.exposure", se_col = "se.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure",  pval_col = "pval.exposure", eaf_col = "eaf.exposure", min_pval = 1e-200, chr_col = "chr.exposure", pos_col = "pos.exposure", log_pval = FALSE)

## Run the loop for analyses, sensitivity analyses and figures
for(i in bugs){
  
  # Read in outcome data
  location <- paste("./Data/reverse_MR_MiBioGen_results/", i, ".txt", sep = "")
  outcome <- read.table(location, sep="\t", header = T) 
  outcome$outcome <- i
  ## create the alt.allele column
  outcome <- outcome %>%
    mutate(alt.allele = ifelse(eff.allele == substr(alleles, 1, 1),
                               substr(alleles, 3, 3),
                               substr(alleles, 1, 1)))
  
  head(outcome)
  
   outcome <- outcome %>% separate(SNP, into = c("Chromosome", "Position"), sep = ":", remove = FALSE)
  assign(paste("out_",i,sep=""), outcome)
  outcome_list[[i]] <- assign(i,outcome)
  
  # Format outcome data for analyses
  outcome_dat <- format_data(outcome, type = "outcome", snps = NULL, header = TRUE, snp_col = "SNP", phenotype_col = "outcome", beta_col = "beta", se_col = "SE", effect_allele_col = "eff.allele", other_allele_col = "alt.allele",  pval_col = "P.wz",  min_pval = 1e-200, chr_col = "chromosome", pos_col = "position", samplesize_col = "N", log_pval = FALSE)
  head(outcome_dat) 
  
  # Harmonize between exposure and outcome data
  prostate_microbiome_dat <- harmonise_data(prostate_data, outcome_dat, action =1)
  
  prostate_microbiome_dat$units.outcome <- "SD"
  prostate_microbiome_dat$ncase.exposure <- 29863
  prostate_microbiome_dat$ncontrol.exposure <- 55586
  prostate_microbiome_dat$prevalence.exposure <- 0.1
  prostate_microbiome_dat$units.exposure <- "log odds"
  
  
  MiBioGen_steiger <- steiger_filtering(prostate_microbiome_dat)
  MiBioGen_steiger$samplesize.exposure <- MiBioGen_steiger$ncase.exposure + MiBioGen_steiger$ncontrol.exposure
  
  # Step 2: Recalculate F for binary exposures
  MiBioGen_steiger$F_statistic <- (MiBioGen_steiger$rsq.exposure * (MiBioGen_steiger$samplesize.exposure - 1 - 1)) / 
    ((1 - MiBioGen_steiger$rsq.exposure) * 1)
  
  
  # Save the updated table
  write.table(MiBioGen_steiger, paste("./Output/reverseMR_",i,"_steiger_filtering.txt", sep=""), sep = "\t", 
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  #keep_snps <- as.vector(Hughes_steiger$steiger_dir == "TRUE")
  keep_snps <- MiBioGen_steiger$SNP[MiBioGen_steiger$steiger_dir == "TRUE"]
  exposure_steiger <- MiBioGen_steiger[MiBioGen_steiger$SNP %in% keep_snps,]
  breast_microbiome_mr_results_steiger <- mr(exposure_steiger, method_list=c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
  
  # Outputting each MR result into results_list
  a <- cbind.data.frame(breast_microbiome_mr_results_steiger$exposure, breast_microbiome_mr_results_steiger$outcome,breast_microbiome_mr_results_steiger$nsnp,breast_microbiome_mr_results_steiger$method,breast_microbiome_mr_results_steiger$b,breast_microbiome_mr_results_steiger$se,breast_microbiome_mr_results_steiger$pval)
  name <- paste("reverseMR_steiger_filter_",i,sep="")
  reverseMR_results_list[[paste("reverseMR_steiger_filter_",i,sep="")]] <- assign(name,a)
  rm(a,name)
  
  # Merge all results files together and input into text files
  all_reverseMR_results_steiger_filter <- do.call(rbind, reverseMR_results_list)
  rownames(all_reverseMR_results_steiger_filter) <- NULL
  write.table(all_reverseMR_results_steiger_filter, "./Output/ReverseMR_MiBioGen_results_steiger_filter.txt", row.names = F, col.names = T, quote = F, sep="\t")
  
  
  
  
  # Run the MR analyses, heterogeneity statistics and pleiotropy tests
  prostate_microbiome_mr_results <- mr(prostate_microbiome_dat, method_list=c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
  write.table(prostate_microbiome_mr_results,paste("./Output/reverseMR_results_", i,".txt", sep=""), row.names = F, col.names = T, quote = F, sep="\t") 
  
  
  # Outputting each MR result into results_list
  a <- cbind.data.frame(prostate_microbiome_mr_results$exposure, prostate_microbiome_mr_results$outcome,prostate_microbiome_mr_results$nsnp,prostate_microbiome_mr_results$method,prostate_microbiome_mr_results$b,prostate_microbiome_mr_results$se,prostate_microbiome_mr_results$pval)
  name <- paste("reverseMR_",i,sep="")
  reverseMR_results_list[[paste("reverseMR_",i,sep="")]] <- assign(name,a)
  rm(a,name)
  
  # Heterogeneity statistics
  prostate_microbiome_het <- mr_heterogeneity(prostate_microbiome_dat)
  b <- cbind.data.frame(prostate_microbiome_het$exposure, prostate_microbiome_het$outcome, prostate_microbiome_het$method, prostate_microbiome_het$Q_pval)
  name <- paste("reverseMR_het_",i,sep="")
  reverseMR_het_results_list[[paste("reverseMR_",i,sep="")]] <- assign(name, b)
  rm(b,name)
  
  # Pleiotropy tests
  prostate_microbiome_pleio <- mr_pleiotropy_test(prostate_microbiome_dat)
  c <- cbind.data.frame(prostate_microbiome_pleio$exposure, prostate_microbiome_pleio$outcome, prostate_microbiome_pleio$egger_intercept, prostate_microbiome_pleio$se, prostate_microbiome_pleio$pval)
  name <- paste("reverseMR_",i,sep="")
  reverseMR_pleio_results_list[[paste("reverseMR_",i,sep="")]] <- assign(name, c)
  rm(c,name) 
  
  # Scatter plot
  png(paste("./Figures/reverseMR_",i,"_scatter.png", sep=""), width = 1600, height = 2400, res = 300)
  print(mr_scatter_plot(prostate_microbiome_mr_results, prostate_microbiome_dat))
  dev.off()
  
  # Forest plot
  reverseMR_single <- mr_singlesnp(prostate_microbiome_dat)
  png(paste("./Figures/reverseMR_",i,"_forest.png", sep=""), width = 1600, height = 3800, res = 300)
  print(mr_forest_plot(reverseMR_single))
  dev.off()
  
  # Funnel plot
  png(paste("./Figures/reverseMR_",i,"_funnel.png", sep=""), width = 1600, height = 2400, res = 300)
  print(mr_funnel_plot(reverseMR_single))
  dev.off()
  
  # Leave one out analysis
  reverseMR_loo <- mr_leaveoneout(prostate_microbiome_dat)
  png(paste("./Figures/reverseMR_",i,"_loo.png", sep=""), width = 1600, height = 3800, res = 300)
  print(mr_leaveoneout_plot(reverseMR_loo))
  dev.off()
  
  # Remove objects before next iteration of loop
  rm(location, outcome, outcome_restrict, outcome_dat, prostate_microbiome_dat, prostate_microbiome_mr_results)
  
}

## Make a data frame for all outcomes
all_outcomes <- do.call(rbind, outcome_list)
head(all_outcomes)

# Merge all results files together and input into text files
all_reverseMR_results <- do.call(rbind, reverseMR_results_list)
rownames(all_reverseMR_results) <- NULL
write.table(all_reverseMR_results, "./Output/ReverseMR_MiBioGen_results.txt", row.names = F, col.names = T, quote = F, sep="\t")
all_reverseMR_het_results <- do.call(rbind, reverseMR_het_results_list)
rownames(all_reverseMR_het_results) <- NULL
write.table(all_reverseMR_het_results, "./Output/ReverseMR_het_MiBioGen_results.txt", row.names = F, col.names = T, quote = F, sep="\t")
all_reverseMR_pleio_results <- do.call(rbind, reverseMR_pleio_results_list)
rownames(all_reverseMR_pleio_results) <- NULL
write.table(all_reverseMR_pleio_results, "./Output/ReverseMR_pleio_MiBioGen_results.txt", row.names = F, col.names = T, quote = F, sep="\t")


#################################################
## ADDING FOREST PLOTS OF ALL RESULTS TOGETHER ##
#################################################
colnames(all_reverseMR_results) <- c("exposure","outcome","nsnp","method","beta","se","pval")
head(all_reverseMR_results)
fres <- all_reverseMR_results[(all_reverseMR_results$method=="Inverse variance weighted"),]
fres$lci <- fres$beta - (1.96*fres$se)
fres$uci <- fres$beta + (1.96*fres$se)
fres$sig <- fres$pval < 0.05
fres$category <-"<0.05"
fres$label <- fres$outcome
fres$observation <- 1:nrow(fres) 
head(fres)

ggplot(data=fres, aes(y=observation, x=beta, xmin=lci, xmax=uci, colour=method)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(fres), labels=fres$label) +
  labs(title="Effect of prostate cancer on the gut microbiome", x="Beta", y = "Microbial trait (MiBioGen)") +
  geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) +
  theme_classic()
ggsave("./Figures/forestplot_reverseMR_MiBioGen.png", width=7, height=8)



