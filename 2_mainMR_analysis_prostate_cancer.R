######################################################################################################################
######## PROJECT: MR analysis of the gut microbiome on prostate cancer using PRACTICAL prostate cancer GWAS data, UK Biobank prostate cancer GWAS data and FinnGen prostate cancer data
######## Script: Run the main MR analysis of the gut microbiome and prostate cancer
######## Date: 25/09/26
######################################################################################################################
## Clear space
rm(list=ls())

## Set working directory

setwd("/path/to/directory/prostate_cancer/")
getwd()
## INSTALL PACKAGES BEFOREHAND
#install.packages("usethis")
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
#devtools::install_github("NightingaleHealth/ggforestplot")
library(calibrate)
library(ggrepel)
library(ggthemes)
library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(openxlsx)
library(plyr) 
library(dplyr)
library(ggplot2)
library(png)
library(ggforestplot) 





##################################################### PART 1 - PRACTICAL GWAS #####################################################


###########################################################################################
## PREPARING GUT MICROBIOME GWAS FILES USING HUGHES ET AL. AND MIBIOGEN WITH PRACTICAL  ##
###########################################################################################
## Read in exposure datasets prepared by previous step (i.e., "1_extraction_main.R" script)
## Hughes
Hughes_exposure_dat <- read_exposure_data("./Data/Hughes_prostate_exposure_data_main.txt", sep="\t", phenotype_col = "mt", snp_col = "rsid", beta_col = "beta", se_col = "se", pval_col = "P_value", eaf_col = "eaf", effect_allele_col = "allele_B", other_allele_col = "allele_A", samplesize_col = "n")
head(Hughes_exposure_dat)

# MiBioGen
MiBioGen_exposure_dat <- read_exposure_data("./Data/MiBioGen_prostate_exposure_data_main.txt", sep="\t", phenotype_col = "exposure", snp_col = "SNP", beta_col = "beta.exposure", se_col = "se.exposure", pval_col = "pval.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure", samplesize_col = "samplesize.exposure")
head(MiBioGen_exposure_dat)


##################################################
## READING IN THE OUTCOME DATA FROM PRACTICAL - script  1_extraction_main ##
##################################################
## Read in outcome datasets prepared by previous step (i.e., "1_extraction_main.R" script)
Hughes_outcome_dat <- read_outcome_data("./Data/Hughes_prostate_cancer_outcome_data_main.txt", sep="\t", snp_col="rsid", beta_col="beta.outcome", se_col="se.outcome", effect_allele_col = "effect_allele.outcome", other_allele_col = "other_allele.outcome", pval_col = "pval.outcome", eaf_col = "eaf.outcome", phenotype_col = "outcome")
head(Hughes_outcome_dat)

MiBioGen_outcome_dat <- read_outcome_data("./Data/MiBioGen_prostate_outcome_data_main.txt", sep="\t", snp_col="SNP", beta_col="beta.outcome", se_col="se.outcome", effect_allele_col = "effect_allele.outcome", other_allele_col = "other_allele.outcome", pval_col = "pval.outcome", eaf_col = "eaf.outcome", phenotype_col = "outcome")
head(MiBioGen_outcome_dat)


####################################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND PROSTATE CANCER (OUTCOME) ##
####################################################################################
## MAIN ANALYSIS INCLUDING THE TOP HITS FROM OUR PUBLISHED META-ANALYSES
Hughes_dat <- harmonise_data(Hughes_exposure_dat, Hughes_outcome_dat, action =2)
Hughes_mr_results <- mr(Hughes_dat, method_list=c("mr_wald_ratio"))
write.table(Hughes_mr_results, "./Output/Hughes_microbiome_prostate_cancer_results.txt", sep = "\t", col.names = T, row.names = F, quote = F)


# Scatter plot
png(paste("./Figures/Hughes_MR_scatter.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_scatter_plot(Hughes_mr_results, Hughes_dat))
dev.off()

# Forest plot
MR_single <- mr_singlesnp(Hughes_dat)
png(paste("./Figures/Hughes_MR_forest.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_forest_plot(MR_single))
dev.off()

# Funnel plot
png(paste("./Figures/Hughes_MR_funnel.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_funnel_plot(MR_single))
dev.off()

# Leave one out analysis
MR_loo <- mr_leaveoneout(Hughes_dat)
png(paste("./Figures/Hughes_MR_loo.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_leaveoneout_plot(MR_loo))
dev.off()





MiBioGen_dat <- harmonise_data(MiBioGen_exposure_dat, MiBioGen_outcome_dat, action =1)
MiBioGen_mr_results <- mr(MiBioGen_dat, method_list=c("mr_wald_ratio", "mr_ivw"))
write.table(MiBioGen_mr_results, "./Output/MiBioGen_microbiome_prostate_cancer_results.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# Scatter plot
png(paste("./Figures/MiBioGen_MR_scatter.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_scatter_plot(MiBioGen_mr_results, MiBioGen_dat))
dev.off()

# Forest plot
MR_single <- mr_singlesnp(MiBioGen_dat)
png(paste("./Figures/MiBioGen_MR_forest.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_forest_plot(MR_single))
dev.off()

# Funnel plot
png(paste("./Figures/MiBioGen_MR_funnel.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_funnel_plot(MR_single))
dev.off()

# Leave one out analysis
MR_loo <- mr_leaveoneout(MiBioGen_dat)
png(paste("./Figures/MiBioGen_MR_loo.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_leaveoneout_plot(MR_loo))
dev.off()


#######################################
## steiger filtering 
######################################

## Using the harmonised data 

## Add the prevalence to the binary phenotypes from FGFP supplementary data (Supplementary Table ) , where the continuous ones will have a prevalence of 1
Hughes_dat$ncase.exposure <- NA
Hughes_dat$ncontrol.exposure <- NA
Hughes_dat$prevalence.exposure <- NA

#C_Gammaproteobacteria_RNT
Hughes_dat$ncase.exposure[Hughes_dat$exposure == "C_Gammaproteobacteria_RNT"] <- 1626
Hughes_dat$ncontrol.exposure[Hughes_dat$exposure == "C_Gammaproteobacteria_RNT"] <- 1626
Hughes_dat$samplesize.exposure[Hughes_dat$exposure == "C_Gammaproteobacteria_RNT"] <- 1626
Hughes_dat$prevalence.exposure[Hughes_dat$exposure == "C_Gammaproteobacteria_RNT"] <- 1

#F_Sutterellaceae_HB
Hughes_dat$ncase.exposure[Hughes_dat$exposure == "F_Sutterellaceae_HB"] <- 2223*0.934730056
Hughes_dat$ncontrol.exposure[Hughes_dat$exposure == "F_Sutterellaceae_HB"] <- 2223*(1-0.934730056)
Hughes_dat$samplesize.exposure[Hughes_dat$exposure == "F_Sutterellaceae_HB"] <- 2223
Hughes_dat$prevalence.exposure[Hughes_dat$exposure == "F_Sutterellaceae_HB"] <- 0.934730056

#G_Bifidobacterium_RNT
Hughes_dat$ncase.exposure[Hughes_dat$exposure == "G_Bifidobacterium_RNT"] <- 1975
Hughes_dat$ncontrol.exposure[Hughes_dat$exposure == "G_Bifidobacterium_RNT"] <- 1975
Hughes_dat$samplesize.exposure[Hughes_dat$exposure == "G_Bifidobacterium_RNT"] <- 1975
Hughes_dat$prevalence.exposure[Hughes_dat$exposure == "G_Bifidobacterium_RNT"] <- 1

#G_Butyricicoccus_RNT
Hughes_dat$ncase.exposure[Hughes_dat$exposure == "G_Butyricicoccus_RNT"] <- 2223
Hughes_dat$ncontrol.exposure[Hughes_dat$exposure == "G_Butyricicoccus_RNT"] <- 2223
Hughes_dat$samplesize.exposure[Hughes_dat$exposure == "G_Butyricicoccus_RNT"] <- 2223
Hughes_dat$prevalence.exposure[Hughes_dat$exposure == "G_Butyricicoccus_RNT"] <- 1

#G_Dialister_HB
Hughes_dat$ncase.exposure[Hughes_dat$exposure == "G_Dialister_HB"] <- 2223*0.502820306
Hughes_dat$ncontrol.exposure[Hughes_dat$exposure == "G_Dialister_HB"] <- 2223*(1-0.502820306)
Hughes_dat$samplesize.exposure[Hughes_dat$exposure == "G_Dialister_HB"] <- 2223
Hughes_dat$prevalence.exposure[Hughes_dat$exposure == "G_Dialister_HB"] <- 0.502820306

#G_Parabacteroides_RNT
Hughes_dat$ncase.exposure[Hughes_dat$exposure == "G_Parabacteroides_RNT"] <- 2223
Hughes_dat$ncontrol.exposure[Hughes_dat$exposure == "G_Parabacteroides_RNT"] <- 2223
Hughes_dat$samplesize.exposure[Hughes_dat$exposure == "G_Parabacteroides_RNT"] <- 2223
Hughes_dat$prevalence.exposure[Hughes_dat$exposure == "G_Parabacteroides_RNT"] <- 1

#G_unclassified_F_Erysipelotrichaceae_HB
Hughes_dat$ncase.exposure[Hughes_dat$exposure == "G_unclassified_F_Erysipelotrichaceae_HB"] <- 2223*0.466962127
Hughes_dat$ncontrol.exposure[Hughes_dat$exposure == "G_unclassified_F_Erysipelotrichaceae_HB"] <- 2223*(1-0.466962127)
Hughes_dat$samplesize.exposure[Hughes_dat$exposure == "G_unclassified_F_Erysipelotrichaceae_HB"] <- 2223
Hughes_dat$prevalence.exposure[Hughes_dat$exposure == "G_unclassified_F_Erysipelotrichaceae_HB"] <- 0.466962127

#G_unclassified_F_Porphyromonadaceae_RNT
Hughes_dat$ncase.exposure[Hughes_dat$exposure == "G_unclassified_F_Porphyromonadaceae_RNT"] <- 1420
Hughes_dat$ncontrol.exposure[Hughes_dat$exposure == "G_unclassified_F_Porphyromonadaceae_RNT"] <- 1420
Hughes_dat$samplesize.exposure[Hughes_dat$exposure == "G_unclassified_F_Porphyromonadaceae_RNT"] <- 1420
Hughes_dat$prevalence.exposure[Hughes_dat$exposure == "G_unclassified_F_Porphyromonadaceae_RNT"] <- 1

#G_unclassified_O_Bacteroidales_HB
Hughes_dat$ncase.exposure[Hughes_dat$exposure == "G_unclassified_O_Bacteroidales_HB"] <- 2223*0.18654311
Hughes_dat$ncontrol.exposure[Hughes_dat$exposure == "G_unclassified_O_Bacteroidales_HB"] <- 2223*(1-0.18654311)
Hughes_dat$samplesize.exposure[Hughes_dat$exposure == "G_unclassified_O_Bacteroidales_HB"] <- 2223
Hughes_dat$prevalence.exposure[Hughes_dat$exposure == "G_unclassified_O_Bacteroidales_HB"] <- 0.18654311

#G_unclassified_P_Firmicutes_HB
Hughes_dat$ncase.exposure[Hughes_dat$exposure == "G_unclassified_P_Firmicutes_HB"] <- 2223*0.856970185
Hughes_dat$ncontrol.exposure[Hughes_dat$exposure == "G_unclassified_P_Firmicutes_HB"] <- 2223*(1-0.856970185)
Hughes_dat$samplesize.exposure[Hughes_dat$exposure == "G_unclassified_P_Firmicutes_HB"] <- 2223
Hughes_dat$prevalence.exposure[Hughes_dat$exposure == "G_unclassified_P_Firmicutes_HB"] <- 0.856970185

#G_unclassified_P_Firmicutes_RNT
Hughes_dat$ncase.exposure[Hughes_dat$exposure == "G_unclassified_P_Firmicutes_RNT"] <- 1904
Hughes_dat$ncontrol.exposure[Hughes_dat$exposure == "G_unclassified_P_Firmicutes_RNT"] <- 1904
Hughes_dat$samplesize.exposure[Hughes_dat$exposure == "G_unclassified_P_Firmicutes_RNT"] <- 1904
Hughes_dat$prevalence.exposure[Hughes_dat$exposure == "G_unclassified_P_Firmicutes_RNT"] <- 1

#G_Veillonella_HB
Hughes_dat$ncase.exposure[Hughes_dat$exposure == "G_Veillonella_HB"] <- 2223*0.278404512
Hughes_dat$ncontrol.exposure[Hughes_dat$exposure == "G_Veillonella_HB"] <- 2223*(1-0.278404512)
Hughes_dat$samplesize.exposure[Hughes_dat$exposure == "G_Veillonella_HB"] <- 2223
Hughes_dat$prevalence.exposure[Hughes_dat$exposure == "G_Veillonella_HB"] <- 0.278404512

#G_Coprococcus_HB
Hughes_dat$ncase.exposure[Hughes_dat$exposure == "G_Coprococcus_HB"] <- 2223*0.941176471
Hughes_dat$ncontrol.exposure[Hughes_dat$exposure == "G_Coprococcus_HB"] <- 2223*(1-0.941176471)
Hughes_dat$samplesize.exposure[Hughes_dat$exposure == "G_Coprococcus_HB"] <- 2223
Hughes_dat$prevalence.exposure[Hughes_dat$exposure == "G_Coprococcus_HB"] <- 0.941176471

#G_Ruminococcus_HB
Hughes_dat$ncase.exposure[Hughes_dat$exposure == "G_Ruminococcus_HB"] <- 2223*0.940370669
Hughes_dat$ncontrol.exposure[Hughes_dat$exposure == "G_Ruminococcus_HB"] <- 2223*(1-0.940370669)
Hughes_dat$samplesize.exposure[Hughes_dat$exposure == "G_Ruminococcus_HB"] <- 2223
Hughes_dat$prevalence.exposure[Hughes_dat$exposure == "G_Ruminococcus_HB"] <- 0.940370669


head(Hughes_dat)

###### Need to find the outcome ncase.outcome ncontrol.outcome from dataset

gwasinfo("ebi-a-GCST006085") # 29863 cases 55586 controls
Hughes_dat$ncase.outcome <- 79148
Hughes_dat$ncontrol.outcome <- 61106
Hughes_dat$prevalence.outcome <- 0.1
Hughes_dat$units.outcome <- "log odds"
Hughes_steiger <- steiger_filtering(Hughes_dat)

binary <- Hughes_steiger[grep("_HB", Hughes_steiger$exposure),]
binary$F<- (binary$rsq.exposure*(binary$samplesize.exposure-1-1))/((1-binary$rsq.exposure)*1)
continuous <- Hughes_steiger[grep("_RNT", Hughes_steiger$exposure),] <- Hughes_steiger[grep("_RNT", Hughes_steiger$exposure),]
continuous$F <- continuous$rsq.exposure*(N-1-k)/((1-continuous$rsq.exposure)*k)

## Combine
combined <- rbind(continuous, binary)
combined

## Summarising 
summary(combined$R2)
summary(combined$F)
combined$F_statistic_exposure <- combined$F
combined$samplesize.outcome <- combined$ncase.outcome + combined$ncontrol.outcome
combined$F_statistic_outcome <- (combined$rsq.outcome*(combined$samplesize.outcome-1-1))/((1-combined$rsq.outcome)*1)
#write.table(combined, "./FStats_Hughes.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(combined, "./Output/Hughes_microbiome_prostate_steiger_filtering_F_stat.txt", sep = "\t", col.names = T, row.names = F, quote = F)




#write.table(Hughes_steiger, "./Output/Hughes_microbiome_prostate_steiger_filtering.txt", sep = "\t", col.names = T, row.names = F, quote = F)


####### MiBioGen 
MiBioGen_dat$units.exposure <- "SD"

MiBioGen_dat$ncase.outcome <- 79148
MiBioGen_dat$ncontrol.outcome <- 61106
MiBioGen_dat$prevalence.outcome <- 0.1
MiBioGen_dat$units.outcome <- "log odds"

MiBioGen_steiger <- steiger_filtering(MiBioGen_dat)
MiBioGen_steiger$F_statistic_exposure <- (MiBioGen_steiger$rsq.exposure*(MiBioGen_steiger$samplesize.exposure-1-1))/((1-MiBioGen_steiger$rsq.exposure)*1)
MiBioGen_steiger$samplesize.outcome <- MiBioGen_steiger$ncase.outcome + MiBioGen_steiger$ncontrol.outcome
MiBioGen_steiger$F_statistic_outcome <- (MiBioGen_steiger$rsq.outcome*(MiBioGen_steiger$samplesize.outcome-1-1))/((1-MiBioGen_steiger$rsq.outcome)*1)

write.table(MiBioGen_steiger, "./Output/MiBioGen_microbiome_prostate_cancer_steiger_filtering_F_stat.txt", sep = "\t", col.names = T, row.names = F, quote = F)

####################################
## ADDING FOREST PLOTS OF RESULTS ##
####################################
## Hughes
Hughes_fres <- Hughes_mr_results
Hughes_fres$sig <- Hughes_fres$pval < 0.05
Hughes_fres$category <-"<0.05"
Hughes_fres$or <- exp(Hughes_fres$b)
Hughes_fres$lci <- exp((Hughes_fres$b - (1.96*Hughes_fres$se)))
Hughes_fres$uci <- exp((Hughes_fres$b + (1.96*Hughes_fres$se)))
Hughes_fres$label <- Hughes_fres$exposure
Hughes_fres$observation <- 1:nrow(Hughes_fres) 
head(Hughes_fres)

ggplot(data=Hughes_fres, aes(y=observation, x=or, xmin=lci, xmax=uci, colour=method)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(Hughes_fres), labels=Hughes_fres$label) +
  labs(title="Effect of the gut microbiome on prostate cancer", x="Odds ratio", y = "Microbial trait (Hughes)") +
  geom_vline(xintercept=1, color="black", linetype="dashed", alpha=.5) +
  theme_classic()
ggsave("./Figures/forestplot_microbiome_prostate_cancer_Hughes.png", width=7, height=8)

# MiBioGen
MiBioGen_fres <- MiBioGen_mr_results
MiBioGen_fres$sig <- MiBioGen_fres$pval < 0.05
MiBioGen_fres$category <-"<0.05"
MiBioGen_fres$or <- exp(MiBioGen_fres$b)
MiBioGen_fres$lci <- exp((MiBioGen_fres$b - (1.96*MiBioGen_fres$se)))
MiBioGen_fres$uci <- exp((MiBioGen_fres$b + (1.96*MiBioGen_fres$se)))
MiBioGen_fres$label <- MiBioGen_fres$exposure
MiBioGen_fres$observation <- 1:nrow(MiBioGen_fres) 
head(MiBioGen_fres)

ggplot(data=MiBioGen_fres, aes(y=observation, x=or, xmin=lci, xmax=uci, colour=method)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(MiBioGen_fres), labels=MiBioGen_fres$label) +
  labs(title="Effect of the gut microbiome on prostate cancer", x="Odds ratio", y = "Microbial trait (MiBioGen)") +
  geom_vline(xintercept=1, color="black", linetype="dashed", alpha=.5) +
  theme_classic()
ggsave("./Figures/forestplot_microbiome_prostate_cancer_MiBioGen.png", width=7, height=8)


##########################################
####### GENERATING A VOLCANO PLOTS #######
##########################################
## Hughes
## Only labelling the ones that meet a multiple testing threshold 
ggplot(Hughes_fres, aes(x=or, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(Hughes_fres, !sig)) +
  geom_point(data=subset(Hughes_fres, sig), aes(colour=category, size=pval < (0.05))) +
  #facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(Hughes_fres, pval<(0.05)), aes(label=label, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3, max.overlaps = Inf) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)
  ) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="Effect of the gut microbiome on prostate cancer", y="P value (-log10)") +
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./Figures/volcanoplot_microbiome_prostate_cancer_Hughes.png", width=10, height=8)

## MiBioGen
## Only labelling the ones that meet a multiple testing threshold 
ggplot(MiBioGen_fres, aes(x=or, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(MiBioGen_fres, !sig)) +
  geom_point(data=subset(MiBioGen_fres, sig), aes(colour=category, size=pval < (0.05))) +
  #facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(MiBioGen_fres, pval<(0.05)), aes(label=label, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3, max.overlaps = Inf) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)
  ) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="Effect of the gut microbiome on prostate cancer", y="P value (-log10)") +
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./Figures/volcanoplot_microbiome_prostate_cancer_MiBioGen.png", width=10, height=8)



##################################################### PART 2 - UK Biobank #####################################################


###############################################
## EXTRACTING OUTCOME DATA FROM IEU OPENGWAS ##
###############################################
## Saving outcome data from IEU OpenGWAS 
# ieu-b-4809 = prostate cancer in Uk Biobank 

## Set up vector containing all outcomes
outcomes <- c("ieu-b-4809") #try extract instruments

## Extract the outcomes
Hughes_UKBB_outcome_dat <- extract_outcome_data(Hughes_exposure_dat$SNP, outcomes = outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, splitsize = 10000, proxy_splitsize = 500)
head(Hughes_UKBB_outcome_dat) 
dim(Hughes_UKBB_outcome_dat) # 14 Hughes SNPs in UK Biobank data
MiBioGen_UKBB_outcome_dat <- extract_outcome_data(MiBioGen_exposure_dat$SNP, outcomes = outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, splitsize = 10000, proxy_splitsize = 500)
head(MiBioGen_UKBB_outcome_dat) 
dim(MiBioGen_UKBB_outcome_dat) # 20 MiBioGen SNPs in UK Biobank data

#########################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND METABS (OUTCOME) ##
#########################################################################
## MAIN ANALYSIS INCLUDING THE TOP HITS FROM OUR PUBLISHED META-ANALYSES
Hughes_UKBB_dat <- harmonise_data(Hughes_exposure_dat, Hughes_UKBB_outcome_dat, action =2)
Hughes_UKBB_mr_results <- mr(Hughes_UKBB_dat, method_list=c("mr_wald_ratio"))
write.table(Hughes_UKBB_mr_results, "./Output/Hughes_microbiome_UKBB_results.txt", sep = "\t", col.names = T, row.names = F, quote = F)

MiBioGen_UKBB_dat <- harmonise_data(MiBioGen_exposure_dat, MiBioGen_UKBB_outcome_dat, action =1)
MiBioGen_UKBB_mr_results <- mr(MiBioGen_UKBB_dat, method_list=c("mr_wald_ratio", "mr_ivw"))
write.table(MiBioGen_UKBB_mr_results, "./Output/MiBioGen_microbiome_UKBB_results.txt", sep = "\t", col.names = T, row.names = F, quote = F)


#######################################
## steiger filtering 
######################################

#######################################################
## Creating exposure phenotype
#######################################################
## Read in the file that has all the instrument information (from David Hughes paper)
#Hughes_exposures <- read.table("./Bug2DiseaseMRformat.txt", header = T, sep = "\t")
#Hughes_exposures$MarkerName <- substr(Hughes_exposures$snpid, 1, nchar(Hughes_exposures$snpid)-4)
#Hughes_exposure_dat <- format_data(Hughes_exposures, phenotype_col = "mt", snp_col = "rsid", beta_col = "beta", se_col = "se", pval_col = "P_value", eaf_col = "eaf", effect_allele_col = "allele_B", other_allele_col = "allele_A", samplesize_col = "n")
#head(Hughes_exposure_dat)

## Add the prevalence to the binary phenotypes from FGFP supplementary data (Supplementary Table ) , where the continuous ones will have a prevalence of 1
Hughes_UKBB_dat$ncase.exposure <- NA
Hughes_UKBB_dat$ncontrol.exposure <- NA
Hughes_UKBB_dat$prevalence.exposure <- NA

#C_Gammaproteobacteria_RNT
Hughes_UKBB_dat$ncase.exposure[Hughes_UKBB_dat$exposure == "C_Gammaproteobacteria_RNT"] <- 1626
Hughes_UKBB_dat$ncontrol.exposure[Hughes_UKBB_dat$exposure == "C_Gammaproteobacteria_RNT"] <- 1626
Hughes_UKBB_dat$samplesize.exposure[Hughes_UKBB_dat$exposure == "C_Gammaproteobacteria_RNT"] <- 1626
Hughes_UKBB_dat$prevalence.exposure[Hughes_UKBB_dat$exposure == "C_Gammaproteobacteria_RNT"] <- 1

#F_Sutterellaceae_HB
Hughes_UKBB_dat$ncase.exposure[Hughes_UKBB_dat$exposure == "F_Sutterellaceae_HB"] <- 2223*0.934730056
Hughes_UKBB_dat$ncontrol.exposure[Hughes_UKBB_dat$exposure == "F_Sutterellaceae_HB"] <- 2223*(1-0.934730056)
Hughes_UKBB_dat$samplesize.exposure[Hughes_UKBB_dat$exposure == "F_Sutterellaceae_HB"] <- 2223
Hughes_UKBB_dat$prevalence.exposure[Hughes_UKBB_dat$exposure == "F_Sutterellaceae_HB"] <- 0.934730056

#G_Bifidobacterium_RNT
Hughes_UKBB_dat$ncase.exposure[Hughes_UKBB_dat$exposure == "G_Bifidobacterium_RNT"] <- 1975
Hughes_UKBB_dat$ncontrol.exposure[Hughes_UKBB_dat$exposure == "G_Bifidobacterium_RNT"] <- 1975
Hughes_UKBB_dat$samplesize.exposure[Hughes_UKBB_dat$exposure == "G_Bifidobacterium_RNT"] <- 1975
Hughes_UKBB_dat$prevalence.exposure[Hughes_UKBB_dat$exposure == "G_Bifidobacterium_RNT"] <- 1

#G_Butyricicoccus_RNT
Hughes_UKBB_dat$ncase.exposure[Hughes_UKBB_dat$exposure == "G_Butyricicoccus_RNT"] <- 2223
Hughes_UKBB_dat$ncontrol.exposure[Hughes_UKBB_dat$exposure == "G_Butyricicoccus_RNT"] <- 2223
Hughes_UKBB_dat$samplesize.exposure[Hughes_UKBB_dat$exposure == "G_Butyricicoccus_RNT"] <- 2223
Hughes_UKBB_dat$prevalence.exposure[Hughes_UKBB_dat$exposure == "G_Butyricicoccus_RNT"] <- 1

#G_Dialister_HB
Hughes_UKBB_dat$ncase.exposure[Hughes_UKBB_dat$exposure == "G_Dialister_HB"] <- 2223*0.502820306
Hughes_UKBB_dat$ncontrol.exposure[Hughes_UKBB_dat$exposure == "G_Dialister_HB"] <- 2223*(1-0.502820306)
Hughes_UKBB_dat$samplesize.exposure[Hughes_UKBB_dat$exposure == "G_Dialister_HB"] <- 2223
Hughes_UKBB_dat$prevalence.exposure[Hughes_UKBB_dat$exposure == "G_Dialister_HB"] <- 0.502820306

#G_Parabacteroides_RNT
Hughes_UKBB_dat$ncase.exposure[Hughes_UKBB_dat$exposure == "G_Parabacteroides_RNT"] <- 2223
Hughes_UKBB_dat$ncontrol.exposure[Hughes_UKBB_dat$exposure == "G_Parabacteroides_RNT"] <- 2223
Hughes_UKBB_dat$samplesize.exposure[Hughes_UKBB_dat$exposure == "G_Parabacteroides_RNT"] <- 2223
Hughes_UKBB_dat$prevalence.exposure[Hughes_UKBB_dat$exposure == "G_Parabacteroides_RNT"] <- 1

#G_unclassified_F_Erysipelotrichaceae_HB
Hughes_UKBB_dat$ncase.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_F_Erysipelotrichaceae_HB"] <- 2223*0.466962127
Hughes_UKBB_dat$ncontrol.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_F_Erysipelotrichaceae_HB"] <- 2223*(1-0.466962127)
Hughes_UKBB_dat$samplesize.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_F_Erysipelotrichaceae_HB"] <- 2223
Hughes_UKBB_dat$prevalence.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_F_Erysipelotrichaceae_HB"] <- 0.466962127

#G_unclassified_F_Porphyromonadaceae_RNT
Hughes_UKBB_dat$ncase.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_F_Porphyromonadaceae_RNT"] <- 1420
Hughes_UKBB_dat$ncontrol.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_F_Porphyromonadaceae_RNT"] <- 1420
Hughes_UKBB_dat$samplesize.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_F_Porphyromonadaceae_RNT"] <- 1420
Hughes_UKBB_dat$prevalence.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_F_Porphyromonadaceae_RNT"] <- 1

#G_unclassified_O_Bacteroidales_HB
Hughes_UKBB_dat$ncase.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_O_Bacteroidales_HB"] <- 2223*0.18654311
Hughes_UKBB_dat$ncontrol.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_O_Bacteroidales_HB"] <- 2223*(1-0.18654311)
Hughes_UKBB_dat$samplesize.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_O_Bacteroidales_HB"] <- 2223
Hughes_UKBB_dat$prevalence.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_O_Bacteroidales_HB"] <- 0.18654311

#G_unclassified_P_Firmicutes_HB
Hughes_UKBB_dat$ncase.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_P_Firmicutes_HB"] <- 2223*0.856970185
Hughes_UKBB_dat$ncontrol.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_P_Firmicutes_HB"] <- 2223*(1-0.856970185)
Hughes_UKBB_dat$samplesize.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_P_Firmicutes_HB"] <- 2223
Hughes_UKBB_dat$prevalence.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_P_Firmicutes_HB"] <- 0.856970185

#G_unclassified_P_Firmicutes_RNT
Hughes_UKBB_dat$ncase.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_P_Firmicutes_RNT"] <- 1904
Hughes_UKBB_dat$ncontrol.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_P_Firmicutes_RNT"] <- 1904
Hughes_UKBB_dat$samplesize.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_P_Firmicutes_RNT"] <- 1904
Hughes_UKBB_dat$prevalence.exposure[Hughes_UKBB_dat$exposure == "G_unclassified_P_Firmicutes_RNT"] <- 1

#G_Veillonella_HB
Hughes_UKBB_dat$ncase.exposure[Hughes_UKBB_dat$exposure == "G_Veillonella_HB"] <- 2223*0.278404512
Hughes_UKBB_dat$ncontrol.exposure[Hughes_UKBB_dat$exposure == "G_Veillonella_HB"] <- 2223*(1-0.278404512)
Hughes_UKBB_dat$samplesize.exposure[Hughes_UKBB_dat$exposure == "G_Veillonella_HB"] <- 2223
Hughes_UKBB_dat$prevalence.exposure[Hughes_UKBB_dat$exposure == "G_Veillonella_HB"] <- 0.278404512

#G_Coprococcus_HB
Hughes_UKBB_dat$ncase.exposure[Hughes_UKBB_dat$exposure == "G_Coprococcus_HB"] <- 2223*0.941176471
Hughes_UKBB_dat$ncontrol.exposure[Hughes_UKBB_dat$exposure == "G_Coprococcus_HB"] <- 2223*(1-0.941176471)
Hughes_UKBB_dat$samplesize.exposure[Hughes_UKBB_dat$exposure == "G_Coprococcus_HB"] <- 2223
Hughes_UKBB_dat$prevalence.exposure[Hughes_UKBB_dat$exposure == "G_Coprococcus_HB"] <- 0.941176471

#G_Ruminococcus_HB
Hughes_UKBB_dat$ncase.exposure[Hughes_UKBB_dat$exposure == "G_Ruminococcus_HB"] <- 2223*0.940370669
Hughes_UKBB_dat$ncontrol.exposure[Hughes_UKBB_dat$exposure == "G_Ruminococcus_HB"] <- 2223*(1-0.940370669)
Hughes_UKBB_dat$samplesize.exposure[Hughes_UKBB_dat$exposure == "G_Ruminococcus_HB"] <- 2223
Hughes_UKBB_dat$prevalence.exposure[Hughes_UKBB_dat$exposure == "G_Ruminococcus_HB"] <- 0.940370669


head(Hughes_UKBB_dat)

###### Need to find the outcome ncase.outcome ncontrol.outcome from dataset

gwasinfo("ieu-b-4809") # 2671 cases 372016 controls
Hughes_UKBB_dat$ncase.outcome <- 9132
Hughes_UKBB_dat$ncontrol.outcome <- 173493
Hughes_UKBB_dat$prevalence.outcome <- 0.1
Hughes_UKBB_dat$units.outcome <- "log odds"
steiger_filtering(Hughes_UKBB_dat)
write.table(steiger_filtering(Hughes_UKBB_dat), "./Output/Hughes_microbiome_UKBB_steiger_filtering.txt", sep = "\t", col.names = T, row.names = F, quote = F)


####### MiBioGen 
MiBioGen_UKBB_dat$units.exposure <- "SD"

MiBioGen_UKBB_dat$ncase.outcome <- 9132
MiBioGen_UKBB_dat$ncontrol.outcome <- 173493
MiBioGen_UKBB_dat$prevalence.outcome <- 0.1
MiBioGen_UKBB_dat$units.outcome <- "log odds"

steiger_filtering(MiBioGen_UKBB_dat)
write.table(steiger_filtering(MiBioGen_UKBB_dat), "./Output/MiBioGen_microbiome_UKBB_steiger_filtering.txt", sep = "\t", col.names = T, row.names = F, quote = F)

####################################
## ADDING FOREST PLOTS OF RESULTS ##
####################################
## Hughes
Hughes_UKBB_fres <- Hughes_UKBB_mr_results
Hughes_UKBB_fres$sig <- Hughes_UKBB_fres$pval < 0.05
Hughes_UKBB_fres$category <-"<0.05"
Hughes_UKBB_fres$or <- exp(Hughes_UKBB_fres$b)
Hughes_UKBB_fres$lci <- exp((Hughes_UKBB_fres$b - (1.96*Hughes_UKBB_fres$se)))
Hughes_UKBB_fres$uci <- exp((Hughes_UKBB_fres$b + (1.96*Hughes_UKBB_fres$se)))
Hughes_UKBB_fres$outcome[Hughes_UKBB_fres$outcome=="prostate cancer || id:ieu-b-4809"] <- "prostate cancer"
Hughes_UKBB_fres$label <- Hughes_UKBB_fres$exposure
Hughes_UKBB_fres$observation <- 1:nrow(Hughes_UKBB_fres) 
head(Hughes_UKBB_fres)

ggplot(data=Hughes_UKBB_fres, aes(y=observation, x=or, xmin=lci, xmax=uci, colour=outcome)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(Hughes_UKBB_fres), labels=Hughes_UKBB_fres$label) +
  labs(title="Effect of the gut microbiome on prostate cancer", x="Odds ratio", y = "Microbial trait (Hughes)") +
  geom_vline(xintercept=1, color="black", linetype="dashed", alpha=.5) +
  theme_classic()
ggsave("./Figures/forestplot_microbiome_UKBB_Hughes.png", width=7, height=8)

# MiBioGen
MiBioGen_UKBB_fres <- MiBioGen_UKBB_mr_results
MiBioGen_UKBB_fres$sig <- MiBioGen_UKBB_fres$pval < 0.05
MiBioGen_UKBB_fres$category <-"<0.05"
MiBioGen_UKBB_fres$or <- exp(MiBioGen_UKBB_fres$b)
MiBioGen_UKBB_fres$lci <- exp((MiBioGen_UKBB_fres$b - (1.96*MiBioGen_UKBB_fres$se)))
MiBioGen_UKBB_fres$uci <- exp((MiBioGen_UKBB_fres$b + (1.96*MiBioGen_UKBB_fres$se)))
MiBioGen_UKBB_fres$outcome[MiBioGen_UKBB_fres$outcome=="prostate cancer || id:ieu-b-4809"] <- "prostate cancer"
MiBioGen_UKBB_fres$label <- MiBioGen_UKBB_fres$exposure
MiBioGen_UKBB_fres$observation <- 1:nrow(MiBioGen_UKBB_fres) 
head(MiBioGen_UKBB_fres)

ggplot(data=MiBioGen_UKBB_fres, aes(y=observation, x=or, xmin=lci, xmax=uci, colour=outcome)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(MiBioGen_UKBB_fres), labels=MiBioGen_UKBB_fres$label) +
  labs(title="Effect of the gut microbiome on prostate cancer", x="Odds ratio", y = "Microbial trait (MiBioGen)") +
  geom_vline(xintercept=1, color="black", linetype="dashed", alpha=.5) +
  theme_classic()
ggsave("./Figures/forestplot_microbiome_UKBB_MiBioGen.png", width=7, height=8)


##########################################
####### GENERATING A VOLCANO PLOTS #######
##########################################
## Hughes
## only labelling the ones that meet a multiple testing threshold 
ggplot(Hughes_UKBB_fres, aes(x=or, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(Hughes_UKBB_fres, !sig)) +
  geom_point(data=subset(Hughes_UKBB_fres, sig), aes(colour=category, size=pval < (0.05))) +
  #facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(Hughes_UKBB_fres, pval<(0.05)), aes(label=label, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3, max.overlaps = Inf) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="Effect of the gut microbiome on prostate cancer", y="P value (-log10)") +
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./Figures/volcanoplot_microbiome_UKBB_prostate_Hughes.png", width=10, height=8)

## MiBioGen
## only labelling the ones that meet a multiple testing threshold 
ggplot(MiBioGen_UKBB_fres, aes(x=or, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(MiBioGen_UKBB_fres, !sig)) +
  geom_point(data=subset(MiBioGen_UKBB_fres, sig), aes(colour=category, size=pval < (0.05))) +
  #facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(MiBioGen_UKBB_fres, pval<(0.05)), aes(label=label, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3, max.overlaps = Inf) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="Effect of the gut microbiome on prostate cancer", y="P value (-log10)") +
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./Figures/volcanoplot_microbiome_UKBB_prostate_MiBioGen.png", width=10, height=8)



######################. Part 3 - Finngen #######################################

## Saving outcome data from IEU OpenGWAS 

outcomes <- c("finn-b-C3_PROSTATE_EXALLC")
## Extract the outcomes
Hughes_finngen_outcome_dat <- extract_outcome_data(Hughes_exposure_dat$SNP, outcomes = outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, splitsize = 10000, proxy_splitsize = 500)
head(Hughes_finngen_outcome_dat) 
MiBioGen_finngen_outcome_dat <- extract_outcome_data(MiBioGen_exposure_dat$SNP, outcomes = outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, splitsize = 10000, proxy_splitsize = 500)
head(MiBioGen_finngen_outcome_dat) 

MiBioGen_finngen_outcome_dat = unique(MiBioGen_finngen_outcome_dat)
#########################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND METABS (OUTCOME) ##
#########################################################################
## MAIN ANALYSIS INCLUDING THE TOP HITS FROM OUR PUBLISHED META-ANALYSES
Hughes_finngen_dat <- harmonise_data(Hughes_exposure_dat, Hughes_finngen_outcome_dat, action =2)
Hughes_finngen_dat <- unique(Hughes_finngen_dat)
Hughes_finngen_mr_results <- mr(Hughes_finngen_dat, method_list=c("mr_wald_ratio"))
write.table(Hughes_finngen_mr_results, "./Output/Hughes_microbiome_FinnGen_results.txt", sep = "\t", col.names = T, row.names = F, quote = F) # 14 mt associations

MiBioGen_finngen_dat <- harmonise_data(MiBioGen_exposure_dat, MiBioGen_finngen_outcome_dat, action =1)
MiBioGen_finngen_mr_results <- mr(MiBioGen_finngen_dat, method_list=c("mr_wald_ratio", "mr_ivw"))
write.table(MiBioGen_finngen_mr_results, "./Output/MiBioGen_microbiome_FinnGen_results.txt", sep = "\t", col.names = T, row.names = F, quote = F)


#######################################
## steiger filtering 
######################################

#######################################################
## Creating exposure phenotype
#######################################################

## Add the prevalence to the binary phenotypes from FGFP supplementary data (Supplementary Table ) , where the continuous ones will have a prevalence of 1
Hughes_finngen_dat$ncase.exposure <- Hughes_finngen_dat$Ncase.exp
Hughes_finngen_dat$ncontrol.exposure <- Hughes_finngen_dat$Ncontrol.exp
Hughes_finngen_dat$prevalence.exposure <- Hughes_finngen_dat$prevalence
head(Hughes_finngen_dat)

###### Need to find the outcome ncase.outcome ncontrol.outcome from dataset

gwasinfo("finn-b-C3_PROSTATE_EXALLC") # 1681 cases 173993 controls
Hughes_finngen_dat$ncase.outcome <- 6311
Hughes_finngen_dat$ncontrol.outcome <- 74685
Hughes_finngen_dat$prevalence.outcome <- 0.1
Hughes_finngen_dat$units.outcome <- "log odds"
steiger_filtering(Hughes_finngen_dat)
write.table(steiger_filtering(Hughes_finngen_dat), "./Output/Hughes_microbiome_finngen_steiger_filtering.txt", sep = "\t", col.names = T, row.names = F, quote = F)


####### MiBioGen 
MiBioGen_finngen_dat$units.exposure <- "SD"

MiBioGen_finngen_dat$ncase.outcome <- 6311
MiBioGen_finngen_dat$ncontrol.outcome <- 74685
MiBioGen_finngen_dat$prevalence.outcome <- 0.1
MiBioGen_finngen_dat$units.outcome <- "log odds"

MiBioGen_finngen_steiger <- steiger_filtering(MiBioGen_finngen_dat)
write.table(MiBioGen_finngen_steiger, "./Output/MiBioGen_microbiome_finngen_steiger_filtering.txt", sep = "\t", col.names = T, row.names = F, quote = F)

####################################
## ADDING FOREST PLOTS OF RESULTS ##
####################################
## Hughes
Hughes_finngen_fres <- Hughes_finngen_mr_results
Hughes_finngen_fres$sig <- Hughes_finngen_fres$pval < 0.05
Hughes_finngen_fres$category <-"<0.05"
Hughes_finngen_fres$or <- exp(Hughes_finngen_fres$b)
Hughes_finngen_fres$lci <- exp((Hughes_finngen_fres$b - (1.96*Hughes_finngen_fres$se)))
Hughes_finngen_fres$uci <- exp((Hughes_finngen_fres$b + (1.96*Hughes_finngen_fres$se)))
Hughes_finngen_fres$outcome[Hughes_finngen_fres$outcome=="Malignant neoplasm of prostate (all cancers excluded) || id:finn-b-C3_PROSTATE_EXALLC"] <- "Malignant neoplasm"
Hughes_finngen_fres$label <- Hughes_finngen_fres$exposure
Hughes_finngen_fres$observation <- 1:nrow(Hughes_finngen_fres) 
head(Hughes_finngen_fres)

ggplot(data=Hughes_finngen_fres, aes(y=observation, x=or, xmin=lci, xmax=uci, colour=outcome)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(Hughes_finngen_fres), labels=Hughes_finngen_fres$label) +
  labs(title="Effect of the gut microbiome on prostate cancer", x="Odds ratio", y = "Microbial trait (Hughes)") +
  geom_vline(xintercept=1, color="black", linetype="dashed", alpha=.5) +
  theme_classic()
ggsave("./Figures/forestplot_microbiome_FinnGen_Hughes.png", width=7, height=8)

# MiBioGen
MiBioGen_finngen_fres <- MiBioGen_finngen_mr_results
MiBioGen_finngen_fres$sig <- MiBioGen_finngen_fres$pval < 0.05 # need to change these to correct number
MiBioGen_finngen_fres$category <-"<0.05"
MiBioGen_finngen_fres$or <- exp(MiBioGen_finngen_fres$b)
MiBioGen_finngen_fres$lci <- exp((MiBioGen_finngen_fres$b - (1.96*MiBioGen_finngen_fres$se)))
MiBioGen_finngen_fres$uci <- exp((MiBioGen_finngen_fres$b + (1.96*MiBioGen_finngen_fres$se)))
MiBioGen_finngen_fres$outcome[MiBioGen_finngen_fres$outcome=="Malignant neoplasm of prostate (all cancers excluded) || id:finn-b-C3_PROSTATE_EXALLC"] <- "Malignant neoplasm"
MiBioGen_finngen_fres$label <- MiBioGen_finngen_fres$exposure
MiBioGen_finngen_fres$observation <- 1:nrow(MiBioGen_finngen_fres) 
head(MiBioGen_finngen_fres)

ggplot(data=MiBioGen_finngen_fres, aes(y=observation, x=or, xmin=lci, xmax=uci, colour=outcome)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(MiBioGen_finngen_fres), labels=MiBioGen_finngen_fres$label) +
  labs(title="Effect of the gut microbiome on prostate cancer", x="Odds ratio", y = "Microbial trait (MiBioGen)") +
  geom_vline(xintercept=1, color="black", linetype="dashed", alpha=.5) +
  theme_classic()
ggsave("./Figures/forestplot_microbiome_FinnGen_MiBioGen.png", width=7, height=8)


##########################################
####### GENERATING A VOLCANO PLOTS #######
##########################################
## Hughes
## only labelling the ones that meet a multiple testing threshold 
ggplot(Hughes_finngen_fres, aes(x=or, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(Hughes_finngen_fres, !sig)) +
  geom_point(data=subset(Hughes_finngen_fres, sig), aes(colour=category, size=pval < (0.05))) +
  #facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(Hughes_finngen_fres, pval<(0.05)), aes(label=label, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3, max.overlaps = Inf) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="Effect of the gut microbiome on prostate cancer", y="P value (-log10)") +
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./Figures/volcanoplot_microbiome_finngen_prostate_Hughes.png", width=10, height=8)

## MiBioGen
## only labelling the ones that meet a multiple testing threshold 
ggplot(MiBioGen_finngen_fres, aes(x=or, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(MiBioGen_finngen_fres, !sig)) +
  geom_point(data=subset(MiBioGen_finngen_fres, sig), aes(colour=category, size=pval < (0.05))) +
  #facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(MiBioGen_finngen_fres, pval<(0.05)), aes(label=label, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3, max.overlaps = Inf) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="Effect of the gut microbiome on prostate cancer", y="P value (-log10)") +
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./Figures/volcanoplot_microbiome_finngen_prostate_MiBioGen.png", width=10, height=8)



##################################################### META-ANALYSING Prostate cancer PRACTICAL, UK Biobank and Finngen MR results ############################################
## Settings
colours <- names(wes_palettes)
discrete_palette <- wes_palette(colours[8], type = "discrete")
psignif <- 0
ci <- 0.95


####################################
## ADDING FOREST PLOTS OF RESULTS ##
####################################
## Set main and sensitivity dataframes
## Hughes
Hughes_main <- Hughes_fres[,c(4,3,6,5,7,8,9)]
colnames(Hughes_main) <- c("exposure", "outcome", "nsnp", "method", "beta", "se", "p")
Hughes_main$exposure <- gsub("RNT", "AB", Hughes_main$exposure)
Hughes_main$exposure <- gsub("HB", "P/A", Hughes_main$exposure)
Hughes_main$dataset <- "PRACTICAL"

Hughes_UKBB_sens <- Hughes_UKBB_fres[,c(4,3,6,5,7,8,9)]
colnames(Hughes_UKBB_sens) <- c("exposure", "outcome", "nsnp", "method", "beta", "se", "p")
Hughes_UKBB_sens$exposure <- gsub("RNT", "AB", Hughes_UKBB_sens$exposure)
Hughes_UKBB_sens$exposure <- gsub("HB", "P/A", Hughes_UKBB_sens$exposure)
Hughes_UKBB_sens$dataset <- "UK Biobank"


Hughes_finngen_sens <- Hughes_finngen_fres[,c(4,3,6,5,7,8,9)]
colnames(Hughes_finngen_sens) <- c("exposure", "outcome", "nsnp", "method", "beta", "se", "p")
Hughes_finngen_sens$exposure <- gsub("RNT", "AB", Hughes_finngen_sens$exposure)
Hughes_finngen_sens$exposure <- gsub("HB", "P/A", Hughes_finngen_sens$exposure)
Hughes_finngen_sens$dataset <- "FinnGen"


## MiBioGen
MiBioGen_main <- MiBioGen_fres[,c(4,3,6,5,7,8,9)]
colnames(MiBioGen_main) <- c("exposure", "outcome", "nsnp", "method", "beta", "se", "p")
MiBioGen_main$dataset <- "PRACTICAL"

MiBioGen_UKBB_sens <- MiBioGen_UKBB_fres[,c(4,3,6,5,7,8,9)]
colnames(MiBioGen_UKBB_sens) <- c("exposure", "outcome", "nsnp", "method", "beta", "se", "p")
MiBioGen_UKBB_sens$dataset <- "UK Biobank"


MiBioGen_finngen_sens <- MiBioGen_finngen_fres[,c(4,3,6,5,7,8,9)]
colnames(MiBioGen_finngen_sens) <- c("exposure", "outcome", "nsnp", "method", "beta", "se", "p")
MiBioGen_finngen_sens$dataset <- "FinnGen"


## Bind together
Hughes_plot_data <- rbind(Hughes_main, Hughes_UKBB_sens, Hughes_finngen_sens)
MiBioGen_plot_data <- rbind(MiBioGen_main, MiBioGen_UKBB_sens, MiBioGen_finngen_sens)

## Now limit to the ones that are consistent across both datasets
Hughes_plot_data <- Hughes_plot_data[ave(1:nrow(Hughes_plot_data), Hughes_plot_data$exposure, FUN=length)>1,]
MiBioGen_plot_data <- MiBioGen_plot_data[ave(1:nrow(MiBioGen_plot_data), MiBioGen_plot_data$exposure, FUN=length)>1,]

## Plot each comparing main and sensitivity analyses (TRICL vs. UKBiobank vs. FinnGen)
## Hughes
ggforestplot::forestplot(df = Hughes_plot_data,
                         name = exposure,
                         estimate = beta,
                         pvalue = p,
                         psignif = psignif,
                         ci = ci,
                         se = se,
                         logodds = T,
                         colour=dataset) +
  ggtitle("Risk of prostate cancer with gut microbiome variation") +
  scale_color_manual(values = c(discrete_palette[1], discrete_palette[2], discrete_palette[3], discrete_palette[4], discrete_palette[5])) +
  xlab("Odds Ratio") +
  theme(legend.title = element_blank())
ggsave("./Figures/Hughes_microbiome_PRACTICAL_UKBB_finngen.png", width=10, height=8)

## MiBioGen
ggforestplot::forestplot(df = MiBioGen_plot_data,
                         name = exposure,
                         estimate = beta,
                         pvalue = p,
                         psignif = psignif,
                         ci = ci,
                         se = se,
                         logodds = T,
                         colour=dataset) +
  ggtitle("Risk of prostate cancer with gut microbiome variation") +
  scale_color_manual(values = c(discrete_palette[1], discrete_palette[2], discrete_palette[3], discrete_palette[4], discrete_palette[5])) +
  xlab("Odds Ratio") +
  theme(legend.title = element_blank())
ggsave("./Figures/MiBioGen_microbiome_PRACTICAL_UKBB_finngen.png", width=10, height=8)

## Simple meta-analysis with Hughes
Hughes_meta <- metagen(Hughes_plot_data$beta, seTE = Hughes_plot_data$se, studlab = Hughes_plot_data$dataset, data = Hughes_plot_data, sm = "OR", fixed = TRUE, random = FALSE)
Hughes_meta_by_exposure <- update(Hughes_meta, subgroup = Hughes_plot_data$exposure, tau.common = FALSE)
summary(Hughes_meta_by_exposure)
pdf("./Figures/Hughes_microbiome_PRACTICAL+UKBB+finngen_meta.pdf", width=10, height=16)
forest(Hughes_meta_by_exposure, studlab = T, common = TRUE, random = FALSE, leftcols = c("dataset"), print.subgroup.name = FALSE, print.subgroup.labels = TRUE, test.subgroup= FALSE, overall = FALSE, overall.hetstat = FALSE, text.common = "Meta-analysis", prediction = FALSE, xlab = "OR of prostate cancer with gut microbiome variation", print.tau2 = FALSE, plotwidth = "15cm", fontsize = 10, spacing = .75)
dev.off()

## Simple meta-analysis with MiBioGen
MiBioGen_meta <- metagen(MiBioGen_plot_data$beta, seTE = MiBioGen_plot_data$se, studlab = MiBioGen_plot_data$dataset, data = MiBioGen_plot_data, sm = "OR", fixed = TRUE, random = FALSE)
MiBioGen_meta_by_exposure <- update(MiBioGen_meta, subgroup = MiBioGen_plot_data$exposure, tau.common = FALSE)
summary(MiBioGen_meta_by_exposure)
pdf("./Figures/MiBioGen_microbiome_PRACTICAL+UKBB+finngen_meta.pdf", width=8, height = 16)
forest(MiBioGen_meta_by_exposure, studlab = T, common = TRUE, random = FALSE, leftcols = c("dataset"), print.subgroup.name = FALSE, print.subgroup.labels = TRUE, test.subgroup= FALSE, overall = FALSE, overall.hetstat = FALSE, text.common = "Meta-analysis", prediction = FALSE, xlab = "OR of prostate cancer with gut microbiome variation", print.tau2 = FALSE, plotwidth = "15cm", fontsize = 6, spacing = .36)
dev.off()

#### only the exposures to follow up in Hughes 
exposures_to_select <- c("G_Parabacteroides_AB", "G_unclassified_P_Firmicutes_AB", "G_unclassified_F_Porphyromonadaceae_AB")
traits_to_follow_up <- Hughes_plot_data[Hughes_plot_data$exposure %in% exposures_to_select, ]
## Simple meta-analysis with MiBioGen
Hughes_meta <- metagen(traits_to_follow_up$beta, seTE = traits_to_follow_up$se, studlab = traits_to_follow_up$dataset, data = traits_to_follow_up, sm = "OR", fixed = TRUE, random = FALSE)
Hughes_meta_by_exposure <- update(Hughes_meta, subgroup = traits_to_follow_up$exposure, tau.common = FALSE)
summary(Hughes_meta_by_exposure)
pdf("./Figures/Hughes_top_MR_traits_microbiome_PRACTICAL+UKBB+finngen_meta.pdf", width=8, height = 3)
forest(Hughes_meta_by_exposure, studlab = T, common = TRUE, random = FALSE, leftcols = c("dataset"), print.subgroup.name = FALSE, print.subgroup.labels = TRUE, test.subgroup= FALSE, overall = FALSE, overall.hetstat = FALSE, text.common = "Meta-analysis", prediction = FALSE, xlab = "OR of prostate cancer with gut microbiome variation", print.tau2 = FALSE, plotwidth = "15cm", fontsize = 6, spacing = .36)
dev.off()


#### only the exposures to follow up in MiBioGen 
exposures_to_select <- c("genus.Allisonella.id.2174")
traits_to_follow_up <- MiBioGen_plot_data[MiBioGen_plot_data$exposure %in% exposures_to_select, ]
## Simple meta-analysis with MiBioGen
MiBioGen_meta <- metagen(traits_to_follow_up$beta, seTE = traits_to_follow_up$se, studlab = traits_to_follow_up$dataset, data = traits_to_follow_up, sm = "OR", fixed = TRUE, random = FALSE)
MiBioGen_meta_by_exposure <- update(MiBioGen_meta, subgroup = traits_to_follow_up$exposure, tau.common = FALSE)
summary(MiBioGen_meta_by_exposure)
pdf("./Figures/MiBioGen_top_MR_traits_microbiome_PRACTICAL+UKBB+finngen_meta.pdf", width=8, height = 3)
forest(MiBioGen_meta_by_exposure, studlab = T, common = TRUE, random = FALSE, leftcols = c("dataset"), print.subgroup.name = FALSE, print.subgroup.labels = TRUE, test.subgroup= FALSE, overall = FALSE, overall.hetstat = FALSE, text.common = "Meta-analysis", prediction = FALSE, xlab = "OR of prostate cancer with gut microbiome variation", print.tau2 = FALSE, plotwidth = "15cm", fontsize = 6, spacing = .36)
dev.off()



