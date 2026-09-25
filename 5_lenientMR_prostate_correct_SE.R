######################################################################################################################
######## PROJECT: MR analysis of the gut microbiome on lung cancer using TRICL, UK Biobank and FinnGen
######## Script: Generate a list of the SNPs associated with the gut microbiome at a more lenient p-value threshold (and directionally consistent)
######## Date: 18/07/24
######################################################################################################################
## Clear space
rm(list=ls())

## Set working directory
setwd("/Users/xh18454/Library/CloudStorage/OneDrive-UniversityofBristol/MR/prostate_cancer/")

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
#install.packages("calibrate")
#install.packages("knitr")
#install.packages("patchwork")
#install.packages("purrr")
#install.packages("ggforestplot")
#install.packages("cowplot")
#install.packages("wesanderson")
#install.packages("meta")
#install.packages("rmeta")
#install.packages("meta")
#install.packages("rmeta")
install.packages("tidyverse")
library(ieugwasr)
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
library(knitr)
library(patchwork)
library(purrr)
library(ggforestplot)
library(cowplot)
library(wesanderson)
library(meta)
library(rmeta)
library(tidyverse)

# Note only need to do for the microbiome traits we want to follow up - G_Parabacteroides_RNT", "family.Streptococcaceae.id.1850", "genus.Intestinibacter.id.11345", "genus.Streptococcus.id.1853"
######################################################################  HUGHES ET AL ######################################################################


################################################
## PREPARING EXPOSURE DATA FROM HUGHES ET AL. ##
################################################
## The following spreadsheet is the Supplementary table document published with the Hughes et al. GWAS
mydata <- readxl::read_excel("./41564_2020_743_MOESM3_ESM.xlsx", sheet = 6 , skip = 4)
head(mydata)
colnames(mydata)

## Identify the tag SNPS
tag_snps <- mydata %>% filter(TagLocus == 1)
dim(tag_snps) # should be 1056

## Identify those that are directionally consistent in all three studies (note: there MUST be an estimate in each of the three studies)
## an estimate can be absent because of MAF and INFO filtering in a particular data set
all_positive <- which(tag_snps$fgfp_beta > 0 & tag_snps$focus_beta > 0 & tag_snps$popgen_beta > 0)
all_negative <- which(tag_snps$fgfp_beta < 0 & tag_snps$focus_beta < 0 & tag_snps$popgen_beta < 0)

## Define a sorted index
directionally_consistent <- sort(c(all_positive, all_negative))

## Define a new data frame and present what the first 20 look like
tag_snps_directionally_consistent <- tag_snps[directionally_consistent, ]
dim(tag_snps_directionally_consistent) # should have 584
tag_snps_directionally_consistent[1:20, c("fgfp_beta","focus_beta","popgen_beta")] %>% 
  knitr::kable() %>% kableExtra::kable_classic()
write.table(tag_snps_directionally_consistent, "./Data/Hughes_tag_snps_directionally_consistent.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## Now identify those that are directionally consistent that also have a p-value less than 1e-05
all_positive_lenient <- which(tag_snps$fgfp_beta > 0 & tag_snps$focus_beta > 0 & tag_snps$popgen_beta > 0 & tag_snps$meta_P < 0.00001)
all_negative_lenient <- which(tag_snps$fgfp_beta < 0 & tag_snps$focus_beta < 0 & tag_snps$popgen_beta < 0 & tag_snps$meta_P < 0.00001)

## Define a sorted index
directionally_consistent_lenient <- sort(c(all_positive_lenient, all_negative_lenient))

## Define a new data frame
tag_snps_directionally_consistent_lenient <- tag_snps[directionally_consistent_lenient, ]
head(tag_snps_directionally_consistent_lenient)
dim(tag_snps_directionally_consistent_lenient) # should have 308
length(unique(tag_snps_directionally_consistent_lenient$mt)) # representing 103 associated microbial traits

## Now define what we need for MR analyses and write out
tag_snps_directionally_consistent_lenient_cols <- tag_snps_directionally_consistent_lenient[,c(1,2,3,4,5,6,7,8,14,15,16,17,18,19,20,21,22)]
names(tag_snps_directionally_consistent_lenient_cols)[names(tag_snps_directionally_consistent_lenient_cols) == 'rsid...4'] <- 'rsid'
head(tag_snps_directionally_consistent_lenient_cols)
write.table(tag_snps_directionally_consistent_lenient_cols, "./Data/tag_snps_directionally_consistent_lenient.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## Define the instruments
tag_snps_directionally_consistent_lenient_cols$MarkerName <- paste(tag_snps_directionally_consistent_lenient_cols$chr,tag_snps_directionally_consistent_lenient_cols$pos,sep=":")
all_lenient_instruments <- as.vector(tag_snps_directionally_consistent_lenient_cols$MarkerName)
unique_lenient_instruments <- unique(tag_snps_directionally_consistent_lenient_cols$MarkerName) # 299 unique instruments out of 308 SNPs

## Check how many of these are associated with the top 14
top14_lenient <- tag_snps_directionally_consistent_lenient_cols[tag_snps_directionally_consistent_lenient_cols$mt %in% c("C_Gammaproteobacteria_RNT","F_Sutterellaceae_HB","G_Bifidobacterium_RNT","G_Butyricicoccus_RNT","G_Dialister_HB","G_Parabacteroides_RNT","G_Ruminococcus_HB","G_unclassified_F_Erysipelotrichaceae_HB","G_unclassified_F_Porphyromonadaceae_RNT","G_unclassified_O_Bacteroidales_HB","G_unclassified_P_Firmicutes_HB","G_unclassified_P_Firmicutes_RNT","G_Veillonella_HB","G_Coprococcus_HB"),]
head(top14_lenient)
dim(top14_lenient) # 76 of these are associated with the 14 original phenotypes. 

# rsid rs565333512 is missing for snpid 20:41277694_T_C - looked up in gnomAD
#rsid rs573805152is missing for snpid 12:9487736_A_G - looked up in gnomAD

top14_lenient$rsid[top14_lenient$snpid=="20:41277694_T_C"] <- "rs565333512" 
top14_lenient$rsid[top14_lenient$snpid=="12:9487736_A_G"] <- "rs573805152"

top14_lenient_snps <- as.vector(top14_lenient$MarkerName)
top14_lenient_rsids <- as.vector(top14_lenient$rsid)
## Check how many SNPs are associated with the one trait that we want to follow up 
buty_lenient_snps <- tag_snps_directionally_consistent_lenient_cols[tag_snps_directionally_consistent_lenient_cols$mt=="G_Butyricicoccus_RNT",] # only increases the number of SNPs by one so not won't be able to undertake pleiotropy robust methods 

## Format exposure data
Hughes_exposure_dat <- format_data(top14_lenient, type = "exposure", snps = NULL, header = TRUE, phenotype_col = "mt", snp_col = "rsid", beta_col = "meta_beta", se_col = "meta_se", effect_allele_col = "effect_allele", other_allele_col = "ref_allele", eaf_col = "EAF", pval_col = "meta_P", samplesize_col = "N", min_pval = 1e-200, chr_col = "chr", pos_col = "pos", log_pval = FALSE)
head(Hughes_exposure_dat)
dim(Hughes_exposure_dat)


###########################################################################
## SEARCHING FOR THE MICROBIOME-SPECIFIC SNPS AND FORMATING OUTCOME DATA ##
###########################################################################
## Search for the SNPs of interest in the TRICL data
#mbg <- read.table("/Users/xh18454/Library/CloudStorage/OneDrive-UniversityofBristol/MR analysis/Data/exposure_MBG.txt", header = T, sep = "\t")


prostate_data <- extract_outcome_data(
  snps = top14_lenient_rsids,
  outcomes = 'ebi-a-GCST006085' 
)
View(prostate_data) # 75 of the 76 SNPs available
prostate_data = unique(prostate_data)
prostate_data$MarkerName <- paste(prostate_data$chr, ":", prostate_data$pos, sep = "")
dim(prostate_data)

prostate_snp_data <- subset(prostate_data, MarkerName %in% top14_lenient_snps)
dim(prostate_snp_data) #75 of the 76 associated with the top 14 phenotypes 
prostate_snp_data$rsid <- prostate_snp_data$SNP
unique(prostate_snp_data$rsid) # only 68 of the 76 SNPs are included

## Check which ones are missing and whether this is due to using the MarkerName rather than rsid
all_lenient_snps <- as.vector(top14_lenient$rsid)
all_lenient_snps # there are 76
unique(all_lenient_snps) # 75 of these are unique, which is weird considering 299 of those using the marker name are unique
prostate_test_snp_data <- subset(prostate_data, SNP %in% all_lenient_snps) 
dim(prostate_test_snp_data) 
rm(prostate_test_snp_data)


## Now checking which SNPs aren't present
prostate_hughes_merge <- merge(top14_lenient, prostate_snp_data, by.x = c("rsid"), by.y = c("rsid"))
head(prostate_hughes_merge)
dim(prostate_hughes_merge) # 
unique(prostate_hughes_merge$mt) # representing 12 individual taxa of the original 14
write.table(prostate_hughes_merge, "./Data/Hughes_prostate_merge_lenient.txt", row.names = F, col.names = T, quote = F, sep = "\t")


## Check how many SNPs are available in the outcome data associated with the one trait that we want to follow up
Parabacteroides_available <- prostate_hughes_merge[prostate_hughes_merge$mt=="G_Parabacteroides_RNT",] # snp is available so don't think there is any need for finding proxies
P_Firmicutes_available <- prostate_hughes_merge[prostate_hughes_merge$mt=="G_unclassified_P_Firmicutes_RNT",] # snp is available so don't think there is any need for finding proxies
F_Porphyromonadaceae <- prostate_hughes_merge[prostate_hughes_merge$mt=="G_unclassified_F_Porphyromonadaceae_RNT",] # snp is available so don't think there is any need for finding proxies


head(prostate_hughes_merge)
## Formating outcome data
outcome_dat <- format_data(prostate_hughes_merge, type = "outcome", snps = NULL, header = TRUE, snp_col = "rsid", beta_col = "beta.outcome", se_col = "se.outcome", effect_allele_col = "effect_allele.outcome", other_allele_col = "other_allele.outcome", eaf_col = "eaf.outcome",  pval_col = "pval.outcome", min_pval = 1e-200, chr_col = "chr.y", pos_col = "pos.y", samplesize_col = "samplesize.outcome", log_pval = FALSE)
outcome_dat$outcome <- "prostate cancer"
head(outcome_dat)
dim(outcome_dat) # 75 SNPs available 
####################################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND LUNG CANCER (OUTCOME) ##
####################################################################################
## MAIN ANALYSIS INCLUDING THE TOP HITS FROM OUR PUBLISHED META-ANALYSES
Hughes_dat <- harmonise_data(Hughes_exposure_dat, outcome_dat, action =2)
Hughes_mr_results <- mr(Hughes_dat, method_list=c("mr_wald_ratio", "mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
write.table(Hughes_mr_results, "./Output/Hughes_microbiome_PRACTICAL_results_lenient.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# Iterate over each exposure in Hughes_mr_results
for (exp in unique(Hughes_mr_results$exposure)) {
  
  # Filter the data for the current exposure
  current_results <- Hughes_mr_results[Hughes_mr_results$exposure == exp, ]
  current_dat <- Hughes_dat[Hughes_dat$exposure == exp, ]
  
  # Scatter plot
  png(paste0("./Figures/Hughes_lenient_scatter_", exp, ".png"), width = 2400, height = 2000, res = 300)
  print(mr_scatter_plot(current_results, current_dat))
  dev.off()
  
  # Forest plot
  MR_single <- mr_singlesnp(current_dat)
  png(paste0("./Figures/Hughes_lenient_forest_", exp, ".png"), width = 2400, height = 2000, res = 300)
  print(mr_forest_plot(MR_single))
  dev.off()
  
  # Funnel plot
  png(paste0("./Figures/Hughes_lenient_funnel_", exp, ".png"), width = 2400, height = 2000, res = 300)
  print(mr_funnel_plot(MR_single))
  dev.off()
  
  # Leave one out analysis
  MR_loo <- mr_leaveoneout(current_dat)
  png(paste0("./Figures/Hughes_lenient_loo_", exp, ".png"), width = 2400, height = 2000, res = 300)
  print(mr_leaveoneout_plot(MR_loo))
  dev.off()
}

################################
## Steiger filtering 
################################


gwasinfo("ebi-a-GCST006085") # 29863 cases 55586 controls
Hughes_dat$ncase.outcome <- 79148
Hughes_dat$ncontrol.outcome <- 61106
Hughes_dat$prevalence.outcome <- 0.1
Hughes_dat$units.outcome <- "log odds"

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


Hughes_steiger <- steiger_filtering(Hughes_dat)


# Add the k values to Hughes_steiger based on each exposure
Hughes_steiger$k <- 1


## Firstly calculate R2 for binary traits (for this we can use the get_r_from_lor function)
table(Hughes_steiger$exposure)
binary <- Hughes_steiger[grep("_HB", Hughes_steiger$exposure),]
binary$R <- NA
binary$R <- get_r_from_lor(lor = binary$beta.exposure, af = binary$eaf.exposure, ncase = binary$ncase.exposure, ncontrol = binary$ncontrol.exposure, prevalence = binary$prevalence.exposure, model = "logit", correction = FALSE)
binary$R2 <- binary$R^2
binary$F<- (binary$R2*(binary$samplesize.exposure-1-1))/((1-binary$R2)*1)
binary$R <- NULL

## Then calculate R2 for continuous traits 
continuous <- Hughes_steiger[grep("_RNT", Hughes_steiger$exposure),]
eaf <- continuous$eaf.exposure     
b <- continuous$beta.exposure
se <- continuous$se.exposure
p <- continuous$pval.exposure
snp <- continuous$SNP
N <- continuous$samplesize.exposure
k <- 1
continuous$VG <- 2*(b^2)*eaf*(1-eaf)
continuous$PV <- 2*(b^2)*eaf*(1-eaf) + ((se^2)*2*N*eaf*(1-eaf))
continuous$R2 <- continuous$VG/continuous$PV
continuous$F <- continuous$R2*(N-1-k)/((1-continuous$R2)*k)
continuous$VG <- NULL
continuous$PV <- NULL

## Combine
combined <- rbind(continuous, binary)
combined

## Summarising 
summary(combined$R2)
summary(combined$F)

# Save the updated table
write.table(combined, "./Output/Hughes_microbiome_PRACTICAL_steiger_filtering_lenient_Fstat_manual.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Separate binary and continuous exposures
binary <- Hughes_steiger[grep("_HB", Hughes_steiger$exposure), ]
continuous <- Hughes_steiger[grep("_RNT", Hughes_steiger$exposure), ]

# Step 2: Recalculate F for binary exposures
binary$F_statistic <- (binary$rsq.exposure * (binary$samplesize.exposure - 1 - binary$k)) / 
  ((1 - binary$rsq.exposure) * binary$k)

# Step 3: Recalculate F for continuous exposures
continuous$F_statistic <- (continuous$rsq.exposure * (continuous$samplesize.exposure - 1 - continuous$k)) / 
  ((1 - continuous$rsq.exposure) * continuous$k)

# Step 4: Combine the binary and continuous data back together
combined <- rbind(continuous, binary)

# Summarise results
summary(combined$F_statistic)

# Save the updated table
write.table(combined, "./Output/Hughes_microbiome_PRACTICAL_steiger_filtering_lenient_Fstat.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)


#write.table(Hughes_steiger, "./Output/Hughes_microbiome_PRACTICAL_steiger_filtering_lenient.txt", sep = "\t", col.names = T, row.names = F, quote = F)


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
ggsave("./Figures/forestplot_microbiome_PRACTICAL_Hughes_lenient.png", width=7, height=8)

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
ggsave("./Figures/volcanoplot_microbiome_PRACTICAL_Hughes_lenient.png", width=10, height=8)


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
dim(Hughes_UKBB_outcome_dat) # 73 Hughes SNPs in UK Biobank data


#########################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND METABS (OUTCOME) ##
#########################################################################
## MAIN ANALYSIS INCLUDING THE TOP HITS FROM OUR PUBLISHED META-ANALYSES
Hughes_UKBB_dat <- harmonise_data(Hughes_exposure_dat, Hughes_UKBB_outcome_dat, action =2)
Hughes_UKBB_mr_results <- mr(Hughes_UKBB_dat, method_list=c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
write.table(Hughes_UKBB_mr_results, "./Output/Hughes_microbiome_UKBB_results_lenient.txt", sep = "\t", col.names = T, row.names = F, quote = F)

gwasinfo("ieu-b-4809") # 9132 cases 173493 controls
Hughes_UKBB_dat$ncase.outcome <- 9132
Hughes_UKBB_dat$ncontrol.outcome <- 173493
Hughes_UKBB_dat$prevalence.outcome <- 0.1
Hughes_UKBB_dat$units.outcome <- "log odds"
steiger_filtering(Hughes_UKBB_dat)
write.table(steiger_filtering(Hughes_UKBB_dat), "./Output/Hughes_microbiome_UKBB_steiger_filtering_lenient.txt", sep = "\t", col.names = T, row.names = F, quote = F)

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
Hughes_UKBB_fres$outcome[Hughes_UKBB_fres$outcome=="prostate Cancer || id:ieu-b-4809"] <- "prostate cancer"
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
ggsave("./Figures/forestplot_microbiome_UKBB_Hughes_lenient.png", width=7, height=8)


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
ggsave("./Figures/volcanoplot_microbiome_UKBB_prostate_Hughes_lenient.png", width=10, height=8)


######################. Part 3 - Finngen #######################################

## Saving outcome data from IEU OpenGWAS 
# finn-b-C3_BRONCHUS_LUNG_EXALLC = Malignant neoplasm of bronchus and lung (all cancers excluded) in Finngen

outcomes <- c("finn-b-C3_PROSTATE_EXALLC")
## Extract the outcomes
Hughes_finngen_outcome_dat <- extract_outcome_data(Hughes_exposure_dat$SNP, outcomes = outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, splitsize = 10000, proxy_splitsize = 500)
head(Hughes_finngen_outcome_dat) 
Hughes_finngen_outcome_dat <- unique(Hughes_finngen_outcome_dat)

#########################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND METABS (OUTCOME) ##
#########################################################################
## MAIN ANALYSIS INCLUDING THE TOP HITS FROM OUR PUBLISHED META-ANALYSES
Hughes_finngen_dat <- harmonise_data(Hughes_exposure_dat, Hughes_finngen_outcome_dat, action =2)
Hughes_finngen_dat <- unique(Hughes_finngen_dat)
Hughes_finngen_mr_results <- mr(Hughes_finngen_dat, method_list=c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
write.table(Hughes_finngen_mr_results, "./Output/Hughes_microbiome_FinnGen_results_lenient.txt", sep = "\t", col.names = T, row.names = F, quote = F) # 14 mt associations


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
Hughes_finngen_dat$ncontrol.outcome <-74685
Hughes_finngen_dat$prevalence.outcome <- 0.1
Hughes_finngen_dat$units.outcome <- "log odds"
steiger_filtering(Hughes_finngen_dat)
write.table(steiger_filtering(Hughes_finngen_dat), "./Output/Hughes_microbiome_finngen_steiger_filtering_lenient.txt", sep = "\t", col.names = T, row.names = F, quote = F)


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
ggsave("./Figures/forestplot_microbiome_FinnGen_prostate_Hughes_lenient.png", width=7, height=8)



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
ggsave("./Figures/volcanoplot_microbiome_finngen_prostate_Hughes_lenient.png", width=10, height=8)


##################################################### META-ANALYSING PRACTICAL, UK Biobank and Finngen #####################################################
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


## Bind together
Hughes_plot_data <- rbind(Hughes_main, Hughes_UKBB_sens, Hughes_finngen_sens)

## Now limit to the ones that are consistent across both datasets
Hughes_plot_data <- Hughes_plot_data[ave(1:nrow(Hughes_plot_data), Hughes_plot_data$exposure, FUN=length)>1,]

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
ggsave("./Figures/Hughes_microbiome_PRACTICAL_UKBB_finngen_lenient.png", width=10, height=8)


## Simple meta-analysis with Hughes
Hughes_meta <- metagen(Hughes_plot_data$beta, seTE = Hughes_plot_data$se, studlab = Hughes_plot_data$dataset, data = Hughes_plot_data, sm = "OR", fixed = TRUE, random = FALSE)
Hughes_meta_by_exposure <- update(Hughes_meta, subgroup = Hughes_plot_data$exposure, tau.common = FALSE)
summary(Hughes_meta_by_exposure)
pdf("./Figures/Hughes_microbiome_PRACTICAL+UKBB+finngen_meta_lenient.pdf", width=10, height=16)
forest(Hughes_meta_by_exposure, studlab = T, common = TRUE, random = FALSE, leftcols = c("dataset"), print.subgroup.name = FALSE, print.subgroup.labels = TRUE, test.subgroup= FALSE, overall = FALSE, overall.hetstat = FALSE, text.common = "Meta-analysis", prediction = FALSE, xlab = "OR of lung cancer with gut microbiome variation", print.tau2 = FALSE, plotwidth = "15cm", fontsize = 10, spacing = .75)
dev.off()


#### creating plots for just the traits we want to follow up

####################################
## ADDING FOREST PLOTS OF RESULTS ##
####################################
## Hughes

Hughes_fres_G_unclassified_P_Firmicutes_RNT <- Hughes_fres[Hughes_fres$exposure=="G_unclassified_P_Firmicutes_RNT",]
Hughes_fres_G_unclassified_F_Porphyromonadaceae_RNT_available <- Hughes_fres[Hughes_fres$exposure=="G_unclassified_F_Porphyromonadaceae_RNT",]
Hughes_fres_G_Parabacteroides_RNT_available <- Hughes_fres[Hughes_fres$exposure=="G_Parabacteroides_RNT",]

Hughes_fres <- rbind(Hughes_fres_G_unclassified_P_Firmicutes_RNT, Hughes_fres_G_unclassified_F_Porphyromonadaceae_RNT_available, Hughes_fres_G_Parabacteroides_RNT_available)


##
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
  labs(title="Effect of the gut microbiome on overall colorectal cancer", x="Odds ratio", y = "Microbial trait (Hughes)") +
  geom_vline(xintercept=1, color="black", linetype="dashed", alpha=.5) +
  theme_classic()
ggsave("./Figures/forestplot_microbiome_prostate_Hughes_traits_to_follow_up_lenient.png", width=7, height=8)

######################################################################  MIBIOGEN ######################################################################
#### Using the European only meta-analysis data, for each microbial trait

# only need to do for the traits we want to follow up  - genus.Allisonella.id.2174)

###### 
genus.Allisonella.id.2174 <- data.table::fread("/Users/xh18454/Library/CloudStorage/OneDrive-UniversityofBristol/MR/lenient_MR/betafilewMAF/genus.Allisonella.id.2174.summaryMAF.txt.gz")
head(genus.Allisonella.id.2174)
dim(genus.Allisonella.id.2174)
Allisonella_hits <- genus.Allisonella.id.2174[genus.Allisonella.id.2174$P.weightedSumZ<=1e-5,]
Allisonella_hits <- as.data.frame(Allisonella_hits)

Allisonella_hits_formated <- format_data(Allisonella_hits, type = "exposure", snps = NULL, header = TRUE, phenotype_col = "bac", snp_col = "rsID", beta_col = "beta", se_col = "SD", eaf_col = "weighted_MAF", effect_allele_col = "eff.allele", other_allele_col = "ref.allele",  pval_col = "P.weightedSumZ", samplesize_col = "N", min_pval = 1e-200, z_col = "Z.weightedSumZ", chr_col = "chr", pos_col = "bp", log_pval = FALSE)
dim(Allisonella_hits_formated)
# clump SNPs
Allisonella_hits_clumped <- clump_data(Allisonella_hits_formated, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, pop = "EUR") # there should be 30
Allisonella_hits_clumped$MarkerName <- paste(Allisonella_hits_clumped$chr.exposure, ":", Allisonella_hits_clumped$pos.exposure, sep ="")
Allisonella_hits_clumped <- Allisonella_hits_clumped[order(Allisonella_hits_clumped$SNP),]
head(Allisonella_hits_clumped)
dim(Allisonella_hits_clumped) # 11 clumped SNPs
write.table(Allisonella_hits_clumped,"./Data/Allisonella_GW_clumped_hits_lenient.txt", col.names = T, row.names = F, quote = F, sep = "\t")



### Find clumped SNPs which are directionally consistent (agreement above 70%)
Allisonella_snps <- Allisonella_hits_clumped$MarkerName

Allisonella_z_agree <- data.table::fread("/Users/xh18454/Library/CloudStorage/OneDrive-UniversityofBristol/MR/lenient_MR/z_agreement/genus.Allisonella.id.2174_cohortz.gz", header = T, sep = "\t")
Allisonella_z_agree_clumped_snps <- subset(Allisonella_z_agree, SNPName %in% Allisonella_snps)
Allisonella_z_agree_clumped_snps <- Allisonella_z_agree_clumped_snps[, c(1,2,3,4,5,7,28,29,30)]
head(Allisonella_z_agree_clumped_snps) # SNP 

Allisonella_merged <- merge(Allisonella_hits_clumped, Allisonella_z_agree_clumped_snps, by.x = "SNP", by.y = "rsID") # 
dim(Allisonella_merged)
Allisonella_new_hits <- Allisonella_merged[Allisonella_merged$percentagree>=1,] # keep the SNPs where the percentage of studies directionally consistent is 100%
Allisonella_new_hits # 4 SNPs with 100% directional consistency

MiBioGen_hits <- Allisonella_new_hits
write.table(MiBioGen_hits,"./Data/Allisonella_GW_clumped_hits_lenient_directionally_consistent_correct_SE.txt", col.names = T, row.names = F, quote = F, sep = "\t")

###########################################################################
## Run the analysis for directionally consistent snps 
###########################################################################
## Format exposure data
MiBioGen_exposure_dat <- format_data(MiBioGen_hits, type = "exposure", snps = NULL, header = TRUE, phenotype_col = "exposure", snp_col = "SNP", beta_col = "beta.exposure", se_col = "se.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure", eaf_col = "eaf.exposure", pval_col = "pval.exposure", samplesize_col = "samplesize.exposure", min_pval = 1e-200, chr_col = "chr.exposure", pos_col = "pos.exposure", log_pval = FALSE)
head(MiBioGen_exposure_dat)
dim(MiBioGen_exposure_dat) # 380 SNP-exposures 

lenient_rsid <- MiBioGen_exposure_dat$SNP 

prostate_data <- extract_outcome_data(
  snps = lenient_rsid,
  outcomes = 'ebi-a-GCST006085' 
)
View(prostate_data)
prostate_data = unique(prostate_data)
prostate_data$MarkerName <- paste(prostate_data$chr, ":", prostate_data$pos, sep = "")
dim(prostate_data)

# Now checking which SNPs aren't present
prostate_MiBioGen_merge <- merge(MiBioGen_exposure_dat, prostate_data, by.x = c("SNP"), by.y = c("SNP"))
head(prostate_MiBioGen_merge)
dim(prostate_MiBioGen_merge)

## Formating outcome data
MiBioGen_outcome_dat <- format_data(prostate_MiBioGen_merge, type = "outcome", snps = NULL, header = TRUE, snp_col = "SNP", beta_col = "beta.outcome", se_col = "se.outcome", effect_allele_col = "effect_allele.outcome", other_allele_col = "other_allele.outcome", eaf_col = "eaf.outcome",  pval_col = "pval.outcome", min_pval = 1e-200, chr_col = "chr", pos_col = "pos", log_pval = FALSE)
MiBioGen_outcome_dat$outcome <- "prostate cancer"
head(MiBioGen_outcome_dat) 

####################################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND LUNG CANCER (OUTCOME) ##
####################################################################################
## MAIN ANALYSIS INCLUDING THE TOP HITS FROM OUR PUBLISHED META-ANALYSES
MiBioGen_dat <- harmonise_data(MiBioGen_exposure_dat, MiBioGen_outcome_dat, action =1)
MiBioGen_mr_results <- mr(MiBioGen_dat, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
write.table(MiBioGen_mr_results, "./Output/MiBioGen_microbiome_PRACTICAL_results_lenient_directionally_consistent_correct_SE.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# Scatter plot
png(paste("./Figures/MiBioGen_lenient_directionally_consistent_scatter.png", sep=""), width = 2400, res = 300)
print(mr_scatter_plot(MiBioGen_mr_results, MiBioGen_dat))
dev.off()

# Forest plot
MR_single <- mr_singlesnp(MiBioGen_dat)
png(paste("./Figures/MiBioGen_lenient_directionally_consistent_forest.png", sep=""), width = 2400, res = 300)
print(mr_forest_plot(MR_single))
dev.off()

# Funnel plot
png(paste("./Figures/MiBioGen_lenient_directionally_consistent_funnel.png", sep=""), width = 2400, res = 300)
print(mr_funnel_plot(MR_single))
dev.off()

# Leave one out analysis
MR_loo <- mr_leaveoneout(MiBioGen_dat)
png(paste("./Figures/MiBioGen_lenient_directionally_consistent_loo.png", sep=""), width = 2400, res = 300)
print(mr_leaveoneout_plot(MR_loo))
dev.off()

################################
## Steiger filtering 
################################

## Need to change Hughes to MiBioGen 
gwasinfo("ebi-a-GCST006085") # 29863 cases 55586 controls
MiBioGen_dat$ncase.outcome <- 79148
MiBioGen_dat$ncontrol.outcome <- 61106
MiBioGen_dat$prevalence.outcome <- 0.1
MiBioGen_dat$units.outcome <- "log odds"
MiBioGen_steiger <- steiger_filtering(MiBioGen_dat)

# Add the k values to Hughes_steiger based on each exposure
MiBioGen_steiger$k <- 1


# Recalculate F_statistic_exposure and F_statistic_outcome
MiBioGen_steiger$F_statistic_exposure <- (MiBioGen_steiger$rsq.exposure * (MiBioGen_steiger$samplesize.exposure - 1 - MiBioGen_steiger$k)) / 
  ((1 - MiBioGen_steiger$rsq.exposure) * MiBioGen_steiger$k)

## Then calculate R2 for continuous traits 
continuous <- MiBioGen_steiger
eaf <- continuous$eaf.exposure     
b <- continuous$beta.exposure
se <- continuous$se.exposure
p <- continuous$pval.exposure
snp <- continuous$SNP
N <- continuous$samplesize.exposure
k <- 1
continuous$VG <- 2*(b^2)*eaf*(1-eaf)
continuous$PV <- 2*(b^2)*eaf*(1-eaf) + ((se^2)*2*N*eaf*(1-eaf))
continuous$R2 <- continuous$VG/continuous$PV
continuous$F_statistic <- continuous$R2*(N-1-k)/((1-continuous$R2)*k)
continuous$VG <- NULL
continuous$PV <- NULL


write.table(continuous, "./Output/MiBioGen_microbiome_PRACTICAL_steiger_filtering_lenient_directionally_consistent_Fstat_correct_SE.txt", sep = "\t", col.names = T, row.names = F, quote = F)


write.table(MiBioGen_steiger, "./Output/MiBioGen_microbiome_PRACTICAL_steiger_filtering_lenient.txt", sep = "\t", col.names = T, row.names = F, quote = F)

MiBioGen_fres <- MiBioGen_mr_results
MiBioGen_fres$sig <- MiBioGen_fres$pval < 0.05
MiBioGen_fres$category <-"<0.05"
MiBioGen_fres$or <- exp(MiBioGen_fres$b)
MiBioGen_fres$lci <- exp((MiBioGen_fres$b - (1.96*MiBioGen_fres$se)))
MiBioGen_fres$uci <- exp((MiBioGen_fres$b + (1.96*MiBioGen_fres$se)))
MiBioGen_fres$label <- MiBioGen_fres$exposure
MiBioGen_fres$observation <- 1:nrow(MiBioGen_fres) 
head(MiBioGen_fres)


#### Forest plot ########
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
ggsave("./Figures/forestplot_microbiome_PRACTICAL_MiBioGen_lenient_sensitivity_directionally_consistent.png", width=12, height=8)

###### Volcano plot

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
ggsave("./Figures/volcanoplot_microbiome_PRACTICAL_MiBioGen_lenient_sensitivity_directionally_consistent.png", width=12, height=8)

############# part 2 UK Biobank #####################################
## Set up vector containing all outcomes
outcomes <- c("ieu-b-4809") #try extract instruments

## Extract the outcomes
MiBioGen_UKBB_outcome_dat <- extract_outcome_data(MiBioGen_exposure_dat$SNP, outcomes = outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, splitsize = 10000, proxy_splitsize = 500)
head(MiBioGen_UKBB_outcome_dat) 
dim(MiBioGen_UKBB_outcome_dat) # 73 Hughes SNPs in UK Biobank data

#########################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND METABS (OUTCOME) ##
#########################################################################
## MAIN ANALYSIS INCLUDING THE TOP HITS FROM OUR PUBLISHED META-ANALYSES

MiBioGen_UKBB_dat <- harmonise_data(MiBioGen_exposure_dat, MiBioGen_UKBB_outcome_dat, action =1)
MiBioGen_UKBB_mr_results <- mr(MiBioGen_UKBB_dat, method_list=c("mr_wald_ratio","mr_ivw", "mr_weighted_median","mr_weighted_mode", "mr_egger_regression"))
write.table(MiBioGen_UKBB_mr_results, "./Output/MiBioGen_EUR_microbiome_UKBB_results_lenient_sensitivity_directionally_consistent.txt", sep = "\t", col.names = T, row.names = F, quote = F)

gwasinfo("ieu-b-4809")
####### MiBioGen 
MiBioGen_UKBB_dat$units.exposure <- "SD"

MiBioGen_UKBB_dat$ncase.outcome <- 9132
MiBioGen_UKBB_dat$ncontrol.outcome <- 173493
MiBioGen_UKBB_dat$prevalence.outcome <- 0.1
MiBioGen_UKBB_dat$units.outcome <- "log odds"

steiger_filtering(MiBioGen_UKBB_dat)
write.table(steiger_filtering(MiBioGen_UKBB_dat), "./Output/MiBioGen_EUR_microbiome_UKBB_steiger_filtering_lenient_sensitivity_directionally_consistent.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# MiBioGen
MiBioGen_UKBB_fres <- MiBioGen_UKBB_mr_results
MiBioGen_UKBB_fres$sig <- MiBioGen_UKBB_fres$pval < 0.05
MiBioGen_UKBB_fres$category <-"<0.05"
MiBioGen_UKBB_fres$or <- exp(MiBioGen_UKBB_fres$b)
MiBioGen_UKBB_fres$lci <- exp((MiBioGen_UKBB_fres$b - (1.96*MiBioGen_UKBB_fres$se)))
MiBioGen_UKBB_fres$uci <- exp((MiBioGen_UKBB_fres$b + (1.96*MiBioGen_UKBB_fres$se)))
MiBioGen_UKBB_fres$outcome[MiBioGen_UKBB_fres$outcome=="Prostate cancer || id:ieu-b-4809"] <- "prostate cancer"
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
ggsave("./Figures/forestplot_microbiome_UKBB_MiBioGen_lenient_sensitivity_directionally_consistent.png", width=10, height=8)

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
ggsave("./Figures/volcanoplot_microbiome_UKBB_pancreatic_MiBioGen_lenient_sensitivity_directionally_consistent.png", width=10, height=8)


##################################################### PART 3 - FINNGEN #####################################################


###############################################
## EXTRACTING OUTCOME DATA FROM IEU OPENGWAS ##
###############################################
## Saving outcome data from IEU OpenGWAS 


## Set up vector containing all outcomes
outcomes <- c("finn-b-C3_PROSTATE_EXALLC")

## Extract the outcomes

MiBioGen_finngen_outcome_dat <- extract_outcome_data(MiBioGen_exposure_dat$SNP, outcomes = outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, splitsize = 10000, proxy_splitsize = 500)
head(MiBioGen_finngen_outcome_dat)


#########################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND METABS (OUTCOME) ##
#########################################################################
## MAIN ANALYSIS INCLUDING THE TOP HITS FROM OUR PUBLISHED META-ANALYSES

MiBioGen_finngen_dat <- harmonise_data(MiBioGen_exposure_dat, MiBioGen_finngen_outcome_dat, action =1)
MiBioGen_finngen_mr_results <- mr(MiBioGen_finngen_dat, method_list=c("mr_wald_ratio","mr_ivw", "mr_weighted_median","mr_weighted_mode", "mr_egger_regression"))
write.table(MiBioGen_finngen_mr_results, "./Output/MiBioGen_microbiome_FinnGen_results_lenient_sensitivity_directionally_consistent.txt", sep = "\t", col.names = T, row.names = F, quote = F)

gwasinfo("finn-b-C3_PROSTATE_EXALLC")
####### MiBioGen 
MiBioGen_finngen_dat$units.exposure <- "SD"

MiBioGen_finngen_dat$ncase.outcome <- 6311
MiBioGen_finngen_dat$ncontrol.outcome <- 74685
MiBioGen_finngen_dat$prevalence.outcome <- 0.1
MiBioGen_finngen_dat$units.outcome <- "log odds"

MiBioGen_finngen_steiger <- steiger_filtering(MiBioGen_finngen_dat)
write.table(MiBioGen_finngen_steiger, "./Output/MiBioGen_microbiome_finngen_steiger_filtering_lenient_sensitivity_directionally_consistent.txt", sep = "\t", col.names = T, row.names = F, quote = F)


# MiBioGen
MiBioGen_finngen_fres <- MiBioGen_finngen_mr_results
MiBioGen_finngen_fres$sig <- MiBioGen_finngen_fres$pval < 0.05
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
ggsave("./Figures/forestplot_microbiome_FinnGen_MiBioGen_lenient_sensitivity_directionally_consistent.png", width=10, height=8)

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
ggsave("./Figures/volcanoplot_microbiome_finngen_pancreatic_MiBioGen_lenient_sensitivity_directionally_consistent.png", width=10, height=8)



###########################################################################
## Run the analysis  - lenient threshold
###########################################################################

MiBioGen_hits <- Allisonella_hits_clumped
dim(MiBioGen_hits)## 39 SNPs 
write.table(MiBioGen_hits,"./Data/Allisonella_exposure_data_clumped_hits_lenient_sensitivity_correct_SE.txt", col.names = T, row.names = F, quote = F, sep = "\t")
lenient_MarkerName <- MiBioGen_hits$MarkerName


## Format exposure data
MiBioGen_exposure_dat <- format_data(MiBioGen_hits, type = "exposure", snps = NULL, header = TRUE, phenotype_col = "exposure", snp_col = "SNP", beta_col = "beta.exposure", se_col = "se.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure", eaf_col = "eaf.exposure", pval_col = "pval.exposure", samplesize_col = "samplesize.exposure", min_pval = 1e-200, chr_col = "chr.exposure", pos_col = "pos.exposure", log_pval = FALSE)
head(MiBioGen_exposure_dat)
dim(MiBioGen_exposure_dat) # 380 SNP-exposures 


lenient_rsid <- MiBioGen_exposure_dat$SNP 

EUR.bim <- read.table("/Users/xh18454/OneDrive - University of Bristol/MR/coloc_locus_zoom/reference_panel/processed/EUR/EUR.bim", header = F)
stats.stats <- read.table("/Users/xh18454/OneDrive - University of Bristol/MR/coloc_locus_zoom/reference_panel/processed/EUR/stats.stats", header = T)
head(EUR.bim)
dim(EUR.bim)
dim(stats.stats)
head(stats.stats)
colnames(stats.stats)[1] <- "rsid"
colnames(EUR.bim)[1] <- "chr"
colnames(EUR.bim)[2] <- "rsid"
colnames(EUR.bim)[3] <- "V3"
colnames(EUR.bim)[4] <- "pos"
colnames(EUR.bim)[5] <- "Allele1"
colnames(EUR.bim)[6] <- "Allele2"

EUR.bim_stats_merged <- merge(EUR.bim, stats.stats, by.x = "rsid", by.y = "rsid")
EUR.bim_stats_merged$MarkerName <- paste(EUR.bim_stats_merged$chr, ":", EUR.bim_stats_merged$pos, sep = "")


# Merge the dataframes based on MiBioGen_coloc_combined$SNP and EUR.bim_stats_merged$MarkerName
merged_data <- merge(MiBioGen_exposure_dat, EUR.bim_stats_merged, by.x = "SNP", by.y = "rsid", all.x = TRUE)


for (j in 1:nrow(merged_data)) {
  if (merged_data$effect_allele.exposure[j] == merged_data$Allele1[j]) {
    merged_data$eaf.exposure[j] <- merged_data$MAF[j]
  } else {
    merged_data$eaf.exposure[j] <- 1 - merged_data$MAF[j]
  }
}

# Update the missing eaf.outcome values with the corresponding MAF values from EUR.bim_stats_merged
MiBioGen_exposure_dat$eaf.exposure <- merged_data$eaf.exposure

# View the updated MiBioGen_coloc_combined dataframe
head(MiBioGen_exposure_dat)
any(is.na(MiBioGen_exposure_dat$eaf.exposure))
sum(is.na(MiBioGen_exposure_dat$eaf.exposure)) # all SNPs have EAF

write.table(MiBioGen_exposure_dat, "./Data/MiBioGen_lenient_directionally_consistent_breast_1000G_eaf.txt", sep = "\t", col.names = T, row.names = F, quote = F)  
dim(MiBioGen_exposure_dat)
MiBioGen_exposure_dat$eaf.exposure <- as.numeric(MiBioGen_exposure_dat$eaf.exposure)
lenient_MarkerName <- MiBioGen_hits$MarkerName





prostate_data <- extract_outcome_data(
  snps = lenient_rsid,
  outcomes = 'ebi-a-GCST006085' 
)
View(prostate_data)
prostate_data = unique(prostate_data)
prostate_data$MarkerName <- paste(prostate_data$chr, ":", prostate_data$pos, sep = "")
dim(prostate_data)

# Now checking which SNPs aren't present
prostate_MiBioGen_merge <- merge(MiBioGen_exposure_dat, prostate_data, by.x = c("SNP"), by.y = c("SNP"))
head(prostate_MiBioGen_merge)
dim(prostate_MiBioGen_merge)

## Formating outcome data
MiBioGen_outcome_dat <- format_data(prostate_MiBioGen_merge, type = "outcome", snps = NULL, header = TRUE, snp_col = "SNP", beta_col = "beta.outcome", se_col = "se.outcome", effect_allele_col = "effect_allele.outcome", other_allele_col = "other_allele.outcome", eaf_col = "eaf.outcome",  pval_col = "pval.outcome", min_pval = 1e-200, chr_col = "chr", pos_col = "pos", log_pval = FALSE)
MiBioGen_outcome_dat$outcome <- "prostate cancer"
head(MiBioGen_outcome_dat) 

####################################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND LUNG CANCER (OUTCOME) ##
####################################################################################
## MAIN ANALYSIS INCLUDING THE TOP HITS FROM OUR PUBLISHED META-ANALYSES
MiBioGen_dat <- harmonise_data(MiBioGen_exposure_dat, MiBioGen_outcome_dat, action =1)
MiBioGen_mr_results <- mr(MiBioGen_dat, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
write.table(MiBioGen_mr_results, "./Output/MiBioGen_microbiome_PRACTICAL_results_lenient_sensitivity_correct_SE.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# Scatter plot
png(paste("./Figures/MiBioGen_lenient_scatter.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_scatter_plot(MiBioGen_mr_results, MiBioGen_dat))
dev.off()

# Forest plot
MR_single <- mr_singlesnp(MiBioGen_dat)
png(paste("./Figures/MiBioGen_lenient_forest.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_forest_plot(MR_single))
dev.off()

# Funnel plot
png(paste("./Figures/MiBioGen_lenient_funnel.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_funnel_plot(MR_single))
dev.off()

# Leave one out analysis
MR_loo <- mr_leaveoneout(MiBioGen_dat)
png(paste("./Figures/MiBioGen_lenient_loo.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_leaveoneout_plot(MR_loo))
dev.off()

################################
## Steiger filtering 
################################


gwasinfo("ebi-a-GCST006085")
MiBioGen_dat$ncase.outcome <- 79148
MiBioGen_dat$ncontrol.outcome <- 61106
MiBioGen_dat$prevalence.outcome <- 0.1
MiBioGen_dat$units.outcome <- "log odds"
MiBioGen_steiger <- steiger_filtering(MiBioGen_dat)


# Add the k values to Hughes_steiger based on each exposure
MiBioGen_steiger$k <- 1

## Then calculate R2 for continuous traits 
continuous <- MiBioGen_steiger
eaf <- continuous$eaf.exposure     
b <- continuous$beta.exposure
se <- continuous$se.exposure
p <- continuous$pval.exposure
snp <- continuous$SNP
N <- continuous$samplesize.exposure
k <- 1
continuous$VG <- 2*(b^2)*eaf*(1-eaf)
continuous$PV <- 2*(b^2)*eaf*(1-eaf) + ((se^2)*2*N*eaf*(1-eaf))
continuous$R2 <- continuous$VG/continuous$PV
continuous$F_statistic <- continuous$R2*(N-1-k)/((1-continuous$R2)*k)
continuous$VG <- NULL
continuous$PV <- NULL


# Recalculate F_statistic_exposure and F_statistic_outcome
#MiBioGen_steiger$F_statistic_exposure <- (MiBioGen_steiger$rsq.exposure * (MiBioGen_steiger$samplesize.exposure - 1 - MiBioGen_steiger$k)) / 
#  ((1 - MiBioGen_steiger$rsq.exposure) * MiBioGen_steiger$k)

#write.table(MiBioGen_steiger, "./Output/MiBioGen_microbiome_PRACTICAL_steiger_filtering_lenient_snps_sensitivity_Fstat.txt", sep = "\t", col.names = T, row.names = F, quote = F)


write.table(continuous, "./Output/MiBioGen_microbiome_PRACTICAL_steiger_filtering_all_lenient_snps_sensitivity_correct_SE.txt", sep = "\t", col.names = T, row.names = F, quote = F)

MiBioGen_fres <- MiBioGen_mr_results
MiBioGen_fres$sig <- MiBioGen_fres$pval < 0.05
MiBioGen_fres$category <-"<0.05"
MiBioGen_fres$or <- exp(MiBioGen_fres$b)
MiBioGen_fres$lci <- exp((MiBioGen_fres$b - (1.96*MiBioGen_fres$se)))
MiBioGen_fres$uci <- exp((MiBioGen_fres$b + (1.96*MiBioGen_fres$se)))
MiBioGen_fres$label <- MiBioGen_fres$exposure
MiBioGen_fres$observation <- 1:nrow(MiBioGen_fres) 
head(MiBioGen_fres)


#### Forest plot ########
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
ggsave("./Figures/forestplot_microbiome_PRACTICAL_MiBioGen_lenient_sensitivity.png", width=12, height=8)

###### Volcano plot

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
ggsave("./Figures/volcanoplot_microbiome_PRACTICAL_MiBioGen_lenient_sensitivity.png", width=12, height=8)

############# part 2 UK Biobank #####################################
## Set up vector containing all outcomes
outcomes <- c("ieu-b-4809") #try extract instruments

## Extract the outcomes
MiBioGen_UKBB_outcome_dat <- extract_outcome_data(MiBioGen_exposure_dat$SNP, outcomes = outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, splitsize = 10000, proxy_splitsize = 500)
head(MiBioGen_UKBB_outcome_dat) 
dim(MiBioGen_UKBB_outcome_dat) # 73 Hughes SNPs in UK Biobank data

#########################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND METABS (OUTCOME) ##
#########################################################################
## MAIN ANALYSIS INCLUDING THE TOP HITS FROM OUR PUBLISHED META-ANALYSES

MiBioGen_UKBB_dat <- harmonise_data(MiBioGen_exposure_dat, MiBioGen_UKBB_outcome_dat, action =1)
MiBioGen_UKBB_mr_results <- mr(MiBioGen_UKBB_dat, method_list=c("mr_wald_ratio","mr_ivw", "mr_weighted_median","mr_weighted_mode", "mr_egger_regression"))
write.table(MiBioGen_UKBB_mr_results, "./Output/MiBioGen_EUR_microbiome_UKBB_results_lenient_sensitivity.txt", sep = "\t", col.names = T, row.names = F, quote = F)

gwasinfo("ieu-b-4809")
####### MiBioGen 
MiBioGen_UKBB_dat$units.exposure <- "SD"

MiBioGen_UKBB_dat$ncase.outcome <- 9132
MiBioGen_UKBB_dat$ncontrol.outcome <- 173493
MiBioGen_UKBB_dat$prevalence.outcome <- 0.1
MiBioGen_UKBB_dat$units.outcome <- "log odds"

steiger_filtering(MiBioGen_UKBB_dat)
write.table(steiger_filtering(MiBioGen_UKBB_dat), "./Output/MiBioGen_EUR_microbiome_UKBB_steiger_filtering_lenient_sensitivity.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# MiBioGen
MiBioGen_UKBB_fres <- MiBioGen_UKBB_mr_results
MiBioGen_UKBB_fres$sig <- MiBioGen_UKBB_fres$pval < 0.05
MiBioGen_UKBB_fres$category <-"<0.05"
MiBioGen_UKBB_fres$or <- exp(MiBioGen_UKBB_fres$b)
MiBioGen_UKBB_fres$lci <- exp((MiBioGen_UKBB_fres$b - (1.96*MiBioGen_UKBB_fres$se)))
MiBioGen_UKBB_fres$uci <- exp((MiBioGen_UKBB_fres$b + (1.96*MiBioGen_UKBB_fres$se)))
MiBioGen_UKBB_fres$outcome[MiBioGen_UKBB_fres$outcome=="Prostate cancer || id:ieu-b-4809"] <- "prostate cancer"
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
ggsave("./Figures/forestplot_microbiome_UKBB_MiBioGen_lenient_sensitivity.png", width=10, height=8)

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
ggsave("./Figures/volcanoplot_microbiome_UKBB_pancreatic_MiBioGen_lenient_sensitivity.png", width=10, height=8)


##################################################### PART 3 - FINNGEN #####################################################


###############################################
## EXTRACTING OUTCOME DATA FROM IEU OPENGWAS ##
###############################################
## Saving outcome data from IEU OpenGWAS 


## Set up vector containing all outcomes
outcomes <- c("finn-b-C3_PROSTATE_EXALLC")

## Extract the outcomes

MiBioGen_finngen_outcome_dat <- extract_outcome_data(MiBioGen_exposure_dat$SNP, outcomes = outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, splitsize = 10000, proxy_splitsize = 500)
head(MiBioGen_finngen_outcome_dat)


#########################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND METABS (OUTCOME) ##
#########################################################################
## MAIN ANALYSIS INCLUDING THE TOP HITS FROM OUR PUBLISHED META-ANALYSES

MiBioGen_finngen_dat <- harmonise_data(MiBioGen_exposure_dat, MiBioGen_finngen_outcome_dat, action =1)
MiBioGen_finngen_mr_results <- mr(MiBioGen_finngen_dat, method_list=c("mr_wald_ratio","mr_ivw", "mr_weighted_median","mr_weighted_mode", "mr_egger_regression"))
write.table(MiBioGen_finngen_mr_results, "./Output/MiBioGen_microbiome_FinnGen_results_lenient.txt", sep = "\t", col.names = T, row.names = F, quote = F)

gwasinfo("finn-b-C3_PROSTATE_EXALLC")
####### MiBioGen 
MiBioGen_finngen_dat$units.exposure <- "SD"

MiBioGen_finngen_dat$ncase.outcome <- 6311
MiBioGen_finngen_dat$ncontrol.outcome <- 74685
MiBioGen_finngen_dat$prevalence.outcome <- 0.1
MiBioGen_finngen_dat$units.outcome <- "log odds"

MiBioGen_finngen_steiger <- steiger_filtering(MiBioGen_finngen_dat)
write.table(MiBioGen_finngen_steiger, "./Output/MiBioGen_microbiome_finngen_steiger_filtering_lenient_sensitivity.txt", sep = "\t", col.names = T, row.names = F, quote = F)


# MiBioGen
MiBioGen_finngen_fres <- MiBioGen_finngen_mr_results
MiBioGen_finngen_fres$sig <- MiBioGen_finngen_fres$pval < 0.05
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
ggsave("./Figures/forestplot_microbiome_FinnGen_MiBioGen_lenient_sensitivity.png", width=10, height=8)

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
ggsave("./Figures/volcanoplot_microbiome_finngen_pancreatic_MiBioGen_lenient.png", width=10, height=8)



