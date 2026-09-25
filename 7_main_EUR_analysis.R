######################################################################################################################
######## PROJECT: MR analysis of the gut microbiome on lung cancer using TRICL lung cancer (ieu-a-987) GWAS data and UK Biobank and FinnGen lung cancer GWAS
######## Script: Run the main MR analysis of the gut microbiome and lung cancer
######## Date: 18/07/24
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
#install.packages("MASS")
#install.packages("calibrate")
#install.packages("knitr")
#install.packages("patchwork")
#install.packages("purrr")
#install.packages("ggforestplot") #Not able to install in this version of R.
#devtools::install_github("NightingaleHealth/ggforestplot")
#install.packages("cowplot")
#install.packages("wesanderson")
#install.packages("meta")
#install.packages("rmeta")
#install.packages("openxlsx")
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
library(openxlsx)




##################################################### PART 1 - PRACTICAL GWAS #####################################################


genus.Allisonella.id.2174 <- data.table::fread("./Data/betas_withMAF_updated/genus.Allisonella.id.2174.summaryMAF.txt.gz") ## European-only cohorts included in the meta-analysis
genus.Allisonella.id.2174_rs602075 <- genus.Allisonella.id.2174[genus.Allisonella.id.2174$rsID == "rs602075", ]

MiBioGen_EUR_top_hits <- as.data.frame(genus.Allisonella.id.2174_rs602075)

dim(MiBioGen_EUR_top_hits) # 27 of the 30 top hits available 
write.table(MiBioGen_EUR_top_hits, "./Data/MiBioGen_gw_top_hits_in_EUR_only_data.txt", sep = "\t", col.names = T, row.names = F, quote = F)

############################################# 



# MiBioGen
MiBioGen_exposure_dat <- format_data(MiBioGen_EUR_top_hits, type = "exposure", snps = NULL, header = TRUE, phenotype_col = "bac", snp_col = "rsID", beta_col = "beta", se_col = "SE", effect_allele_col = "eff.allele", other_allele_col = "ref.allele",  pval_col = "P.weightedSumZ", samplesize_col = "N", min_pval = 1e-200, z_col = "Z.weightedSumZ", chr_col = "chr", pos_col = "bp", log_pval = FALSE, eaf_col = "weighted_MAF")

head(MiBioGen_exposure_dat)
MiBioGen_exposure_dat <- MiBioGen_exposure_dat %>% arrange(SNP)

MiBioGen_outcome_dat <- read_outcome_data("./Data/MiBioGen_prostate_outcome_data_main.txt", sep="\t", snp_col="SNP", beta_col="beta.outcome", se_col="se.outcome", effect_allele_col = "effect_allele.outcome", other_allele_col = "other_allele.outcome", pval_col = "pval.outcome", eaf_col = "eaf.outcome", phenotype_col = "outcome")
head(MiBioGen_outcome_dat)

MiBioGen_dat <- harmonise_data(MiBioGen_exposure_dat, MiBioGen_outcome_dat, action =1)
MiBioGen_mr_results <- mr(MiBioGen_dat, method_list=c("mr_wald_ratio", "mr_ivw"))
write.table(MiBioGen_mr_results, "./Output/MiBioGen_EUR_only_microbiome_prostate_cancer_results.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# Scatter plot
png(paste("./Figures/MiBioGen_EUR_only_MR_scatter.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_scatter_plot(MiBioGen_mr_results, MiBioGen_dat))
dev.off()

# Forest plot
MR_single <- mr_singlesnp(MiBioGen_dat)
png(paste("./Figures/MiBioGen_EUR_only_MR_forest.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_forest_plot(MR_single))
dev.off()

# Funnel plot
png(paste("./Figures/MiBioGen_EUR_only_MR_funnel.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_funnel_plot(MR_single))
dev.off()

# Leave one out analysis
MR_loo <- mr_leaveoneout(MiBioGen_dat)
png(paste("./Figures/MiBioGen_EUR_only_MR_loo.png", sep=""), width = 2400, height = 2000, res = 300)
print(mr_leaveoneout_plot(MR_loo))
dev.off()


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

write.table(MiBioGen_steiger, "./Output/MiBioGen_EUR_only_microbiome_prostate_cancer_steiger_filtering_F_stat.txt", sep = "\t", col.names = T, row.names = F, quote = F)

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
ggsave("./Figures/forestplot_microbiome_prostate_cancer_MiBioGen_EUR_only.png", width=7, height=8)


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
ggsave("./Figures/volcanoplot_microbiome_prostate_cancer_MiBioGen_EUR_only.png", width=10, height=8)


##################################################### PART 2 - UK Biobank #####################################################


###############################################
## EXTRACTING OUTCOME DATA FROM IEU OPENGWAS ##
###############################################
## Saving outcome data from IEU OpenGWAS 
# ieu-b-4809 = prostate cancer in Uk Biobank 

## Set up vector containing all outcomes
outcomes <- c("ieu-b-4809") #try extract instruments

## Extract the outcomes

MiBioGen_UKBB_outcome_dat <- extract_outcome_data(MiBioGen_exposure_dat$SNP, outcomes = outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, splitsize = 10000, proxy_splitsize = 500)
head(MiBioGen_UKBB_outcome_dat) 
dim(MiBioGen_UKBB_outcome_dat) # 20 MiBioGen SNPs in UK Biobank data

#########################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND METABS (OUTCOME) ##
#########################################################################
## MAIN ANALYSIS INCLUDING THE TOP HITS FROM OUR PUBLISHED META-ANALYSES

MiBioGen_UKBB_dat <- harmonise_data(MiBioGen_exposure_dat, MiBioGen_UKBB_outcome_dat, action =1)
MiBioGen_UKBB_mr_results <- mr(MiBioGen_UKBB_dat, method_list=c("mr_wald_ratio", "mr_ivw"))
write.table(MiBioGen_UKBB_mr_results, "./Output/MiBioGen_microbiome_UKBB_results.txt", sep = "\t", col.names = T, row.names = F, quote = F)

####### MiBioGen 
MiBioGen_UKBB_dat$units.exposure <- "SD"

MiBioGen_UKBB_dat$ncase.outcome <- 9132
MiBioGen_UKBB_dat$ncontrol.outcome <- 173493
MiBioGen_UKBB_dat$prevalence.outcome <- 0.1
MiBioGen_UKBB_dat$units.outcome <- "log odds"

steiger_filtering(MiBioGen_UKBB_dat)
write.table(steiger_filtering(MiBioGen_UKBB_dat), "./Output/MiBioGen_EUR_Only_microbiome_UKBB_steiger_filtering.txt", sep = "\t", col.names = T, row.names = F, quote = F)


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
ggsave("./Figures/forestplot_microbiome_UKBB_MiBioGen_EUR_only.png", width=7, height=8)


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
ggsave("./Figures/volcanoplot_microbiome_UKBB_prostate_MiBioGen_EUR_only.png", width=10, height=8)


######################. Part 3 - Finngen #######################################

## Saving outcome data from IEU OpenGWAS 
# finn-b-C3_BRONCHUS_LUNG_EXALLC = Malignant neoplasm of bronchus and lung (all cancers excluded) in Finngen

outcomes <- c("finn-b-C3_PROSTATE_EXALLC")
## Extract the outcomes

MiBioGen_finngen_outcome_dat <- extract_outcome_data(MiBioGen_exposure_dat$SNP, outcomes = outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, splitsize = 10000, proxy_splitsize = 500)
head(MiBioGen_finngen_outcome_dat) 

MiBioGen_finngen_outcome_dat = unique(MiBioGen_finngen_outcome_dat)
#########################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND METABS (OUTCOME) ##
#########################################################################
## MAIN ANALYSIS INCLUDING THE TOP HITS FROM OUR PUBLISHED META-ANALYSES

MiBioGen_finngen_dat <- harmonise_data(MiBioGen_exposure_dat, MiBioGen_finngen_outcome_dat, action =1)
MiBioGen_finngen_mr_results <- mr(MiBioGen_finngen_dat, method_list=c("mr_wald_ratio", "mr_ivw"))
write.table(MiBioGen_finngen_mr_results, "./Output/MiBioGen_EUR_only_microbiome_FinnGen_results.txt", sep = "\t", col.names = T, row.names = F, quote = F)

####### MiBioGen 
MiBioGen_finngen_dat$units.exposure <- "SD"

MiBioGen_finngen_dat$ncase.outcome <- 6311
MiBioGen_finngen_dat$ncontrol.outcome <- 74685
MiBioGen_finngen_dat$prevalence.outcome <- 0.1
MiBioGen_finngen_dat$units.outcome <- "log odds"

MiBioGen_finngen_steiger <- steiger_filtering(MiBioGen_finngen_dat)
write.table(MiBioGen_finngen_steiger, "./Output/MiBioGen_EUR_only_microbiome_finngen_steiger_filtering.txt", sep = "\t", col.names = T, row.names = F, quote = F)

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
ggsave("./Figures/forestplot_microbiome_FinnGen_MiBioGen_EUR_only.png", width=7, height=8)

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
ggsave("./Figures/volcanoplot_microbiome_finngen_prostate_MiBioGen_EUR_only.png", width=10, height=8)

##################################################### META-ANALYSING TRICL, UK Biobank and Finngen #####################################################
## Settings
colours <- names(wes_palettes)
discrete_palette <- wes_palette(colours[8], type = "discrete")
psignif <- 0
ci <- 0.95

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
MiBioGen_plot_data <- rbind(MiBioGen_main, MiBioGen_UKBB_sens, MiBioGen_finngen_sens)

## Now limit to the ones that are consistent across both datasets
MiBioGen_plot_data <- MiBioGen_plot_data[ave(1:nrow(MiBioGen_plot_data), MiBioGen_plot_data$exposure, FUN=length)>1,]


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
ggsave("./Figures/MiBioGen_EUR_only_microbiome_PRACTICAL_UKBB_finngen.png", width=10, height=8)



## Simple meta-analysis with MiBioGen
MiBioGen_meta <- metagen(MiBioGen_plot_data$beta, seTE = MiBioGen_plot_data$se, studlab = MiBioGen_plot_data$dataset, data = MiBioGen_plot_data, sm = "OR", fixed = TRUE, random = FALSE)
MiBioGen_meta_by_exposure <- update(MiBioGen_meta, subgroup = MiBioGen_plot_data$exposure, tau.common = FALSE)
summary(MiBioGen_meta_by_exposure)
pdf("./Figures/MiBioGen_EUR_only_microbiome_PRACTICAL+UKBB+finngen_meta.pdf", width=8, height = 16)
forest(MiBioGen_meta_by_exposure, studlab = T, common = TRUE, random = FALSE, leftcols = c("dataset"), print.subgroup.name = FALSE, print.subgroup.labels = TRUE, test.subgroup= FALSE, overall = FALSE, overall.hetstat = FALSE, text.common = "Meta-analysis", prediction = FALSE, xlab = "OR of prostate cancer with gut microbiome variation", print.tau2 = FALSE, plotwidth = "15cm", fontsize = 6, spacing = .36)
dev.off()



#### only the exposures to follow up 
exposures_to_select <- c("genus.Allisonella.id.2174")
traits_to_follow_up <- MiBioGen_plot_data[MiBioGen_plot_data$exposure %in% exposures_to_select, ]
## Simple meta-analysis with MiBioGen
MiBioGen_meta <- metagen(traits_to_follow_up$beta, seTE = traits_to_follow_up$se, studlab = traits_to_follow_up$dataset, data = traits_to_follow_up, sm = "OR", fixed = TRUE, random = FALSE)
MiBioGen_meta_by_exposure <- update(MiBioGen_meta, subgroup = traits_to_follow_up$exposure, tau.common = FALSE)
summary(MiBioGen_meta_by_exposure)
pdf("./Figures/MiBioGen_top_MR_traits_EUR_only_microbiome_PRACTICAL+UKBB+finngen_meta.pdf", width=8, height = 3)
forest(MiBioGen_meta_by_exposure, studlab = T, common = TRUE, random = FALSE, leftcols = c("dataset"), print.subgroup.name = FALSE, print.subgroup.labels = TRUE, test.subgroup= FALSE, overall = FALSE, overall.hetstat = FALSE, text.common = "Meta-analysis", prediction = FALSE, xlab = "OR of prostate cancer with gut microbiome variation", print.tau2 = FALSE, plotwidth = "15cm", fontsize = 6, spacing = .36)
dev.off()


#### creating plots for just the traits we want to follow up

####################################
## ADDING FOREST PLOTS OF RESULTS ##
####################################


MiBioGen_fres <- MiBioGen_fres[MiBioGen_fres$exposure=="genus.Allisonella.id.2174",]

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
  labs(title="Effect of the gut microbiome on endometrial cancer", x="Odds ratio", y = "Microbial trait (MiBioGen)") +
  geom_vline(xintercept=1, color="black", linetype="dashed", alpha=.5) +
  theme_classic()
ggsave("./Figures/forestplot_microbiome_prostate_cancer_MiBioGen_EUR_only_traits_to_follow_up.png", width=7, height=8)


