##############################################################################################################################
## PROJECT: MR analysis of the gut microbiome on prostate cancer using PRACTICAL, UK Biobank and FinnGen
## Script: Extraction of FGFP microbiome GWAS results based on regions around defined set of SNPs (for coloc) on HPC
## Date: 25/09/26
##############################################################################################################################
## Clear space
rm(list=ls())

## Set working directory
setwd("/path/to/directory/prostate_cancer")

## Load the library data.table
library(data.table)
#install.packages("R.utils")
library(R.utils)
library(tidyr)
#################################################### MIBIOGEN ####################################################

## Read in the flat text file with the extract details (these include the SNPs around which we want to extract the summary-level data)
## This was copied over to the working directory
d <- read.table("./MiBioGen_coloc_snps.txt", header = T, fill = T, sep = "\t")
d <- as.data.frame(d)

## Rename columns
names(d) <- c("MT", "rsid", "chr", "base_pair")

d$MT <- as.character(d$MT)
#d$MT <- gsub("G_u_","G_unclassified_",d$MT)

## Directory with the data
# change directory of the data
data_dir <- "/path/to/directory/full_MiBioGen_all_ancestries_summary_statistics/"

## file names in data_dir
f <- list.files(data_dir)

## Iterate over taxa to extract data from
for(i in 1:nrow(d)){
  cat(paste0("Now processing mt number ", i,"\n"))
  ## parameters
  mt <- as.character( d[i,"MT"] )
  rsid <- as.character( d[i,"rsid"] )
  chr <- d[i,"chr"]
  bp <- d[i, "base_pair"]
  ### find the file I need
  w <- grep(mt, f)
  ## if 1 file found
  ## read in the data
  if(length(w) == 1){
    gwas_file_name <- paste0(data_dir, f[w])
    gwas_data <- data.table::fread(gwas_file_name)

    ## create chromosome and position columns 
gwas_data <- separate(gwas_data, col = "SNP", into = c("chromosome", "position"), sep = ":", remove = FALSE)

## Convert position to numeric if it's not already
gwas_data$position <- as.numeric(gwas_data$position)

### subset the data
out <- gwas_data[gwas_data$chromosome == chr & gwas_data$position >= (bp - 1e6) & gwas_data$position <= (bp + 1e6), ]

## write to file
out_name <- paste0("./results/", mt, "_", rsid, "_plusminus_1e6.txt")
cat(paste0("Now writing data for ", mt," to file.\n"))
write.table(out, file = out_name, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  }
  rm(gwas_data)
  rm(out)
}
    ### subset the data
  #  out <- gwas_data[chromosome=chr & position >= bp-1e6 & position <= bp+1e6]
    ## write to file
  #  out_name <- paste0("./results/", mt, "_", rsid, "_plusminus_1e6.txt")
  #  cat(paste0("Now writing data for ", mt," to file.\n"))
  #  write.table(out, file = out_name, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
 # }
 # rm(gwas_data)
 # rm(out)



