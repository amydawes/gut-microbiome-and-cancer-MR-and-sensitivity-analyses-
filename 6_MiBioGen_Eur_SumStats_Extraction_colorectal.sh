#!/bin/bash
#SBATCH --job-name=extraction
#SBATCH --output=extraction_output
#SBATCH --error=extraction_error
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --ntasks-per-node=1
#SBATCH --account=xxxxxxxx
#--------------------------

######################################################################################################################
## PROJECT: MR analysis of the gut microbiome on prostate cancer using PRACTICAL, UK Biobank, FinnGen
## Script: Extraction of MiBioGen microbiome GWAS results based on a defined set of SNPs
## Date: 25/09/26
######################################################################################################################
## Working on Blue Pebble
## Set working directory
cd /path/to/directory/prostate_cancer

## Made a flat text file called "prostate_cancer_snps.txt" with the rsids to extract and copied over to Blue Pebble
dos2unix prostate_cancer_snps.txt
head prostate_cancer_snps.txt

## Define SNPLIST that is in your working directory
SNPLIST=/path/to/directory/prostate_cancer/prostate_cancer_snps.txt

## Define the directory that holds the FGFP GWAS Sum.Stats. files
DATADIR=/directory/to/mibiogen_data/mibiogen/summary/mbqtl/raw/meta_alec/EUR/betafilewMAF

## Make a directory in the result folder to hold the data
mkdir results

## Define the directory to place the results
OUTDIR=/path/to/directory/prostate_cancer/results/

	
############################
## TEST the FIRST 3 files ##
############################
## STEP 1
ls -1 ${DATADIR} | grep txt.gz | head -n 3 | while read filename; do echo ${filename}; zcat ${DATADIR}${filename} | head -n 1  > ${OUTDIR}${filename%_allchr.txt.gz}.txt; done

## STEP 2
ls -1 ${DATADIR} | grep txt.gz | head -n 3 | while read filename; do echo ${filename}; TAXA=${filename%_allchr.txt.gz}; zcat ${DATADIR}${filename} | awk 'FNR==NR {a[$1]; next} FNR> 1 && $1 in a' ${SNPLIST} - >> ${OUTDIR}${filename%_allchr.txt.gz}.txt; done

## Check these output files
wc -l ${OUTDIR}/*.txt

## Remove these just in case
rm ${OUTDIR}/*.txt


########################
##  Run all the files ##
########################
## STEP 1
ls -1 ${DATADIR} | grep txt.gz | while read filename; do echo ${filename}; zcat ${DATADIR}${filename} | head -n 1  > ${OUTDIR}${filename%_allchr.txt.gz}.txt; done

## STEP 2
ls -1 ${DATADIR} | grep txt.gz | while read filename; do echo ${filename}; TAXA=${filename%_allchr.txt.gz}; zcat ${DATADIR}${filename} | awk 'FNR==NR {a[$1]; next} FNR> 1 && $1 in a' ${SNPLIST} - >> ${OUTDIR}${filename%_allchr.txt.gz}.txt; done

## Check these output files
wc -l ${OUTDIR}/*.txt

###################################################### END ######################################################
