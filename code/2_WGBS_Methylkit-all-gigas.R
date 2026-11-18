# Title: WGBS analysis.R
# Author: Matthew George; mattgeorgephd@gmail.com
# Date: 03/2021

rm(list=ls()) # clear all

## First time? - install these ####################################################################
# install.packages("readr")
# install.packages("tidyverse")
# install.packages("devtools")
# install.packages("vegan")
# install.packages("pheatmap")
# install.packages("gplots")
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install(version = "3.12")
# BiocManager::install(c("GenomicFeatures", "AnnotationDbi","methylKit"))
# BiocManager::install("GenomeInfoDb", force = TRUE)
# browseVignettes("methylKit") #methylKit manual

## Setup work space #############################################################################

# Load R packages
library(rstudioapi)
require(tidyverse)
require(devtools)
require(vegan)
require(pheatmap) 
require(gplots)
require(methylKit)
require(readr)
library(readxl)
library(ggplot2)
library(RColorBrewer)

# Grab the WD from the file location and set it
library(rstudioapi)
home <- getActiveDocumentContext()$path
setwd(dirname(home)); getwd(); setwd('..')
home <- getwd(); getwd()

## set theme

my_theme <- theme(line              = element_line(size=1.5),
                  rect              = element_rect(size=1.5),
                  text              = element_text(size=16,color="black"),
                  panel.background  = element_blank(),
                  panel.grid.major  = element_blank(), 
                  panel.grid.minor  = element_blank(),
                  axis.text.x       = element_text(size=16,color="black"),
                  axis.text.y       = element_text(size=16,color="black"),
                  axis.title.x      = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                  axis.title.y      = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.ticks.x      = element_line(color="black"),
                  axis.ticks.y      = element_line(color="black"),
                  # axis.line         = element_line(color = "black", size = 0.1),
                  panel.border      = element_rect(color = "black", fill=NA, size=1.5),
                  legend.key        = element_blank()) # removes background of legend bullets

## Load previous session data #################################################################################################

setwd("E:/WGBS-gigas-ploidy-desiccation")
# save.image("20211005_controls-only.RData") #Save R Data in case R crashes
# load("20211003_DML-ploidy.RData") #Load R Data
# save.image("20210923_methylKit_post-process.RData")
# load("20210923_methylKit_post-process.RData")

## Download files from Gannet ##################################################################################################

# ```{bash}
# wget -r -l1 --no-parent --no-directories -A.merged_CpG_evidence.cov https://gannet.fish.washington.edu/spartina/project-gigas-oa-meth/output/bismark-roslin/ #Download files from gannet to the current directory (where script is saved). Exclude directory structure
#   ```
# 
# ```{bash}
# mkdir ../output/06-methylKit #Make methylkit output directory
# ```
# 
# ```{bash}
# mv *.merged_CpG_evidence.cov ../output/06-methylKit #Move files from current directory to analyses folder
# ```


## Find files, generate list, methRead, filter, calculateDiffMeth #######################

# Put all .cov files into a list for analysis
setwd(home); getwd() # reset working directory to home dir
meta <- read_excel("metadata.xlsx", sheet = "list", col_names = TRUE)
analysisFiles <- as.list(meta$path)

# Create data.frame with sample info
sampleMetadata <- data.frame("sampleID"        =  meta$sampleID,
                             "ploidy"          =  meta$ploidy,
                             "heat"            =  meta$trt)
sampleMetadata$sampleID <- as.character(sampleMetadata$sampleID)

head(sampleMetadata) #Confirm data.frame creation


# STEP 1 - methRead: create a methylation object from the coverage files, and include sample ID and treatment information.
processedFiles_ploidy <- methylKit::methRead(analysisFiles,
                                      sample.id = as.list(sampleMetadata$sampleID),
                                      assembly = "Roslin",
                                      treatment = sampleMetadata$ploidy,
                                      pipeline = "bismarkCoverage",
                                      mincov = 2) #Process files. 0 = diploid, 1 = triploid.

processedFiles_heat  <- methylKit::methRead(analysisFiles,
                                             sample.id = as.list(sampleMetadata$sampleID),
                                             assembly = "Roslin",
                                             treatment = sampleMetadata$heat,
                                             pipeline = "bismarkCoverage",
                                             mincov = 2) #Process files. 0 = control, 1 = heat.

# STEP 2 - filter by 10x coverage
processedFilteredFilesCov10_ploidy <- methylKit::filterByCoverage(processedFiles_ploidy,
                                                          lo.count = 10, lo.perc = NULL,
                                                          high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)

processedFilteredFilesCov10_heat   <- methylKit::filterByCoverage(processedFiles_heat,
                                                                  lo.count = 10, lo.perc = NULL,
                                                                  high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)

# STEP 3 - unite/destrand(if applicable)
methylationInformationFilteredCov10_ploidy <- methylKit::unite(processedFilteredFilesCov10_ploidy,
                                                       destrand = FALSE,
                                                       mc.cores = 8)

methylationInformationFilteredCov10_heat <- methylKit::unite(processedFilteredFilesCov10_heat,
                                                        destrand = FALSE,
                                                        mc.cores = 8)

# head(methylationInformationFilteredCov10_ploidy)
nrow(methylationInformationFilteredCov10_ploidy) # 480579 total CpG

# head(methylationInformationFilteredCov10_heat)
nrow(methylationInformationFilteredCov10_heat) # 480579 total CpG

write.csv(methylationInformationFilteredCov10_ploidy, "DML/CpG/ALL-CpG-Cov10-all-gigas-ploidy.csv") #Save table as .csv
write.csv(methylationInformationFilteredCov10_heat,   "DML/CpG/ALL-CpG-Cov10-all-gigas-heat.csv") #Save table as .csv

# STEP 4 - Calculate differential methylation statistics (calculateDiffMeth)

covariatePLOIDY <- data.frame("PLOIDY" = c(rep("D", times = 5),
                                           rep("T", times = 5),
                                           rep("D", times = 6),
                                           rep("T", times = 6))) #Create dataframe with PLOIDY covariate information


covariateHEAT  <- data.frame("HEAT"    = c(rep("yes", times = 10),
                                           rep("no", times = 12))) #Create dataframe with HEAT covariate information


# calculateDiffMeth() - PLOIDY
differentialMethylationStatsPLOIDYCov10 <- methylKit::calculateDiffMeth(methylationInformationFilteredCov10_ploidy,
                                                                        covariates = covariateHEAT,
                                                                        mc.cores = 1) #Calculate differential methylation statistics based on treatment indication from methRead. Include HEAT as a covariate.
length(differentialMethylationStatsPLOIDYCov10$chr) # 480579

# calculateDiffMeth() - HEAT
differentialMethylationStatsHEATCov10 <- methylKit::calculateDiffMeth(methylationInformationFilteredCov10_heat,
                                                                      covariates = covariatePLOIDY,
                                                                      mc.cores = 1) #Calculate differential methylation statistics based on treatment indication from methRead. Include HEAT as a covariate.
length(differentialMethylationStatsHEATCov10$chr) # 480579

write.csv(differentialMethylationStatsPLOIDYCov10, "DML/ALL_gigas/DML-stats-PLOIDY-Cov10.csv") #Save table as .csv
write.csv(differentialMethylationStatsHEATCov10, "DML/ALL_gigas/DML-stats-HEAT-Cov10.csv") #Save table as .csv
write.csv(differentialMethylationStatsINTERACTIONCov10, "DML/ALL_gigas/DML-stats-INTERACTION-Cov10.csv") #Save table as .csv

## Calculate DML using % difference thresholds (getMethylDiff) #################

# PLOIDY
diffMethStatsPLOIDYCov10Diff5 <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 5, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff5$chr) #Count the number of DML --> 8319

diffMethStatsPLOIDYCov10Diff5Hyper <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 5, qvalue = 0.01, type = "hyper") #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff5Hyper$chr) #Count the number of DML --> 3496

diffMethStatsPLOIDYCov10Diff5Hypo <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 5, qvalue = 0.01,type = "hypo") #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff5Hypo$chr) #Count the number of DML --> 4823

diffMethStatsPLOIDYCov10Diff10 <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 10, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff10$chr) #Count the number of DML --> 5352

diffMethStatsPLOIDYCov10Diff10Hyper <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 10, qvalue = 0.01, type = "hyper") #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff10Hyper$chr) #Count the number of DML --> 2204

diffMethStatsPLOIDYCov10Diff10Hypo <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 10, qvalue = 0.01,type = "hypo") #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff10Hypo$chr) #Count the number of DML --> 3148

diffMethStatsPLOIDYCov10Diff15 <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 15, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff15$chr) #Count the number of DML --> 2872

diffMethStatsPLOIDYCov10Diff15Hyper <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 15, qvalue = 0.01, type = "hyper") #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff15Hyper$chr) #Count the number of DML --> 1228

diffMethStatsPLOIDYCov10Diff15Hypo <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 15, qvalue = 0.01,type = "hypo") #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff15Hypo$chr) #Count the number of DML --> 1644

diffMethStatsPLOIDYCov10Diff20 <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 20, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff20$chr) #Count the number of DML --> 1038

diffMethStatsPLOIDYCov10Diff20Hyper <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 20, qvalue = 0.01, type = "hyper") #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff20Hyper$chr) #Count the number of DML --> 442

diffMethStatsPLOIDYCov10Diff20Hypo <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 20, qvalue = 0.01,type = "hypo") #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff20Hypo$chr) #Count the number of DML --> 596

diffMethStatsPLOIDYCov10Diff25 <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 25, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff25$chr) #Count the number of DML --> 372

diffMethStatsPLOIDYCov10Diff25Hyper <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 25, qvalue = 0.01, type = "hyper") #Identify loci that are at least 25% different
length(diffMethStatsPLOIDYCov10Diff25Hyper$chr) #Count the number of DML --> 146

diffMethStatsPLOIDYCov10Diff25Hypo <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 25, qvalue = 0.01,type = "hypo") #Identify loci that are at least 25% different
length(diffMethStatsPLOIDYCov10Diff25Hypo$chr) #Count the number of DML --> 226

diffMethStatsPLOIDYCov10Diff30 <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 30, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff30$chr) #Count the number of DML --> 131

diffMethStatsPLOIDYCov10Diff30Hyper <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 30, qvalue = 0.01, type = "hyper") #Identify loci that are at least 30% different
length(diffMethStatsPLOIDYCov10Diff30Hyper$chr) #Count the number of DML --> 45

diffMethStatsPLOIDYCov10Diff30Hypo <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 30, qvalue = 0.01,type = "hypo") #Identify loci that are at least 30% different
length(diffMethStatsPLOIDYCov10Diff30Hypo$chr) #Count the number of DML --> 86

diffMethStatsPLOIDYCov10Diff35 <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 35, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff35$chr) #Count the number of DML --> 29

diffMethStatsPLOIDYCov10Diff35Hyper <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 35, qvalue = 0.01, type = "hyper") #Identify loci that are at least 35% different
length(diffMethStatsPLOIDYCov10Diff35Hyper$chr) #Count the number of DML --> 12

diffMethStatsPLOIDYCov10Diff35Hypo <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 35, qvalue = 0.01,type = "hypo") #Identify loci that are at least 35% different
length(diffMethStatsPLOIDYCov10Diff35Hypo$chr) #Count the number of DML --> 17

diffMethStatsPLOIDYCov10Diff40 <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 40, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff40$chr) #Count the number of DML --> 11

diffMethStatsPLOIDYCov10Diff40Hyper <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 40, qvalue = 0.01, type = "hyper") #Identify loci that are at least 40% different
length(diffMethStatsPLOIDYCov10Diff40Hyper$chr) #Count the number of DML --> 3

diffMethStatsPLOIDYCov10Diff40Hypo <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 40, qvalue = 0.01,type = "hypo") #Identify loci that are at least 40% different
length(diffMethStatsPLOIDYCov10Diff40Hypo$chr) #Count the number of DML --> 8

diffMethStatsPLOIDYCov10Diff45 <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 45, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsPLOIDYCov10Diff45$chr) #Count the number of DML --> 2

diffMethStatsPLOIDYCov10Diff45Hyper <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 45, qvalue = 0.01, type = "hyper") #Identify loci that are at least 45% different
length(diffMethStatsPLOIDYCov10Diff45Hyper$chr) #Count the number of DML --> 0

diffMethStatsPLOIDYCov10Diff45Hypo <- methylKit::getMethylDiff(differentialMethylationStatsPLOIDYCov10, difference = 45, qvalue = 0.01,type = "hypo") #Identify loci that are at least 45% different
length(diffMethStatsPLOIDYCov10Diff45Hypo$chr) #Count the number of DML --> 2

write.csv(diffMethStatsPLOIDYCov10Diff20,      "DML/ALL_gigas/DML-getMethylDiff-PLOIDY-Cov10-20.csv") #Save table as .csv
write.csv(diffMethStatsPLOIDYCov10Diff20Hyper, "DML/ALL_gigas/DML-getMethylDiff-PLOIDY-Cov10-20-Hyper.csv") #Save table as .csv
write.csv(diffMethStatsPLOIDYCov10Diff20Hypo,  "DML/ALL_gigas/DML-getMethylDiff-PLOIDY-Cov10-20-Hypo.csv") #Save table as .csv


# HEAT
diffMethStatsHEATCov10Diff5 <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 5, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff5$chr) #Count the number of DML --> 16125

diffMethStatsHEATCov10Diff5Hyper <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 5, qvalue = 0.01, type = "hyper") #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff5Hyper$chr) #Count the number of DML --> 8618

diffMethStatsHEATCov10Diff5Hypo <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 5, qvalue = 0.01,type = "hypo") #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff5Hypo$chr) #Count the number of DML --> 7507

diffMethStatsHEATCov10Diff10 <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 10, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff10$chr) #Count the number of DML --> 10747

diffMethStatsHEATCov10Diff10Hyper <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 10, qvalue = 0.01, type = "hyper") #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff10Hyper$chr) #Count the number of DML --> 5591

diffMethStatsHEATCov10Diff10Hypo <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 10, qvalue = 0.01,type = "hypo") #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff10Hypo$chr) #Count the number of DML --> 5156

diffMethStatsHEATCov10Diff15 <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 15, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff15$chr) #Count the number of DML --> 5910

diffMethStatsHEATCov10Diff15Hyper <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 15, qvalue = 0.01, type = "hyper") #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff15Hyper$chr) #Count the number of DML --> 2938

diffMethStatsHEATCov10Diff15Hypo <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 15, qvalue = 0.01,type = "hypo") #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff15Hypo$chr) #Count the number of DML --> 2972

diffMethStatsHEATCov10Diff20 <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 20, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff20$chr) #Count the number of DML --> 2351

diffMethStatsHEATCov10Diff20Hyper <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 20, qvalue = 0.01, type = "hyper") #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff20Hyper$chr) #Count the number of DML --> 1150

diffMethStatsHEATCov10Diff20Hypo <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 20, qvalue = 0.01,type = "hypo") #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff20Hypo$chr) #Count the number of DML --> 5156

diffMethStatsHEATCov10Diff25 <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 25, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff25$chr) #Count the number of DML --> 928

diffMethStatsHEATCov10Diff25Hyper <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 25, qvalue = 0.01, type = "hyper") #Identify loci that are at least 25% different
length(diffMethStatsHEATCov10Diff25Hyper$chr) #Count the number of DML --> 451

diffMethStatsHEATCov10Diff25Hypo <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 25, qvalue = 0.01,type = "hypo") #Identify loci that are at least 25% different
length(diffMethStatsHEATCov10Diff25Hypo$chr) #Count the number of DML --> 477

diffMethStatsHEATCov10Diff30 <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 30, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff30$chr) #Count the number of DML --> 359

diffMethStatsHEATCov10Diff30Hyper <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 30, qvalue = 0.01, type = "hyper") #Identify loci that are at least 30% different
length(diffMethStatsHEATCov10Diff30Hyper$chr) #Count the number of DML --> 176

diffMethStatsHEATCov10Diff30Hypo <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 30, qvalue = 0.01,type = "hypo") #Identify loci that are at least 30% different
length(diffMethStatsHEATCov10Diff30Hypo$chr) #Count the number of DML --> 183

diffMethStatsHEATCov10Diff35 <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 35, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff35$chr) #Count the number of DML --> 139

diffMethStatsHEATCov10Diff35Hyper <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 35, qvalue = 0.01, type = "hyper") #Identify loci that are at least 35% different
length(diffMethStatsHEATCov10Diff35Hyper$chr) #Count the number of DML --> 63

diffMethStatsHEATCov10Diff35Hypo <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 35, qvalue = 0.01,type = "hypo") #Identify loci that are at least 35% different
length(diffMethStatsHEATCov10Diff35Hypo$chr) #Count the number of DML --> 76

diffMethStatsHEATCov10Diff40 <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 40, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff40$chr) #Count the number of DML --> 48

diffMethStatsHEATCov10Diff40Hyper <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 40, qvalue = 0.01, type = "hyper") #Identify loci that are at least 40% different
length(diffMethStatsHEATCov10Diff40Hyper$chr) #Count the number of DML --> 24

diffMethStatsHEATCov10Diff40Hypo <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 40, qvalue = 0.01,type = "hypo") #Identify loci that are at least 40% different
length(diffMethStatsHEATCov10Diff40Hypo$chr) #Count the number of DML --> 24

diffMethStatsHEATCov10Diff45 <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 45, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff45$chr) #Count the number of DML --> 19

diffMethStatsHEATCov10Diff45Hyper <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 45, qvalue = 0.01, type = "hyper") #Identify loci that are at least 45% different
length(diffMethStatsHEATCov10Diff45Hyper$chr) #Count the number of DML --> 9

diffMethStatsHEATCov10Diff45Hypo <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 45, qvalue = 0.01,type = "hypo") #Identify loci that are at least 45% different
length(diffMethStatsHEATCov10Diff45Hypo$chr) #Count the number of DML --> 10

diffMethStatsHEATCov10Diff50 <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 50, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff50$chr) #Count the number of DML --> 4

diffMethStatsHEATCov10Diff50Hyper <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 50, qvalue = 0.01, type = "hyper") #Identify loci that are at least 50% different
length(diffMethStatsHEATCov10Diff50Hyper$chr) #Count the number of DML --> 3

diffMethStatsHEATCov10Diff50Hypo <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 50, qvalue = 0.01,type = "hypo") #Identify loci that are at least 50% different
length(diffMethStatsHEATCov10Diff50Hypo$chr) #Count the number of DML --> 1

write.csv(diffMethStatsHEATCov10Diff20, "DML/ALL_gigas/DML-getMethylDiff-HEAT-Cov10-20.csv") #Save table as .csv
write.csv(diffMethStatsHEATCov10Diff20Hyper, "DML/ALL_gigas/DML-getMethylDiff-HEAT-Cov10-20-Hyper.csv") #Save table as .csv
write.csv(diffMethStatsHEATCov10Diff20Hypo, "DML/ALL_gigas/DML-getMethylDiff-HEAT-Cov10-20-Hypo.csv") #Save table as .csv

# Interaction

# cov 10; 5% diff 
cov10diff5_HEAT   <- getData(diffMethStatsHEATCov10Diff5) # convert to df
cov10diff5_HEAT$row <- 1:nrow(cov10diff5_HEAT) # add a row counter

DML_INTERACTIONCov10Diff5 <- merge(x=getData(diffMethStatsPLOIDYCov10Diff5),
                                    y=cov10diff5_HEAT,
                                    by="row.names") # grab rows where start is shared by x and y
length(DML_INTERACTIONCov10Diff5$row) #Count the number of DML --> 2890

diffMethStatsINTERACTIONCov10Diff5 <- select(diffMethStatsHEATCov10Diff5, 
                                              DML_INTERACTIONCov10Diff5$row) # # have to use 'select' to subset a getMethylDiff() obj; select rows in diffMethStatsHEATCov10Diff5 by row counter
length(diffMethStatsINTERACTIONCov10Diff5$chr) # confirm length
sum(diffMethStatsINTERACTIONCov10Diff5$meth.diff > 0, na.rm=TRUE) #hyper
sum(diffMethStatsINTERACTIONCov10Diff5$meth.diff < 0, na.rm=TRUE) #hypo

# cov 10; 10% diff 
cov10diff10_HEAT   <- getData(diffMethStatsHEATCov10Diff10) # convert to df
cov10diff10_HEAT$row <- 1:nrow(cov10diff10_HEAT) # add a row counter

DML_INTERACTIONCov10Diff10 <- merge(x=getData(diffMethStatsPLOIDYCov10Diff10),
                                   y=cov10diff10_HEAT,
                                   by="row.names") # grab rows where start is shared by x and y
length(DML_INTERACTIONCov10Diff10$row) #Count the number of DML --> 1623

diffMethStatsINTERACTIONCov10Diff10 <- select(diffMethStatsHEATCov10Diff10, 
                                             DML_INTERACTIONCov10Diff10$row) # # have to use 'select' to subset a getMethylDiff() obj; select rows in diffMethStatsHEATCov10Diff5 by row counter
length(diffMethStatsINTERACTIONCov10Diff10$chr) # confirm length
sum(diffMethStatsINTERACTIONCov10Diff10$meth.diff > 0, na.rm=TRUE) #hyper
sum(diffMethStatsINTERACTIONCov10Diff10$meth.diff < 0, na.rm=TRUE) #hypo

# cov 10; 15% diff 
cov10diff15_HEAT   <- getData(diffMethStatsHEATCov10Diff15) # convert to df
cov10diff15_HEAT$row <- 1:nrow(cov10diff15_HEAT) # add a row counter

DML_INTERACTIONCov10Diff15 <- merge(x=getData(diffMethStatsPLOIDYCov10Diff15),
                                    y=cov10diff15_HEAT,
                                    by="row.names") # grab rows where start is shared by x and y
length(DML_INTERACTIONCov10Diff15$row) #Count the number of DML --> 738

diffMethStatsINTERACTIONCov10Diff15 <- select(diffMethStatsHEATCov10Diff15, 
                                              DML_INTERACTIONCov10Diff15$row) # # have to use 'select' to subset a getMethylDiff() obj; select rows in diffMethStatsHEATCov10Diff5 by row counter
length(diffMethStatsINTERACTIONCov10Diff15$chr) # confirm length
sum(diffMethStatsINTERACTIONCov10Diff15$meth.diff > 0, na.rm=TRUE) #hyper
sum(diffMethStatsINTERACTIONCov10Diff15$meth.diff < 0, na.rm=TRUE) #hypo

# cov 10; 20% diff 
cov10diff20_HEAT   <- getData(diffMethStatsHEATCov10Diff20) # convert to df
cov10diff20_HEAT$row <- 1:nrow(cov10diff20_HEAT) # add a row counter

DML_INTERACTIONCov10Diff20 <- merge(x=getData(diffMethStatsPLOIDYCov10Diff20),
                                    y=cov10diff20_HEAT,
                                    by="row.names") # grab rows where start is shared by x and y
length(DML_INTERACTIONCov10Diff20$row) #Count the number of DML --> 200

diffMethStatsINTERACTIONCov10Diff20 <- select(diffMethStatsHEATCov10Diff20, 
                                              DML_INTERACTIONCov10Diff20$row) # # have to use 'select' to subset a getMethylDiff() obj; select rows in diffMethStatsHEATCov10Diff20 by row counter
length(diffMethStatsINTERACTIONCov10Diff20$chr) # confirm length
sum(diffMethStatsINTERACTIONCov10Diff20$meth.diff > 0, na.rm=TRUE) #hyper
sum(diffMethStatsINTERACTIONCov10Diff20$meth.diff < 0, na.rm=TRUE) #hypo

# cov 10; 25% diff 
cov10diff25_HEAT   <- getData(diffMethStatsHEATCov10Diff25) # convert to df
cov10diff25_HEAT$row <- 1:nrow(cov10diff25_HEAT) # add a row counter

DML_INTERACTIONCov10Diff25 <- merge(x=getData(diffMethStatsPLOIDYCov10Diff25),
                                    y=cov10diff25_HEAT,
                                    by="row.names") # grab rows where start is shared by x and y
length(DML_INTERACTIONCov10Diff25$row) #Count the number of DML --> 50

diffMethStatsINTERACTIONCov10Diff25 <- select(diffMethStatsHEATCov10Diff25, 
                                              DML_INTERACTIONCov10Diff25$row) # # have to use 'select' to subset a getMethylDiff() obj; select rows in diffMethStatsHEATCov10Diff5 by row counter
length(diffMethStatsINTERACTIONCov10Diff25$chr) # confirm length
sum(diffMethStatsINTERACTIONCov10Diff25$meth.diff > 0, na.rm=TRUE) #hyper
sum(diffMethStatsINTERACTIONCov10Diff25$meth.diff < 0, na.rm=TRUE) #hypo

# cov 10; 30% diff 
cov10diff30_HEAT   <- getData(diffMethStatsHEATCov10Diff30) # convert to df
cov10diff30_HEAT$row <- 1:nrow(cov10diff30_HEAT) # add a row counter

DML_INTERACTIONCov10Diff30 <- merge(x=getData(diffMethStatsPLOIDYCov10Diff30),
                                    y=cov10diff30_HEAT,
                                    by="row.names") # grab rows where start is shared by x and y
length(DML_INTERACTIONCov10Diff30$row) #Count the number of DML --> 4

diffMethStatsINTERACTIONCov10Diff30 <- select(diffMethStatsHEATCov10Diff30, 
                                              DML_INTERACTIONCov10Diff30$row) # # have to use 'select' to subset a getMethylDiff() obj; select rows in diffMethStatsHEATCov10Diff5 by row counter
length(diffMethStatsINTERACTIONCov10Diff30$chr) # confirm length
sum(diffMethStatsINTERACTIONCov10Diff30$meth.diff > 0, na.rm=TRUE) #hyper
sum(diffMethStatsINTERACTIONCov10Diff30$meth.diff < 0, na.rm=TRUE) #hypo

# cov 10; 35% diff 
cov10diff35_HEAT   <- getData(diffMethStatsHEATCov10Diff35) # convert to df
cov10diff35_HEAT$row <- 1:nrow(cov10diff35_HEAT) # add a row counter

DML_INTERACTIONCov10Diff35 <- merge(x=getData(diffMethStatsPLOIDYCov10Diff35),
                                    y=cov10diff35_HEAT,
                                    by="row.names") # grab rows where start is shared by x and y
length(DML_INTERACTIONCov10Diff35$row) #Count the number of DML --> 1

diffMethStatsINTERACTIONCov10Diff35 <- select(diffMethStatsHEATCov10Diff35, 
                                              DML_INTERACTIONCov10Diff35$row) # # have to use 'select' to subset a getMethylDiff() obj; select rows in diffMethStatsHEATCov10Diff5 by row counter
length(diffMethStatsINTERACTIONCov10Diff35$chr) # confirm length
sum(diffMethStatsINTERACTIONCov10Diff35$meth.diff > 0, na.rm=TRUE) #hyper
sum(diffMethStatsINTERACTIONCov10Diff35$meth.diff < 0, na.rm=TRUE) #hypo

# cov 10; 40% diff 
cov10diff40_HEAT   <- getData(diffMethStatsHEATCov10Diff40) # convert to df
cov10diff40_HEAT$row <- 1:nrow(cov10diff40_HEAT) # add a row counter

DML_INTERACTIONCov10Diff40 <- merge(x=getData(diffMethStatsPLOIDYCov10Diff40),
                                    y=cov10diff40_HEAT,
                                    by="row.names") # grab rows where start is shared by x and y
length(DML_INTERACTIONCov10Diff40$row) #Count the number of DML --> 0

diffMethStatsINTERACTIONCov10Diff40 <- select(diffMethStatsHEATCov10Diff40, 
                                              DML_INTERACTIONCov10Diff40$row) # # have to use 'select' to subset a getMethylDiff() obj; select rows in diffMethStatsHEATCov10Diff5 by row counter
length(diffMethStatsINTERACTIONCov10Diff40$chr) # confirm length
sum(diffMethStatsINTERACTIONCov10Diff40$meth.diff > 0, na.rm=TRUE) #hyper
sum(diffMethStatsINTERACTIONCov10Diff40$meth.diff < 0, na.rm=TRUE) #hypo

# calculate hyper

diffMethStatsINTERACTIONCov10Diff20Hyper <- methylKit::getMethylDiff(diffMethStatsINTERACTIONCov10Diff20,
                                                                     difference = 20, 
                                                                     qvalue = 0.01,
                                                                     type = "hyper") #Identify loci that are at least 20% different
length(diffMethStatsINTERACTIONCov10Diff20Hyper$chr) #Count the number of DML --> 106
head(diffMethStatsINTERACTIONCov10Diff20Hyper) #Confirm creation

# calculate hypo
diffMethStatsINTERACTIONCov10Diff20Hypo <- methylKit::getMethylDiff(diffMethStatsINTERACTIONCov10Diff20,
                                                                    difference = 20, 
                                                                    qvalue = 0.01, 
                                                                    type = "hypo") #Identify loci that are at least 20% different
length(diffMethStatsINTERACTIONCov10Diff20Hypo$chr) #Count the number of DML --> 94
head(diffMethStatsINTERACTIONCov10Diff20Hypo) #Confirm creation

write.csv(diffMethStatsINTERACTIONCov10Diff20, "DML/ALL_gigas/DML-getMethylDiff-INTERACTION-Cov10-20.csv") #Save table as .csv
write.csv(diffMethStatsINTERACTIONCov10Diff20Hyper, "DML/ALL_gigas/DML-getMethylDiff-INTERACTION-Cov10-20-Hyper.csv") #Save table as .csv
write.csv(diffMethStatsINTERACTIONCov10Diff20Hypo, "DML/ALL_gigas/DML-getMethylDiff-INTERACTION-Cov10-20-Hypo.csv") #Save table as .csv



DMLPositions <- rep(0, times = length(diffMethStatsHEATCov10Diff20$chr)) #Create an empty vector
for (i in 1:length(DMLPositions)) {
  DMLPositions[i] <- which(cov10diff20_HEAT == DML_INTERACTIONCov10Diff20$chr)
}

diffMethPerChr(diffMethStatsPLOIDYCov10Diff20,
               plot=TRUE,
               qvalue.cutoff = 0.01,
               meth.cutoff = 20)

diffMethPerChr(diffMethStatsHEATCov10Diff20,
               plot=TRUE,
               qvalue.cutoff = 0.01,
               meth.cutoff = 20)

diffMethPerChr(diffMethStatsINTERACTIONCov10Diff20,
               plot=TRUE,
               qvalue.cutoff = 0.01,
               meth.cutoff = 20)


## Plot % difference ###################################

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path )); setwd(".."); 

diffplot     <- read_excel("analyses/DML_difference_threshold_comparison.xlsx", 
                           sheet = "plot", col_names = TRUE)
diffplot$trt <- factor(diffplot$trt,
                       levels=c("ploidy","desiccation","interaction"), ordered=TRUE)

p1 <- ggplot(diffplot, aes(x=difference,y=DML,color=trt)) +
      geom_point(size = 3, aes(shape = factor(trt))) +
      geom_line(size = 1.2) +
      scale_color_manual(values=c("royalblue2","orangered2","seagreen2")) +
      # geom_point() +
      scale_y_continuous(breaks = seq(0, 15000, 2500), 
                         limits = c(0, 17000),
                         labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
      # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) + 
      my_theme

p1

ggsave("plots/DML_diffthreshold_comparison.png",
       plot   = p1,
       dpi    = 1200,
       device = "png",
       width  = 6,
       height = 5,
       units  = "in")


p2 <- ggplot(diffplot, aes(x=difference,y=ratio_hyper,color=trt)) +
      geom_point(size = 3, aes(shape = factor(trt))) +
      geom_line(size = 1.2) +
      scale_color_manual(values=c("royalblue2","orangered2","seagreen2")) +
      # geom_point() +
      scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
      # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) + 
      my_theme

p2

ggsave("plots/DML_cov10_Hyper_ratio.png",
       plot   = p2,
       dpi    = 1200,
       device = "png",
       width  = 6,
       height = 5,
       units  = "in")

## SAVE POINT #########################################

save.image("R_sessions/2_WGBS_Methylkit-all-gigas-sp1.RData")
# load("20210923_methylKit_post-process.RData")


## DML distribution figure ############################

# Calculate distribution of DML across chromosomes
# Ploidy
DMLchrCounts_PLOIDY <-as.data.frame(table(diffMethStatsPLOIDYCov10Diff20$chr)) #Count the number of DML/chromosome
DMLchrCounts_PLOIDY <- DMLchrCounts_PLOIDY[-c(11:length(DMLchrCounts_PLOIDY$Var1)),]
prefix <- "Chr"
suffix <- seq(1:length(DMLchrCounts_PLOIDY$Var1))
my_names <- paste(prefix,suffix,sep="")
DMLchrCounts_PLOIDY$chr <- my_names
DMLchrCounts_PLOIDY <- DMLchrCounts_PLOIDY[,-1] #Remove column with chromosome RefSeq ID
DMLchrCounts_PLOIDY <- DMLchrCounts_PLOIDY[,c(2,1)] #Reorganize columns
colnames(DMLchrCounts_PLOIDY) <- c("chr", "DMLCount") #Rename columns
head(DMLchrCounts_PLOIDY) #Confirm formatting
write.csv(DMLchrCounts_PLOIDY, "DML/ALL_gigas/PLOIDY-DML-per-Chromosome.csv", row.names = FALSE, quote = FALSE) #Save file

# Ploidy - Hyper
DMLchrCounts_PLOIDYhyper <-as.data.frame(table(diffMethStatsPLOIDYCov10Diff20Hyper$chr)) #Count the number of DML/chromosome
DMLchrCounts_PLOIDYhyper <- DMLchrCounts_PLOIDYhyper[-c(11:length(DMLchrCounts_PLOIDYhyper$Var1)),]
prefix <- "Chr"
suffix <- seq(1:length(DMLchrCounts_PLOIDYhyper$Var1))
my_names <- paste(prefix,suffix,sep="")
DMLchrCounts_PLOIDYhyper$chr <- my_names
DMLchrCounts_PLOIDYhyper <- DMLchrCounts_PLOIDYhyper[,-1] #Remove column with chromosome RefSeq ID
DMLchrCounts_PLOIDYhyper <- DMLchrCounts_PLOIDYhyper[,c(2,1)] #Reorganize columns
colnames(DMLchrCounts_PLOIDYhyper) <- c("chr", "DMLCount") #Rename columns
head(DMLchrCounts_PLOIDYhyper) #Confirm formatting
write.csv(DMLchrCounts_PLOIDYhyper, "DML/ALL_gigas/PLOIDY-DML-per-Chromosome-Hyper.csv", row.names = FALSE, quote = FALSE) #Save file

# Ploidy - Hypo
DMLchrCounts_PLOIDYhypo <-as.data.frame(table(diffMethStatsPLOIDYCov10Diff20Hypo$chr)) #Count the number of DML/chromosome
DMLchrCounts_PLOIDYhypo <- DMLchrCounts_PLOIDYhypo[-c(11:length(DMLchrCounts_PLOIDYhypo$Var1)),]
prefix <- "Chr"
suffix <- seq(1:length(DMLchrCounts_PLOIDYhypo$Var1))
my_names <- paste(prefix,suffix,sep="")
DMLchrCounts_PLOIDYhypo$chr <- my_names
DMLchrCounts_PLOIDYhypo <- DMLchrCounts_PLOIDYhypo[,-1] #Remove column with chromosome RefSeq ID
DMLchrCounts_PLOIDYhypo <- DMLchrCounts_PLOIDYhypo[,c(2,1)] #Reorganize columns
colnames(DMLchrCounts_PLOIDYhypo) <- c("chr", "DMLCount") #Rename columns
head(DMLchrCounts_PLOIDYhypo) #Confirm formatting
write.csv(DMLchrCounts_PLOIDYhypo, "DML/ALL_gigas/PLOIDY-DML-per-Chromosome-Hypo.csv", row.names = FALSE, quote = FALSE) #Save file

# HEAT
DMLchrCounts_HEAT <-as.data.frame(table(diffMethStatsHEATCov10Diff20$chr)) #Count the number of DML/chromosome
DMLchrCounts_HEAT <- DMLchrCounts_HEAT[-c(11:length(DMLchrCounts_HEAT$Var1)),]
prefix <- "Chr"
suffix <- seq(1:length(DMLchrCounts_HEAT$Var1))
my_names <- paste(prefix,suffix,sep="")
DMLchrCounts_HEAT$chr <- my_names
DMLchrCounts_HEAT <- DMLchrCounts_HEAT[,-1] #Remove column with chromosome RefSeq ID
DMLchrCounts_HEAT <- DMLchrCounts_HEAT[,c(2,1)] #Reorganize columns
colnames(DMLchrCounts_HEAT) <- c("chr", "DMLCount") #Rename columns
head(DMLchrCounts_HEAT) #Confirm formatting
write.csv(DMLchrCounts_HEAT, "DML/ALL_gigas/HEAT-DML-per-Chromosome.csv", row.names = FALSE, quote = FALSE) #Save file

# HEAT - Hyper
DMLchrCounts_HEAThyper <-as.data.frame(table(diffMethStatsHEATCov10Diff20Hyper$chr)) #Count the number of DML/chromosome
DMLchrCounts_HEAThyper <- DMLchrCounts_HEAThyper[-c(11:length(DMLchrCounts_HEAThyper$Var1)),]
prefix <- "Chr"
suffix <- seq(1:length(DMLchrCounts_HEAThyper$Var1))
my_names <- paste(prefix,suffix,sep="")
DMLchrCounts_HEAThyper$chr <- my_names
DMLchrCounts_HEAThyper <- DMLchrCounts_HEAThyper[,-1] #Remove column with chromosome RefSeq ID
DMLchrCounts_HEAThyper <- DMLchrCounts_HEAThyper[,c(2,1)] #Reorganize columns
colnames(DMLchrCounts_HEAThyper) <- c("chr", "DMLCount") #Rename columns
head(DMLchrCounts_HEAThyper) #Confirm formatting
write.csv(DMLchrCounts_HEAThyper, "DML/ALL_gigas/HEAT-DML-per-Chromosome-Hyper.csv", row.names = FALSE, quote = FALSE) #Save file

# HEAT - Hypo
DMLchrCounts_HEAThypo <-as.data.frame(table(diffMethStatsHEATCov10Diff20Hypo$chr)) #Count the number of DML/chromosome
DMLchrCounts_HEAThypo <- DMLchrCounts_HEAThypo[-c(11:length(DMLchrCounts_HEAThypo$Var1)),]
prefix <- "Chr"
suffix <- seq(1:length(DMLchrCounts_HEAThypo$Var1))
my_names <- paste(prefix,suffix,sep="")
DMLchrCounts_HEAThypo$chr <- my_names
DMLchrCounts_HEAThypo <- DMLchrCounts_HEAThypo[,-1] #Remove column with chromosome RefSeq ID
DMLchrCounts_HEAThypo <- DMLchrCounts_HEAThypo[,c(2,1)] #Reorganize columns
colnames(DMLchrCounts_HEAThypo) <- c("chr", "DMLCount") #Rename columns
head(DMLchrCounts_HEAThypo) #Confirm formatting
write.csv(DMLchrCounts_HEAThypo, "DML/ALL_gigas/HEAT-DML-per-Chromosome-Hypo.csv", row.names = FALSE, quote = FALSE) #Save file

# INTERACTION
DMLchrCounts_INTERACTION <-as.data.frame(table(diffMethStatsINTERACTIONCov10Diff20$chr)) #Count the number of DML/chromosome
DMLchrCounts_INTERACTION <- DMLchrCounts_INTERACTION[-c(11:length(DMLchrCounts_INTERACTION$Var1)),]
prefix <- "Chr"
suffix <- seq(1:length(DMLchrCounts_INTERACTION$Var1))
my_names <- paste(prefix,suffix,sep="")
DMLchrCounts_INTERACTION$chr <- my_names
DMLchrCounts_INTERACTION <- DMLchrCounts_INTERACTION[,-1] #Remove column with chromosome RefSeq ID
DMLchrCounts_INTERACTION <- DMLchrCounts_INTERACTION[,c(2,1)] #Reorganize columns
colnames(DMLchrCounts_INTERACTION) <- c("chr", "DMLCount") #Rename columns
head(DMLchrCounts_INTERACTION) #Confirm formatting
write.csv(DMLchrCounts_INTERACTION, "DML/ALL_gigas/INTERACTION-DML-per-Chromosome.csv", row.names = FALSE, quote = FALSE) #Save file

# INTERACTION - Hyper
DMLchrCounts_INTERACTIONHyper <-as.data.frame(table(diffMethStatsINTERACTIONCov10Diff20Hyper$chr)) #Count the number of DML/chromosome
DMLchrCounts_INTERACTIONHyper <- DMLchrCounts_INTERACTIONHyper[-c(11:length(DMLchrCounts_INTERACTIONHyper$Var1)),]
prefix <- "Chr"
suffix <- seq(1:length(DMLchrCounts_INTERACTIONHyper$Var1))
my_names <- paste(prefix,suffix,sep="")
DMLchrCounts_INTERACTIONHyper$chr <- my_names
DMLchrCounts_INTERACTIONHyper <- DMLchrCounts_INTERACTIONHyper[,-1] #Remove column with chromosome RefSeq ID
DMLchrCounts_INTERACTIONHyper <- DMLchrCounts_INTERACTIONHyper[,c(2,1)] #Reorganize columns
colnames(DMLchrCounts_INTERACTIONHyper) <- c("chr", "DMLCount") #Rename columns
head(DMLchrCounts_INTERACTIONHyper) #Confirm formatting
write.csv(DMLchrCounts_INTERACTIONHyper, "DML/ALL_gigas/INTERACTIONHyper-DML-per-Chromosome.csv", row.names = FALSE, quote = FALSE) #Save file

# INTERACTION - Hypo
DMLchrCounts_INTERACTIONHypo <-as.data.frame(table(diffMethStatsINTERACTIONCov10Diff20Hypo$chr)) #Count the number of DML/chromosome
DMLchrCounts_INTERACTIONHypo <- DMLchrCounts_INTERACTIONHypo[-c(11:length(DMLchrCounts_INTERACTIONHypo$Var1)),]
prefix <- "Chr"
suffix <- seq(1:length(DMLchrCounts_INTERACTIONHypo$Var1))
my_names <- paste(prefix,suffix,sep="")
DMLchrCounts_INTERACTIONHypo$chr <- my_names
DMLchrCounts_INTERACTIONHypo <- DMLchrCounts_INTERACTIONHypo[,-1] #Remove column with chromosome RefSeq ID
DMLchrCounts_INTERACTIONHypo <- DMLchrCounts_INTERACTIONHypo[,c(2,1)] #Reorganize columns
colnames(DMLchrCounts_INTERACTIONHypo) <- c("chr", "DMLCount") #Rename columns
head(DMLchrCounts_INTERACTIONHypo) #Confirm formatting
write.csv(DMLchrCounts_INTERACTIONHypo, "DML/ALL_gigas/INTERACTIONHypo-DML-per-Chromosome.csv", row.names = FALSE, quote = FALSE) #Save file

# Plot stacked bar plot
library(ggplot2)

DMLdist            <- read_excel("DML/ALL_gigas/DML-gene-distribution/DML-gene-distribution.xlsx", sheet = "plot", col_names = TRUE)
DMLdist$chr    <- factor(DMLdist$chr_num, levels=c("1","2","3","4","5","6","7","8","9","10"),ordered=TRUE)
DMLdist$DML        <- factor(DMLdist$DML, levels=c("INTERACTION_hyper","INTERACTION_hypo",
                                               "HEAT_hyper","HEAT_hypo",
                                               "PLOIDY_hyper","PLOIDY_hypo"),ordered=TRUE)

## Set ggplot theme
my_theme <- theme(line              = element_line(size=1.5),
                  rect              = element_rect(size=1.5),
                  text              = element_text(size=14,color="black"),
                  panel.background  = element_blank(),
                  panel.grid.major  = element_blank(), 
                  panel.grid.minor  = element_blank(),
                  axis.text.x       = element_text(size=16,color="black"),
                  axis.text.y       = element_text(size=16,color="black"),
                  axis.title.x      = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                  axis.title.y      = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.ticks.x      = element_line(color="black"),
                  axis.ticks.y      = element_line(color="black"),
                  # axis.line         = element_line(color = "black", size = 0.1),
                  panel.border      = element_rect(color = "black", fill=NA, size=2),
                  legend.key        = element_blank()) # removes background of legend bullets


p1 <- ggplot(DMLdist) +
      geom_bar(aes(x=chr, y=count, fill=DML),color="black",position="stack",stat="identity") +
      scale_fill_manual(values=c("white","gray","orangered3","orangered1","royalblue3","royalblue1")) +
      # scale_color_manual(values="black") +
      geom_line(aes(x=chr_num, y=genes/5),size=1.2) +
      geom_point(aes(x=chr_num, y=genes/5),size=4) +
      # facet_wrap(~ID) +
      scale_y_continuous(breaks = seq(0, 1000, 250), limits = c(0, 850),sec.axis = sec_axis(~.*5)) +
      # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
      # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
      my_theme

p1

ggsave("plots/DML-gene-distribution.png",
       plot   = p1,
       dpi    = 1200,
       device = "png",
       width  = 9,
       height = 6,
       units  = "in")

## Principal Component Analyses #####################################

# All data
sample.IDs <- as.list(sampleMetadata$sampleID) #Create list of sample IDs
treatmentSpecification <- sampleMetadata$ploidy 

plotCustomization <- data.frame(sample = 1:22,
                                treatmentSpecification) #Create dataframe with sample treatment information
head(plotCustomization) #Confirm dataframe creation

allDataPCA <- PCASamples(methylationInformationFilteredCov10_ploidy, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA) #Look at summary statistics. The first PC explains 15.6% of variation, the second PC explains 12.1% of variation


RColorBrewer::display.brewer.all() #Show all RColorBrewer palettes. I will choose greens.
plotColors <- c("royalblue1","royalblue4","orangered1","orangered4","seagreen1","seagreen4")

pdf("DML/ALL_gigas/PCA-plots/All-gigas-ploidy-PCA.pdf", width = 11, height = 8.5)
par(mar = c(5, 5, 1, 1)) #Specify inner and outer margins
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.DMLDataPCA_ploidy, "sites", col = c(rep(plotColors[1], times = 5), 
                                               rep(plotColors[5], times = 5),
                                               rep(plotColors[1], times = 6),
                                               rep(plotColors[5], times = 6)), 
                                       pch = c(rep(16, times = 5), 
                                               rep(17, times = 5),
                                               rep(16, times = 6), 
                                               rep(17, times = 6)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2
#Add multiple white boxes on top of the default black box to manually change the color
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
ordiellipse(allDataPCA, plotCustomization$treatment, show.groups = "1", col = plotColors[4]) #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, plotCustomization$treatment, show.groups = "0", col = plotColors[2]) #Add confidence ellipse around the samples in ambient pCO2
axis(side =  1, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1 (14.25%)", line = 3, cex = 1.5) #Add x-axis label
axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2 (25.4%)", line = 3, cex = 1.5) #Add y-axis label
legend("topleft", 
       pch = c(16, 17), 
       legend = c("diploid", "triploid"), 
       col = c(plotColors[2], plotColors[4]), 
       cex = 1.7, bty = "n") #Add a legend with information about ambient and elevated samples
dev.off()

# Ploidy - DML

sample.IDs <- as.list(sampleMetadata$sampleID) #Create list of sample IDs
treatmentSpecification <- sampleMetadata$ploidyTreatment

plotCustomization <- data.frame(sample = 1:22,
                                sampleMetadata$ploidy) #Create dataframe with sample treatment information
head(plotCustomization) #Confirm dataframe creation

allDataPCA_ploidy <- PCASamples(methylationInformationFilteredCov10_ploidy, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA_ploidy) #Look at summary statistics. The first PC explains 9.5% of variation, the second PC explains 6.4% of variation

DMLPositions <- rep(0, times = length(diffMethStatsPLOIDYCov10Diff20$chr)) #Create an empty vector with 598 places to store row numbers
for (i in 1:length(DMLPositions)) {
  DMLPositions[i] <- which(getData(diffMethStatsPLOIDYCov10Diff20)$start[i] == getData(methylationInformationFilteredCov10_ploidy)$start)
} #For each DML, save the row number where that DML is found in methylationInformationFilteredCov5
tail(DMLPositions) #Confirm vector was created

DMLMatrix <- methylationInformationFilteredCov10_ploidy[DMLPositions,] #Subset methylationInformationFilteredCov5Destrand to only include DML and save as a new methylBase object
sum((DMLMatrix$start) == (diffMethStatsPLOIDYCov10Diff20$start)) == length(diffMethStatsPLOIDYCov10Diff20$start) #Confirm that start columns are identical. If they are identical, the sum of all TRUE statements should equal the length of the original methylBase object
tail(DMLMatrix) #Confirm methylBase object creation

DMLDataPCA_ploidy <- PCASamples(DMLMatrix, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(DMLDataPCA_ploidy) #Look at summary statistics. The first PC explains 24.0% of variation, the second PC explains 32.8% of variation


pdf("figures/DML_figure/All-gigas-ploidy-PCA.pdf", width = 8, height = 10)
par(mar = c(5, 5, 1, 1)) #Specify inner and outer margins
fig.DMLDataPCA_ploidy <- ordiplot(DMLDataPCA_ploidy, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.DMLDataPCA_ploidy, "sites", col = c(rep(plotColors[1], times = 5), 
                                               rep(plotColors[3], times = 5),
                                               rep(plotColors[1], times = 6),
                                               rep(plotColors[3], times = 6)), 
                                pch = c(rep(16, times = 5), 
                                        rep(17, times = 5),
                                        rep(16, times = 6), 
                                        rep(17, times = 6)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2
#Add multiple white boxes on top of the default black box to manually change the color
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
ordiellipse(DMLDataPCA_ploidy, plotCustomization$sampleMetadata.ploidy, show.groups = "1", col = plotColors[3]) #Add confidence ellipse around the samples in heated
ordiellipse(DMLDataPCA_ploidy, plotCustomization$sampleMetadata.ploidy, show.groups = "0", col = plotColors[1]) #Add confidence ellipse around the samples in control
axis(side =  1, labels = TRUE, col = "black", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1 (24.0%)", line = 3, cex = 1.5) #Add x-axis label
axis(side =  2, labels = TRUE, col = "black", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2 (8.8%)", line = 3, cex = 1.5) #Add y-axis label
legend("topright", 
       pch = c(16, 17), 
       legend = c("diploid", "triploid"), 
       col = c(plotColors[1], plotColors[3]), 
       cex = 1.7, bty = "n") #Add a legend with information about ambient and elevated samples
dev.off()

# HEAT - DML

sample.IDs <- as.list(sampleMetadata$sampleID) #Create list of sample IDs
treatmentSpecification <- sampleMetadata$heatTreatment 

plotCustomization <- data.frame(sample = 1:22,
                                sampleMetadata$heat) #Create dataframe with sample treatment information
head(plotCustomization) #Confirm dataframe creation

allDataPCA <- PCASamples(methylationInformationFilteredCov10_heat, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA) #Look at summary statistics. The first PC explains 15.6% of variation, the second PC explains 12.1% of variation

DMLPositions <- rep(0, times = length(diffMethStatsHEATCov10Diff20$chr)) #Create an empty vector with 598 places to store row numbers
for (i in 1:length(DMLPositions)) {
  DMLPositions[i] <- which(getData(diffMethStatsHEATCov10Diff20)$start[i] == getData(methylationInformationFilteredCov10_heat)$start)
} #For each DML, save the row number where that DML is found in methylationInformationFilteredCov5
tail(DMLPositions) #Confirm vector was created

DMLMatrix <- methylationInformationFilteredCov10_heat[DMLPositions,] #Subset methylationInformationFilteredCov5Destrand to only include DML and save as a new methylBase object
sum((DMLMatrix$start) == (diffMethStatsHEATCov10Diff20$start)) == length(diffMethStatsHEATCov10Diff20$start) #Confirm that start columns are identical. If they are identical, the sum of all TRUE statements should equal the length of the original methylBase object
tail(DMLMatrix) #Confirm methylBase object creation

DMLDataPCA_HEAT <- PCASamples(DMLMatrix, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(DMLDataPCA_HEAT) #Look at summary statistics. The first PC explains 25.4% of variation, the second PC explains 7.8% of variation


pdf("figures/DML_figure/All-gigas-HEAT-PCA.pdf", width = 8, height = 10)
par(mar = c(5, 5, 1, 1)) #Specify inner and outer margins
fig.DMLDataPCA_HEAT <- ordiplot(DMLDataPCA_HEAT, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.DMLDataPCA_HEAT, "sites", col = c(rep(plotColors[1], times = 10), 
                                             rep(plotColors[3], times = 12)), 
                                             pch = c(rep(16, times = 10), 
                                                     rep(17, times = 12)), 
                                             cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2
#Add multiple white boxes on top of the default black box to manually change the color
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
ordiellipse(DMLDataPCA_HEAT, plotCustomization$sampleMetadata.heat, show.groups = "1", col = plotColors[1]) #Add confidence ellipse around the samples in triploid
ordiellipse(DMLDataPCA_HEAT, plotCustomization$sampleMetadata.heat, show.groups = "0", col = plotColors[3]) #Add confidence ellipse around the samples in diploid
axis(side =  1, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1 (25.4%)", line = 3, cex = 1.5) #Add x-axis label
axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2 (7.8%)", line = 3, cex = 1.5) #Add y-axis label
legend("topright", 
       pch = c(16, 17), 
       legend = c("control", "desiccation"), 
       col = c(plotColors[2], plotColors[3]), 
       cex = 1.7, bty = "n") #Add a legend with information about ambient and elevated samples
dev.off()

## Heat Map ######################################
# Ploidy

DMLPositions <- rep(0, times = length(diffMethStatsPLOIDYCov10Diff20$chr)) #Create an empty vector with 598 places to store row numbers
for (i in 1:length(DMLPositions)) {
  DMLPositions[i] <- which(getData(diffMethStatsPLOIDYCov10Diff20)$start[i] == getData(methylationInformationFilteredCov10_ploidy)$start)
} #For each DML, save the row number where that DML is found in methylationInformationFilteredCov5Destrand
tail(DMLPositions) #Confirm vector was created

DMLMatrix <- methylationInformationFilteredCov10_heat[DMLPositions,] #Subset methylationInformationFilteredCov5Destrand to only include DML and save as a new methylBase object
sum((DMLMatrix$start) == (diffMethStatsPLOIDYCov10Diff20$start)) == length(diffMethStatsPLOIDYCov10Diff20$start) #Confirm that start columns are identical. If they are identical, the sum of all TRUE statements should equal the length of the original methylBase object
tail(DMLMatrix) #Confirm methylBase object creation

percMethDML <- percMethylation(DMLMatrix, rowids = TRUE) #Get percent methylation for all samples at DML. Include row IDs (chr, start, end) information
head(percMethDML) #Confirm percent methylation matrix was created

p1 <- pheatmap(percMethDML, color = rev(brewer.pal(8,'RdYlBu')),
         cluster_rows = TRUE, 
         clustering_distance_rows = "manhattan", # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
         clustering_method = "complete",
         treeheight_row = 100, 
         show_rownames = FALSE,
         cluster_cols = TRUE, 
         clustering_distance_cols = "manhattan", 
         treeheight_col = 40, 
         show_colnames = FALSE,
         # annotation_col = data.frame(pCO2 = factor(rep(c("Ambient","Treatment"), each = 5))),
         # annotation_colors = list(pCO2 = c(Ambient = "black", Treatment = "black")), 
         annotation_legend = FALSE, annotation_names_col = FALSE,
         legend = TRUE) #Create heatmap using pheatmap using percMethDML and plotColors color scheme. Cluster rows and columns using euclidean distances. Adjust the dendogram tree heights and do not show any row or column names. Create a dataframe with treatment information using annotation_col. Use annotation_colors to indicate colors for treatment ("grey90") and ambient ("grey10") samples. Do not include an annotation_legend or name for annotations (annotatino_names_col). Include a legend.

p1

ggsave("figures/heatmap_ploidy.png",
       plot   = p1,
       dpi    = 1200,
       device = "png",
       width  = 10,
       height = 10,
       units  = "in")

# HEAT
DMLPositions <- rep(0, times = length(diffMethStatsHEATCov10Diff20$chr)) #Create an empty vector with 598 places to store row numbers
for (i in 1:length(DMLPositions)) {
  DMLPositions[i] <- which(getData(diffMethStatsHEATCov10Diff20)$start[i] == getData(methylationInformationFilteredCov10_heat)$start)
} #For each DML, save the row number where that DML is found in methylationInformationFilteredCov5Destrand
tail(DMLPositions) #Confirm vector was created

DMLMatrix <- methylationInformationFilteredCov10_heat[DMLPositions,] #Subset methylationInformationFilteredCov5Destrand to only include DML and save as a new methylBase object
sum((DMLMatrix$start) == (diffMethStatsHEATCov10Diff20$start)) == length(diffMethStatsHEATCov10Diff20$start) #Confirm that start columns are identical. If they are identical, the sum of all TRUE statements should equal the length of the original methylBase object
tail(DMLMatrix) #Confirm methylBase object creation

percMethDML <- percMethylation(DMLMatrix, rowids = TRUE) #Get percent methylation for all samples at DML. Include row IDs (chr, start, end) information
head(percMethDML) #Confirm percent methylation matrix was created

p2 <- pheatmap(percMethDML, color = rev(brewer.pal(8,'RdYlBu')),
               cluster_rows = TRUE, 
               clustering_distance_rows = "manhattan",
               clustering_method = "complete",
               treeheight_row = 100, 
               show_rownames = FALSE,
               cluster_cols = TRUE, 
               clustering_distance_cols = "manhattan", 
               treeheight_col = 40, 
               show_colnames = FALSE,
               # annotation_col = data.frame(pCO2 = factor(rep(c("Ambient","Treatment"), each = 5))),
               # annotation_colors = list(pCO2 = c(Ambient = "black", Treatment = "black")), 
               annotation_legend = FALSE, annotation_names_col = FALSE,
               legend = TRUE) #Create heatmap using pheatmap using percMethDML and plotColors color scheme. Cluster rows and columns using euclidean distances. Adjust the dendogram tree heights and do not show any row or column names. Create a dataframe with treatment information using annotation_col. Use annotation_colors to indicate colors for treatment ("grey90") and ambient ("grey10") samples. Do not include an annotation_legend or name for annotations (annotatino_names_col). Include a legend.


p2

ggsave("figures/heatmap_HEAT.png",
       plot   = p2,
       dpi    = 1200,
       device = "png",
       width  = 10,
       height = 10,
       units  = "in")


# INTERACTION
DMLPositions <- rep(0, times = length(diffMethStatsINTERACTIONCov10Diff20$chr)) #Create an empty vector with 598 places to store row numbers
for (i in 1:length(DMLPositions)) {
  DMLPositions[i] <- which(getData(diffMethStatsINTERACTIONCov10Diff20)$start[i] == getData(methylationInformationFilteredCov10_heat)$start)
} #For each DML, save the row number where that DML is found in methylationInformationFilteredCov5Destrand
tail(DMLPositions) #Confirm vector was created

DMLMatrix <- methylationInformationFilteredCov10_heat[DMLPositions,] #Subset methylationInformationFilteredCov5Destrand to only include DML and save as a new methylBase object
sum((DMLMatrix$start) == (diffMethStatsINTERACTIONCov10Diff20$start)) == length(diffMethStatsINTERACTIONCov10Diff20$start) #Confirm that start columns are identical. If they are identical, the sum of all TRUE statements should equal the length of the original methylBase object
tail(DMLMatrix) #Confirm methylBase object creation

percMethDML <- percMethylation(DMLMatrix, rowids = TRUE) #Get percent methylation for all samples at DML. Include row IDs (chr, start, end) information
head(percMethDML) #Confirm percent methylation matrix was created

p3 <- pheatmap(percMethDML, color = rev(brewer.pal(8,'RdYlBu')),
               cluster_rows = TRUE, 
               clustering_distance_rows = "manhattan",
               clustering_method = "complete",
               treeheight_row = 100, 
               show_rownames = FALSE,
               cluster_cols = TRUE, 
               clustering_distance_cols = "manhattan", 
               treeheight_col = 40, 
               show_colnames = FALSE,
               # annotation_col = data.frame(pCO2 = factor(rep(c("Ambient","Treatment"), each = 5))),
               # annotation_colors = list(pCO2 = c(Ambient = "black", Treatment = "black")), 
               annotation_legend = FALSE, annotation_names_col = FALSE,
               legend = TRUE) #Create heatmap using pheatmap using percMethDML and plotColors color scheme. Cluster rows and columns using euclidean distances. Adjust the dendogram tree heights and do not show any row or column names. Create a dataframe with treatment information using annotation_col. Use annotation_colors to indicate colors for treatment ("grey90") and ambient ("grey10") samples. Do not include an annotation_legend or name for annotations (annotatino_names_col). Include a legend.


p3

ggsave("figures/heatmap_INTERACTION.png",
       plot   = p3,
       dpi    = 1200,
       device = "png",
       width  = 10,
       height = 10,
       units  = "in")


## Convert DML list to Bed files #######################################

DML <- data.frame(diffMethStatsPLOIDYCov10Diff20$chr,
                  diffMethStatsPLOIDYCov10Diff20$start,
                  diffMethStatsPLOIDYCov10Diff20$end,
                  diffMethStatsPLOIDYCov10Diff20$meth.diff) # arrange(chr, start) %>% #Join + and - strand information to be saved as a BED file, and avoid writing information in scientific notation
write_delim(DML, "bed_files/DML-allgigas-Cov10-Diff20-PLOIDY.bed",  delim = '\t', col_names = FALSE) #Save data as a BED file

DML <- data.frame(diffMethStatsPLOIDYCov10Diff20Hyper$chr,
                  diffMethStatsPLOIDYCov10Diff20Hyper$start,
                  diffMethStatsPLOIDYCov10Diff20Hyper$end,
                  diffMethStatsPLOIDYCov10Diff20Hyper$meth.diff) # arrange(chr, start) %>% #Join + and - strand information to be saved as a BED file, and avoid writing information in scientific notation
write_delim(DML, "bed_files/DML-allgigas-Cov10-Diff20-PLOIDYHyper.bed",  delim = '\t', col_names = FALSE) #Save data as a BED file

DML <- data.frame(diffMethStatsPLOIDYCov10Diff20Hypo$chr,
                  diffMethStatsPLOIDYCov10Diff20Hypo$start,
                  diffMethStatsPLOIDYCov10Diff20Hypo$end,
                  diffMethStatsPLOIDYCov10Diff20Hypo$meth.diff) # arrange(chr, start) %>% #Join + and - strand information to be saved as a BED file, and avoid writing information in scientific notation
write_delim(DML, "bed_files/DML-allgigas-Cov10-Diff20-PLOIDYHypo.bed",  delim = '\t', col_names = FALSE) #Save data as a BED file

DML <- data.frame(diffMethStatsHEATCov10Diff20$chr,
                  diffMethStatsHEATCov10Diff20$start,
                  diffMethStatsHEATCov10Diff20$end,
                  diffMethStatsHEATCov10Diff20$meth.diff) # arrange(chr, start) %>% #Join + and - strand information to be saved as a BED file, and avoid writing information in scientific notation
write_delim(DML, "bed_files/DML-allgigas-Cov10-Diff20-HEAT.bed",  delim = '\t', col_names = FALSE) #Save data as a BED file

DML <- data.frame(diffMethStatsHEATCov10Diff20Hyper$chr,
                  diffMethStatsHEATCov10Diff20Hyper$start,
                  diffMethStatsHEATCov10Diff20Hyper$end,
                  diffMethStatsHEATCov10Diff20Hyper$meth.diff) # arrange(chr, start) %>% #Join + and - strand information to be saved as a BED file, and avoid writing information in scientific notation
write_delim(DML, "bed_files/DML-allgigas-Cov10-Diff20-HEATHyper.bed",  delim = '\t', col_names = FALSE) #Save data as a BED file

DML <- data.frame(diffMethStatsHEATCov10Diff20Hypo$chr,
                  diffMethStatsHEATCov10Diff20Hypo$start,
                  diffMethStatsHEATCov10Diff20Hypo$end,
                  diffMethStatsHEATCov10Diff20Hypo$meth.diff) # arrange(chr, start) %>% #Join + and - strand information to be saved as a BED file, and avoid writing information in scientific notation
write_delim(DML, "bed_files/DML-allgigas-Cov10-Diff20-HEATHypo.bed",  delim = '\t', col_names = FALSE) #Save data as a BED file

DML <- data.frame(diffMethStatsINTERACTIONCov10Diff20$chr,
                  diffMethStatsINTERACTIONCov10Diff20$start,
                  diffMethStatsINTERACTIONCov10Diff20$end,
                  diffMethStatsINTERACTIONCov10Diff20$meth.diff) # arrange(chr, start) %>% #Join + and - strand information to be saved as a BED file, and avoid writing information in scientific notation
write_delim(DML, "bed_files/DML-allgigas-Cov10-Diff20-INTERACTION.bed",  delim = '\t', col_names = FALSE) #Save data as a BED file

DML <- data.frame(diffMethStatsINTERACTIONCov10Diff20Hyper$chr,
                  diffMethStatsINTERACTIONCov10Diff20Hyper$start,
                  diffMethStatsINTERACTIONCov10Diff20Hyper$end,
                  diffMethStatsINTERACTIONCov10Diff20Hyper$meth.diff) # arrange(chr, start) %>% #Join + and - strand information to be saved as a BED file, and avoid writing information in scientific notation
write_delim(DML, "bed_files/DML-allgigas-Cov10-Diff20-INTERACTIONHyper.bed",  delim = '\t', col_names = FALSE) #Save data as a BED file

DML <- data.frame(diffMethStatsINTERACTIONCov10Diff20Hypo$chr,
                  diffMethStatsINTERACTIONCov10Diff20Hypo$start,
                  diffMethStatsINTERACTIONCov10Diff20Hypo$end,
                  diffMethStatsINTERACTIONCov10Diff20Hypo$meth.diff) # arrange(chr, start) %>% #Join + and - strand information to be saved as a BED file, and avoid writing information in scientific notation
write_delim(DML, "bed_files/DML-allgigas-Cov10-Diff20-INTERACTIONHypo.bed",  delim = '\t', col_names = FALSE) #Save data as a BED file


## DML-genome-feature-location Stacked bar plot #####################################

data         <- read_excel("bisulfide_analysis/WGBS/analyses/DML_locations-all-datasets.xlsx", sheet = "plot", col_names = TRUE)
data$feature <- factor(data$feature, levels = c("exon_UTR",
                                                "CDS",
                                                "intron",
                                                "upstream_flank",
                                                "downstream_flank",
                                                "lncRNA",
                                                "TE",
                                                "intergenic"), ordered=TRUE) 

# Plot
p1 <- ggplot(data = data, aes(y = B_combo_contrast, x = 'DML', fill = feature)) + 
  geom_bar(stat  = "identity",
           color = "black",
           #fill  = rep(v,length.out=18),
           size  = 1) +
  scale_fill_brewer(palette = "YlGnBu") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20)) +
  theme( line              = element_line(size=1.5),
         rect              = element_rect(size=1.5),
         text              = element_text(size=22,color="black"),
         panel.background  = element_blank(),
         panel.grid.major  = element_blank(), 
         panel.grid.minor  = element_blank(),
         axis.text.x       = element_text(size=22,color="black"),
         axis.text.y       = element_text(size=22,color="black"),
         # axis.ticks.y      = element_blank(),
         axis.title.x      = element_blank(),
         axis.title.y      = element_blank(),
         axis.line         = element_line(color = "black", size = 0.5),
         panel.border      = element_rect(color = "black", fill=NA, size=2))

p1

ggsave("bisulfide_analysis/WGBS/plots/B_combo_contrast.png",
       plot   = p1,
       dpi    = 1200,
       device = "png",
       width  = 4.8,
       height = 6,
       units  = "in")

## DML-genome-feature-location Stacked bar plot ########################################

library(hrbrthemes)

dd <- read_excel("DML_summary.xlsx", sheet = "all", col_names = TRUE)

a <- dd %>%
        group_by(ID) %>%
        filter(n() == 2)

data      <- read_excel("DML_summary.xlsx", sheet = "plot", col_names = TRUE)
data$type <- factor(data$type, levels = c("++","+-","-+","--"))

p <- ggplot(data, aes(x=num, y=ploidy_diff)) +
      geom_segment(
        aes(x=num, xend=num, y=0, yend=heat_diff), 
        color="orange", 
        size=ifelse(data$num %in% c("A","D"), 1.3, 0.7),
        position=position_dodge()
      ) +
      geom_point(aes(x=num,y=heat_diff),
        color="orange", 
        size=ifelse(data$num %in% c("A","D"), 5, 2),
        position=position_dodge()
      ) +
      geom_segment(
        aes(x=num, xend=num, y=0, yend=ploidy_diff), 
        color="blue", 
        size=ifelse(data$num %in% c("A","D"), 1.3, 0.7)
      ) +
      geom_point(aes(x=num,y=ploidy_diff),
                 color="blue", 
                 size=ifelse(data$num %in% c("A","D"), 5, 2)
      ) +
  # theme_ipsum() +
  coord_flip() +
  # scale_color_manual(values = c("red","green","blue","black")) +
  theme(line              = element_line(size=2),
        rect              = element_rect(size=2),
        text              = element_text(size=22,color="black"),
        panel.background  = element_blank(),
        panel.grid.major  = element_blank(),
        axis.ticks.x      = element_blank(),
        axis.ticks.y      = element_line(color = "black"),
        axis.text         = element_text(size=22,color="black"),
        panel.grid.minor  = element_blank(),
        axis.title.x      = element_blank(),
        axis.line         = element_blank(),
        panel.border      = element_rect(color = "black", fill=NA, size=2),
        legend.position   = 'none',
        legend.key        = element_rect(fill = 'white'))



p

ggsave("p.png",
       plot   = p,
       dpi    = 600,
       device = "png",
       width  = 7,
       height = 9,
       units  = "in")

##########################################################################################

setwd('E:/WGBS-gigas-ploidy-desiccation/DML_lists/genome_feature')

a <- read.csv('ALL-CpG-Cov10-Ronit.csv')
b <- read.csv('ALL-CpG-Cov10-Yaamini.csv')

joined_tibble <- left_join(b, a, 
                           by = "X")

c <- joined_tibble$start.x - joined_tibble$start.y

