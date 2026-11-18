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

##Setup work space #############################################################################

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
home <- getActiveDocumentContext()$path
setwd(dirname(home)); getwd(); setwd('..')
home <- getwd(); getwd()

## set theme

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


## Find files and generate list ##################################################################################################

# Put all .cov files into a list for analysis
setwd(home); getwd() # reset working directory to home dir
meta <- read_excel("metadata.xlsx", sheet = "list", col_names = TRUE)

analysisFiles <- as.list(meta$path)
analysisFiles_diploid     <- as.list(c(meta$path[1:5],meta$path[11:16]))
analysisFiles_diploid_ID  <- as.list(as.character(c(meta$sampleID[1:5],meta$sampleID[11:16])))
analysisFiles_diploid_trt <- c(meta$trt[1:5],meta$trt[11:16])
analysisFiles_triploid     <- as.list(c(meta$path[6:10],meta$path[17:22]))
analysisFiles_triploid_ID  <- as.list(as.character(c(meta$sampleID[6:10],meta$sampleID[17:22])))
analysisFiles_triploid_trt <- c(meta$trt[6:10],meta$trt[17:22])

analysisFiles_control     <- as.list(c(meta$path[11:22]))
analysisFiles_control_ID  <- as.list(as.character(c(meta$sampleID[11:22])))
analysisFiles_control_trt <- c(meta$trt[11:22])
analysisFiles_heat     <- as.list(c(meta$path[1:10]))
analysisFiles_heat_ID  <- as.list(as.character(c(meta$sampleID[1:10])))
analysisFiles_heat_trt <- c(meta$trt[1:10])


##### Explore coverage and plot 

## STEP 1 - methRead
# all data
processedFiles <- methylKit::methRead(analysisFiles,
                                              sample.id = as.list(meta$sampleID),
                                              assembly = "Roslin",
                                              treatment = meta$ploidy,
                                              pipeline = "bismarkCoverage",
                                              mincov = 1) #Process files. x1 coverage
# diploid
processedFiles_diploid <- methylKit::methRead(analysisFiles_diploid,
                                      sample.id = analysisFiles_diploid_ID,
                                      assembly = "Roslin",
                                      treatment = analysisFiles_diploid_trt,
                                      pipeline = "bismarkCoverage",
                                      mincov = 1) #Process files. x1 coverage
# triploid
processedFiles_triploid <- methylKit::methRead(analysisFiles_triploid,
                                              sample.id = analysisFiles_triploid_ID,
                                              assembly = "Roslin",
                                              treatment = analysisFiles_triploid_trt,
                                              pipeline = "bismarkCoverage",
                                              mincov = 1) #Process files. x1 coverage
# control
processedFiles_control <- methylKit::methRead(analysisFiles_control,
                                              sample.id = analysisFiles_control_ID,
                                              assembly = "Roslin",
                                              treatment = analysisFiles_control_trt,
                                              pipeline = "bismarkCoverage",
                                              mincov = 1) #Process files. x1 coverage
# heat
processedFiles_heat <- methylKit::methRead(analysisFiles_heat,
                                               sample.id = analysisFiles_heat_ID,
                                               assembly = "Roslin",
                                               treatment = analysisFiles_heat_trt,
                                               pipeline = "bismarkCoverage",
                                               mincov = 1) #Process files. x1 coverage

## STEP 2 - filterByCoverage
# all data
processedFilteredFilesCov1 <- methylKit::filterByCoverage(processedFiles,
                                                          lo.count = 1, lo.perc = NULL,
                                                          high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov3 <- methylKit::filterByCoverage(processedFiles,
                                                          lo.count = 3, lo.perc = NULL,
                                                          high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov5 <- methylKit::filterByCoverage(processedFiles,
                                                          lo.count = 5, lo.perc = NULL,
                                                          high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov10 <- methylKit::filterByCoverage(processedFiles,
                                                           lo.count = 10, lo.perc = NULL,
                                                           high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov20 <- methylKit::filterByCoverage(processedFiles,
                                                           lo.count = 20, lo.perc = NULL,
                                                           high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)

# diploid
processedFilteredFilesCov1_diploid <- methylKit::filterByCoverage(processedFiles_diploid,
                                                          lo.count = 1, lo.perc = NULL,
                                                          high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov3_diploid <- methylKit::filterByCoverage(processedFiles_diploid,
                                                          lo.count = 3, lo.perc = NULL,
                                                          high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov5_diploid <- methylKit::filterByCoverage(processedFiles_diploid,
                                                          lo.count = 5, lo.perc = NULL,
                                                          high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov10_diploid <- methylKit::filterByCoverage(processedFiles_diploid,
                                                                   lo.count = 10, lo.perc = NULL,
                                                                   high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov20_diploid <- methylKit::filterByCoverage(processedFiles_diploid,
                                                                   lo.count = 20, lo.perc = NULL,
                                                                   high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
# triploid
processedFilteredFilesCov1_triploid <- methylKit::filterByCoverage(processedFiles_triploid,
                                                                  lo.count = 1, lo.perc = NULL,
                                                                  high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov3_triploid <- methylKit::filterByCoverage(processedFiles_triploid,
                                                                  lo.count = 3, lo.perc = NULL,
                                                                  high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov5_triploid <- methylKit::filterByCoverage(processedFiles_triploid,
                                                                  lo.count = 5, lo.perc = NULL,
                                                                  high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov10_triploid <- methylKit::filterByCoverage(processedFiles_triploid,
                                                                   lo.count = 10, lo.perc = NULL,
                                                                   high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov20_triploid <- methylKit::filterByCoverage(processedFiles_triploid,
                                                                   lo.count = 20, lo.perc = NULL,
                                                                   high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
# control
processedFilteredFilesCov1_control <- methylKit::filterByCoverage(processedFiles_control,
                                                                  lo.count = 1, lo.perc = NULL,
                                                                  high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov3_control <- methylKit::filterByCoverage(processedFiles_control,
                                                                  lo.count = 3, lo.perc = NULL,
                                                                  high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov5_control <- methylKit::filterByCoverage(processedFiles_control,
                                                                  lo.count = 5, lo.perc = NULL,
                                                                  high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov10_control <- methylKit::filterByCoverage(processedFiles_control,
                                                                   lo.count = 10, lo.perc = NULL,
                                                                   high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov20_control <- methylKit::filterByCoverage(processedFiles_control,
                                                                   lo.count = 20, lo.perc = NULL,
                                                                   high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
# heat
processedFilteredFilesCov1_heat <- methylKit::filterByCoverage(processedFiles_heat,
                                                                   lo.count = 1, lo.perc = NULL,
                                                                   high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov3_heat <- methylKit::filterByCoverage(processedFiles_heat,
                                                                   lo.count = 3, lo.perc = NULL,
                                                                   high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov5_heat <- methylKit::filterByCoverage(processedFiles_heat,
                                                                   lo.count = 5, lo.perc = NULL,
                                                                   high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov10_heat <- methylKit::filterByCoverage(processedFiles_heat,
                                                                    lo.count = 10, lo.perc = NULL,
                                                                    high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
processedFilteredFilesCov20_heat <- methylKit::filterByCoverage(processedFiles_heat,
                                                                    lo.count = 20, lo.perc = NULL,
                                                                    high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)
## STEP 3 - unite/destrand(if necessary)
# all
methylationInformationFiltered <- methylKit::unite(processedFiles,
                                                   destrand = FALSE,
                                                   mc.cores = 8)
# all data
methylationInformationFilteredCov1 <- methylKit::unite(processedFilteredFilesCov1,
                                                       destrand = FALSE,
                                                       mc.cores = 8)
methylationInformationFilteredCov3 <- methylKit::unite(processedFilteredFilesCov3,
                                                       destrand = FALSE,
                                                       mc.cores = 8)
methylationInformationFilteredCov5 <- methylKit::unite(processedFilteredFilesCov5,
                                                       destrand = FALSE,
                                                       mc.cores = 8)
methylationInformationFilteredCov10 <- methylKit::unite(processedFilteredFilesCov10,
                                                        destrand = FALSE,
                                                        mc.cores = 8)
methylationInformationFilteredCov20 <- methylKit::unite(processedFilteredFilesCov20,
                                                        destrand = FALSE,
                                                        mc.cores = 8)


# diploid
methylationInformationFilteredCov1_diploid <- methylKit::unite(processedFilteredFilesCov1_diploid,
                                                       destrand = FALSE,
                                                       mc.cores = 8)
methylationInformationFilteredCov3_diploid <- methylKit::unite(processedFilteredFilesCov3_diploid,
                                                       destrand = FALSE,
                                                       mc.cores = 8)
methylationInformationFilteredCov5_diploid <- methylKit::unite(processedFilteredFilesCov5_diploid,
                                                       destrand = FALSE,
                                                       mc.cores = 8)
methylationInformationFilteredCov10_diploid <- methylKit::unite(processedFilteredFilesCov10_diploid,
                                                                destrand = FALSE,
                                                                mc.cores = 8)
methylationInformationFilteredCov20_diploid <- methylKit::unite(processedFilteredFilesCov20_diploid,
                                                                destrand = FALSE,
                                                                mc.cores = 8)

# triploid
methylationInformationFilteredCov1_triploid <- methylKit::unite(processedFilteredFilesCov1_triploid,
                                                               destrand = FALSE,
                                                               mc.cores = 8)
methylationInformationFilteredCov3_triploid <- methylKit::unite(processedFilteredFilesCov3_triploid,
                                                               destrand = FALSE,
                                                               mc.cores = 8)
methylationInformationFilteredCov5_triploid <- methylKit::unite(processedFilteredFilesCov5_triploid,
                                                               destrand = FALSE,
                                                               mc.cores = 8)
methylationInformationFilteredCov10_triploid <- methylKit::unite(processedFilteredFilesCov10_triploid,
                                                                destrand = FALSE,
                                                                mc.cores = 8)
methylationInformationFilteredCov20_triploid <- methylKit::unite(processedFilteredFilesCov20_triploid,
                                                                destrand = FALSE,
                                                                mc.cores = 8)

# control
methylationInformationFilteredCov1_control <- methylKit::unite(processedFilteredFilesCov1_control,
                                                               destrand = FALSE,
                                                               mc.cores = 8)
methylationInformationFilteredCov3_control <- methylKit::unite(processedFilteredFilesCov3_control,
                                                               destrand = FALSE,
                                                               mc.cores = 8)
methylationInformationFilteredCov5_control <- methylKit::unite(processedFilteredFilesCov5_control,
                                                               destrand = FALSE,
                                                               mc.cores = 8)
methylationInformationFilteredCov10_control <- methylKit::unite(processedFilteredFilesCov10_control,
                                                                destrand = FALSE,
                                                                mc.cores = 8)
methylationInformationFilteredCov20_control <- methylKit::unite(processedFilteredFilesCov20_control,
                                                                destrand = FALSE,
                                                                mc.cores = 8)
# heat
methylationInformationFilteredCov1_heat <- methylKit::unite(processedFilteredFilesCov1_heat,
                                                                destrand = FALSE,
                                                                mc.cores = 8)
methylationInformationFilteredCov3_heat <- methylKit::unite(processedFilteredFilesCov3_heat,
                                                                destrand = FALSE,
                                                                mc.cores = 8)
methylationInformationFilteredCov5_heat <- methylKit::unite(processedFilteredFilesCov5_heat,
                                                                destrand = FALSE,
                                                                mc.cores = 8)
methylationInformationFilteredCov10_heat <- methylKit::unite(processedFilteredFilesCov10_heat,
                                                                 destrand = FALSE,
                                                                 mc.cores = 8)
methylationInformationFilteredCov20_heat <- methylKit::unite(processedFilteredFilesCov20_heat,
                                                                 destrand = FALSE,
                                                                 mc.cores = 8)
## STEP 4 - count CpG
nrow(methylationInformationFilteredCov1) # 5852939 total CpG
nrow(methylationInformationFilteredCov3) #  3714898 total CpG
nrow(methylationInformationFilteredCov5) #  2297236 total CpG
nrow(methylationInformationFilteredCov10) #  480579 total CpG
nrow(methylationInformationFilteredCov20) #  29884 total CpG

nrow(methylationInformationFilteredCov1_diploid) #  6583553 total CpG
nrow(methylationInformationFilteredCov3_diploid) #  4410404 total CpG
nrow(methylationInformationFilteredCov5_diploid) #  2928292 total CpG
nrow(methylationInformationFilteredCov10_diploid) #  766054 total CpG
nrow(methylationInformationFilteredCov20_diploid) #  47820 total CpG

nrow(methylationInformationFilteredCov1_triploid) #  861750 total CpG
nrow(methylationInformationFilteredCov3_triploid) #  4527090 total CpG
nrow(methylationInformationFilteredCov5_triploid) #  3023664 total CpG
nrow(methylationInformationFilteredCov10_triploid) #  861750 total CpG
nrow(methylationInformationFilteredCov20_triploid) #  54998 total CpG

nrow(methylationInformationFilteredCov1_control) #   6337867 total CpG
nrow(methylationInformationFilteredCov3_control) #   3993030 total CpG
nrow(methylationInformationFilteredCov5_control) #   2465533 total CpG
nrow(methylationInformationFilteredCov10_control) #   529068 total CpG
nrow(methylationInformationFilteredCov20_control) #   35027 total CpG

nrow(methylationInformationFilteredCov1_heat) #   7261925 total CpG
nrow(methylationInformationFilteredCov3_heat) #   5416388 total CpG
nrow(methylationInformationFilteredCov5_heat) #   4194066 total CpG
nrow(methylationInformationFilteredCov10_heat) #   1921856 total CpG
nrow(methylationInformationFilteredCov20_heat) #   193521 total CpG

## Sample-specific coverage and methylation statistics - by ploidy
# diploid, 1x
nFiles <- length(analysisFiles_diploid) #Count number of samples
fileName <- data.frame("nameBase"  = rep("plots/percent_CpG_coverage_plots/diploid/percent-CpG-methylation", times = nFiles),
                       "nameBase2" = rep("plots/percent_CpG_coverage_plots/diploid/percent-CpG-coverage", times = nFiles),
                       "sample.ID" = 1:11) #Create new dataframe for filenames
head(fileName) #Confirm dataframe creation

fileName$actualFileName1 <- paste(fileName$nameBase, "-Filtered", "-1x_coverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG methylation plot
fileName$actualFileName2 <- paste(fileName$nameBase2, "-Filtered", "-1x_coverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG coverage plot

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName1[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getMethylationStats(processedFilteredFilesCov1_diploid[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
  dev.off() #Turn off plotting device
} #Plot and save %CpG methylation information

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName2[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getCoverageStats(processedFilteredFilesCov1_diploid[[i]], plot = TRUE, both.strands = FALSE) #Get CpG coverage information
  dev.off() #Turn off plotting device
} #Plot and save CpG coverage information

# diploid, 10x
nFiles <- length(analysisFiles_diploid) #Count number of samples
fileName <- data.frame("nameBase"  = rep("plots/percent_CpG_coverage_plots/diploid/10x/percent-CpG-methylation", times = nFiles),
                       "nameBase2" = rep("plots/percent_CpG_coverage_plots/diploid/10x/percent-CpG-coverage", times = nFiles),
                       "sample.ID" = 1:11) #Create new dataframe for filenames
head(fileName) #Confirm dataframe creation

fileName$actualFileName1 <- paste(fileName$nameBase, "-Filtered", "-10x_coverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG methylation plot
fileName$actualFileName2 <- paste(fileName$nameBase2, "-Filtered", "-10x_coverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG coverage plot

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName1[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getMethylationStats(processedFilteredFilesCov10_diploid[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
  dev.off() #Turn off plotting device
} #Plot and save %CpG methylation information

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName2[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getCoverageStats(processedFilteredFilesCov10_diploid[[i]], plot = TRUE, both.strands = FALSE) #Get CpG coverage information
  dev.off() #Turn off plotting device
} #Plot and save CpG coverage information

# triploid,1x
nFiles <- length(analysisFiles_triploid) #Count number of samples
fileName <- data.frame("nameBase"  = rep("plots/percent_CpG_coverage_plots/triploid/percent-CpG-methylation", times = nFiles),
                       "nameBase2" = rep("plots/percent_CpG_coverage_plots/triploid/percent-CpG-coverage", times = nFiles),
                       "sample.ID" = 1:11) #Create new dataframe for filenames
head(fileName) #Confirm dataframe creation

fileName$actualFileName1 <- paste(fileName$nameBase, "-Filtered", "-1x_coverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG methylation plot
fileName$actualFileName2 <- paste(fileName$nameBase2, "-Filtered", "-1x_coverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG coverage plot

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName1[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getMethylationStats(processedFilteredFilesCov1_triploid[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
  dev.off() #Turn off plotting device
} #Plot and save %CpG methylation information

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName2[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getCoverageStats(processedFilteredFilesCov1_triploid[[i]], plot = TRUE, both.strands = FALSE) #Get CpG coverage information
  dev.off() #Turn off plotting device
} #Plot and save CpG coverage information

# control, 1x
nFiles <- length(analysisFiles_control) #Count number of samples
fileName <- data.frame("nameBase"  = rep("plots/percent_CpG_coverage_plots/control/percent-CpG-methylation", times = nFiles),
                       "nameBase2" = rep("plots/percent_CpG_coverage_plots/control/percent-CpG-coverage", times = nFiles),
                       "sample.ID" = 1:12) #Create new dataframe for filenames
head(fileName) #Confirm dataframe creation

fileName$actualFileName1 <- paste(fileName$nameBase, "-Filtered", "-1x_coverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG methylation plot
fileName$actualFileName2 <- paste(fileName$nameBase2, "-Filtered", "-1x_coverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG coverage plot

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName1[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getMethylationStats(processedFilteredFilesCov1_control[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
  dev.off() #Turn off plotting device
} #Plot and save %CpG methylation information

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName2[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getCoverageStats(processedFilteredFilesCov1_control[[i]], plot = TRUE, both.strands = FALSE) #Get CpG coverage information
  dev.off() #Turn off plotting device
} #Plot and save CpG coverage information

# heat, 1x
nFiles <- length(analysisFiles_heat) #Count number of samples
fileName <- data.frame("nameBase"  = rep("plots/percent_CpG_coverage_plots/heat/percent-CpG-methylation", times = nFiles),
                       "nameBase2" = rep("plots/percent_CpG_coverage_plots/heat/percent-CpG-coverage", times = nFiles),
                       "sample.ID" = 1:10) #Create new dataframe for filenames
head(fileName) #Confirm dataframe creation

fileName$actualFileName1 <- paste(fileName$nameBase, "-Filtered", "-1x_coverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG methylation plot
fileName$actualFileName2 <- paste(fileName$nameBase2, "-Filtered", "-1x_coverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG coverage plot

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName1[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getMethylationStats(processedFilteredFilesCov1_heat[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
  dev.off() #Turn off plotting device
} #Plot and save %CpG methylation information

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName2[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getCoverageStats(processedFilteredFilesCov1_heat[[i]], plot = TRUE, both.strands = FALSE) #Get CpG coverage information
  dev.off() #Turn off plotting device
} #Plot and save CpG coverage information

##################

nFiles <- length(analysisFiles) #Count number of samples
fileName <- data.frame("nameBase"  = rep("plots/percent_CpG_coverage_plots/percent-CpG-methylation", times = nFiles),
                       "nameBase2" = rep("plots/percent_CpG_coverage_plots/percent-CpG-coverage", times = nFiles),
                       "sample.ID" = 1:22) #Create new dataframe for filenames
head(fileName) #Confirm dataframe creation

fileName$actualFileName1 <- paste(fileName$nameBase, "-Filtered", "-1x_coverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG methylation plot
fileName$actualFileName2 <- paste(fileName$nameBase2, "-Filtered", "-1x_coverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG coverage plot

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName1[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getMethylationStats(processedFilteredFilesCov1[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
  dev.off() #Turn off plotting device
} #Plot and save %CpG methylation information

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName2[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getCoverageStats(processedFilteredFilesCov1[[i]], plot = TRUE, both.strands = FALSE) #Get CpG coverage information
  dev.off() #Turn off plotting device
} #Plot and save CpG coverage information


clusteringInformationFilteredCov1 <- methylKit::clusterSamples(methylationInformationFilteredCov1, dist = "correlation", method = "ward", plot = FALSE) #Save cluster information as a new object

# jpeg(filename = "plots/percent_CpG_coverage_plots/Full-Sample-Pearson-Correlation-Plot-FilteredCov1Destrand.jpeg", height = 1000, width = 1000) #Save file with designated name
# methylKit::getCorrelation(methylationInformationFilteredCov1, plot = TRUE) #Understand correlation between methylation patterns in different samples
# dev.off()



### Plot

setwd(home); getwd() # reset working directory to home dir
total_cpg             <- read_excel("analyses/coverage_summary.xlsx", sheet = "total_cpg", col_names = TRUE)
total_cpg$ploidy      <- factor(total_cpg$ploidy,levels=c("diploid","triploid","both","control","heat"),ordered=TRUE)
total_meth            <- read_excel("analyses/coverage_summary.xlsx", sheet = "ploidy_plot", col_names = TRUE)
total_meth$ploidy     <- factor(total_meth$ploidy,levels=c("diploid","triploid","both"),ordered=TRUE)
total_meth$bin        <- factor(total_meth$bin,levels=c('0-10',
                                                        '11-20',
                                                        '21-30',
                                                        '31-40',
                                                        '41-50',
                                                        '51-60',
                                                        '61-70',
                                                        '71-80',
                                                        '81-90',
                                                        '91-100'), ordered=TRUE)

p1 <- ggplot(total_cpg, aes(x=coverage,y=cpg,group=ploidy,color=ploidy)) +
      geom_smooth(size = 1.2,
                  method = lm, formula = y ~ splines::bs(x, 4), se = FALSE) +    
      geom_point(size = 5, pch = c(rep(16, times = 5), 
                                   rep(16, times = 5),
                                   rep(17, times = 5),
                                   rep(15, times = 5),
                                   rep(15, times = 5))) +
      scale_color_manual(values=c("royalblue1",
                                  "orangered1",
                                  "chartreuse4",
                                  "royalblue1",
                                  "orangered1")) +
      # geom_line(size = 1.2, color="chartreuse4") +
      # geom_point() +
      # scale_y_continuous(breaks = seq(0, 15000, 5000), limits = c(0, 17000)) +
      # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) + 
      my_theme

p1

p2 <- ggplot(total_meth, aes(x=bin, y=percent,fill=ploidy)) +
              geom_errorbar(aes(ymin = percent-se, ymax = percent+se),position = position_dodge(width = 0.8), width = 0) +
              geom_bar(stat="identity", line="black",position = position_dodge(width = 0.8)) +
              scale_fill_manual(values=c("royalblue1","orangered1","chartreuse4")) +
              # geom_point() +
              scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 90)) +
              scale_x_discrete(expand=c(0,1)) + 
              my_theme

p2

total_meth            <- read_excel("analyses/coverage_summary.xlsx", sheet = "trt_plot", col_names = TRUE)
total_meth$trt     <- factor(total_meth$trt,levels=c("control","heat","both"),ordered=TRUE)
total_meth$bin        <- factor(total_meth$bin,levels=c('0-10',
                                                        '11-20',
                                                        '21-30',
                                                        '31-40',
                                                        '41-50',
                                                        '51-60',
                                                        '61-70',
                                                        '71-80',
                                                        '81-90',
                                                        '91-100'), ordered=TRUE)


p3 <- ggplot(total_meth, aes(x=bin, y=percent,fill=trt)) +
      geom_errorbar(aes(ymin = percent-se, ymax = percent+se),position = position_dodge(width = 0.8), width = 0) +
      geom_bar(stat="identity", line="black",position = position_dodge(width = 0.8)) +
      scale_fill_manual(values=c("royalblue1","orangered1","chartreuse4")) +
      # geom_point() +
      scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 90)) +
      scale_x_discrete(expand=c(0,1)) + 
      my_theme

p3

ggsave("plots/coverage_totalCpG.png",
       plot   = p1,
       dpi    = 1200,
       device = "png",
       width  = 6,
       height = 5,
       units  = "in")

ggsave("plots/methylated_percentage_ploidy.png",
       plot   = p2,
       dpi    = 1200,
       device = "png",
       width  = 6,
       height = 5.3,
       units  = "in")

ggsave("plots/methylated_percentage_trt.jpeg",
       plot   = p3,
       dpi    = 1200,
       device = "jpeg",
       width  = 6,
       height = 5.3,
       units  = "in")


### PCA

allDataPCA <- PCASamples(methylationInformationFilteredCov1, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA) #Look at summary statistics. The first PC explains 11.0% of variation, the second PC explains 6.63% of variation

pdf("plots/all-sample-PCA.pdf", width = 10, height = 10)

par(mar = c(5, 5, 1, 1)) #Specify inner and outer margins

fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), 
                           type = "none", display = "sites", cex = 0.5, 
                           xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points

points(fig.allDataPCA, "sites", col = c(rep(alpha("royalblue1", 0.8), times = 5), 
                                        rep(alpha("orangered1", 0.8), times = 5),
                                        rep(alpha("royalblue1", 0.8), times = 6),
                                        rep(alpha("orangered1", 0.8), times = 6)), 
       pch = c(rep(17, times = 11), 
               rep(16, times = 11)), cex = 3) #Add each sample. Purple = diploid, pink = triploid, 16 = circle, 17 = triangle

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

ordiellipse(allDataPCA, sampleMetadata$ploidy, show.groups = "0", col = "royalblue1") #Add confidence ellipse - diploid
ordiellipse(allDataPCA, sampleMetadata$ploidy, show.groups = "1", col = "orangered1") #Add confidence ellipse - triploid

ordiellipse(allDataPCA, sampleMetadata$heat, show.groups = "0", col = "royalblue1", lty = 2) #Add confidence ellipse - control
ordiellipse(allDataPCA, sampleMetadata$heat, show.groups = "1", col = "orangered1", lty = 2) #Add confidence ellipse - desiccation

axis(side =  1, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1 (8.4%)", line = 3, cex = 1.5) #Add x-axis label

axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2 (6.9%)", line = 3, cex = 1.5) #Add y-axis label

legend("bottomright", 
       pch = c(15, 15, 16, 17), 
       legend = c("diploid", "triploid", "control", "desiccation"), 
       col = c("royalblue1","orangered1", "grey20", "grey20"), 
       cex = 1.7, bty = "n") #Add a legend with information about ambient and elevated samples

dev.off()


