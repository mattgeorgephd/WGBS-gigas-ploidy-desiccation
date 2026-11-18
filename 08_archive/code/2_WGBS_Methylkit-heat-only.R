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
setwd(dirname(home)); getwd()
setwd('..')
home <- getwd(); getwd()

## Load previous session data #################################################################################################

setwd("E:/WGBS-gigas-ploidy-desiccation")
save.image("R_sessions/20211005_Ronit-SHOCK.RData") #Save R Data in case R crashes
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

# Put all .cov files into a list for analysis.
setwd("E:/WGBS-gigas-ploidy-desiccation/") # set working directory to location of .cov files
# meta <- read_excel("metadata.xlsx", sheet = "list", col_names = TRUE)
meta <- read_excel("metadata.xlsx", sheet = "list3", col_names = TRUE)

analysisFiles <- as.list(meta$path)

# setwd(dirname(home)); getwd() # reset working directory to home dir

####################################################################################################
## Run Methylkit

# Create data.frame with sample info
sampleMetadata <- data.frame("sampleID" = meta$sampleID,
                             "ploidyTreatment" =  meta$ploidy,
                             "heatTreatment"   =  meta$trt) # Specify which treatment the samples were from. All animals were subjected to desiccation. 0 = diploid, 1 = triploid
sampleMetadata$sampleID <- as.character(sampleMetadata$sampleID)

head(sampleMetadata) #Confirm data.frame creation
# rm(meta)

# Use `methRead` to create a methylation object from the coverage files, and include sample ID and treatment information.
processedFiles <- methylKit::methRead(analysisFiles,
                                      sample.id = as.list(sampleMetadata$sampleID),
                                      assembly = "Roslin",
                                      treatment = sampleMetadata$heatTreatment,
                                      pipeline = "bismarkCoverage",
                                      mincov = 2) #Process files. Treatment specified based on ploidy status. Use mincov = 2 to quickly process reads.

processedFilteredFilesCov10 <- methylKit::filterByCoverage(processedFiles,
                                                          lo.count = 10, lo.perc = NULL,
                                                          high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)

methylationInformationFilteredCov10 <- methylKit::unite(processedFilteredFilesCov10,
                                                       destrand = FALSE,
                                                       mc.cores = 8)

head(methylationInformationFilteredCov10)
nrow(methylationInformationFilteredCov10)

write.csv(methylationInformationFilteredCov10, "DML_lists/ALL-CpG-Cov10-Ronit.csv") #Save table as .csv


differentialMethylationStatsHEATCov10 <- methylKit::calculateDiffMeth(methylationInformationFilteredCov10,
                                                                      mc.cores = 1) #Calculate differential methylation statistics based on treatment indication from methRead. Include HEAT as a covariate.

head(differentialMethylationStatsHEATCov10) #Look at differential methylation statistics

write.csv(differentialMethylationStatsHEATCov10, "DML_lists/DML-stats-SHOCK-Cov10-Ronit.csv") #Save table as .csv

# diffMethStatsPloidyCov10Diff20 <- methylKit::getMethylDiff(differentialMethylationStatsPloidyCov10, difference = 20, qvalue = 0.01) #Identify loci that are at least 20% different
# length(diffMethStatsPloidyCov10Diff20$chr) #Count the number of DML --> 1038
# head(diffMethStatsPloidyCov10Diff20) #Confirm creation

diffMethStatsHEATCov10Diff20 <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 20, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsHEATCov10Diff20$chr) #Count the number of DML --> 1038
head(diffMethStatsHEATCov10Diff20) #Confirm creation

write.csv(diffMethStatsHEATCov10Diff20, "DML_lists/DML-getMethylDiff-SHOCK-Cov10-20-Ronit.csv") #Save table as .csv

#### Hypermethylated DML
diffMethStats20FilteredCov10Hyper <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 20, qvalue = 0.01, type = "hyper") #Identify hypermethylated loci that are at least 20% different
head(diffMethStats20FilteredCov10Hyper) #Confirm creation
nrow(diffMethStats20FilteredCov10Hyper)
write.csv(diffMethStats20FilteredCov10Hyper, "DML_lists/DML-getMethylDiff-Cov10-20-Ronit-SHOCK-hyper.csv") #Save table as .csv

##### Hypomethylated DML
diffMethStats20FilteredCov10Hypo <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 20, qvalue = 0.01, type = "hypo") #Identify hypomethylated loci that are at least 50% different
head(diffMethStats20FilteredCov10Hypo) #Confirm creation
nrow(diffMethStats20FilteredCov10Hypo)
write.csv(diffMethStats20FilteredCov10Hypo, "DML_lists/DML-getMethylDiff-Cov10-20-Ronit-SHOCK-hypo.csv") #Save table as .csv

dd <- read_excel("DML_lists/DML_summary.xlsx", sheet = "ronit", col_names = TRUE)

a <- dd %>%
  group_by(ID) %>%
  filter(n() == 2)

write.csv(a, "DML_lists/Ronit_overlap.csv")

####################################################################################################
## Run Methylkit - just Yaamani's controls

# Create data.frame with sample info
sampleMetadata <- data.frame("sampleID" = meta$sampleID,
                             "ploidyTreatment" =  meta$ploidy,
                             "heatTreatment"   =  meta$trt) # Specify which treatment the samples were from. All animals were subjected to desiccation. 0 = diploid, 1 = triploid
sampleMetadata$sampleID <- as.character(sampleMetadata$sampleID)

head(sampleMetadata) #Confirm data.frame creation
# rm(meta)

processedFiles <- methylKit::methRead(analysisFiles,
                                      sample.id = as.list(sampleMetadata$sampleID),
                                      assembly = "Roslin",
                                      treatment = sampleMetadata$ploidyTreatment,
                                      pipeline = "bismarkCoverage",
                                      mincov = 2) #Process files. Treatment specified based on ploidy status. Use mincov = 2 to quickly process reads.

processedFilteredFilesCov10 <- methylKit::filterByCoverage(processedFiles,
                                                          lo.count = 10, lo.perc = NULL,
                                                          high.count = NULL, high.perc = 99.9) %>% methylKit::normalizeCoverage(.)

methylationInformationFilteredCov10 <- methylKit::unite(processedFilteredFilesCov10,
                                                       destrand = FALSE,
                                                       mc.cores = 8)

head(methylationInformationFilteredCov10)
nrow(methylationInformationFilteredCov10)

write.csv(methylationInformationFilteredCov10, "ALL-CpG-Cov10-Yaamini-controls.csv") #Save table as .csv


differentialMethylationStatsPloidyCov10 <- methylKit::calculateDiffMeth(methylationInformationFilteredCov10,
                                                                        mc.cores = 1) #Calculate differential methylation statistics based on treatment indication from methRead. Include HEAT as a covariate.

head(differentialMethylationStatsHEATCov10) #Look at differential methylation statistics

write.csv(differentialMethylationStatsHEATCov10, "DML-stats-ploidy-Cov10-Yaamini.csv") #Save table as .csv

# diffMethStatsPloidyCov10Diff20 <- methylKit::getMethylDiff(differentialMethylationStatsPloidyCov10, difference = 20, qvalue = 0.01) #Identify loci that are at least 20% different
# length(diffMethStatsPloidyCov10Diff20$chr) #Count the number of DML --> 1038
# head(diffMethStatsPloidyCov10Diff20) #Confirm creation

diffMethStatsCov10Diff20 <- methylKit::getMethylDiff(differentialMethylationStatsPloidyCov10, difference = 20, qvalue = 0.01) #Identify loci that are at least 20% different
length(diffMethStatsCov10Diff20$chr) #Count the number of DML --> 1038
head(diffMethStatsCov10Diff20) #Confirm creation


write.csv(diffMethStatsCov10Diff20, "DML-getMethylDiff-Yaamini-Cov10-20.csv") #Save table as .csv

#### Hypermethylated DML
diffMethStats20FilteredCov10Hyper <- methylKit::getMethylDiff(differentialMethylationStatsPloidyCov10, difference = 20, qvalue = 0.01, type = "hyper") #Identify hypermethylated loci that are at least 20% different
head(diffMethStats20FilteredCov10Hyper) #Confirm creation
write.csv(diffMethStats20FilteredCov10Hyper, "DML-getMethylDiff-Yaamini-Cov10-20-hyper.csv") #Save table as .csv

##### Hypomethylated DML
diffMethStats20FilteredCov10Hypo <- methylKit::getMethylDiff(differentialMethylationStatsPloidyCov10, difference = 20, qvalue = 0.01, type = "hypo") #Identify hypomethylated loci that are at least 50% different
head(diffMethStats20FilteredCov10Hypo) #Confirm creation
write.csv(diffMethStats20FilteredCov10Hypo, "DML-getMethylDiff-Yaamini-Cov10-20-hypo.csv") #Save table as .csv













# Identify DML - HEAT

# # Reorganize samples based on HEAT since ploidy was used as treatment variable
# 
# methylationInformationFilteredCov10T <- methylKit::reorganize(processedFilteredFilesCov10,
#                                                              sample.id = sampleMetadata$sampleID,
#                                                              treatment = sampleMetadata$heatTreatment)
# 
#   
# methylationInformationFilteredCov10 <- methylKit::unite(methylationInformationFilteredCov10T, destrand = FALSE, mc.cores = 8, min.per.group = 8L) #Re-assign treatment information to filtered and normalized 5x CpG loci, then unite 5x loci with coverage in at least 8/12 samples/treatment
# 
# 
# covariatePLOIDY <- data.frame("PLOIDY" = c(rep("H", times = 10),
#                                            rep("C", times = 12))) #Create dataframe with HEAT covariate information



# #### Hypermethylated DML
# diffMethStats20FilteredCov10Hyper <- methylKit::getMethylDiff(differentialMethylationStatsPloidyCov10, difference = 20, qvalue = 0.01, type = "hyper") #Identify hypermethylated loci that are at least 20% different
# head(diffMethStats20FilteredCov10Hyper) #Confirm creation
# write.csv(diffMethStats20FilteredCov10Hyper, "DML-getMethylDiff-ploidy-Cov10-20-hyper.csv") #Save table as .csv
# 
# ##### Hypomethylated DML
# diffMethStats20FilteredCov10Hypo <- methylKit::getMethylDiff(differentialMethylationStatsPloidyCov10, difference = 20, qvalue = 0.01, type = "hypo") #Identify hypermethylated loci that are at least 50% different
# head(diffMethStats20FilteredCov10Hypo) #Confirm creation
# write.csv(diffMethStats20FilteredCov10Hypo, "DML-getMethylDiff-ploidy-Cov10-20-hypo.csv") #Save table as .csv
# 

##### Hypermethylated DML
diffMethStats20FilteredCov10Hyper <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 20, qvalue = 0.01, type = "hyper") #Identify hypermethylated loci that are at least 20% different
head(diffMethStats20FilteredCov10Hyper) #Confirm creation
write.csv(diffMethStats20FilteredCov10Hyper, "DML-getMethylDiff-HEAT-Cov10-20-hyper.csv") #Save table as .csv

##### Hypomethylated DML
diffMethStats20FilteredCov10Hypo <- methylKit::getMethylDiff(differentialMethylationStatsHEATCov10, difference = 20, qvalue = 0.01, type = "hypo") #Identify hypermethylated loci that are at least 50% different
head(diffMethStats20FilteredCov10Hypo) #Confirm creation
write.csv(diffMethStats20FilteredCov10Hypo, "DML-getMethylDiff-HEAT-Cov10-20-hypo.csv") #Save table as .csv



#### DML distribution figure

##### Calculate distribution of DML across chromosomes
DMLchrCounts <-as.data.frame(table(diffMethStatsHEATCov10Diff20$chr)) #Count the number of DML/chromosome
DMLchrCounts <- DMLchrCounts[-c(11:length(DMLchrCounts$Var1)),]
prefix <- "Chr"
suffix <- seq(1:length(DMLchrCounts$Var1))
my_names <- paste(prefix,suffix,sep="")
DMLchrCounts$chr <- my_names
DMLchrCounts <- DMLchrCounts[,-1] #Remove column with chromosome RefSeq ID
DMLchrCounts <- DMLchrCounts[,c(2,1)] #Reorganize columns
colnames(DMLchrCounts) <- c("chr", "DMLCount") #Rename columns
head(DMLchrCounts) #Confirm formatting
write.csv(DMLchrCounts, "DML-per-Chromosome.csv", row.names = FALSE, quote = FALSE) #Save file

##### [*C. gigas* genome information from NCBI](https://www.ncbi.nlm.nih.gov/genome/?term=pacific+oyster)

DMLchrCounts$geneCount <- c(3425, 3526, 3796, 3124, 3915, 4129, 4097, 3617, 1829, 3713) #Create column with number of gene sequences in each chromosome
DMLchrCounts$DMLbyGenes <- (DMLchrCounts$DMLCount / DMLchrCounts$geneCount)*100
head(DMLchrCounts) #Confirm column creation

DMLchrCounts$chrLengthMb <- c(55.79, 73.22, 58.32, 53.13, 73.55, 60.15, 62.11, 58.46, 37.09, 57.54) #Create column with number Mb in each chromosome
DMLchrCounts$DMLbyChrLength <- (DMLchrCounts$DMLCount / DMLchrCounts$chrLengthMb) #Normalize DML counts by chr length in Mb. Axes should indicate that these counts are multiplied by 10e-6
head(DMLchrCounts) #Confirm column creation

DMLchrCounts$chrCpGCounts <- c(1394962, 1237469, 1568195, 1200917, 2008680, 1121671, 1251184, 1633946, 2300335, 740913) #Create column with number of CpGs in each chromosome
DMLchrCounts$DMLbyChrCpG <- (DMLchrCounts$DMLCount / DMLchrCounts$chrCpGCounts) #Normalize DML counts by number of CpGs in each chromosome
range(DMLchrCounts$DMLbyChrCpG) #Look at range of values
head(DMLchrCounts) #Confirm column creation

### Need to update these numbers for c.gigas
DMLchrCounts$chrCpGCounts <- c(1394962, 1237469, 1568195, 1200917, 2008680, 1121671, 1251184, 1633946, 2300335, 740913) #Create column with number of CpGs in each chromosome
DMLchrCounts$DMLbyChrCpG <- (DMLchrCounts$DMLCount / DMLchrCounts$chrCpGCounts) #Normalize DML counts by number of CpGs in each chromosome
range(DMLchrCounts$DMLbyChrCpG) #Look at range of values
head(DMLchrCounts) #Confirm column creation

### Create Figure

pdf("2021-10-04-DML-and-Gene-Distribution.pdf", height = 8.5, width = 12) #Save figure as a pdf
par(mar = c(5,7,2,10)) #Change figure boundaries
DMLbarplot <- barplot(as.matrix(t(DMLchrCounts$DMLbyChrCpG)),
                      axes = FALSE, 
                      names.arg = DMLchrCounts$chr,
                      cex.names = 1.5,
                      xlim = c(0.7,11.5),
                      ylim = c(0,60e-05),
                      col = "grey80") #Create a barplot and save as a new object. Use axes = FALSE to remove the y-axis and names.arg to set labels on the x-axis. The object contains x coordinates for bars, so xlim is set at 12 to compensate for maximum value of 11.5
mtext(side = 1, "Chromosome", las=1, line = 3, cex = 1.5) #Add x-axis label
axis(side = 2, line = 1.5, at = seq(0, 60e-05, by = 10e-05), las = 2, col = "grey80", cex.axis = 1.2) #Add y-axis for DML counts
mtext(side = 2, "Number DML per 10,000 CpGs", line = 5, cex = 1.5) #Add y-axis label for DML counts
par(new = TRUE) #Create a new plot
plot(x = DMLbarplot,
     y = DMLchrCounts$geneCount,
     type = "b",
     axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i",
     pch = 16, col = "grey20",
     xlim = c(0,12), ylim = c(0,6500)) #Plot points and lines (type = "b") for gene count by chromosome. Use the coordinates from DMLbarplot (x = DMLbarplot) and set xlim = (0,12) so plots are lined up. Use axes = FALSE to remove both axes. Remove x and y lables (xlab = ""; ylab = ""). Set ylim = (0,6500) to account for max y-values. Use xaxs and yxaxs to remove space between axes.
axis(side = 4, line = 1.5, at = seq(0, 6500, by = 500), las = 2, col = "grey80", cex.axis = 1.2) #Add y-axis for gene sequence counts
mtext(side = 4, "Number of Genes", line = 6, cex = 1.5) #Add y-axis label for gene sequence counts
dev.off() #Turn off plotting device


###### Principal Component Analyses

sample.IDs <- as.list(sampleMetadata$sampleID) #Create list of sample IDs
treatmentSpecification <- sampleMetadata$ploidyTreatment 

plotCustomization <- data.frame(sample = 1:12,
                                treatmentSpecification) #Create dataframe with sample treatment information
head(plotCustomization) #Confirm dataframe creation

allDataPCA <- PCASamples(methylationInformationFilteredCov10, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA) #Look at summary statistics. The first PC explains 15.6% of variation, the second PC explains 12.1% of variation


RColorBrewer::display.brewer.all() #Show all RColorBrewer palettes. I will choose greens.
plotColors <- rev(RColorBrewer::brewer.pal(5, "GnBu")) #Create a color palette for the barplots. Use 5 green-blue shades from RColorBrewer. Reverse theorder so the darkest shade is used first.

pdf("2021-10-04-All-Yaamini-PCA.pdf", width = 11, height = 8.5)
par(mar = c(5, 5, 1, 1)) #Specify inner and outer margins
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep(plotColors[2], times = 5), rep(plotColors[4], times = 5)), pch = c(rep(16, times = 5), rep(17, times = 5)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2
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

####### DML

# The first thing I need to do is subset the data so it only includes DML. Then, I can use similar code to what I used above to create PCA plots.

DMLPositions <- rep(0, times = length(diffMethStatsCov10Diff20$chr)) #Create an empty vector with 598 places to store row numbers
for (i in 1:length(DMLPositions)) {
  DMLPositions[i] <- which(getData(diffMethStatsCov10Diff20)$start[i] == getData(methylationInformationFilteredCov10)$start)
} #For each DML, save the row number where that DML is found in methylationInformationFilteredCov5
tail(DMLPositions) #Confirm vector was created

DMLMatrix <- methylationInformationFilteredCov10[DMLPositions,] #Subset methylationInformationFilteredCov5Destrand to only include DML and save as a new methylBase object
sum((DMLMatrix$start) == (diffMethStatsCov10Diff20$start)) == length(diffMethStatsCov10Diff20$start) #Confirm that start columns are identical. If they are identical, the sum of all TRUE statements should equal the length of the original methylBase object
tail(DMLMatrix) #Confirm methylBase object creation

DMLDataPCA <- PCASamples(DMLMatrix, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(DMLDataPCA) #Look at summary statistics. The first PC explains 24.0% of variation, the second PC explains 32.8% of variation


pdf("2021-10-04-DML-Only-Yaamini-PCA.pdf", width = 11, height = 8.5)
par(mar = c(5, 5, 1, 1)) #Specify inner and outer margins
fig.DMLDataPCA <- ordiplot(DMLDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.DMLDataPCA, "sites", col = c(rep(plotColors[2], times = 5), rep(plotColors[4], times = 5)), pch = c(rep(16, times = 5), rep(17, times = 5)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2
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
ordiellipse(DMLDataPCA, plotCustomization$treatment, show.groups = "1", col = plotColors[4]) #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(DMLDataPCA, plotCustomization$treatment, show.groups = "0", col = plotColors[2]) #Add confidence ellipse around the samples in ambient pCO2
axis(side =  1, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1 (36.8%)", line = 3, cex = 1.5) #Add x-axis label
axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2 (11.9%)", line = 3, cex = 1.5) #Add y-axis label
legend("topright", 
       pch = c(16, 17), 
       legend = c("control", "heated"), 
       col = c(plotColors[2], plotColors[4]), 
       cex = 1.7, bty = "n") #Add a legend with information about ambient and elevated samples
dev.off()

#### Multipanel plot

pdf("output/2021-05-26-PCA-Multpanel.pdf", height = 8.5, width = 11) #Save plot
par(mfrow = c(1, 2), oma = c(5, 2, 2, 0), mar = c(0, 3, 0, 5)) #Set up parameters for multipanel plot
#All CpG Loci
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep(plotColors[2], times = 5), rep(plotColors[4], times = 5)), pch = c(rep(16, times = 5), rep(17, times = 5)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2
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
ordiellipse(allDataPCA, plotCustomization$treatment, show.groups = "1", col = plotColors[4]) #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, plotCustomization$treatment, show.groups = "0", col = plotColors[2]) #Add confidence ellipse around the samples in ambient pCO2
axis(side =  1, at = seq(-400, 200, 200), col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1 (15.6%)", line = 3, cex = 1.5) #Add x-axis label
axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2 (18.1%)", line = 3, cex = 1.5) #Add y-axis label
mtext(side = 3, line = -5, adj = c(-100,0), text = "    a. All CpG Loci", cex = 1.5)

fig.DMLDataPCA <- ordiplot(DMLDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-20,20), ylim = c(-10,20), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.DMLDataPCA, "sites", col = c(rep(plotColors[2], times = 5), rep(plotColors[4], times = 5)), pch = c(rep(16, times = 5), rep(17, times = 5)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2
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
ordiellipse(DMLDataPCA, plotCustomization$treatment, show.groups = "1", col = plotColors[4]) #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(DMLDataPCA, plotCustomization$treatment, show.groups = "0", col = plotColors[2]) #Add confidence ellipse around the samples in ambient pCO2
axis(side =  1, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1 (61.8%)", line = 3, cex = 1.5) #Add x-axis label
axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2 (7.5%)", line = 3, cex = 1.5) #Add y-axis label
mtext(side = 3, line = -5, adj = c(-100, 0), text = "    b. DML", cex = 1.5) #Add test category
legend(x = 0, y = -22, 
       pch = c(16, 17), 
       legend = c("Control", "Elevated"), 
       col = c(plotColors[2], plotColors[4]),
       y.intersp = 1, x.intersp = 1,
       cex = 1.7, bty = "n") #Add a legend with information about ambient and elevated samples
dev.off()

######################################################################################################
### Heat Map

percMethDML <- percMethylation(DMLMatrix, rowids = TRUE) #Get percent methylation for all samples at DML. Include row IDs (chr, start, end) information
head(percMethDML) #Confirm percent methylation matrix was created

pheatmap(percMethDML, color = rev(plotColors),
         cluster_rows = TRUE, clustering_distance_rows = "euclidean", treeheight_row = 70, show_rownames = FALSE,
         cluster_cols = TRUE, clustering_distance_cols = "euclidean", treeheight_col = 40, show_colnames = FALSE,
         annotation_col = data.frame(pCO2 = factor(rep(c("Ambient","Treatment"), each = 5))),
         annotation_colors = list(pCO2 = c(Ambient = "grey90", Treatment = "grey10")), 
         annotation_legend = FALSE, annotation_names_col = FALSE,
         legend = TRUE) #Create heatmap using pheatmap using percMethDML and plotColors color scheme. Cluster rows and columns using euclidean distances. Adjust the dendogram tree heights and do not show any row or column names. Create a dataframe with treatment information using annotation_col. Use annotation_colors to indicate colors for treatment ("grey90") and ambient ("grey10") samples. Do not include an annotation_legend or name for annotations (annotatino_names_col). Include a legend.

pdf("2021-05-27-DML-Only-Heatmap.pdf", height = 8.5, width = 11)
par(oma = c(0, 1, 0, 0)) #Adjust outer margins
heatmap.2(percMethDML, col = rev(plotColors), scale = "none", margins = c(1,1),
          trace = "none", tracecol = "black",
          labRow = FALSE, labCol = FALSE, 
          ColSideColors = c(rep("grey90", times = 5), rep("grey10", times = 5)),
          key = TRUE, keysize = 1.8, density.info = "density", key.title = "", key.xlab = "% Methylation", key.ylab = "",
          key.par = list(cex.lab = 2.0, cex.axis = 1.5)) #Create heatmap using heatmap.2 from gplots package using percMethDML data. Use plotColors but do not scale data, label rows, or label columns. Use ColSideColors to indicate colors for treatment and ambient samples. Add a legend using key, and adjust keysize. Have key display density data with density.info. Do not add a key title or y-axis label, and label x axis with key.xlab.
mtext("Density", cex = 1.6, las = 3, adj = 0.8, padj = -29) #Manually add y-axis label for key since heatmap.2 doesn't let you change font size
dev.off()


############################################################################################################
### Convert DML list to Bed files

DML <- data.frame(diffMethStatsHEATCov10Diff20$chr,
                  diffMethStatsHEATCov10Diff20$start,
                  diffMethStatsHEATCov10Diff20$end,
                  diffMethStatsHEATCov10Diff20$meth.diff) # arrange(chr, start) %>% #Join + and - strand information to be saved as a BED file, and avoid writing information in scientific notation
write_delim(DML, "bed_files/2021-10-04-DML-Locations-SHOCK-20.bed",  delim = '\t', col_names = FALSE) #Save data as a BED file


############################################################################################################
### DML-genome-feature-location Stacked bar plot

data         <- read_excel("analyses/DML_locations_genome_feature.xlsx", sheet = "plot", col_names = TRUE)
data$feature <- factor(data$feature, levels = c("exon","intron","intergenic","TE","UTR"), ordered=TRUE) 

v = RColorBrewer::brewer.pal()

# Plot
p1 <- ggplot(data = data, aes(y = all_CpG_percent, x = 'Loci (10x coverage)', fill = feature)) + 
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

ggsave("plots/All_CpGs.png",
       plot   = p1,
       dpi    = 1200,
       device = "png",
       width  = 4,
       height = 6,
       units  = "in")

# Plot
p1 <- ggplot(data = data, aes(y = x10_20_percent, x = 'DML', fill = feature)) + 
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

ggsave("plots/DML.png",
       plot   = p1,
       dpi    = 1200,
       device = "png",
       width  = 4,
       height = 6,
       units  = "in")

############################################################################################################
### DML-genome-feature-location Stacked bar plot

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
