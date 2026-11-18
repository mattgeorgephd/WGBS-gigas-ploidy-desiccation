# Title: Functional enrichment analysis w/ topGO - BOTH DATASETS
# Author: Matthew George; mattgeorgephd@gmail.com
# Date: 06/24/2021

## clear
rm(list=ls())

# # First Time? - install these
# install.packages("BiocManager")
# install.packages("SparseM")
# install.packages("magick")
# install.packages("gclus")
# BiocManager::install(version = "3.13")
# BiocManager::install("BiocGenerics")
# BiocManager::install("graph")
# BiocManager::install("Biobase")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("GO.db")
# BiocManager::install("topGO", force = TRUE)
# BiocManager::install("simplifyEnrichment") #Needs newest R Version (4.1.0)
# BiocManager::install("rrvgo", force = TRUE)

options("install.lock"=FALSE)

## Load R packages
library(rstudioapi)
library(readxl)
require(tidyverse)
library(BiocManager)
library(topGO)
library(simplifyEnrichment)
library(magick)
library(gclus)
library(rrvgo)

## Grab the WD from the file location and set it
home <- getActiveDocumentContext()$path
setwd(".."); getwd()

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
                  panel.border      = element_rect(color = "black", fill=NA, size=1.5),
                  legend.key        = element_blank()) # removes background of legend bullets



### NOTES ####################################################################################################

# Chandra Rajan et al. (2020) Methods:
# https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fgcb.15675&file=gcb15675-sup-0002-FileS2.txt

# Bioconductor SimplifyEnrichment Methods:
# http://www.bioconductor.org/packages/devel/bioc/vignettes/simplifyEnrichment/inst/doc/simplifyEnrichment.html

### LOAD ####################################################################################################
# Before launching into enrichment, I want to create lists of GOterms associated with each DML list for descriptive purposes.

setwd('E:/WGBS-gigas-ploidy-desiccation/'); getwd()

PLOIDYDMLGOterms <- read.table("data/gene_annotations/BOTH/DML-Cov10-20-BOTH-PLOIDY.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", 
                               header = FALSE, sep = "\t", 
                               col.names = c("transcript", "geneID", "chr", "start", "end","meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for ploidy DML
head(PLOIDYDMLGOterms)

HEATDMLGOterms <- read.table("data/gene_annotations/BOTH/DML-Cov10-20-BOTH-HEAT.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", 
                               header = FALSE, sep = "\t", 
                               col.names = c("transcript", "geneID", "chr", "start", "end","meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for ploidy DML
head(HEATDMLGOterms)

COMBODMLGOterms <- read.table("data/gene_annotations/BOTH/DML-Cov10-20-BOTH-COMBO.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", 
                               header = FALSE, sep = "\t", 
                               col.names = c("transcript", "geneID", "chr", "start", "end","meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for ploidy DML
head(COMBODMLGOterms)

# # Rearrange table
# ploidyDMLGOtermCounts <- as.data.frame(table(ploidyDMLGOterms$GOID)) #Create table with GOID and counts
# colnames(ploidyDMLGOtermCounts) <- c("GOID", "Freq") #Add column names
# range(ploidyDMLGOtermCounts$Freq) #Get range of GOID frequency
# nrow(ploidyDMLGOtermCounts) #Number of unique GOIDs associated with ploidy-DML


####### Sumarize and filter GO Terms  ###########
# Filter by BP, CC, or MF
PLOIDYDMLGOtermsBP <- PLOIDYDMLGOterms %>% filter(., PLOIDYDMLGOterms$GOcat == "P")
PLOIDYDMLGOtermsCC <- PLOIDYDMLGOterms %>% filter(., PLOIDYDMLGOterms$GOcat == "C")
PLOIDYDMLGOtermsMF <- PLOIDYDMLGOterms %>% filter(., PLOIDYDMLGOterms$GOcat == "F")

HEATDMLGOtermsBP <- HEATDMLGOterms %>% filter(., HEATDMLGOterms$GOcat == "P")
HEATDMLGOtermsCC <- HEATDMLGOterms %>% filter(., HEATDMLGOterms$GOcat == "C")
HEATDMLGOtermsMF <- HEATDMLGOterms %>% filter(., HEATDMLGOterms$GOcat == "F")

COMBODMLGOtermsBP <- COMBODMLGOterms %>% filter(., COMBODMLGOterms$GOcat == "P")
COMBODMLGOtermsCC <- COMBODMLGOterms %>% filter(., COMBODMLGOterms$GOcat == "C")
COMBODMLGOtermsMF <- COMBODMLGOterms %>% filter(., COMBODMLGOterms$GOcat == "F")

# Count GOIDs
PLOIDYDMLGOtermCountsBP <- as.data.frame(table(PLOIDYDMLGOtermsBP$GOID))
PLOIDYDMLGOtermCountsCC <- as.data.frame(table(PLOIDYDMLGOtermsCC$GOID))
PLOIDYDMLGOtermCountsMF <- as.data.frame(table(PLOIDYDMLGOtermsMF$GOID))

HEATDMLGOtermCountsBP <- as.data.frame(table(HEATDMLGOtermsBP$GOID))
HEATDMLGOtermCountsCC <- as.data.frame(table(HEATDMLGOtermsCC$GOID))
HEATDMLGOtermCountsMF <- as.data.frame(table(HEATDMLGOtermsMF$GOID))

COMBODMLGOtermCountsBP <- as.data.frame(table(COMBODMLGOtermsBP$GOID))
COMBODMLGOtermCountsCC <- as.data.frame(table(COMBODMLGOtermsCC$GOID))
COMBODMLGOtermCountsMF <- as.data.frame(table(COMBODMLGOtermsMF$GOID))

nrow(PLOIDYDMLGOtermCountsBP) #662 of unique GOIDs associated with PLOIDY-DML
nrow(PLOIDYDMLGOtermCountsCC) #178 of unique GOIDs associated with PLOIDY-DML
nrow(PLOIDYDMLGOtermCountsMF) #214 of unique GOIDs associated with PLOIDY-DML

nrow(HEATDMLGOtermCountsBP) #2224 of unique GOIDs associated with HEAT-DML
nrow(HEATDMLGOtermCountsCC) #460 of unique GOIDs associated with HEAT-DML
nrow(HEATDMLGOtermCountsMF) #703 of unique GOIDs associated with HEAT-DML

nrow(COMBODMLGOtermCountsBP) #517 of unique GOIDs associated with COMBO-DML
nrow(COMBODMLGOtermCountsCC) #151 of unique GOIDs associated with COMBO-DML
nrow(COMBODMLGOtermCountsMF) #158 of unique GOIDs associated with COMBO-DML

# data <- data.frame(
#   group=c("ploidy","heat","interaction"),
#   value=c(662,2224,517)
# )
# 
# ggplot(data, aes(x="", y=value, fill=group)) +
#   geom_bar(stat="identity", width=1, color="White") +
#   scale_fill_manual(values = c("orangered1","seagreen3","royalblue1")) +
#   coord_polar("y", start=0) +
#   my_theme
# 

# Add column names
colnames(PLOIDYDMLGOtermCountsBP) <- c("GOID", "Freq")
colnames(PLOIDYDMLGOtermCountsCC) <- c("GOID", "Freq")
colnames(PLOIDYDMLGOtermCountsMF) <- c("GOID", "Freq")

colnames(HEATDMLGOtermCountsBP) <- c("GOID", "Freq")
colnames(HEATDMLGOtermCountsCC) <- c("GOID", "Freq")
colnames(HEATDMLGOtermCountsMF) <- c("GOID", "Freq")

colnames(COMBODMLGOtermCountsBP) <- c("GOID", "Freq")
colnames(COMBODMLGOtermCountsCC) <- c("GOID", "Freq")
colnames(COMBODMLGOtermCountsMF) <- c("GOID", "Freq")

# Filter by frequency
PLOIDYDMLGOtermCountsBP <- PLOIDYDMLGOtermCountsBP %>% filter(PLOIDYDMLGOtermCountsBP$Freq != 0)
PLOIDYDMLGOtermCountsCC <- PLOIDYDMLGOtermCountsCC %>% filter(PLOIDYDMLGOtermCountsCC$Freq != 0)
PLOIDYDMLGOtermCountsMF <- PLOIDYDMLGOtermCountsMF %>% filter(PLOIDYDMLGOtermCountsMF$Freq != 0)

HEATDMLGOtermCountsBP <- HEATDMLGOtermCountsBP %>% filter(HEATDMLGOtermCountsBP$Freq != 0)
HEATDMLGOtermCountsCC <- HEATDMLGOtermCountsCC %>% filter(HEATDMLGOtermCountsCC$Freq != 0)
HEATDMLGOtermCountsMF <- HEATDMLGOtermCountsMF %>% filter(HEATDMLGOtermCountsMF$Freq != 0)

COMBODMLGOtermCountsBP <- COMBODMLGOtermCountsBP %>% filter(COMBODMLGOtermCountsBP$Freq != 0)
COMBODMLGOtermCountsCC <- COMBODMLGOtermCountsCC %>% filter(COMBODMLGOtermCountsCC$Freq != 0)
COMBODMLGOtermCountsMF <- COMBODMLGOtermCountsMF %>% filter(COMBODMLGOtermCountsMF$Freq != 0)

# Save GOterms associated with BP, CC, or MF
write.csv(PLOIDYDMLGOtermCountsBP, "GOterms/BOTH-PLOIDY/PLOIDY-BP-GOterms.csv", quote = FALSE, row.names = FALSE)
write.csv(PLOIDYDMLGOtermCountsCC, "GOterms/BOTH-PLOIDY/PLOIDY-CC-GOterms.csv", quote = FALSE, row.names = FALSE)
write.csv(PLOIDYDMLGOtermCountsMF, "GOterms/BOTH-PLOIDY/PLOIDY-MF-GOterms.csv", quote = FALSE, row.names = FALSE)

write.csv(HEATDMLGOtermCountsBP, "GOterms/BOTH-HEAT/HEAT-BP-GOterms.csv", quote = FALSE, row.names = FALSE)
write.csv(HEATDMLGOtermCountsCC, "GOterms/BOTH-HEAT/HEAT-CC-GOterms.csv", quote = FALSE, row.names = FALSE)
write.csv(HEATDMLGOtermCountsMF, "GOterms/BOTH-HEAT/HEAT-MF-GOterms.csv", quote = FALSE, row.names = FALSE)

write.csv(COMBODMLGOtermCountsBP, "GOterms/BOTH-COMBO/COMBO-BP-GOterms.csv", quote = FALSE, row.names = FALSE)
write.csv(COMBODMLGOtermCountsCC, "GOterms/BOTH-COMBO/COMBO-CC-GOterms.csv", quote = FALSE, row.names = FALSE)
write.csv(COMBODMLGOtermCountsMF, "GOterms/BOTH-COMBO/COMBO-MF-GOterms.csv", quote = FALSE, row.names = FALSE)

########### Load genetic background #################
# Load data for gene enrichment
geneID2GO <- topGO::readMappings(file = "data/gene_annotations/BOTH/geneid2go-union1x-BOTH.tab") #Loading the GO annotations and GeneIDs. Each line has one transcript ID and all associated GOterms
str(head(geneID2GO)) #Confirm file structure
geneNames <- names(geneID2GO) #Extract names to use as gene universe
head(geneNames)

########## Link DML to gene transcripts ##################
# I'm using the list of genes with DML as my list of interest. 
PLOIDYGenes <- PLOIDYDMLGOterms$transcript #Extract transcript ID
PLOIDYGeneList <- factor(as.integer(geneNames %in% PLOIDYGenes))
names(PLOIDYGeneList) <- geneNames #PLOIDYGenes
str(PLOIDYGeneList)

HEATGenes   <- HEATDMLGOterms$transcript #Extract transcript ID
HEATGeneList <- factor(as.integer(geneNames %in% HEATGenes))
names(HEATGeneList) <- geneNames #HEATGenes
str(HEATGeneList)

COMBOGenes  <- COMBODMLGOterms$transcript #Extract transcript ID
COMBOGeneList <- factor(as.integer(geneNames %in% COMBOGenes))
names(COMBOGeneList) <- geneNames #COMBOGenes
str(COMBOGeneList)

########## Run TopGo #########################
# Create topGOdata - Biological processes
PLOIDYGOdataBP <- new("topGOdata", ontology = "BP", allGenes = PLOIDYGeneList,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
PLOIDYGOdataBP #Get summary of object

HEATGOdataBP <- new("topGOdata", ontology = "BP", allGenes = HEATGeneList,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
HEATGOdataBP #Get summary of object

COMBOGOdataBP <- new("topGOdata", ontology = "BP", allGenes = COMBOGeneList,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
COMBOGOdataBP #Get summary of object

# Create topGOdata - Cellular Components
PLOIDYGOdataCC <- new("topGOdata", ontology = "CC", allGenes = PLOIDYGeneList,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
PLOIDYGOdataCC #Get summary of object

HEATGOdataCC <- new("topGOdata", ontology = "CC", allGenes = HEATGeneList,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
HEATGOdataCC #Get summary of object

COMBOGOdataCC <- new("topGOdata", ontology = "CC", allGenes = COMBOGeneList,
                     annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
COMBOGOdataCC #Get summary of object

# Create topGOdata - Molecular Functions
PLOIDYGOdataMF <- new("topGOdata", ontology = "MF", allGenes = PLOIDYGeneList,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
PLOIDYGOdataMF #Get summary of object

HEATGOdataMF <- new("topGOdata", ontology = "MF", allGenes = HEATGeneList,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
HEATGOdataMF #Get summary of object

COMBOGOdataMF <- new("topGOdata", ontology = "MF", allGenes = COMBOGeneList,
                     annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
COMBOGOdataMF #Get summary of object


############# Save Point 1 ################
# setwd('E:/WGBS-gigas-ploidy-desiccation/data/gene_annotations/BOTH')
# save.image("20211124_save_point_1.RData") #Save R Data in case R crashes
# # load("R_sessions/20210626_Functional-enrichment-analysis_1.RData") #Save R Data in case R crashes

############# Statistical Tests ################
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# PLOIDY
resultFisher.PLOIDYBP <- getSigGroups(PLOIDYGOdataBP, test.stat)
pvalFis.PLOIDYBP <- score(resultFisher.PLOIDYBP) #Extract p-values
allRes.PLOIDYBP <- GenTable(PLOIDYGOdataBP, classic = resultFisher.PLOIDYBP, ranksOf = "classic", orderBy = "classic", topNodes = length(pvalFis.PLOIDYBP)) #Create a statistical results table with statistical test results. Order by p-value (classic), and include all results (topNodes)
head(allRes.PLOIDYBP)

resultFisher.PLOIDYCC <- getSigGroups(PLOIDYGOdataCC, test.stat)
pvalFis.PLOIDYCC <- score(resultFisher.PLOIDYCC) #Extract p-values
allRes.PLOIDYCC <- GenTable(PLOIDYGOdataCC, classic = resultFisher.PLOIDYCC, ranksOf = "classic", orderBy = "classic", topNodes = length(pvalFis.PLOIDYCC)) #Create a statistical results table with statistical test results. Order by p-value (classic), and include all results (topNodes)
head(allRes.PLOIDYCC)

resultFisher.PLOIDYMF <- getSigGroups(PLOIDYGOdataMF, test.stat)
pvalFis.PLOIDYMF <- score(resultFisher.PLOIDYMF) #Extract p-values
allRes.PLOIDYMF <- GenTable(PLOIDYGOdataMF, classic = resultFisher.PLOIDYMF, ranksOf = "classic", orderBy = "classic", topNodes = length(pvalFis.PLOIDYMF)) #Create a statistical results table with statistical test results. Order by p-value (classic), and include all results (topNodes)
head(allRes.PLOIDYMF)

# HEAT
resultFisher.HEATBP <- getSigGroups(HEATGOdataBP, test.stat)
pvalFis.HEATBP <- score(resultFisher.HEATBP) #Extract p-values
allRes.HEATBP <- GenTable(HEATGOdataBP, classic = resultFisher.HEATBP, ranksOf = "classic", orderBy = "classic", topNodes = length(pvalFis.HEATBP)) #Create a statistical results table with statistical test results. Order by p-value (classic), and include all results (topNodes)
head(allRes.HEATBP)

resultFisher.HEATCC <- getSigGroups(HEATGOdataCC, test.stat)
pvalFis.HEATCC <- score(resultFisher.HEATCC) #Extract p-values
allRes.HEATCC <- GenTable(HEATGOdataCC, classic = resultFisher.HEATCC, ranksOf = "classic", orderBy = "classic", topNodes = length(pvalFis.HEATCC)) #Create a statistical results table with statistical test results. Order by p-value (classic), and include all results (topNodes)
head(allRes.HEATCC)

resultFisher.HEATMF <- getSigGroups(HEATGOdataMF, test.stat)
pvalFis.HEATMF <- score(resultFisher.HEATMF) #Extract p-values
allRes.HEATMF <- GenTable(HEATGOdataMF, classic = resultFisher.HEATMF, ranksOf = "classic", orderBy = "classic", topNodes = length(pvalFis.HEATMF)) #Create a statistical results table with statistical test results. Order by p-value (classic), and include all results (topNodes)
head(allRes.HEATMF)

# COMBO
resultFisher.COMBOBP <- getSigGroups(COMBOGOdataBP, test.stat)
pvalFis.COMBOBP <- score(resultFisher.COMBOBP) #Extract p-values
# resultFisher.COMBOBP
allRes.COMBOBP <- GenTable(COMBOGOdataBP, classic = resultFisher.COMBOBP, ranksOf = "classic", orderBy = "classic", topNodes = length(pvalFis.COMBOBP)) #Create a statistical results table with statistical test results. Order by p-value (classic), and include all results (topNodes)
head(allRes.COMBOBP)

resultFisher.COMBOCC <- getSigGroups(COMBOGOdataCC, test.stat)
pvalFis.COMBOCC <- score(resultFisher.COMBOCC) #Extract p-values
allRes.COMBOCC <- GenTable(COMBOGOdataCC, classic = resultFisher.COMBOCC, ranksOf = "classic", orderBy = "classic", topNodes = length(pvalFis.COMBOCC)) #Create a statistical results table with statistical test results. Order by p-value (classic), and include all results (topNodes)
head(allRes.COMBOCC)

resultFisher.COMBOMF <- getSigGroups(COMBOGOdataMF, test.stat)
pvalFis.COMBOMF <- score(resultFisher.COMBOMF) #Extract p-values
allRes.COMBOMF <- GenTable(COMBOGOdataMF, classic = resultFisher.COMBOMF, ranksOf = "classic", orderBy = "classic", topNodes = length(pvalFis.COMBOMF)) #Create a statistical results table with statistical test results. Order by p-value (classic), and include all results (topNodes)
head(allRes.COMBOMF)

# Save results
setwd('E:/WGBS-gigas-ploidy-desiccation/')
write.csv(allRes.PLOIDYBP, "GOterms/BOTH-PLOIDY/PLOIDY-BP-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe
write.csv(allRes.HEATBP,   "GOterms/BOTH-HEAT/HEAT-BP-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe
write.csv(allRes.COMBOBP,  "GOterms/BOTH-COMBO/COMBO-BP-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe

write.csv(allRes.PLOIDYCC, "GOterms/BOTH-PLOIDY/PLOIDY-CC-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe
write.csv(allRes.HEATCC,   "GOterms/BOTH-HEAT/HEAT-CC-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe
write.csv(allRes.COMBOCC,  "GOterms/BOTH-COMBO/COMBO-CC-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe

write.csv(allRes.PLOIDYMF, "GOterms/BOTH-PLOIDY/PLOIDY-MF-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe
write.csv(allRes.HEATMF,   "GOterms/BOTH-HEAT/HEAT-MF-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe
write.csv(allRes.COMBOMF,  "GOterms/BOTH-COMBO/COMBO-MF-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe

############# Match enriched GOterms with general annotation information and product ID #################
# Filter significantly enriched GOterms, only keep GOID, GOterm, and p-value
sigRes.PLOIDYBP <- allRes.PLOIDYBP[1:10,c(1, 2, 6)] 
colnames(sigRes.PLOIDYBP) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.HEATBP <- allRes.HEATBP[1:10,c(1, 2, 6)] 
colnames(sigRes.HEATBP) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.COMBOBP <- allRes.COMBOBP[1:10,c(1, 2, 6)] 
colnames(sigRes.COMBOBP) <- c("GOID", "GOterm", "p.value") #Change column names

sigRes.PLOIDYCC <- allRes.PLOIDYCC[1:10,c(1, 2, 6)] 
colnames(sigRes.PLOIDYCC) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.HEATCC <- allRes.HEATCC[1:10,c(1, 2, 6)] 
colnames(sigRes.HEATCC) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.COMBOCC <- allRes.COMBOCC[1:10,c(1, 2, 6)] 
colnames(sigRes.COMBOCC) <- c("GOID", "GOterm", "p.value") #Change column names

sigRes.PLOIDYMF <- allRes.PLOIDYMF[1:10,c(1, 2, 6)] 
colnames(sigRes.PLOIDYMF) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.HEATMF <- allRes.HEATMF[1:10,c(1, 2, 6)] 
colnames(sigRes.HEATMF) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.COMBOMF <- allRes.COMBOMF[1:10,c(1, 2, 6)] 
colnames(sigRes.COMBOMF) <- c("GOID", "GOterm", "p.value") #Change column names

############## Extract genes associated with each GOterm from topGO's annotations
## Analysis BP ###########
# PLOIDY
allDMLGOterms <- read.table("data/gene_annotations/BOTH/DML-Cov10-20-BOTH-PLOIDY.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", header = FALSE, sep = "\t", col.names = c("transcript", "geneID", "chr", "start", "end", "meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for all DML
head(allDMLGOterms)
nrow(allDMLGOterms) # 10979
allDMLGOtermsBP <- allDMLGOterms %>% filter(., allDMLGOterms$GOcat == "P") #Filter biological process terms
nrow(allDMLGOtermsBP) # 4379
sigRes.allBP <- allRes.PLOIDYBP[1:10,c(1, 2, 6)] #Filter significantly enriched GOterms, only keep GOID, GOterm, and p-value
colnames(sigRes.allBP) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.allBP

sigRes.allBPGenes <- genesInTerm(PLOIDYGOdataBP, whichGO = allRes.PLOIDYBP$GO.ID[1:10]) #Extract genes associated with each GOterm from topGO's annotations
sigRes.allBPGenes <- as.data.frame(vapply(sigRes.allBPGenes, paste, collapse = ",", character(1L))) #Collapse all transcript IDs into one column and convert into a dataframe. GO IDs will be the rownames
head(sigRes.allBPGenes)
colnames(sigRes.allBPGenes) <- c("transcript")
sigRes.allBPGenes$GOID <- rownames(sigRes.allBPGenes) #Convert row names to a separate column
rownames(sigRes.allBPGenes) <- 1:nrow(sigRes.allBPGenes) #Convert row names to numeric values
sigRes.allBPGenes <- separate_rows(sigRes.allBPGenes, transcript, sep = ",") #Separate transcripts so each transcript is on a separate row
head(sigRes.allBPGenes) #Confirm formatting

sigRes.PLOIDYBPtopGO <- left_join(sigRes.allBP, sigRes.allBPGenes, by = "GOID") #Add transcript information from topGO to GOIDs
sum(is.na(sigRes.PLOIDYBPtopGO$transcript)) #sum = 0, there are no N/As and every GOID has a transcript
nrow(sigRes.PLOIDYBPtopGO) # 199
head(sigRes.PLOIDYBPtopGO) #Confirm formatting

sigRes.allBPAnnot <- left_join(sigRes.PLOIDYBPtopGO, allDMLGOtermsBP, by = c("GOID", "transcript")) #Additional annotations. Left join (merge) based on matching transcript IDs. There are some rows with transcript IDs without corresponding annotations because there was no matching GOterm from my annotations.
colnames(sigRes.allBPAnnot)[2] <- "topGOterm" #Rename to indicate GOterms from topGO
colnames(sigRes.allBPAnnot)[10] <- "GOterm" #Rename to indicate a GOterms from my annotation
sigRes.allBPAnnot <- unique(sigRes.allBPAnnot[,-12]) #Drop GOcat column and retain only unique lines 
nrow(sigRes.allBPAnnot) #Count the number of rows
length(na.omit(unique(sigRes.allBPAnnot$geneID))) #1 unique genes
length(na.omit(unique(sigRes.allBPAnnot$transcript))) #128 unique transcripts
tail(sigRes.allBPAnnot) #Confirm formatting 

write.csv(sigRes.allBPAnnot, "GOterms/BOTH-PLOIDY/all-BP-EnrichedGO-DML-withTranscript.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allBPAnnotnoT <- unique(sigRes.allBPAnnot[,-4]) #Drop transcript column and keep only unique rows
head(sigRes.allBPAnnotnoT) #Confirm formatting
# nrow(sigRes.allBPAnnotnoT)

write.csv(sigRes.allBPAnnotnoT, "GOterms/BOTH-PLOIDY/all-BP-EnrichedGO-DML.csv", quote = FALSE, row.names = FALSE) #Save dataframe

geneProductIDs <- read.delim("GFF/cgigas_uk_roslin_v1_rna_from_genomic_annot.transcript.tab", sep = "\t", header = FALSE, col.names = c("transcript", "product")) #Import transcript and gene product information
head(geneProductIDs) #Confirm import

sigRes.allBPProduct <- left_join(sigRes.allBPAnnot, geneProductIDs, by = "transcript") #Add product annotations and retain rows that may not have annotations
nrow(sigRes.allBPProduct) == nrow(sigRes.allBPAnnot) #Confirm all transcripts have an annotation
head(sigRes.allBPProduct) #Confirm formatting
nrow(sigRes.allBPProduct)

write.csv(sigRes.allBPProduct, "GOterms/BOTH-PLOIDY/all-BP-EnrichedGO-DML-Product.csv", quote = FALSE, row.names = FALSE) #Save dataframe

## Count the number of hyper- or hypomethylated DML in genes with enriched GOterms

sigRes.allBPMethDiff <- unique(sigRes.allBPProduct %>%
                                 dplyr::select(., transcript) %>%
                                 left_join(., PLOIDYDMLGOterms, by = "transcript") %>%
                                 filter(., GOcat == "P") %>%
                                 dplyr::select(., geneID, chr, start, end, meth.diff)) #Isolate transcript column; join with DML list; remove all rows that aren't BP GOterms, select gene ID/chr/start/end/meth.diff columns
head(sigRes.allBPMethDiff) #Confirm formatting


nrow(sigRes.allBPMethDiff) #Number of unqiue DML in genes with enriched GOterms
sum(sigRes.allBPMethDiff$meth.diff >= 0) #Number of unique DML hypermethylated - 5
sum(sigRes.allBPMethDiff$meth.diff <= 0) #Number of unique DML hypomethylated - 5

write.csv(sigRes.allBPMethDiff, "GOterms/BOTH-PLOIDY/GENELIST-BP-EnrichedGO-DMLOnly.csv", quote = FALSE, row.names = FALSE) #Save dataframe as csv

# HEAT
allDMLGOterms <- read.table("data/gene_annotations/BOTH/DML-Cov10-20-BOTH-HEAT.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", header = FALSE, sep = "\t", col.names = c("transcript", "geneID", "chr", "start", "end", "meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for all DML
head(allDMLGOterms)
nrow(allDMLGOterms) # 81109
allDMLGOtermsBP <- allDMLGOterms %>% filter(., allDMLGOterms$GOcat == "P") #Filter biological process terms
nrow(allDMLGOtermsBP) # 39688
sigRes.allBP <- allRes.HEATBP[1:10,c(1, 2, 6)] #Filter significantly enriched GOterms, only keep GOID, GOterm, and p-value
colnames(sigRes.allBP) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.allBP

sigRes.allBPGenes <- genesInTerm(HEATGOdataBP, whichGO = allRes.HEATBP$GO.ID[1:10]) #Extract genes associated with each GOterm from topGO's annotations
sigRes.allBPGenes <- as.data.frame(vapply(sigRes.allBPGenes, paste, collapse = ",", character(1L))) #Collapse all transcript IDs into one column and convert into a dataframe. GO IDs will be the rownames
# head(sigRes.allBPGenes)
colnames(sigRes.allBPGenes) <- c("transcript")
sigRes.allBPGenes$GOID <- rownames(sigRes.allBPGenes) #Convert row names to a separate column
rownames(sigRes.allBPGenes) <- 1:nrow(sigRes.allBPGenes) #Convert row names to numeric values
sigRes.allBPGenes <- separate_rows(sigRes.allBPGenes, transcript, sep = ",") #Separate transcripts so each transcript is on a separate row
head(sigRes.allBPGenes) #Confirm formatting

sigRes.HEATBPtopGO <- left_join(sigRes.allBP, sigRes.allBPGenes, by = "GOID") #Add transcript information from topGO to GOIDs
sum(is.na(sigRes.HEATBPtopGO$transcript)) #sum = 0, there are no N/As and every GOID has a transcript
nrow(sigRes.HEATBPtopGO) #Count the number of rows
head(sigRes.HEATBPtopGO) #Confirm formatting

sigRes.allBPAnnot <- left_join(sigRes.HEATBPtopGO, allDMLGOtermsBP, by = c("GOID", "transcript")) #Additional annotations. Left join (merge) based on matching transcript IDs. There are some rows with transcript IDs without corresponding annotations because there was no matching GOterm from my annotations.
colnames(sigRes.allBPAnnot)[2] <- "topGOterm" #Rename to indicate GOterms from topGO
colnames(sigRes.allBPAnnot)[10] <- "GOterm" #Rename to indicate a GOterms from my annotation
sigRes.allBPAnnot <- unique(sigRes.allBPAnnot[,-12]) #Drop GOcat column and retain only unique lines 
nrow(sigRes.allBPAnnot) #Count the number of rows
length(na.omit(unique(sigRes.allBPAnnot$geneID))) #3 unique genes
length(na.omit(unique(sigRes.allBPAnnot$transcript))) #93 unique transcripts
tail(sigRes.allBPAnnot) #Confirm formatting 

write.csv(sigRes.allBPAnnot, "GOterms/BOTH-HEAT/all-BP-EnrichedGO-DML-withTranscript.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allBPAnnotnoT <- unique(sigRes.allBPAnnot[,-4]) #Drop transcript column and keep only unique rows
nrow(sigRes.allBPAnnotnoT)
head(sigRes.allBPAnnotnoT) #Confirm formatting

write.csv(sigRes.allBPAnnotnoT, "GOterms/BOTH-HEAT/all-BP-EnrichedGO-DML.csv", quote = FALSE, row.names = FALSE) #Save dataframe

geneProductIDs <- read.delim("GFF/cgigas_uk_roslin_v1_rna_from_genomic_annot.transcript.tab", sep = "\t", header = FALSE, col.names = c("transcript", "product")) #Import transcript and gene product information
head(geneProductIDs) #Confirm import

sigRes.allBPProduct <- left_join(sigRes.allBPAnnot, geneProductIDs, by = "transcript") #Add product annotations and retain rows that may not have annotations
nrow(sigRes.allBPProduct) == nrow(sigRes.allBPAnnot) #Confirm all transcripts have an annotation
nrow(sigRes.allBPProduct)
head(sigRes.allBPProduct) #Confirm formatting

write.csv(sigRes.allBPProduct, "GOterms/BOTH-HEAT/all-BP-EnrichedGO-DML-Product.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allBPMethDiff <- unique(sigRes.allBPProduct %>%
                                 dplyr::select(., transcript) %>%
                                 left_join(., HEATDMLGOterms, by = "transcript") %>%
                                 filter(., GOcat == "P") %>%
                                 dplyr::select(., geneID, chr, start, end, meth.diff)) #Isolate transcript column; join with DML list; remove all rows that aren't BP GOterms, select gene ID/chr/start/end/meth.diff columns
head(sigRes.allBPMethDiff) #Confirm formatting


nrow(sigRes.allBPMethDiff) #Number of unqiue DML in genes with enriched GOterms - 7
sum(sigRes.allBPMethDiff$meth.diff >= 0) #Number of unique DML hypermethylated - 2
sum(sigRes.allBPMethDiff$meth.diff <= 0) #Number of unique DML hypomethylated - 5

write.csv(sigRes.allBPMethDiff, "GOterms/BOTH-HEAT/GENELIST-BP-EnrichedGO-DMLOnly.csv", quote = FALSE, row.names = FALSE) #Save dataframe as csv

# COMBO
allDMLGOterms <- read.table("data/gene_annotations/BOTH/DML-Cov10-20-BOTH-COMBO.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", header = FALSE, sep = "\t", col.names = c("transcript", "geneID", "chr", "start", "end", "meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for all DML
head(allDMLGOterms)
nrow(allDMLGOterms) # 6941
allDMLGOtermsBP <- allDMLGOterms %>% filter(., allDMLGOterms$GOcat == "P") #Filter biological process terms
nrow(allDMLGOtermsBP) # 3332
sigRes.allBP <- allRes.COMBOBP[1:10,c(1, 2, 6)] #Filter significantly enriched GOterms, only keep GOID, GOterm, and p-value
colnames(sigRes.allBP) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.allBP

sigRes.allBPGenes <- genesInTerm(COMBOGOdataBP, whichGO = allRes.COMBOBP$GO.ID[1:10]) #Extract genes associated with each GOterm from topGO's annotations
sigRes.allBPGenes <- as.data.frame(vapply(sigRes.allBPGenes, paste, collapse = ",", character(1L))) #Collapse all transcript IDs into one column and convert into a dataframe. GO IDs will be the rownames
# head(sigRes.allBPGenes)
colnames(sigRes.allBPGenes) <- c("transcript")
sigRes.allBPGenes$GOID <- rownames(sigRes.allBPGenes) #Convert row names to a separate column
rownames(sigRes.allBPGenes) <- 1:nrow(sigRes.allBPGenes) #Convert row names to numeric values
sigRes.allBPGenes <- separate_rows(sigRes.allBPGenes, transcript, sep = ",") #Separate transcripts so each transcript is on a separate row
head(sigRes.allBPGenes) #Confirm formatting

sigRes.COMBOBPtopGO <- left_join(sigRes.allBP, sigRes.allBPGenes, by = "GOID") #Add transcript information from topGO to GOIDs
sum(is.na(sigRes.COMBOBPtopGO$transcript)) #sum = 0, there are no N/As and every GOID has a transcript
nrow(sigRes.COMBOBPtopGO) #Count the number of rows
head(sigRes.COMBOBPtopGO) #Confirm formatting

sigRes.allBPAnnot <- left_join(sigRes.COMBOBPtopGO, allDMLGOtermsBP, by = c("GOID", "transcript")) #Additional annotations. Left join (merge) based on matching transcript IDs. There are some rows with transcript IDs without corresponding annotations because there was no matching GOterm from my annotations.
colnames(sigRes.allBPAnnot)[2] <- "topGOterm" #Rename to indicate GOterms from topGO
colnames(sigRes.allBPAnnot)[10] <- "GOterm" #Rename to indicate a GOterms from my annotation
sigRes.allBPAnnot <- unique(sigRes.allBPAnnot[,-12]) #Drop GOcat column and retain only unique lines 
nrow(sigRes.allBPAnnot) #Count the number of rows
length(na.omit(unique(sigRes.allBPAnnot$geneID))) #1 unique genes
length(na.omit(unique(sigRes.allBPAnnot$transcript))) #55 unique transcripts
tail(sigRes.allBPAnnot) #Confirm formatting 

write.csv(sigRes.allBPAnnot, "GOterms/BOTH-COMBO/all-BP-EnrichedGO-DML-withTranscript.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allBPAnnotnoT <- unique(sigRes.allBPAnnot[,-4]) #Drop transcript column and keep only unique rows
head(sigRes.allBPAnnotnoT) #Confirm formatting

write.csv(sigRes.allBPAnnotnoT, "GOterms/BOTH-COMBO/all-BP-EnrichedGO-DML.csv", quote = FALSE, row.names = FALSE) #Save dataframe

geneProductIDs <- read.delim("GFF/cgigas_uk_roslin_v1_rna_from_genomic_annot.transcript.tab", sep = "\t", header = FALSE, col.names = c("transcript", "product")) #Import transcript and gene product information
head(geneProductIDs) #Confirm import

sigRes.allBPProduct <- left_join(sigRes.allBPAnnot, geneProductIDs, by = "transcript") #Add product annotations and retain rows that may not have annotations
nrow(sigRes.allBPProduct) == nrow(sigRes.allBPAnnot) #Confirm all transcripts have an annotation
nrow(sigRes.allBPProduct)
head(sigRes.allBPProduct) #Confirm formatting

write.csv(sigRes.allBPProduct, "GOterms/BOTH-COMBO/all-BP-EnrichedGO-DML-Product.csv", quote = FALSE, row.names = FALSE) #Save dataframe

## Count the number of hyper- or hypomethylated DML in genes with enriched GOterms

sigRes.allBPMethDiff <- unique(sigRes.allBPProduct %>%
                                 dplyr::select(., transcript) %>%
                                 left_join(., COMBODMLGOterms, by = "transcript") %>%
                                 filter(., GOcat == "P") %>%
                                 dplyr::select(., geneID, chr, start, end, meth.diff)) #Isolate transcript column; join with DML list; remove all rows that aren't BP GOterms, select gene ID/chr/start/end/meth.diff columns
head(sigRes.allBPMethDiff) #Confirm formatting


nrow(sigRes.allBPMethDiff) #Number of unqiue DML in genes with enriched GOterms - 3
sum(sigRes.allBPMethDiff$meth.diff >= 0) #Number of unique DML hypermethylated - 2
sum(sigRes.allBPMethDiff$meth.diff <= 0) #Number of unique DML hypomethylated - 1

write.csv(sigRes.allBPMethDiff, "GOterms/BOTH-COMBO/GENELIST-BP-EnrichedGO-DMLOnly.csv", quote = FALSE, row.names = FALSE) #Save dataframe as csv

## Analysis CC ###########
# PLOIDY
allDMLGOterms <- read.table("data/gene_annotations/BOTH/DML-Cov10-20-BOTH-PLOIDY.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", header = FALSE, sep = "\t", col.names = c("transcript", "geneID", "chr", "start", "end", "meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for all DML
head(allDMLGOterms)
nrow(allDMLGOterms) # 10979
allDMLGOtermsCC <- allDMLGOterms %>% filter(., allDMLGOterms$GOcat == "C") #Filter biological process terms
nrow(allDMLGOtermsCC) # 4145
sigRes.allCC <- allRes.PLOIDYCC[1:10,c(1, 2, 6)] #Filter significantly enriched GOterms, only keep GOID, GOterm, and p-value
colnames(sigRes.allCC) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.allCC

sigRes.allCCGenes <- genesInTerm(PLOIDYGOdataCC, whichGO = allRes.PLOIDYCC$GO.ID[1:10]) #Extract genes associated with each GOterm from topGO's annotations
sigRes.allCCGenes <- as.data.frame(vapply(sigRes.allCCGenes, paste, collapse = ",", character(1L))) #Collapse all transcript IDs into one column and convert into a dataframe. GO IDs will be the rownames
# head(sigRes.allCCGenes)
colnames(sigRes.allCCGenes) <- c("transcript")
sigRes.allCCGenes$GOID <- rownames(sigRes.allCCGenes) #Convert row names to a separate column
rownames(sigRes.allCCGenes) <- 1:nrow(sigRes.allCCGenes) #Convert row names to numeric values
sigRes.allCCGenes <- separate_rows(sigRes.allCCGenes, transcript, sep = ",") #Separate transcripts so each transcript is on a separate row
head(sigRes.allCCGenes) #Confirm formatting

sigRes.PLOIDYCCtopGO <- left_join(sigRes.allCC, sigRes.allCCGenes, by = "GOID") #Add transcript information from topGO to GOIDs
sum(is.na(sigRes.PLOIDYCCtopGO$transcript)) #sum = 0, there are no N/As and every GOID has a transcript
nrow(sigRes.PLOIDYCCtopGO) #Count the number of rows - 1083
head(sigRes.PLOIDYCCtopGO) #Confirm formatting

sigRes.allCCAnnot <- left_join(sigRes.PLOIDYCCtopGO, allDMLGOtermsCC, by = c("GOID", "transcript")) #Additional annotations. Left join (merge) based on matching transcript IDs. There are some rows with transcript IDs without corresponding annotations because there was no matching GOterm from my annotations.
colnames(sigRes.allCCAnnot)[2] <- "topGOterm" #Rename to indicate GOterms from topGO
colnames(sigRes.allCCAnnot)[10] <- "GOterm" #Rename to indicate a GOterms from my annotation
sigRes.allCCAnnot <- unique(sigRes.allCCAnnot[,-12]) #Drop GOcat column and retain only unique lines 
nrow(sigRes.allCCAnnot) #Count the number of rows - 1420
length(na.omit(unique(sigRes.allCCAnnot$geneID))) #12 unique genes
length(na.omit(unique(sigRes.allCCAnnot$transcript))) #601 unique transcripts
# tail(sigRes.allCCAnnot) #Confirm formatting 

write.csv(sigRes.allCCAnnot, "GOterms/BOTH-PLOIDY/all-CC-EnrichedGO-DML-withTranscript.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allCCAnnotnoT <- unique(sigRes.allCCAnnot[,-4]) #Drop transcript column and keep only unique rows
head(sigRes.allCCAnnotnoT) #Confirm formatting
# nrow(sigRes.allCCAnnotnoT)

write.csv(sigRes.allCCAnnotnoT, "GOterms/BOTH-PLOIDY/all-CC-EnrichedGO-DML.csv", quote = FALSE, row.names = FALSE) #Save dataframe

geneProductIDs <- read.delim("GFF/cgigas_uk_roslin_v1_rna_from_genomic_annot.transcript.tab", sep = "\t", header = FALSE, col.names = c("transcript", "product")) #Import transcript and gene product information
head(geneProductIDs) #Confirm import

sigRes.allCCProduct <- left_join(sigRes.allCCAnnot, geneProductIDs, by = "transcript") #Add product annotations and retain rows that may not have annotations
nrow(sigRes.allCCProduct) == nrow(sigRes.allCCAnnot) #Confirm all transcripts have an annotation
head(sigRes.allCCProduct) #Confirm formatting
nrow(sigRes.allCCProduct)

write.csv(sigRes.allCCProduct, "GOterms/BOTH-PLOIDY/all-CC-EnrichedGO-DML-Product.csv", quote = FALSE, row.names = FALSE) #Save dataframe

## Count the number of hyper- or hypomethylated DML in genes with enriched GOterms

sigRes.allCCMethDiff <- unique(sigRes.allCCProduct %>%
                                 dplyr::select(., transcript) %>%
                                 left_join(., PLOIDYDMLGOterms, by = "transcript") %>%
                                 filter(., GOcat == "C") %>%
                                 dplyr::select(., geneID, chr, start, end, meth.diff)) #Isolate transcript column; join with DML list; remove all rows that aren't CC GOterms, select gene ID/chr/start/end/meth.diff columns
head(sigRes.allCCMethDiff) #Confirm formatting


nrow(sigRes.allCCMethDiff) #Number of unqiue DML in genes with enriched GOterms - 24
sum(sigRes.allCCMethDiff$meth.diff >= 0) #Number of unique DML hypermethylated - 12
sum(sigRes.allCCMethDiff$meth.diff <= 0) #Number of unique DML hypomethylated - 12

write.csv(sigRes.allCCMethDiff, "GOterms/BOTH-PLOIDY/GENELIST-CC-EnrichedGO-DMLOnly.csv", quote = FALSE, row.names = FALSE) #Save dataframe as csv

# HEAT
allDMLGOterms <- read.table("data/gene_annotations/BOTH/DML-Cov10-20-BOTH-HEAT.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", header = FALSE, sep = "\t", col.names = c("transcript", "geneID", "chr", "start", "end", "meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for all DML
head(allDMLGOterms)
nrow(allDMLGOterms) # 81109
allDMLGOtermsCC <- allDMLGOterms %>% filter(., allDMLGOterms$GOcat == "C") #Filter biological process terms
nrow(allDMLGOtermsCC) # 26084
sigRes.allCC <- allRes.HEATCC[1:10,c(1, 2, 6)] #Filter significantly enriched GOterms, only keep GOID, GOterm, and p-value
colnames(sigRes.allCC) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.allCC

sigRes.allCCGenes <- genesInTerm(HEATGOdataCC, whichGO = allRes.HEATCC$GO.ID[1:10]) #Extract genes associated with each GOterm from topGO's annotations
sigRes.allCCGenes <- as.data.frame(vapply(sigRes.allCCGenes, paste, collapse = ",", character(1L))) #Collapse all transcript IDs into one column and convert into a dataframe. GO IDs will be the rownames
# head(sigRes.allCCGenes)
colnames(sigRes.allCCGenes) <- c("transcript")
sigRes.allCCGenes$GOID <- rownames(sigRes.allCCGenes) #Convert row names to a separate column
rownames(sigRes.allCCGenes) <- 1:nrow(sigRes.allCCGenes) #Convert row names to numeric values
sigRes.allCCGenes <- separate_rows(sigRes.allCCGenes, transcript, sep = ",") #Separate transcripts so each transcript is on a separate row
head(sigRes.allCCGenes) #Confirm formatting

sigRes.HEATCCtopGO <- left_join(sigRes.allCC, sigRes.allCCGenes, by = "GOID") #Add transcript information from topGO to GOIDs
sum(is.na(sigRes.HEATCCtopGO$transcript)) #sum = 0, there are no N/As and every GOID has a transcript
nrow(sigRes.HEATCCtopGO) #Count the number of rows
head(sigRes.HEATCCtopGO) #Confirm formatting

sigRes.allCCAnnot <- left_join(sigRes.HEATCCtopGO, allDMLGOtermsCC, by = c("GOID", "transcript")) #Additional annotations. Left join (merge) based on matching transcript IDs. There are some rows with transcript IDs without corresponding annotations because there was no matching GOterm from my annotations.
colnames(sigRes.allCCAnnot)[2] <- "topGOterm" #Rename to indicate GOterms from topGO
colnames(sigRes.allCCAnnot)[10] <- "GOterm" #Rename to indicate a GOterms from my annotation
sigRes.allCCAnnot <- unique(sigRes.allCCAnnot[,-12]) #Drop GOcat column and retain only unique lines 
nrow(sigRes.allCCAnnot) #1698 Count the number of rows
length(na.omit(unique(sigRes.allCCAnnot$geneID))) #17 unique genes
length(na.omit(unique(sigRes.allCCAnnot$transcript))) #727 unique transcripts
# tail(sigRes.allCCAnnot) #Confirm formatting 

write.csv(sigRes.allCCAnnot, "GOterms/BOTH-HEAT/all-CC-EnrichedGO-DML-withTranscript.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allCCAnnotnoT <- unique(sigRes.allCCAnnot[,-4]) #Drop transcript column and keep only unique rows
nrow(sigRes.allCCAnnotnoT)
head(sigRes.allCCAnnotnoT) #Confirm formatting

write.csv(sigRes.allCCAnnotnoT, "GOterms/BOTH-HEAT/all-CC-EnrichedGO-DML.csv", quote = FALSE, row.names = FALSE) #Save dataframe

geneProductIDs <- read.delim("GFF/cgigas_uk_roslin_v1_rna_from_genomic_annot.transcript.tab", sep = "\t", header = FALSE, col.names = c("transcript", "product")) #Import transcript and gene product information
head(geneProductIDs) #Confirm import

sigRes.allCCProduct <- left_join(sigRes.allCCAnnot, geneProductIDs, by = "transcript") #Add product annotations and retain rows that may not have annotations
nrow(sigRes.allCCProduct) == nrow(sigRes.allCCAnnot) #Confirm all transcripts have an annotation
nrow(sigRes.allCCProduct)
head(sigRes.allCCProduct) #Confirm formatting

write.csv(sigRes.allCCProduct, "GOterms/BOTH-HEAT/all-CC-EnrichedGO-DML-Product.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allCCMethDiff <- unique(sigRes.allCCProduct %>%
                                 dplyr::select(., transcript) %>%
                                 left_join(., HEATDMLGOterms, by = "transcript") %>%
                                 filter(., GOcat == "C") %>%
                                 dplyr::select(., geneID, chr, start, end, meth.diff)) #Isolate transcript column; join with DML list; remove all rows that aren't CC GOterms, select gene ID/chr/start/end/meth.diff columns
head(sigRes.allCCMethDiff) #Confirm formatting


nrow(sigRes.allCCMethDiff) #Number of unqiue DML in genes with enriched GOterms - 54
sum(sigRes.allCCMethDiff$meth.diff >= 0) #Number of unique DML hypermethylated - 30
sum(sigRes.allCCMethDiff$meth.diff <= 0) #Number of unique DML hypomethylated - 24

write.csv(sigRes.allCCMethDiff, "GOterms/BOTH-HEAT/GENELIST-CC-EnrichedGO-DMLOnly.csv", quote = FALSE, row.names = FALSE) #Save dataframe as csv

# COMBO
allDMLGOterms <- read.table("data/gene_annotations/BOTH/DML-Cov10-20-BOTH-COMBO.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", header = FALSE, sep = "\t", col.names = c("transcript", "geneID", "chr", "start", "end", "meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for all DML
head(allDMLGOterms)
nrow(allDMLGOterms) # 6941
allDMLGOtermsCC <- allDMLGOterms %>% filter(., allDMLGOterms$GOcat == "C") #Filter biological process terms
nrow(allDMLGOtermsCC) # 2337
sigRes.allCC <- allRes.COMBOCC[1:10,c(1, 2, 6)] #Filter significantly enriched GOterms, only keep GOID, GOterm, and p-value
colnames(sigRes.allCC) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.allCC

sigRes.allCCGenes <- genesInTerm(COMBOGOdataCC, whichGO = allRes.COMBOCC$GO.ID[1:10]) #Extract genes associated with each GOterm from topGO's annotations
sigRes.allCCGenes <- as.data.frame(vapply(sigRes.allCCGenes, paste, collapse = ",", character(1L))) #Collapse all transcript IDs into one column and convert into a dataframe. GO IDs will be the rownames
# head(sigRes.allCCGenes)
colnames(sigRes.allCCGenes) <- c("transcript")
sigRes.allCCGenes$GOID <- rownames(sigRes.allCCGenes) #Convert row names to a separate column
rownames(sigRes.allCCGenes) <- 1:nrow(sigRes.allCCGenes) #Convert row names to numeric values
sigRes.allCCGenes <- separate_rows(sigRes.allCCGenes, transcript, sep = ",") #Separate transcripts so each transcript is on a separate row
head(sigRes.allCCGenes) #Confirm formatting

sigRes.COMBOCCtopGO <- left_join(sigRes.allCC, sigRes.allCCGenes, by = "GOID") #Add transcript information from topGO to GOIDs
sum(is.na(sigRes.COMBOCCtopGO$transcript)) #sum = 0, there are no N/As and every GOID has a transcript
nrow(sigRes.COMBOCCtopGO) #Count the number of rows
head(sigRes.COMBOCCtopGO) #Confirm formatting

sigRes.allCCAnnot <- left_join(sigRes.COMBOCCtopGO, allDMLGOtermsCC, by = c("GOID", "transcript")) #Additional annotations. Left join (merge) based on matching transcript IDs. There are some rows with transcript IDs without corresponding annotations because there was no matching GOterm from my annotations.
colnames(sigRes.allCCAnnot)[2] <- "topGOterm" #Rename to indicate GOterms from topGO
colnames(sigRes.allCCAnnot)[10] <- "GOterm" #Rename to indicate a GOterms from my annotation
sigRes.allCCAnnot <- unique(sigRes.allCCAnnot[,-12]) #Drop GOcat column and retain only unique lines 
nrow(sigRes.allCCAnnot) #Count the number of rows
length(na.omit(unique(sigRes.allCCAnnot$geneID))) #1 unique genes
length(na.omit(unique(sigRes.allCCAnnot$transcript))) #55 unique transcripts
# tail(sigRes.allCCAnnot) #Confirm formatting 

write.csv(sigRes.allCCAnnot, "GOterms/BOTH-COMBO/all-CC-EnrichedGO-DML-withTranscript.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allCCAnnotnoT <- unique(sigRes.allCCAnnot[,-4]) #Drop transcript column and keep only unique rows
head(sigRes.allCCAnnotnoT) #Confirm formatting

write.csv(sigRes.allCCAnnotnoT, "GOterms/BOTH-COMBO/all-CC-EnrichedGO-DML.csv", quote = FALSE, row.names = FALSE) #Save dataframe

geneProductIDs <- read.delim("GFF/cgigas_uk_roslin_v1_rna_from_genomic_annot.transcript.tab", sep = "\t", header = FALSE, col.names = c("transcript", "product")) #Import transcript and gene product information
head(geneProductIDs) #Confirm import

sigRes.allCCProduct <- left_join(sigRes.allCCAnnot, geneProductIDs, by = "transcript") #Add product annotations and retain rows that may not have annotations
nrow(sigRes.allCCProduct) == nrow(sigRes.allCCAnnot) #Confirm all transcripts have an annotation
nrow(sigRes.allCCProduct)
head(sigRes.allCCProduct) #Confirm formatting

write.csv(sigRes.allCCProduct, "GOterms/BOTH-COMBO/all-CC-EnrichedGO-DML-Product.csv", quote = FALSE, row.names = FALSE) #Save dataframe

## Count the number of hyper- or hypomethylated DML in genes with enriched GOterms

sigRes.allCCMethDiff <- unique(sigRes.allCCProduct %>%
                                 dplyr::select(., transcript) %>%
                                 left_join(., COMBODMLGOterms, by = "transcript") %>%
                                 filter(., GOcat == "C") %>%
                                 dplyr::select(., geneID, chr, start, end, meth.diff)) #Isolate transcript column; join with DML list; remove all rows that aren't CC GOterms, select gene ID/chr/start/end/meth.diff columns
head(sigRes.allCCMethDiff) #Confirm formatting


nrow(sigRes.allCCMethDiff) #Number of unqiue DML in genes with enriched GOterms - 16
sum(sigRes.allCCMethDiff$meth.diff >= 0) #Number of unique DML hypermethylated - 8
sum(sigRes.allCCMethDiff$meth.diff <= 0) #Number of unique DML hypomethylated - 8

write.csv(sigRes.allCCMethDiff, "GOterms/BOTH-COMBO/GENELIST-CC-EnrichedGO-DMLOnly.csv", quote = FALSE, row.names = FALSE) #Save dataframe as csv

## Analysis MF ###########
# PLOIDY
allDMLGOterms <- read.table("data/gene_annotations/BOTH/DML-Cov10-20-BOTH-PLOIDY.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", header = FALSE, sep = "\t", col.names = c("transcript", "geneID", "chr", "start", "end", "meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for all DML
head(allDMLGOterms)
nrow(allDMLGOterms) # 10979
allDMLGOtermsMF <- allDMLGOterms %>% filter(., allDMLGOterms$GOcat == "F") #Filter biological process terms
nrow(allDMLGOtermsMF) # 2443
sigRes.allMF <- allRes.PLOIDYMF[1:10,c(1, 2, 6)] #Filter significantly enriched GOterms, only keep GOID, GOterm, and p-value
colnames(sigRes.allMF) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.allMF

sigRes.allMFGenes <- genesInTerm(PLOIDYGOdataMF, whichGO = allRes.PLOIDYMF$GO.ID[1:10]) #Extract genes associated with each GOterm from topGO's annotations
sigRes.allMFGenes <- as.data.frame(vapply(sigRes.allMFGenes, paste, collapse = ",", character(1L))) #Collapse all transcript IDs into one column and convert into a dataframe. GO IDs will be the rownames
# head(sigRes.allMFGenes)
colnames(sigRes.allMFGenes) <- c("transcript")
sigRes.allMFGenes$GOID <- rownames(sigRes.allMFGenes) #Convert row names to a separate column
rownames(sigRes.allMFGenes) <- 1:nrow(sigRes.allMFGenes) #Convert row names to numeric values
sigRes.allMFGenes <- separate_rows(sigRes.allMFGenes, transcript, sep = ",") #Separate transcripts so each transcript is on a separate row
head(sigRes.allMFGenes) #Confirm formatting

sigRes.PLOIDYMFtopGO <- left_join(sigRes.allMF, sigRes.allMFGenes, by = "GOID") #Add transcript information from topGO to GOIDs
sum(is.na(sigRes.PLOIDYMFtopGO$transcript)) #sum = 0, there are no N/As and every GOID has a transcript
nrow(sigRes.PLOIDYMFtopGO) #Count the number of rows
head(sigRes.PLOIDYMFtopGO) #Confirm formatting

sigRes.allMFAnnot <- left_join(sigRes.PLOIDYMFtopGO, allDMLGOtermsMF, by = c("GOID", "transcript")) #Additional annotations. Left join (merge) based on matching transcript IDs. There are some rows with transcript IDs without corresponding annotations because there was no matching GOterm from my annotations.
colnames(sigRes.allMFAnnot)[2] <- "topGOterm" #Rename to indicate GOterms from topGO
colnames(sigRes.allMFAnnot)[10] <- "GOterm" #Rename to indicate a GOterms from my annotation
sigRes.allMFAnnot <- unique(sigRes.allMFAnnot[,-12]) #Drop GOcat column and retain only unique lines 
nrow(sigRes.allMFAnnot) # 392 Count the number of rows
length(na.omit(unique(sigRes.allMFAnnot$geneID))) #6 unique genes
length(na.omit(unique(sigRes.allMFAnnot$transcript))) #178 unique transcripts
tail(sigRes.allMFAnnot) #Confirm formatting 

write.csv(sigRes.allMFAnnot, "GOterms/BOTH-PLOIDY/all-MF-EnrichedGO-DML-withTranscript.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allMFAnnotnoT <- unique(sigRes.allMFAnnot[,-4]) #Drop transcript column and keep only unique rows
head(sigRes.allMFAnnotnoT) #Confirm formatting
# nrow(sigRes.allMFAnnotnoT)

write.csv(sigRes.allMFAnnotnoT, "GOterms/BOTH-PLOIDY/all-MF-EnrichedGO-DML.csv", quote = FALSE, row.names = FALSE) #Save dataframe

geneProductIDs <- read.delim("GFF/cgigas_uk_roslin_v1_rna_from_genomic_annot.transcript.tab", sep = "\t", header = FALSE, col.names = c("transcript", "product")) #Import transcript and gene product information
head(geneProductIDs) #Confirm import

sigRes.allMFProduct <- left_join(sigRes.allMFAnnot, geneProductIDs, by = "transcript") #Add product annotations and retain rows that may not have annotations
nrow(sigRes.allMFProduct) == nrow(sigRes.allMFAnnot) #Confirm all transcripts have an annotation
head(sigRes.allMFProduct) #Confirm formatting
nrow(sigRes.allMFProduct)

write.csv(sigRes.allMFProduct, "GOterms/BOTH-PLOIDY/all-MF-EnrichedGO-DML-Product.csv", quote = FALSE, row.names = FALSE) #Save dataframe

## Count the number of hyper- or hypomethylated DML in genes with enriched GOterms

sigRes.allMFMethDiff <- unique(sigRes.allMFProduct %>%
                                 dplyr::select(., transcript) %>%
                                 left_join(., PLOIDYDMLGOterms, by = "transcript") %>%
                                 filter(., GOcat == "F") %>%
                                 dplyr::select(., geneID, chr, start, end, meth.diff)) #Isolate transcript column; join with DML list; remove all rows that aren't MF GOterms, select gene ID/chr/start/end/meth.diff columns
head(sigRes.allMFMethDiff) #Confirm formatting


nrow(sigRes.allMFMethDiff) #Number of unqiue DML in genes with enriched GOterms - 15
sum(sigRes.allMFMethDiff$meth.diff >= 0) #Number of unique DML hypermethylated - 5
sum(sigRes.allMFMethDiff$meth.diff <= 0) #Number of unique DML hypomethylated - 10

write.csv(sigRes.allMFMethDiff, "GOterms/BOTH-PLOIDY/GENELIST-MF-EnrichedGO-DMLOnly.csv", quote = FALSE, row.names = FALSE) #Save dataframe as csv

# HEAT
allDMLGOterms <- read.table("data/gene_annotations/BOTH/DML-Cov10-20-BOTH-HEAT.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", header = FALSE, sep = "\t", col.names = c("transcript", "geneID", "chr", "start", "end", "meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for all DML
head(allDMLGOterms)
nrow(allDMLGOterms) # 81109
allDMLGOtermsMF <- allDMLGOterms %>% filter(., allDMLGOterms$GOcat == "F") #Filter biological process terms
nrow(allDMLGOtermsMF) # 15213
sigRes.allMF <- allRes.HEATMF[1:10,c(1, 2, 6)] #Filter significantly enriched GOterms, only keep GOID, GOterm, and p-value
colnames(sigRes.allMF) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.allMF

sigRes.allMFGenes <- genesInTerm(HEATGOdataMF, whichGO = allRes.HEATMF$GO.ID[1:10]) #Extract genes associated with each GOterm from topGO's annotations
sigRes.allMFGenes <- as.data.frame(vapply(sigRes.allMFGenes, paste, collapse = ",", character(1L))) #Collapse all transcript IDs into one column and convert into a dataframe. GO IDs will be the rownames
# head(sigRes.allMFGenes)
colnames(sigRes.allMFGenes) <- c("transcript")
sigRes.allMFGenes$GOID <- rownames(sigRes.allMFGenes) #Convert row names to a separate column
rownames(sigRes.allMFGenes) <- 1:nrow(sigRes.allMFGenes) #Convert row names to numeric values
sigRes.allMFGenes <- separate_rows(sigRes.allMFGenes, transcript, sep = ",") #Separate transcripts so each transcript is on a separate row
head(sigRes.allMFGenes) #Confirm formatting

sigRes.HEATMFtopGO <- left_join(sigRes.allMF, sigRes.allMFGenes, by = "GOID") #Add transcript information from topGO to GOIDs
sum(is.na(sigRes.HEATMFtopGO$transcript)) #sum = 0, there are no N/As and every GOID has a transcript
nrow(sigRes.HEATMFtopGO) #Count the number of rows
head(sigRes.HEATMFtopGO) #Confirm formatting

sigRes.allMFAnnot <- left_join(sigRes.HEATMFtopGO, allDMLGOtermsMF, by = c("GOID", "transcript")) #Additional annotations. Left join (merge) based on matching transcript IDs. There are some rows with transcript IDs without corresponding annotations because there was no matching GOterm from my annotations.
colnames(sigRes.allMFAnnot)[2] <- "topGOterm" #Rename to indicate GOterms from topGO
colnames(sigRes.allMFAnnot)[10] <- "GOterm" #Rename to indicate a GOterms from my annotation
sigRes.allMFAnnot <- unique(sigRes.allMFAnnot[,-12]) #Drop GOcat column and retain only unique lines 
nrow(sigRes.allMFAnnot) # 5801 Count the number of rows
length(na.omit(unique(sigRes.allMFAnnot$geneID))) #49 unique genes
length(na.omit(unique(sigRes.allMFAnnot$transcript))) #1841 unique transcripts
tail(sigRes.allMFAnnot) #Confirm formatting 

write.csv(sigRes.allMFAnnot, "GOterms/BOTH-HEAT/all-MF-EnrichedGO-DML-withTranscript.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allMFAnnotnoT <- unique(sigRes.allMFAnnot[,-4]) #Drop transcript column and keep only unique rows
nrow(sigRes.allMFAnnotnoT)
head(sigRes.allMFAnnotnoT) #Confirm formatting

write.csv(sigRes.allMFAnnotnoT, "GOterms/BOTH-HEAT/all-MF-EnrichedGO-DML.csv", quote = FALSE, row.names = FALSE) #Save dataframe

geneProductIDs <- read.delim("GFF/cgigas_uk_roslin_v1_rna_from_genomic_annot.transcript.tab", sep = "\t", header = FALSE, col.names = c("transcript", "product")) #Import transcript and gene product information
head(geneProductIDs) #Confirm import

sigRes.allMFProduct <- left_join(sigRes.allMFAnnot, geneProductIDs, by = "transcript") #Add product annotations and retain rows that may not have annotations
nrow(sigRes.allMFProduct) == nrow(sigRes.allMFAnnot) #Confirm all transcripts have an annotation
nrow(sigRes.allMFProduct)
head(sigRes.allMFProduct) #Confirm formatting

write.csv(sigRes.allMFProduct, "GOterms/BOTH-HEAT/all-MF-EnrichedGO-DML-Product.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allMFMethDiff <- unique(sigRes.allMFProduct %>%
                                 dplyr::select(., transcript) %>%
                                 left_join(., HEATDMLGOterms, by = "transcript") %>%
                                 filter(., GOcat == "F") %>%
                                 dplyr::select(., geneID, chr, start, end, meth.diff)) #Isolate transcript column; join with DML list; remove all rows that aren't MF GOterms, select gene ID/chr/start/end/meth.diff columns
head(sigRes.allMFMethDiff) #Confirm formatting


nrow(sigRes.allMFMethDiff) #Number of unqiue DML in genes with enriched GOterms - 134
sum(sigRes.allMFMethDiff$meth.diff >= 0) #Number of unique DML hypermethylated - 87
sum(sigRes.allMFMethDiff$meth.diff <= 0) #Number of unique DML hypomethylated - 47

write.csv(sigRes.allMFMethDiff, "GOterms/BOTH-HEAT/GENELIST-MF-EnrichedGO-DMLOnly.csv", quote = FALSE, row.names = FALSE) #Save dataframe as csv

# COMBO
allDMLGOterms <- read.table("data/gene_annotations/BOTH/DML-Cov10-20-BOTH-COMBO.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", header = FALSE, sep = "\t", col.names = c("transcript", "geneID", "chr", "start", "end", "meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for all DML
head(allDMLGOterms)
nrow(allDMLGOterms) # 6941
allDMLGOtermsMF <- allDMLGOterms %>% filter(., allDMLGOterms$GOcat == "F") #Filter biological process terms
nrow(allDMLGOtermsMF) # 1262
sigRes.allMF <- allRes.COMBOMF[1:10,c(1, 2, 6)] #Filter significantly enriched GOterms, only keep GOID, GOterm, and p-value
colnames(sigRes.allMF) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.allMF

sigRes.allMFGenes <- genesInTerm(COMBOGOdataMF, whichGO = allRes.COMBOMF$GO.ID[1:10]) #Extract genes associated with each GOterm from topGO's annotations
sigRes.allMFGenes <- as.data.frame(vapply(sigRes.allMFGenes, paste, collapse = ",", character(1L))) #Collapse all transcript IDs into one column and convert into a dataframe. GO IDs will be the rownames
# head(sigRes.allMFGenes)
colnames(sigRes.allMFGenes) <- c("transcript")
sigRes.allMFGenes$GOID <- rownames(sigRes.allMFGenes) #Convert row names to a separate column
rownames(sigRes.allMFGenes) <- 1:nrow(sigRes.allMFGenes) #Convert row names to numeric values
sigRes.allMFGenes <- separate_rows(sigRes.allMFGenes, transcript, sep = ",") #Separate transcripts so each transcript is on a separate row
head(sigRes.allMFGenes) #Confirm formatting

sigRes.COMBOMFtopGO <- left_join(sigRes.allMF, sigRes.allMFGenes, by = "GOID") #Add transcript information from topGO to GOIDs
sum(is.na(sigRes.COMBOMFtopGO$transcript)) #sum = 0, there are no N/As and every GOID has a transcript
nrow(sigRes.COMBOMFtopGO) #Count the number of rows
head(sigRes.COMBOMFtopGO) #Confirm formatting

sigRes.allMFAnnot <- left_join(sigRes.COMBOMFtopGO, allDMLGOtermsMF, by = c("GOID", "transcript")) #Additional annotations. Left join (merge) based on matching transcript IDs. There are some rows with transcript IDs without corresponding annotations because there was no matching GOterm from my annotations.
colnames(sigRes.allMFAnnot)[2] <- "topGOterm" #Rename to indicate GOterms from topGO
colnames(sigRes.allMFAnnot)[10] <- "GOterm" #Rename to indicate a GOterms from my annotation
sigRes.allMFAnnot <- unique(sigRes.allMFAnnot[,-12]) #Drop GOcat column and retain only unique lines 
nrow(sigRes.allMFAnnot) #Count the number of rows
length(na.omit(unique(sigRes.allMFAnnot$geneID))) #4 unique genes
length(na.omit(unique(sigRes.allMFAnnot$transcript))) #207 unique transcripts
# tail(sigRes.allMFAnnot) #Confirm formatting 

write.csv(sigRes.allMFAnnot, "GOterms/BOTH-COMBO/all-MF-EnrichedGO-DML-withTranscript.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allMFAnnotnoT <- unique(sigRes.allMFAnnot[,-4]) #Drop transcript column and keep only unique rows
head(sigRes.allMFAnnotnoT) #Confirm formatting

write.csv(sigRes.allMFAnnotnoT, "GOterms/BOTH-COMBO/all-MF-EnrichedGO-DML.csv", quote = FALSE, row.names = FALSE) #Save dataframe

geneProductIDs <- read.delim("GFF/cgigas_uk_roslin_v1_rna_from_genomic_annot.transcript.tab", sep = "\t", header = FALSE, col.names = c("transcript", "product")) #Import transcript and gene product information
head(geneProductIDs) #Confirm import

sigRes.allMFProduct <- left_join(sigRes.allMFAnnot, geneProductIDs, by = "transcript") #Add product annotations and retain rows that may not have annotations
nrow(sigRes.allMFProduct) == nrow(sigRes.allMFAnnot) #Confirm all transcripts have an annotation
nrow(sigRes.allMFProduct)
head(sigRes.allMFProduct) #Confirm formatting

write.csv(sigRes.allMFProduct, "GOterms/BOTH-COMBO/all-MF-EnrichedGO-DML-Product.csv", quote = FALSE, row.names = FALSE) #Save dataframe

## Count the number of hyper- or hypomethylated DML in genes with enriched GOterms

sigRes.allMFMethDiff <- unique(sigRes.allMFProduct %>%
                                 dplyr::select(., transcript) %>%
                                 left_join(., COMBODMLGOterms, by = "transcript") %>%
                                 filter(., GOcat == "F") %>%
                                 dplyr::select(., geneID, chr, start, end, meth.diff)) #Isolate transcript column; join with DML list; remove all rows that aren't MF GOterms, select gene ID/chr/start/end/meth.diff columns
head(sigRes.allMFMethDiff) #Confirm formatting


nrow(sigRes.allMFMethDiff) #Number of unqiue DML in genes with enriched GOterms - 8
sum(sigRes.allMFMethDiff$meth.diff >= 0) #Number of unique DML hypermethylated - 5
sum(sigRes.allMFMethDiff$meth.diff <= 0) #Number of unique DML hypomethylated - 3

write.csv(sigRes.allMFMethDiff, "GOterms/BOTH-COMBO/GENELIST-MF-EnrichedGO-DMLOnly.csv", quote = FALSE, row.names = FALSE) #Save dataframe as csv

# RONIT ONLY - HEAT - BP

# HEAT
allDMLGOterms <- read.table("data/gene_annotations/RONIT/DML-Cov10-20-RONIT-HEAT.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", header = FALSE, sep = "\t", col.names = c("transcript", "geneID", "chr", "start", "end", "meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for all DML
head(allDMLGOterms)
nrow(allDMLGOterms) # 114226
allDMLGOtermsBP <- allDMLGOterms %>% filter(., allDMLGOterms$GOcat == "P") #Filter biological process terms
nrow(allDMLGOtermsBP) # 15213
sigRes.allBP <- allRes.HEATBP[1:10,c(1, 2, 6)] #Filter significantly enriched GOterms, only keep GOID, GOterm, and p-value
colnames(sigRes.allBP) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.allBP

sigRes.allBPGenes <- genesInTerm(HEATGOdataBP, whichGO = allRes.HEATBP$GO.ID[1:10]) #Extract genes associated with each GOterm from topGO's annotations
sigRes.allBPGenes <- as.data.frame(vapply(sigRes.allBPGenes, paste, collapse = ",", character(1L))) #Collapse all transcript IDs into one column and convert into a dataframe. GO IDs will be the rownames
# head(sigRes.allBPGenes)
colnames(sigRes.allBPGenes) <- c("transcript")
sigRes.allBPGenes$GOID <- rownames(sigRes.allBPGenes) #Convert row names to a separate column
rownames(sigRes.allBPGenes) <- 1:nrow(sigRes.allBPGenes) #Convert row names to numeric values
sigRes.allBPGenes <- separate_rows(sigRes.allBPGenes, transcript, sep = ",") #Separate transcripts so each transcript is on a separate row
head(sigRes.allBPGenes) #Confirm formatting

sigRes.HEATBPtopGO <- left_join(sigRes.allBP, sigRes.allBPGenes, by = "GOID") #Add transcript information from topGO to GOIDs
sum(is.na(sigRes.HEATBPtopGO$transcript)) #sum = 0, there are no N/As and every GOID has a transcript
nrow(sigRes.HEATBPtopGO) #Count the number of rows
head(sigRes.HEATBPtopGO) #Confirm formatting

sigRes.allBPAnnot <- left_join(sigRes.HEATBPtopGO, allDMLGOtermsBP, by = c("GOID", "transcript")) #Additional annotations. Left join (merge) based on matching transcript IDs. There are some rows with transcript IDs without corresponding annotations because there was no matching GOterm from my annotations.
colnames(sigRes.allBPAnnot)[2] <- "topGOterm" #Rename to indicate GOterms from topGO
colnames(sigRes.allBPAnnot)[10] <- "GOterm" #Rename to indicate a GOterms from my annotation
sigRes.allBPAnnot <- unique(sigRes.allBPAnnot[,-12]) #Drop GOcat column and retain only unique lines 
nrow(sigRes.allBPAnnot) # 5801 Count the number of rows
length(na.omit(unique(sigRes.allBPAnnot$geneID))) #49 unique genes
length(na.omit(unique(sigRes.allBPAnnot$transcript))) #1841 unique transcripts
tail(sigRes.allBPAnnot) #Confirm formatting 

write.csv(sigRes.allBPAnnot, "GOterms/RONIT-HEAT/all-BP-EnrichedGO-DML-withTranscript.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allBPAnnotnoT <- unique(sigRes.allBPAnnot[,-4]) #Drop transcript column and keep only unique rows
nrow(sigRes.allBPAnnotnoT)
head(sigRes.allBPAnnotnoT) #Confirm formatting

write.csv(sigRes.allBPAnnotnoT, "GOterms/RONIT-HEAT/all-BP-EnrichedGO-DML.csv", quote = FALSE, row.names = FALSE) #Save dataframe

geneProductIDs <- read.delim("GFF/cgigas_uk_roslin_v1_rna_from_genomic_annot.transcript.tab", sep = "\t", header = FALSE, col.names = c("transcript", "product")) #Import transcript and gene product information
head(geneProductIDs) #Confirm import

sigRes.allBPProduct <- left_join(sigRes.allBPAnnot, geneProductIDs, by = "transcript") #Add product annotations and retain rows that may not have annotations
nrow(sigRes.allBPProduct) == nrow(sigRes.allBPAnnot) #Confirm all transcripts have an annotation
nrow(sigRes.allBPProduct)
head(sigRes.allBPProduct) #Confirm formatting

write.csv(sigRes.allBPProduct, "GOterms/RONIT-HEAT/all-BP-EnrichedGO-DML-Product.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allBPMethDiff <- unique(sigRes.allBPProduct %>%
                                 dplyr::select(., transcript) %>%
                                 left_join(., HEATDMLGOterms, by = "transcript") %>%
                                 filter(., GOcat == "F") %>%
                                 dplyr::select(., geneID, chr, start, end, meth.diff)) #Isolate transcript column; join with DML list; remove all rows that aren't BP GOterms, select gene ID/chr/start/end/meth.diff columns
head(sigRes.allBPMethDiff) #Confirm formatting


nrow(sigRes.allBPMethDiff) #Number of unqiue DML in genes with enriched GOterms - 7
sum(sigRes.allBPMethDiff$meth.diff >= 0) #Number of unique DML hypermethylated - 2
sum(sigRes.allBPMethDiff$meth.diff <= 0) #Number of unique DML hypomethylated - 5

write.csv(sigRes.allBPMethDiff, "GOterms/RONIT-HEAT/GENELIST-BP-EnrichedGO-DMLOnly.csv", quote = FALSE, row.names = FALSE) #Save dataframe as csv

## Generate GO_similarity matrices #########################################################################################################
# PLOIDY
sigRes.PLOIDYBP <- allRes.PLOIDYBP[,c(1, 6)] #Filter significantly enriched GOterms, only keep GOID and p-value
sigRes.PLOIDYBP <- filter(sigRes.PLOIDYBP,classic<0.01)
colnames(sigRes.PLOIDYBP) <- c("GOID", "p.value") #Change column names
head(sigRes.PLOIDYBP)
matPLOIDYBP <- GO_similarity(go_id = sigRes.PLOIDYBP$GOID, ont = "BP") #Calculate the semantic similarity matrix using the Rel method (default)

sigRes.PLOIDYCC <- allRes.PLOIDYCC[,c(1, 6)] #Filter significantly enriched GOterms, only keep GOID and p-value
sigRes.PLOIDYCC <- filter(sigRes.PLOIDYCC,classic<0.01)
colnames(sigRes.PLOIDYCC) <- c("GOID", "p.value") #Change column names
head(sigRes.PLOIDYCC)
matPLOIDYCC <- GO_similarity(go_id = sigRes.PLOIDYCC$GOID, ont = "CC") #Calculate the semantic similarity matrix using the Rel method (default)

sigRes.PLOIDYMF <- allRes.PLOIDYMF[,c(1, 6)] #Filter significantly enriched GOterms, only keep GOID and p-value
sigRes.PLOIDYMF <- filter(sigRes.PLOIDYMF,classic<0.01)
colnames(sigRes.PLOIDYMF) <- c("GOID", "p.value") #Change column names
head(sigRes.PLOIDYMF)
matPLOIDYMF <- GO_similarity(go_id = sigRes.PLOIDYMF$GOID, ont = "MF") #Calculate the semantic similarity matrix using the Rel method (default)

# HEAT
sigRes.HEATBP <- allRes.HEATBP[,c(1, 6)] #Filter significantly enriched GOterms, only keep GOID and p-value
sigRes.HEATBP <- filter(sigRes.HEATBP,classic<0.01)
colnames(sigRes.HEATBP) <- c("GOID", "p.value") #Change column names
head(sigRes.HEATBP)
matHEATBP <- GO_similarity(go_id = sigRes.HEATBP$GOID, ont = "BP") #Calculate the semantic similarity matrix using the Rel method (default)

sigRes.HEATCC <- allRes.HEATCC[,c(1, 6)] #Filter significantly enriched GOterms, only keep GOID and p-value
sigRes.HEATCC <- filter(sigRes.HEATCC,classic<0.01)
colnames(sigRes.HEATCC) <- c("GOID", "p.value") #Change column names
head(sigRes.HEATCC)
matHEATCC <- GO_similarity(go_id = sigRes.HEATCC$GOID, ont = "CC") #Calculate the semantic similarity matrix using the Rel method (default)

sigRes.HEATMF <- allRes.HEATMF[,c(1, 6)] #Filter significantly enriched GOterms, only keep GOID and p-value
sigRes.HEATMF <- filter(sigRes.HEATMF,classic<0.01)
colnames(sigRes.HEATMF) <- c("GOID", "p.value") #Change column names
head(sigRes.HEATMF)
matHEATMF <- GO_similarity(go_id = sigRes.HEATMF$GOID, ont = "MF") #Calculate the semantic similarity matrix using the Rel method (default)

# PLOIDY
sigRes.COMBOBP <- allRes.COMBOBP[,c(1, 6)] #Filter significantly enriched GOterms, only keep GOID and p-value
sigRes.COMBOBP <- filter(sigRes.COMBOBP,classic<0.01)
colnames(sigRes.COMBOBP) <- c("GOID", "p.value") #Change column names
head(sigRes.COMBOBP)
matCOMBOBP <- GO_similarity(go_id = sigRes.COMBOBP$GOID, ont = "BP") #Calculate the semantic similarity matrix using the Rel method (default)

sigRes.COMBOCC <- allRes.COMBOCC[,c(1, 6)] #Filter significantly enriched GOterms, only keep GOID and p-value
sigRes.COMBOCC <- filter(sigRes.COMBOCC,classic<0.01)
colnames(sigRes.COMBOCC) <- c("GOID", "p.value") #Change column names
head(sigRes.COMBOCC)
matCOMBOCC <- GO_similarity(go_id = sigRes.COMBOCC$GOID, ont = "CC") #Calculate the semantic similarity matrix using the Rel method (default)

sigRes.COMBOMF <- allRes.COMBOMF[,c(1, 6)] #Filter significantly enriched GOterms, only keep GOID and p-value
sigRes.COMBOMF <- filter(sigRes.COMBOMF,classic<0.01)
colnames(sigRes.COMBOMF) <- c("GOID", "p.value") #Change column names
head(sigRes.COMBOMF)
matCOMBOMF <- GO_similarity(go_id = sigRes.COMBOMF$GOID, ont = "MF") #Calculate the semantic similarity matrix using the Rel method (default)

### PLOT ##################
# PLOIDY
pdf("figures/Enrichment_figure/allgigas-PLOIDY-simplifyEnrichment-BP.pdf", width = 13, height = 9) #Save figure
simplifyGO(matPLOIDYBP,
           method = "binary_cut",
           # column_title = "", 
           # col = rev(plotColors), 
           fontsize_range = c(10,25),
           word_cloud_grob_param = list(max_width = 100)) #Plot GOterms based on semantic similarity. Do not include a column title. Set colors to be plot colors, and set fontsize to range from 10 to 40. Pass arguments to word_cloud_grob_param to dictate the colors of the words and maximum width
dev.off()

dfPLOIDYBP <- simplifyGO(matPLOIDYBP, plot = FALSE) #Cluster GOterms by semantic similarity
head(dfPLOIDYBP) #Dataframe with GOIDs, GOterms, and cluster
sort(table(dfPLOIDYBP$cluster)) #Obtain size of clusters

pdf("figures/Enrichment_figure/allgigas-PLOIDY-simplifyEnrichment-CC.pdf", width = 13, height = 9) #Save figure
simplifyGO(matPLOIDYCC,
           method = "binary_cut",
           # column_title = "", 
           # col = rev(plotColors), 
           fontsize_range = c(10,25),
           word_cloud_grob_param = list(max_width = 100)) #Plot GOterms based on semantic similarity. Do not include a column title. Set colors to be plot colors, and set fontsize to range from 10 to 40. Pass arguments to word_cloud_grob_param to dictate the colors of the words and maximum width
dev.off()

dfPLOIDYCC <- simplifyGO(matPLOIDYCC, plot = FALSE) #Cluster GOterms by semantic similarity
head(dfPLOIDYCC) #Dataframe with GOIDs, GOterms, and cluster
sort(table(dfPLOIDYCC$cluster)) #Obtain size of clusters

pdf("figures/Enrichment_figure/allgigas-PLOIDY-simplifyEnrichment-MF.pdf", width = 13, height = 9) #Save figure
simplifyGO(matPLOIDYMF,
           method = "binary_cut",
           # column_title = "", 
           # col = rev(plotColors), 
           fontsize_range = c(10,25),
           word_cloud_grob_param = list(max_width = 100)) #Plot GOterms based on semantic similarity. Do not include a column title. Set colors to be plot colors, and set fontsize to range from 10 to 40. Pass arguments to word_cloud_grob_param to dictate the colors of the words and maximum width
dev.off()

dfPLOIDYMF <- simplifyGO(matPLOIDYMF, plot = FALSE) #Cluster GOterms by semantic similarity
head(dfPLOIDYMF) #Dataframe with GOIDs, GOterms, and cluster
sort(table(dfPLOIDYMF$cluster)) #Obtain size of clusters

## HEAT
pdf("figures/Enrichment_figure/allgigas-HEAT-simplifyEnrichment-BP.pdf", width = 13, height = 9) #Save figure
simplifyGO(matHEATBP,
           method = "binary_cut",
           # column_title = "", 
           # col = rev(plotColors), 
           fontsize_range = c(10,25),
           word_cloud_grob_param = list(max_width = 100)) #Plot GOterms based on semantic similarity. Do not include a column title. Set colors to be plot colors, and set fontsize to range from 10 to 40. Pass arguments to word_cloud_grob_param to dictate the colors of the words and maximum width
dev.off()

dfHEATBP <- simplifyGO(matHEATBP, plot = FALSE) #Cluster GOterms by semantic similarity
head(dfHEATBP) #Dataframe with GOIDs, GOterms, and cluster
sort(table(dfHEATBP$cluster)) #Obtain size of clusters

pdf("figures/Enrichment_figure/allgigas-HEAT-simplifyEnrichment-CC.pdf", width = 13, height = 9) #Save figure
simplifyGO(matHEATCC,
           method = "binary_cut",
           # column_title = "", 
           # col = rev(plotColors), 
           fontsize_range = c(10,25),
           word_cloud_grob_param = list(max_width = 100)) #Plot GOterms based on semantic similarity. Do not include a column title. Set colors to be plot colors, and set fontsize to range from 10 to 40. Pass arguments to word_cloud_grob_param to dictate the colors of the words and maximum width
dev.off()

dfHEATCC <- simplifyGO(matHEATCC, plot = FALSE) #Cluster GOterms by semantic similarity
head(dfHEATCC) #Dataframe with GOIDs, GOterms, and cluster
sort(table(dfHEATCC$cluster)) #Obtain size of clusters

pdf("figures/Enrichment_figure/allgigas-HEAT-simplifyEnrichment-MF.pdf", width = 13, height = 9) #Save figure
simplifyGO(matHEATMF,
           method = "binary_cut",
           # column_title = "", 
           # col = rev(plotColors), 
           fontsize_range = c(10,25),
           word_cloud_grob_param = list(max_width = 100)) #Plot GOterms based on semantic similarity. Do not include a column title. Set colors to be plot colors, and set fontsize to range from 10 to 40. Pass arguments to word_cloud_grob_param to dictate the colors of the words and maximum width
dev.off()

dfHEATMF <- simplifyGO(matHEATMF, plot = FALSE) #Cluster GOterms by semantic similarity
head(dfHEATMF) #Dataframe with GOIDs, GOterms, and cluster
sort(table(dfHEATMF$cluster)) #Obtain size of clusters

## COMBO
pdf("figures/Enrichment_figure/allgigas-INTERACTION-simplifyEnrichment-BP.pdf", width = 13, height = 9) #Save figure
simplifyGO(matCOMBOBP,
           method = "binary_cut",
           # column_title = "", 
           # col = rev(plotColors), 
           fontsize_range = c(10,25),
           word_cloud_grob_param = list(max_width = 100)) #Plot GOterms based on semantic similarity. Do not include a column title. Set colors to be plot colors, and set fontsize to range from 10 to 40. Pass arguments to word_cloud_grob_param to dictate the colors of the words and maximum width
dev.off()

dfCOMBOBP <- simplifyGO(matCOMBOBP, plot = FALSE) #Cluster GOterms by semantic similarity
head(dfCOMBOBP) #Dataframe with GOIDs, GOterms, and cluster
sort(table(dfCOMBOBP$cluster)) #Obtain size of clusters

pdf("figures/Enrichment_figure/allgigas-INTERACTION-simplifyEnrichment-CC.pdf", width = 13, height = 9) #Save figure
simplifyGO(matCOMBOCC,
           method = "binary_cut",
           # column_title = "", 
           # col = rev(plotColors), 
           fontsize_range = c(10,25),
           word_cloud_grob_param = list(max_width = 100)) #Plot GOterms based on semantic similarity. Do not include a column title. Set colors to be plot colors, and set fontsize to range from 10 to 40. Pass arguments to word_cloud_grob_param to dictate the colors of the words and maximum width
dev.off()

dfCOMBOCC <- simplifyGO(matCOMBOCC, plot = FALSE) #Cluster GOterms by semantic similarity
head(dfCOMBOCC) #Dataframe with GOIDs, GOterms, and cluster
sort(table(dfCOMBOCC$cluster)) #Obtain size of clusters

pdf("figures/Enrichment_figure/allgigas-INTERACTION-simplifyEnrichment-MF.pdf", width = 13, height = 9) #Save figure
simplifyGO(matCOMBOMF,
           method = "binary_cut",
           # column_title = "", 
           # col = rev(plotColors), 
           fontsize_range = c(10,25),
           word_cloud_grob_param = list(max_width = 100)) #Plot GOterms based on semantic similarity. Do not include a column title. Set colors to be plot colors, and set fontsize to range from 10 to 40. Pass arguments to word_cloud_grob_param to dictate the colors of the words and maximum width
dev.off()

dfCOMBOMF <- simplifyGO(matCOMBOMF, plot = FALSE) #Cluster GOterms by semantic similarity
head(dfCOMBOMF) #Dataframe with GOIDs, GOterms, and cluster
sort(table(dfCOMBOMF$cluster)) #Obtain size of clusters

### SAVE POINT - pt. 2 ####################################################################################################

# setwd("..")
# save.image("GOterms/BOTH-COMBO/2021121_Functional-enrichment-analysis_2.RData") #Save R Data in case R crashes
load("GOterms/BOTH-COMBO/2021121_Functional-enrichment-analysis_2.RData") #Save R Data in case R crashes

### Test best clustering method ###################################################################################################

pdf("GOterms/BOTH-COMBO/clustering_method_comparison_BP.pdf", width = 13, height = 13) #Save figure
compare_clustering_methods(matBP, plot_type = "heatmap")
dev.off()

pdf("GOterms/BOTH-COMBO/clustering_method_comparison_CC.pdf", width = 13, height = 13) #Save figure
compare_clustering_methods(matCC, plot_type = "heatmap")
dev.off()

pdf("GOterms/BOTH-COMBO/clustering_method_comparison_MF.pdf", width = 13, height = 13) #Save figure
compare_clustering_methods(matMF, plot_type = "heatmap")
dev.off()



