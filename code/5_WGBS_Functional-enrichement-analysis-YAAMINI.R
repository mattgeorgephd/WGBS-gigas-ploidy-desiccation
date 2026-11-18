# Title: Functional enrichment analysis w/ topGO - YAAMINI's DATASET
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
# home <- getActiveDocumentContext()$path

### NOTES ####################################################################################################

# Chandra Rajan et al. (2020) Methods:
# https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fgcb.15675&file=gcb15675-sup-0002-FileS2.txt

# Bioconductor SimplifyEnrichment Methods:
# http://www.bioconductor.org/packages/devel/bioc/vignettes/simplifyEnrichment/inst/doc/simplifyEnrichment.html

### LOAD ####################################################################################################
# Before launching into enrichment, I want to create lists of GOterms associated with each DML list for descriptive purposes.

setwd('E:/WGBS-gigas-ploidy-desiccation/'); getwd()

PLOIDYDMLGOterms <- read.table("data/gene_annotations/YAAMINI/DML-Cov10-20-YAAMINI-PLOIDY.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", 
                               header = FALSE, sep = "\t", 
                               col.names = c("transcript", "geneID", "chr", "start", "end","meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for ploidy DML
head(PLOIDYDMLGOterms)
nrow(PLOIDYDMLGOterms)

####### Sumarize and filter GO Terms  ###########
# Filter by BP, CC, or MF
PLOIDYDMLGOtermsBP <- PLOIDYDMLGOterms %>% filter(., PLOIDYDMLGOterms$GOcat == "P")
PLOIDYDMLGOtermsCC <- PLOIDYDMLGOterms %>% filter(., PLOIDYDMLGOterms$GOcat == "C")
PLOIDYDMLGOtermsMF <- PLOIDYDMLGOterms %>% filter(., PLOIDYDMLGOterms$GOcat == "F")

# Count GOIDs
PLOIDYDMLGOtermCountsBP <- as.data.frame(table(PLOIDYDMLGOtermsBP$GOID))
PLOIDYDMLGOtermCountsCC <- as.data.frame(table(PLOIDYDMLGOtermsCC$GOID))
PLOIDYDMLGOtermCountsMF <- as.data.frame(table(PLOIDYDMLGOtermsMF$GOID))

nrow(PLOIDYDMLGOtermCountsBP) #2166 of unique GOIDs associated with PLOIDY-DML
nrow(PLOIDYDMLGOtermCountsCC) #431 of unique GOIDs associated with PLOIDY-DML
nrow(PLOIDYDMLGOtermCountsMF) #650 of unique GOIDs associated with PLOIDY-DML

# Add column names
colnames(PLOIDYDMLGOtermCountsBP) <- c("GOID", "Freq")
colnames(PLOIDYDMLGOtermCountsCC) <- c("GOID", "Freq")
colnames(PLOIDYDMLGOtermCountsMF) <- c("GOID", "Freq")

# Filter by frequency
PLOIDYDMLGOtermCountsBP <- PLOIDYDMLGOtermCountsBP %>% filter(PLOIDYDMLGOtermCountsBP$Freq != 0)
PLOIDYDMLGOtermCountsCC <- PLOIDYDMLGOtermCountsCC %>% filter(PLOIDYDMLGOtermCountsCC$Freq != 0)
PLOIDYDMLGOtermCountsMF <- PLOIDYDMLGOtermCountsMF %>% filter(PLOIDYDMLGOtermCountsMF$Freq != 0)

# Save GOterms associated with BP, CC, or MF
write.csv(PLOIDYDMLGOtermCountsBP, "GOterms/YAAMINI-PLOIDY/PLOIDY-BP-GOterms.csv", quote = FALSE, row.names = FALSE)
write.csv(PLOIDYDMLGOtermCountsCC, "GOterms/YAAMINI-PLOIDY/PLOIDY-CC-GOterms.csv", quote = FALSE, row.names = FALSE)
write.csv(PLOIDYDMLGOtermCountsMF, "GOterms/YAAMINI-PLOIDY/PLOIDY-MF-GOterms.csv", quote = FALSE, row.names = FALSE)

########### Load genetic background #################
# Load data for gene enrichment
geneID2GO <- topGO::readMappings(file = "data/gene_annotations/YAAMINI/geneid2go-union1x-YAAMINI.tab") #Loading the GO annotations and GeneIDs. Each line has one transcript ID and all associated GOterms
str(head(geneID2GO)) #Confirm file structure
geneNames <- names(geneID2GO) #Extract names to use as gene universe
head(geneNames)

########## Link DML to gene transcripts ##################
# I'm using the list of genes with DML as my list of interest. 
PLOIDYGenes <- PLOIDYDMLGOterms$transcript #Extract transcript ID
PLOIDYGeneList <- factor(as.integer(geneNames %in% PLOIDYGenes))
names(PLOIDYGeneList) <- geneNames #PLOIDYGenes
str(PLOIDYGeneList)

########## Run TopGo #########################
# Create topGOdata - Biological processes
PLOIDYGOdataBP <- new("topGOdata", ontology = "BP", allGenes = PLOIDYGeneList,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
PLOIDYGOdataBP #Get summary of object

# Create topGOdata - Cellular Components
PLOIDYGOdataCC <- new("topGOdata", ontology = "CC", allGenes = PLOIDYGeneList,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
PLOIDYGOdataCC #Get summary of object

# Create topGOdata - Molecular Functions
PLOIDYGOdataMF <- new("topGOdata", ontology = "MF", allGenes = PLOIDYGeneList,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
PLOIDYGOdataMF #Get summary of object

############# Save Point 1 ################
# setwd('E:/WGBS-gigas-ploidy-desiccation/data/gene_annotations/YAAMINI')
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

# Save results
setwd('E:/WGBS-gigas-ploidy-desiccation/')
write.csv(allRes.PLOIDYBP, "GOterms/YAAMINI-PLOIDY/PLOIDY-BP-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe
write.csv(allRes.PLOIDYCC, "GOterms/YAAMINI-PLOIDY/PLOIDY-CC-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe
write.csv(allRes.PLOIDYMF, "GOterms/YAAMINI-PLOIDY/PLOIDY-MF-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe

############# Match enriched GOterms with general annotation information and product ID #################
# Filter significantly enriched GOterms, only keep GOID, GOterm, and p-value
sigRes.PLOIDYBP <- allRes.PLOIDYBP[1:10,c(1, 2, 6)] 
colnames(sigRes.PLOIDYBP) <- c("GOID", "GOterm", "p.value") #Change column names

sigRes.PLOIDYCC <- allRes.PLOIDYCC[1:10,c(1, 2, 6)] 
colnames(sigRes.PLOIDYCC) <- c("GOID", "GOterm", "p.value") #Change column names

sigRes.PLOIDYMF <- allRes.PLOIDYMF[1:10,c(1, 2, 6)] 
colnames(sigRes.PLOIDYMF) <- c("GOID", "GOterm", "p.value") #Change column names

############## Extract genes associated with each GOterm from topGO's annotations
## Analysis BP ###########
# PLOIDY
allDMLGOterms <- read.table("data/gene_annotations/YAAMINI/DML-Cov10-20-YAAMINI-PLOIDY.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", header = FALSE, sep = "\t", col.names = c("transcript", "geneID", "chr", "start", "end", "meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for all DML
head(allDMLGOterms)
nrow(allDMLGOterms) # 70762
allDMLGOtermsBP <- allDMLGOterms %>% filter(., allDMLGOterms$GOcat == "P") #Filter biological process terms
nrow(allDMLGOtermsBP)
sigRes.allBP <- allRes.PLOIDYBP[1:10,c(1, 2, 6)] #Filter significantly enriched GOterms, only keep GOID, GOterm, and p-value
colnames(sigRes.allBP) <- c("GOID", "GOterm", "p.value") #Change column names
sigRes.allBP

sigRes.allBPGenes <- genesInTerm(PLOIDYGOdataBP, whichGO = allRes.PLOIDYBP$GO.ID[1:10]) #Extract genes associated with each GOterm from topGO's annotations
sigRes.allBPGenes <- as.data.frame(vapply(sigRes.allBPGenes, paste, collapse = ",", character(1L))) #Collapse all transcript IDs into one column and convert into a dataframe. GO IDs will be the rownames
# head(sigRes.allBPGenes)
colnames(sigRes.allBPGenes) <- c("transcript")
sigRes.allBPGenes$GOID <- rownames(sigRes.allBPGenes) #Convert row names to a separate column
rownames(sigRes.allBPGenes) <- 1:nrow(sigRes.allBPGenes) #Convert row names to numeric values
sigRes.allBPGenes <- separate_rows(sigRes.allBPGenes, transcript, sep = ",") #Separate transcripts so each transcript is on a separate row
head(sigRes.allBPGenes) #Confirm formatting

sigRes.PLOIDYBPtopGO <- left_join(sigRes.allBP, sigRes.allBPGenes, by = "GOID") #Add transcript information from topGO to GOIDs
sum(is.na(sigRes.PLOIDYBPtopGO$transcript)) #sum = 0, there are no N/As and every GOID has a transcript
nrow(sigRes.PLOIDYBPtopGO) #Count the number of rows
head(sigRes.PLOIDYBPtopGO) #Confirm formatting

sigRes.allBPAnnot <- left_join(sigRes.PLOIDYBPtopGO, allDMLGOtermsBP, by = c("GOID", "transcript")) #Additional annotations. Left join (merge) based on matching transcript IDs. There are some rows with transcript IDs without corresponding annotations because there was no matching GOterm from my annotations.
colnames(sigRes.allBPAnnot)[2] <- "topGOterm" #Rename to indicate GOterms from topGO
colnames(sigRes.allBPAnnot)[10] <- "GOterm" #Rename to indicate a GOterms from my annotation
sigRes.allBPAnnot <- unique(sigRes.allBPAnnot[,-12]) #Drop GOcat column and retain only unique lines 
nrow(sigRes.allBPAnnot) #Count the number of rows
length(na.omit(unique(sigRes.allBPAnnot$geneID))) #29 unique genes
length(na.omit(unique(sigRes.allBPAnnot$transcript))) #950 unique transcripts
tail(sigRes.allBPAnnot) #Confirm formatting 

write.csv(sigRes.allBPAnnot, "GOterms/YAAMINI-PLOIDY/all-BP-EnrichedGO-DML-withTranscript.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allBPAnnotnoT <- unique(sigRes.allBPAnnot[,-4]) #Drop transcript column and keep only unique rows
head(sigRes.allBPAnnotnoT) #Confirm formatting
# nrow(sigRes.allBPAnnotnoT)

write.csv(sigRes.allBPAnnotnoT, "GOterms/YAAMINI-PLOIDY/all-BP-EnrichedGO-DML.csv", quote = FALSE, row.names = FALSE) #Save dataframe

geneProductIDs <- read.delim("GFF/cgigas_uk_roslin_v1_rna_from_genomic_annot.transcript.tab", sep = "\t", header = FALSE, col.names = c("transcript", "product")) #Import transcript and gene product information
head(geneProductIDs) #Confirm import

sigRes.allBPProduct <- left_join(sigRes.allBPAnnot, geneProductIDs, by = "transcript") #Add product annotations and retain rows that may not have annotations
nrow(sigRes.allBPProduct) == nrow(sigRes.allBPAnnot) #Confirm all transcripts have an annotation
head(sigRes.allBPProduct) #Confirm formatting
nrow(sigRes.allBPProduct)

write.csv(sigRes.allBPProduct, "GOterms/YAAMINI-PLOIDY/all-BP-EnrichedGO-DML-Product.csv", quote = FALSE, row.names = FALSE) #Save dataframe

## Count the number of hyper- or hypomethylated DML in genes with enriched GOterms

sigRes.allBPMethDiff <- unique(sigRes.allBPProduct %>%
                                 dplyr::select(., transcript) %>%
                                 left_join(., PLOIDYDMLGOterms, by = "transcript") %>%
                                 filter(., GOcat == "P") %>%
                                 dplyr::select(., geneID, chr, start, end, meth.diff)) #Isolate transcript column; join with DML list; remove all rows that aren't BP GOterms, select gene ID/chr/start/end/meth.diff columns
head(sigRes.allBPMethDiff) #Confirm formatting


nrow(sigRes.allBPMethDiff) #Number of unqiue DML in genes with enriched GOterms
sum(sigRes.allBPMethDiff$meth.diff >= 0) #Number of unique DML hypermethylated
sum(sigRes.allBPMethDiff$meth.diff <= 0) #Number of unique DML hypomethylated

write.csv(sigRes.allBPMethDiff, "GOterms/YAAMINI-PLOIDY/GENELIST-BP-EnrichedGO-DMLOnly.csv", quote = FALSE, row.names = FALSE) #Save dataframe as csv

## Analysis CC ###########
# PLOIDY
allDMLGOterms <- read.table("data/gene_annotations/YAAMINI/DML-Cov10-20-YAAMINI-PLOIDY.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", header = FALSE, sep = "\t", col.names = c("transcript", "geneID", "chr", "start", "end", "meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for all DML
head(allDMLGOterms)
nrow(allDMLGOterms) # 10979
allDMLGOtermsCC <- allDMLGOterms %>% filter(., allDMLGOterms$GOcat == "C") #Filter biological process terms
nrow(allDMLGOtermsCC)
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
nrow(sigRes.PLOIDYCCtopGO) #Count the number of rows
head(sigRes.PLOIDYCCtopGO) #Confirm formatting

sigRes.allCCAnnot <- left_join(sigRes.PLOIDYCCtopGO, allDMLGOtermsCC, by = c("GOID", "transcript")) #Additional annotations. Left join (merge) based on matching transcript IDs. There are some rows with transcript IDs without corresponding annotations because there was no matching GOterm from my annotations.
colnames(sigRes.allCCAnnot)[2] <- "topGOterm" #Rename to indicate GOterms from topGO
colnames(sigRes.allCCAnnot)[10] <- "GOterm" #Rename to indicate a GOterms from my annotation
sigRes.allCCAnnot <- unique(sigRes.allCCAnnot[,-12]) #Drop GOcat column and retain only unique lines 
nrow(sigRes.allCCAnnot) #Count the number of rows
length(na.omit(unique(sigRes.allCCAnnot$geneID))) #21 unique genes
length(na.omit(unique(sigRes.allCCAnnot$transcript))) #2705 unique transcripts
# tail(sigRes.allCCAnnot) #Confirm formatting 

write.csv(sigRes.allCCAnnot, "GOterms/YAAMINI-PLOIDY/all-CC-EnrichedGO-DML-withTranscript.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allCCAnnotnoT <- unique(sigRes.allCCAnnot[,-4]) #Drop transcript column and keep only unique rows
head(sigRes.allCCAnnotnoT) #Confirm formatting
# nrow(sigRes.allCCAnnotnoT)

write.csv(sigRes.allCCAnnotnoT, "GOterms/YAAMINI-PLOIDY/all-CC-EnrichedGO-DML.csv", quote = FALSE, row.names = FALSE) #Save dataframe

geneProductIDs <- read.delim("GFF/cgigas_uk_roslin_v1_rna_from_genomic_annot.transcript.tab", sep = "\t", header = FALSE, col.names = c("transcript", "product")) #Import transcript and gene product information
head(geneProductIDs) #Confirm import

sigRes.allCCProduct <- left_join(sigRes.allCCAnnot, geneProductIDs, by = "transcript") #Add product annotations and retain rows that may not have annotations
nrow(sigRes.allCCProduct) == nrow(sigRes.allCCAnnot) #Confirm all transcripts have an annotation
head(sigRes.allCCProduct) #Confirm formatting
nrow(sigRes.allCCProduct)

write.csv(sigRes.allCCProduct, "GOterms/YAAMINI-PLOIDY/all-CC-EnrichedGO-DML-Product.csv", quote = FALSE, row.names = FALSE) #Save dataframe

## Count the number of hyper- or hypomethylated DML in genes with enriched GOterms

sigRes.allCCMethDiff <- unique(sigRes.allCCProduct %>%
                                 dplyr::select(., transcript) %>%
                                 left_join(., PLOIDYDMLGOterms, by = "transcript") %>%
                                 filter(., GOcat == "C") %>%
                                 dplyr::select(., geneID, chr, start, end, meth.diff)) #Isolate transcript column; join with DML list; remove all rows that aren't CC GOterms, select gene ID/chr/start/end/meth.diff columns
head(sigRes.allCCMethDiff) #Confirm formatting


nrow(sigRes.allCCMethDiff) #Number of unqiue DML in genes with enriched GOterms
sum(sigRes.allCCMethDiff$meth.diff >= 0) #Number of unique DML hypermethylated
sum(sigRes.allCCMethDiff$meth.diff <= 0) #Number of unique DML hypomethylated

write.csv(sigRes.allCCMethDiff, "GOterms/YAAMINI-PLOIDY/GENELIST-CC-EnrichedGO-DMLOnly.csv", quote = FALSE, row.names = FALSE) #Save dataframe as csv

## Analysis MF ###########
# PLOIDY
allDMLGOterms <- read.table("data/gene_annotations/YAAMINI/DML-Cov10-20-YAAMINI-PLOIDY.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", header = FALSE, sep = "\t", col.names = c("transcript", "geneID", "chr", "start", "end", "meth.diff","GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for all DML
head(allDMLGOterms)
nrow(allDMLGOterms) # 10979
allDMLGOtermsMF <- allDMLGOterms %>% filter(., allDMLGOterms$GOcat == "F") #Filter biological process terms
nrow(allDMLGOtermsMF)
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
nrow(sigRes.allMFAnnot) #Count the number of rows
length(na.omit(unique(sigRes.allMFAnnot$geneID))) #1 unique genes
length(na.omit(unique(sigRes.allMFAnnot$transcript))) #128 unique transcripts
tail(sigRes.allMFAnnot) #Confirm formatting 

write.csv(sigRes.allMFAnnot, "GOterms/YAAMINI-PLOIDY/all-MF-EnrichedGO-DML-withTranscript.csv", quote = FALSE, row.names = FALSE) #Save dataframe

sigRes.allMFAnnotnoT <- unique(sigRes.allMFAnnot[,-4]) #Drop transcript column and keep only unique rows
head(sigRes.allMFAnnotnoT) #Confirm formatting
# nrow(sigRes.allMFAnnotnoT)

write.csv(sigRes.allMFAnnotnoT, "GOterms/YAAMINI-PLOIDY/all-MF-EnrichedGO-DML.csv", quote = FALSE, row.names = FALSE) #Save dataframe

geneProductIDs <- read.delim("GFF/cgigas_uk_roslin_v1_rna_from_genomic_annot.transcript.tab", sep = "\t", header = FALSE, col.names = c("transcript", "product")) #Import transcript and gene product information
head(geneProductIDs) #Confirm import

sigRes.allMFProduct <- left_join(sigRes.allMFAnnot, geneProductIDs, by = "transcript") #Add product annotations and retain rows that may not have annotations
nrow(sigRes.allMFProduct) == nrow(sigRes.allMFAnnot) #Confirm all transcripts have an annotation
head(sigRes.allMFProduct) #Confirm formatting
nrow(sigRes.allMFProduct)

write.csv(sigRes.allMFProduct, "GOterms/YAAMINI-PLOIDY/all-MF-EnrichedGO-DML-Product.csv", quote = FALSE, row.names = FALSE) #Save dataframe

## Count the number of hyper- or hypomethylated DML in genes with enriched GOterms

sigRes.allMFMethDiff <- unique(sigRes.allMFProduct %>%
                                 dplyr::select(., transcript) %>%
                                 left_join(., PLOIDYDMLGOterms, by = "transcript") %>%
                                 filter(., GOcat == "F") %>%
                                 dplyr::select(., geneID, chr, start, end, meth.diff)) #Isolate transcript column; join with DML list; remove all rows that aren't MF GOterms, select gene ID/chr/start/end/meth.diff columns
head(sigRes.allMFMethDiff) #Confirm formatting


nrow(sigRes.allMFMethDiff) #Number of unqiue DML in genes with enriched GOterms
sum(sigRes.allMFMethDiff$meth.diff >= 0) #Number of unique DML hypermethylated
sum(sigRes.allMFMethDiff$meth.diff <= 0) #Number of unique DML hypomethylated

write.csv(sigRes.allMFMethDiff, "GOterms/YAAMINI-PLOIDY/GENELIST-MF-EnrichedGO-DMLOnly.csv", quote = FALSE, row.names = FALSE) #Save dataframe as csv

### Summarize GOterms #########################################################################################################

allDMLGOtermCounts <- as.data.frame(table(PLOIDYDMLGOterms$GOID)) #Create table with GOID and counts
colnames(allDMLGOtermCounts) <- c("GOID", "Freq") #Add column names
range(allDMLGOtermCounts$Freq) #Get range of GOID frequency
nrow(allDMLGOtermCounts) #Number of unique GOIDs associated with all-DML

allDMLGOtermsBP <- PLOIDYDMLGOterms %>% filter(., PLOIDYDMLGOterms$GOcat == "P") #Filter biological process terms
allDMLGOtermCountsBP <- as.data.frame(table(allDMLGOtermsBP$GOID)) #Count GOIDs
colnames(allDMLGOtermCountsBP) <- c("GOID", "Freq") #Add column names
allDMLGOtermCountsBP <- allDMLGOtermCountsBP %>% filter(allDMLGOtermCountsBP$Freq != 0)
range(allDMLGOtermCountsBP$Freq) #Get range of GOID frequency
nrow(allDMLGOtermCountsBP) #Number of unique GOIDs associated with all-DML

# write.csv(allDMLGOtermCountsBP, "GOterms/YAAMINI-PLOIDY/PLOIDY-BP-GOterms.csv", quote = FALSE, row.names = FALSE)

write.csv(allRes.PLOIDYBP, "GOterms/YAAMINI-PLOIDY/PLOIDY-BP-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe
write.csv(allRes.PLOIDYCC, "GOterms/YAAMINI-PLOIDY/PLOIDY-CC-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe
write.csv(allRes.PLOIDYMF, "GOterms/YAAMINI-PLOIDY/PLOIDY-MF-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe

#### EXTRA ########################################################################################################




# # Create all geneList
# allGenes <- PLOIDYDMLGOterms$transcript #Extract transcript ID from dataframe that does not include SNP-DML overlaps
# allGeneList <- factor(as.integer(geneNames %in% PLOIDYDMLGOterms))
# names(allGeneList) <- allGenes
# str(allGeneList)
# head(allGeneList)
# 
# # BP geneList
# GenesBP <- PLOIDYDMLGOterms$transcript #Extract transcript ID from dataframe that does not include SNP-DML overlaps
# allGeneList <- factor(as.integer(geneNames %in% PLOIDYDMLGOterms))
# names(allGeneList) <- allGenes
# str(allGeneList)
# head(allGeneList)

# Create topGo object
# allGOdataBP <- new("topGOdata", ontology = "BP", allGenes = allGeneList,
#                    annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
# allGOdataBP #Get summary of object


### Generate GO_similarity matrixes #########################################################################################################
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


### PLOT ##################
# PLOIDY
pdf("GOterms/YAAMINI-PLOIDY/simplifyEnrichment-BP.pdf", width = 13, height = 9) #Save figure
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

pdf("GOterms/YAAMINI-PLOIDY/simplifyEnrichment-CC.pdf", width = 13, height = 9) #Save figure
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

pdf("GOterms/YAAMINI-PLOIDY/simplifyEnrichment-MF.pdf", width = 13, height = 9) #Save figure
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

### SAVE POINT - pt. 2 ####################################################################################################

# setwd("..")
# save.image("GOterms/YAAMINI-COMBO/2021121_Functional-enrichment-analysis_2.RData") #Save R Data in case R crashes
load("GOterms/YAAMINI-COMBO/2021121_Functional-enrichment-analysis_2.RData") #Save R Data in case R crashes

### Test best clustering method ###################################################################################################

pdf("GOterms/YAAMINI-PLOIDY/clustering_method_comparison_BP.pdf", width = 13, height = 13) #Save figure
compare_clustering_methods(matPLOIDYBP, plot_type = "heatmap")
dev.off()

pdf("GOterms/YAAMINI-PLOIDY/clustering_method_comparison_CC.pdf", width = 13, height = 13) #Save figure
compare_clustering_methods(matPLOIDYCC, plot_type = "heatmap")
dev.off()

pdf("GOterms/YAAMINI-PLOIDY/clustering_method_comparison_MF.pdf", width = 13, height = 13) #Save figure
compare_clustering_methods(matPLOIDYMF, plot_type = "heatmap")
dev.off()






