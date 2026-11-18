# Title: Functional enrichment analysis w/ topGO
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
# BiocManager::install("rrvgo")

options("install.lock"=FALSE)

## Load R packages
library(rstudioapi)
library(readxl)
library(tidyverse)
library(BiocManager)
library(topGO)
library(simplifyEnrichment)
library(magick)
library(gclus)
library(rrvgo)

## Grab the WD from the file location and set it
home <- getActiveDocumentContext()$path
setwd(dirname(home)); getwd()

#######################################################################################################

# Chandra Rajan et al. (2020) Methods:
# https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fgcb.15675&file=gcb15675-sup-0002-FileS2.txt

# Bioconductor SimplifyEnrichment Methods:
# http://www.bioconductor.org/packages/devel/bioc/vignettes/simplifyEnrichment/inst/doc/simplifyEnrichment.html

#######################################################################################################
# Before launching into enrichment, I want to create lists of GOterms associated with each DML list for descriptive purposes.

setwd(dirname(home)); setwd('..'); getwd()

ploidyDMLGOterms <- read.table("data/gene_annotations/DML-ploidy-Cov10-20.GeneIDs.geneOverlap.transcriptIDs.GOAnnot", 
                               header = FALSE, sep = "\t", 
                               col.names = c("transcript", "geneID", "chr", "start", "end", "GOID", "GOterm", "GOSlim", "GOcat")) #Import list of GO term annotations for ploidy DML
head(ploidyDMLGOterms)

######################################################################################################

# Rearrange table
ploidyDMLGOtermCounts <- as.data.frame(table(ploidyDMLGOterms$GOID)) #Create table with GOID and counts
colnames(ploidyDMLGOtermCounts) <- c("GOID", "Freq") #Add column names
range(ploidyDMLGOtermCounts$Freq) #Get range of GOID frequency
nrow(ploidyDMLGOtermCounts) #Number of unique GOIDs associated with ploidy-DML


# Biological processes
ploidyDMLGOtermsBP <- ploidyDMLGOterms %>% filter(., ploidyDMLGOterms$GOcat == "P") #Filter biological process terms
ploidyDMLGOtermCountsBP <- as.data.frame(table(ploidyDMLGOtermsBP$GOID)) #Count GOIDs
colnames(ploidyDMLGOtermCountsBP) <- c("GOID", "Freq") #Add column names
ploidyDMLGOtermCountsBP <- ploidyDMLGOtermCountsBP %>% filter(ploidyDMLGOtermCountsBP$Freq != 0)
range(ploidyDMLGOtermCountsBP$Freq) #Get range of GOID frequency
nrow(ploidyDMLGOtermCountsBP) #Number of unique GOIDs associated with ploidy-DML

write.csv(ploidyDMLGOtermCountsBP, "data/ploidy-BP-GOterms.csv", quote = FALSE, row.names = FALSE)

# Load data for gene enrichment
geneID2GO <- topGO::readMappings(file = "data/geneid2go-union1x.tab") #Loading the GO annotations and GeneIDs. Each line has one transcript ID and all associated GOterms
str(head(geneID2GO)) #Confirm file structure
geneNames <- names(geneID2GO) #Extract names to use as gene universe
head(geneNames)

# I'm using the list of genes with DML as my list of interest. 
ploidyGenes <- ploidyDMLGOterms$transcript #Extract transcript ID
ploidyGeneList <- factor(as.integer(geneNames %in% ploidyGenes))
names(ploidyGeneList) <- geneNames #ploidyGenes
str(ploidyGeneList)


### Biological processes
ploidyGOdataBP <- new("topGOdata", ontology = "BP", allGenes = ploidyGeneList,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
ploidyGOdataBP #Get summary of object

test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher.ploidyBP <- getSigGroups(ploidyGOdataBP, test.stat)
resultFisher.ploidyBP

pvalFis.ploidyBP <- score(resultFisher.ploidyBP) #Extract p-values
head(pvalFis.ploidyBP)
hist(pvalFis.ploidyBP, 50, xlab = "p-values") #Plot histogram of p-values

allRes.ploidyBP <- GenTable(ploidyGOdataBP, classic = resultFisher.ploidyBP, ranksOf = "classic", orderBy = "classic", topNodes = length(pvalFis.ploidyBP)) #Create a statistical results table with statistical test results. Order by p-value (classic), and include all results (topNodes)
head(allRes.ploidyBP)

write.csv(allRes.ploidyBP, "tables/ploidy-BP-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe

### Cellular components

ploidyGOdataCC <- new("topGOdata", ontology = "CC", allGenes = ploidyGeneList,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
ploidyGOdataCC #Get summary of object

resultFisher.ploidyCC <- getSigGroups(ploidyGOdataCC, test.stat)
resultFisher.ploidyCC

pvalFis.ploidyCC <- score(resultFisher.ploidyCC) #Extract p-values
head(pvalFis.ploidyCC)
hist(pvalFis.ploidyCC, 50, xlab = "p-values") #Plot histogram of p-values

allRes.ploidyCC <- GenTable(ploidyGOdataCC, classic = resultFisher.ploidyCC, ranksOf = "classic", orderBy = "classic", topNodes = length(pvalFis.ploidyCC)) #Create a statistical results table with statistical test results. Order by p-value (classic), and include all results (topNodes)
head(allRes.ploidyCC)

write.csv(allRes.ploidyCC, "tables/ploidy-CC-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe


### Molecular function
ploidyGOdataMF <- new("topGOdata", ontology = "MF", allGenes = ploidyGeneList,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
ploidyGOdataMF #Get summary of object

resultFisher.ploidyMF <- getSigGroups(ploidyGOdataMF, test.stat) #Extract significant GOterms using Fisher exact test
resultFisher.ploidyMF #Get results summary

pvalFis.ploidyMF <- score(resultFisher.ploidyMF) #Extract p-values
head(pvalFis.ploidyMF)
hist(pvalFis.ploidyMF, 50, xlab = "p-values") #Plot histogram of p-values

allRes.ploidyMF <- GenTable(ploidyGOdataMF, classic = resultFisher.ploidyMF, ranksOf = "classic", orderBy = "classic", topNodes = length(pvalFis.ploidyMF)) #Create a statistical results table with statistical test results. Order by p-value (classic), and include all results (topNodes)
head(allRes.ploidyMF)

write.csv(allRes.ploidyMF, "tables/ploidy-MF-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe


#######################################################################################################
## SAVE POINT - pt. 1
setwd("..")
# save.image("R_sessions/20210626_Functional-enrichment-analysis_1.RData") #Save R Data in case R crashes
# load("R_sessions/20210626_Functional-enrichment-analysis_1.RData") #Save R Data in case R crashes

############################################################################################################
## Remove SNPs? - unclear if this is important/how it interacts with Bismark outputs

# uniqueSNPAll <- read.table("../10_DML-characterization/DML-pH-50-Cov5-All-unique-CT-SNPs.bed", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "meth.diff")) #Import overlaps between unique C/T SNPs and all DML produced with bedtools
# uniqueSNPAll$chrStartEnd <- paste(uniqueSNPAll$chr, "-", uniqueSNPAll$start, "-", uniqueSNPAll$end, sep = "") #Create a column with chr, start, and end information separated by "-"
# head(uniqueSNPAll) #Confirm import

# Filter by unique SNPs
# allDMLGOterms <- left_join(x = allDMLGOterms, y = uniqueSNPAll, by = "chrStartEnd") #Add all columns from uniqueSNPAll to allDMLGOterms
# allDMLGOtermsFiltered <- filter(allDMLGOterms, is.na(chr.y) == "TRUE")  %>%
#   select(!c(chrStartEnd, chr.y, start.y, end.y, meth.diff.y)) #Keep all rows where chr.y has an "NA", since this indicates there was no matching entry in the DML-SNP overlap dataframe. Then drop chrStartEnd, chr.y, start.y, end.y, and meth.diff.y
# length(allDMLGOtermsFiltered[,1]) #17810 rows without SNP overlaps
# colnames(allDMLGOtermsFiltered) <- c("transcript", "geneID", "chr", "start", "end", "meth.diff","GOID", "GOterm", "GOSlim", "GOcat") #Rename columns for consistency
# head(allDMLGOtermsFiltered) #Confirm formatting

############################################################################################################
## Summarize GOterms

allDMLGOtermCounts <- as.data.frame(table(ploidyDMLGOterms$GOID)) #Create table with GOID and counts
colnames(allDMLGOtermCounts) <- c("GOID", "Freq") #Add column names
range(allDMLGOtermCounts$Freq) #Get range of GOID frequency
nrow(allDMLGOtermCounts) #Number of unique GOIDs associated with all-DML

allDMLGOtermsBP <- ploidyDMLGOterms %>% filter(., ploidyDMLGOterms$GOcat == "P") #Filter biological process terms
allDMLGOtermCountsBP <- as.data.frame(table(allDMLGOtermsBP$GOID)) #Count GOIDs
colnames(allDMLGOtermCountsBP) <- c("GOID", "Freq") #Add column names
allDMLGOtermCountsBP <- allDMLGOtermCountsBP %>% filter(allDMLGOtermCountsBP$Freq != 0)
range(allDMLGOtermCountsBP$Freq) #Get range of GOID frequency
nrow(allDMLGOtermCountsBP) #Number of unique GOIDs associated with all-DML

write.csv(allDMLGOtermCountsBP, "Ploidy-BP-GOterms.csv", quote = FALSE, row.names = FALSE)

# Create all geneList
allGenes <- ploidyDMLGOterms$transcript #Extract transcript ID from dataframe that does not include SNP-DML overlaps
allGeneList <- factor(as.integer(geneNames %in% ploidyDMLGOterms))
names(allGeneList) <- allGenes
str(allGeneList)

# Create topGo object
# allGOdataBP <- new("topGOdata", ontology = "BP", allGenes = allGeneList,
#                    annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
# allGOdataBP #Get summary of object

# write.csv(allRes.ploidyBP, "tables/ploidy-BP-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe
# write.csv(allRes.ploidyCC, "tables/ploidy-CC-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe
# write.csv(allRes.ploidyMF, "tables/ploidy-MF-FisherTestResults.csv", quote = FALSE, row.names = FALSE) #Save dataframe

############################################################################################################
## Generate GO_similarity matrixes

# Biological Processes
sigRes.allBP <- allRes.ploidyBP[,c(1, 6)] #Filter significantly enriched GOterms, only keep GOID and p-value
sigRes.allBP <- filter(sigRes.allBP,classic<0.01)
colnames(sigRes.allBP) <- c("GOID", "p.value") #Change column names
head(sigRes.allBP)
matBP <- GO_similarity(go_id = sigRes.allBP$GOID, ont = "BP") #Calculate the semantic similarity matrix using the Rel method (default)

# Cellular Components
sigRes.allCC <- allRes.ploidyCC[,c(1, 6)] #Filter significantly enriched GOterms, only keep GOID and p-value
sigRes.allCC <- filter(sigRes.allCC,classic<0.01)
colnames(sigRes.allCC) <- c("GOID", "p.value") #Change column names
head(sigRes.allCC)
matCC <- GO_similarity(go_id = sigRes.allCC$GOID, ont = "CC") #Calculate the semantic similarity matrix using the Rel method (default)

# Molecular Functions
sigRes.allMF <- allRes.ploidyMF[,c(1, 6)] #Filter significantly enriched GOterms, only keep GOID and p-value
sigRes.allMF <- filter(sigRes.allMF,classic<0.01)
colnames(sigRes.allMF) <- c("GOID", "p.value") #Change column names
head(sigRes.allMF)
matMF <- GO_similarity(go_id = sigRes.allMF$GOID, ont = "MF") #Calculate the semantic similarity matrix using the Rel method (default)



#######################################################################################################
## SAVE POINT - pt. 2
setwd("..")
# save.image("R_sessions/20210626_Functional-enrichment-analysis_2.RData") #Save R Data in case R crashes
load("R_sessions/20210626_Functional-enrichment-analysis_2.RData") #Save R Data in case R crashes

######################################################################################################
## Test best clustering method

pdf("plots/clustering_method_comparison_BP.pdf", width = 13, height = 13) #Save figure
compare_clustering_methods(matBP, plot_type = "heatmap")
dev.off()

pdf("plots/clustering_method_comparison_CC.pdf", width = 13, height = 13) #Save figure
compare_clustering_methods(matCC, plot_type = "heatmap")
dev.off()

pdf("plots/clustering_method_comparison_MF.pdf", width = 13, height = 13) #Save figure
compare_clustering_methods(matMF, plot_type = "heatmap")
dev.off()

######################################################################################################
## Figures

## EXAMPLE/TEST
# set.seed(888)
# go_id = random_GO(500)
# mat = GO_similarity(go_id)
# GO_similarity(go_id, measure = "Wang")
# df = simplifyGO(mat)
# head(df) #Dataframe with GOIDs, GOterms, and cluster
# sort(table(df$cluster)) #Obtain size of clusters

# Biological Processes
pdf("plots/simplifyEnrichment-BP.pdf", width = 13, height = 9) #Save figure
simplifyGO(matBP,
           method = "binary_cut",
           # column_title = "", 
           # col = rev(plotColors), 
           fontsize_range = c(10,25),
           word_cloud_grob_param = list(max_width = 100)) #Plot GOterms based on semantic similarity. Do not include a column title. Set colors to be plot colors, and set fontsize to range from 10 to 40. Pass arguments to word_cloud_grob_param to dictate the colors of the words and maximum width
dev.off()

dfBP <- simplifyGO(matBP, plot = FALSE) #Cluster GOterms by semantic similarity
head(dfBP) #Dataframe with GOIDs, GOterms, and cluster
sort(table(dfBP$cluster)) #Obtain size of clusters

# Cellular Components
pdf("plots/simplifyEnrichment-CC.pdf", width = 13, height = 9) #Save figure
simplifyGO(matCC,
           method = "binary_cut",
           # column_title = "", 
           # col = rev(plotColors), 
           fontsize_range = c(10,25),
           word_cloud_grob_param = list(max_width = 100)) #Plot GOterms based on semantic similarity. Do not include a column title. Set colors to be plot colors, and set fontsize to range from 10 to 40. Pass arguments to word_cloud_grob_param to dictate the colors of the words and maximum width
dev.off()

dfCC <- simplifyGO(matCC, plot = FALSE) #Cluster GOterms by semantic similarity
head(dfCC) #Dataframe with GOIDs, GOterms, and cluster
sort(table(dfCC$cluster)) #Obtain size of clusters

# Molecular Functions
pdf("plots/simplifyEnrichment-MF.pdf", width = 13, height = 9) #Save figure
simplifyGO(matMF,
           method = "binary_cut",
           # column_title = "", 
           # col = rev(plotColors), 
           fontsize_range = c(10,25),
           word_cloud_grob_param = list(max_width = 100)) #Plot GOterms based on semantic similarity. Do not include a column title. Set colors to be plot colors, and set fontsize to range from 10 to 40. Pass arguments to word_cloud_grob_param to dictate the colors of the words and maximum width
dev.off()

dfMF <- simplifyGO(matMF, plot = FALSE) #Cluster GOterms by semantic similarity
head(dfMF) #Dataframe with GOIDs, GOterms, and cluster
sort(table(dfMF$cluster)) #Obtain size of clusters

######################################################################################################
## Simplot - need to make oyster genome

BiocManager::install("rtracklayer")
library(rtracklayer)

txdbOYSTER <- makeTxDbFromUCSC(genome="mm9", tablename="knownGene",
                              transcript_ids=geneID2GO)

## Using that, we can call our function to promote it to an OrgDb object:
odb <- makeOrganismDbFromTxDb(txdb=txdbMouse)

columns(odb)



simMatrix <- calculateSimMatrix(sigRes.allBP$GOID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(sigRes.allBP$p.value, sigRes.allBP$GOID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)

scatterPlot(simMatrix, reducedTerms)

treemapPlot(reducedTerms)

wordcloudPlot(reducedTerms, min.freq=1, colors="black")

# Install Shiny before use
# rrvgo::shiny_rrvgo()

######################################################################################################
## KEGG pathway analysis

# # BiocManager::install(version = "3.12") # Needs R version 4.0
# BiocManager::install("KEGG.db", force=TRUE)
# BiocManager::install("clusterProfiler", force=TRUE, dependencies = TRUE, INSTALL_opts = '--no-lock')
# BiocManager::install("KEGGREST")
# 
# library(KEGGREST)
# library(clusterProfiler)
# library(KEGG.db)
# 
# term_similarity_from_KEGG(matBP)


write.csv(sigRes.allBP,file="output.csv")



