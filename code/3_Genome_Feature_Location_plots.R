# Title: Genome feature location plots
# Author: Matthew George; mattgeorgephd@gmail.com
# Date: 01/2022

###################################################################################################
## Setup work space

rm(list=ls()) # clear all

## First time? - install these
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


###########################################################################################################
### DML-genome-feature-location Stacked bar plot

data         <- read_excel("analyses/DML_locations-all-datasets.xlsx", sheet = "plot", col_names = TRUE)
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





