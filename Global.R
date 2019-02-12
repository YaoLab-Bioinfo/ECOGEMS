
options(warn=-1)

library(IRanges)
library(plotly)
library(LDheatmap)
library(chopsticks)
library(foreach)
library(ape)
library(pegas)
library(plyr)
library(dplyr)
library(ggmap)
library(tidyr)
library(gridExtra)
library(ggtree)
library(grid)
library(snpStats)
library(htmlwidgets)
library(shinycssloaders)
library(shinysky)

source("fetchSnp.R")
source("ld.heatmap.R")
source("phylo.R")
source("hapNet.R")
source("hapConten.R")
source("nucDiv.R")
source("writeHap.R")
source("hapGeo.R")
source("GBrowser.R")
source("anaReg.R")
source("geneStru.R")
source("hapGeoStatic.R")
source("snpInfo.R")
source("validReg.R")
source("alleleFreq.R")

acc.info <- read.table("./data/all.acc.txt", head=T, as.is=T, sep="\t", quote="")  
load("./data/gff.msu.v7.RData")
snp.lst <- read.table("./data/snp.RData.lst", head=T, as.is=T, sep="\t")
load("./data/gene.info.RData")

source("chooser.R")
#all.acc.cho <- acc.info$ID[!is.na(acc.info$Latitude)]
all.acc.cho <- paste(acc.info$ID, acc.info$Name, acc.info$Ecotype, sep=", ")
all.acc.cho <- c("Aus", "Indica", "IndicaI", "IndicaII", "Japonica", "TeJ", 
                 "TrJ", "Or-I", "Or-II", "Or-III", all.acc.cho)

chrInfo <- read.table("./data/chrInfo.txt", head=T, as.is=T, sep="\t")

acc.tree <- read.table("./data/acc.tree.txt", 
                       head=T, as.is=T, sep="\t", row.names = 1)
