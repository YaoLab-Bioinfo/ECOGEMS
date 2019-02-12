
# A function to show the haplotype geographic distribution using SNP data in a specified genomic region.
# Change to the directory of ECOGEMS using the setwd function of R.
# Usage: type the next three lines in R Console without the leading #
# source("Global.R")
# snp.data <- fetchSnp(chr="chr02", start=26550915, end=26552218, mutType=NULL)
# hapGeoStatic(haplotype = hapConten(data = snp.data[[1]], min.freq=50, max.freq=2508, snpSites = NULL))
# Then the result plot would be displayed in a plotting device.
# For more info, please check the Haplotype menu of the ECOGEMS database.

hapGeoStatic <- function(haplotype=NULL) {
  hap.con <- haplotype
  dat <- acc.info
  
  dat <- dat[!is.na(dat$Latitude), ]
  dat$hap <- ""
  
  for (i in 1:length(hap.con[[2]])) {
    hap.id <- names(hap.con[[2]])[i]
    dat$hap[dat$ID %in% hap.con[[2]][[i]]] <- hap.id
  }
  
  dat <- dat[dat$hap!="", ]
  dat$text <- paste(dat[, 1], dat[,2], dat[,3], dat[,4], dat[,7], sep=", ")
  
  load("./data/worldmap.RData")
  
  mp <- mp + geom_point(aes(x=Longitude, y=Latitude, color=hap), size=0.5, data=dat) + 
    scale_x_continuous("", breaks=NULL) + scale_y_continuous("", breaks=NULL) 
  mp <- mp + guides(color=guide_legend(title=NULL))
  mp
}

