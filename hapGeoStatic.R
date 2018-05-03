

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

