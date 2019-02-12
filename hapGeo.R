
# A function to show the haplotype geographic distribution using SNP data in a specified genomic region.
# Change to the directory of ECOGEMS using the setwd function of R.
# Usage: type the next three lines in R Console without the leading #
# source("Global.R")
# snp.data <- fetchSnp(chr="chr02", start=26550915, end=26552218, mutType=NULL)
# hapGeo(haplotype = hapConten(data = snp.data[[1]], min.freq=50, max.freq=2508, snpSites = NULL))
# Then the result plot would be displayed in a plotting device.
# For more info, please check the Haplotype menu of the ECOGEMS database.

hapGeo <- function(haplotype=NULL) {
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
  
  g <- list(
    scope = 'world',
    projection = list(type = 'Equirectangular'),
    showland = TRUE,
    showocean=TRUE,
    showcountries = TRUE,
    showsubunits = TRUE,
    landcolor = "white",
    oceancolor = toRGB("gray90"),
    subunitwidth = 1,
    countrywidth = 0.5,
    subunitcolor = "blue",
    countrycolor = "gray85"
  )
  
  plot_geo(dat, lat = ~Latitude, lon = ~Longitude, color = ~hap) %>%
    add_markers(
      hovertext = ~text
    ) %>% layout(title = 'Geographic distribution of rice accessions with different haplotypes', geo = g)
  
}

