

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

