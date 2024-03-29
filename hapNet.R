
# A function to perform haplotype network analysis using SNP data in a specified genomic region.
# Change to the directory of ECOGEMS using the setwd function of R.
# Usage: type the next three lines in R Console without the leading #
# source("Global.R")
# snp.data <- fetchSnp(chr="chr02", start=26550915, end=26552218, mutType=NULL)
# hapnet.plot <- hapNet(data=snp.data[[1]], legend.x="bottomleft", legend.y=NULL, min.freq=20, max.freq=3000, pop.list=1, snpSites = NULL, scale.ratio=500)
# Then the result plot would be displayed in a plotting device.
# For more info, please check the Haplotype menu of the ECOGEMS database.

hapNet <- function(data=NULL, legend.x="topleft", legend.y=NULL, 
                   min.freq=10, max.freq=400, pop.list=1, snpSites = NULL, ...) {
  if (!is.null(snpSites) && length(snpSites)>=1) {
    data <- data[rownames(data) %in% snpSites, ]
  }
  dat.mat <- data
  dat.mat[is.na(dat.mat)] <- "-"
  dat.mat <- t(dat.mat)
  dat.mat.bin <- ape::as.DNAbin(dat.mat)
  
  h <- pegas::haplotype(dat.mat.bin)
  h <- pegas:::sort.haplotype(h)
  h.sub <- pegas:::subset.haplotype(h, minfreq = min.freq, maxfreq = max.freq)
  dimnames(h.sub)[[1]] <- as.character(as.roman(1:nrow(h.sub)))
  net <- pegas::haploNet(h.sub)
  
  if (pop.list == 3) {
	acc.info$Ecotype[acc.info$Ecotype %in% c("Ind_Int", "IndI", "indica", "IndII")] <- "Ind"
	acc.info$Ecotype[acc.info$Ecotype %in% c("Jap_Int", "TeJ", "TrJ")] <- "Jap"
	acc.info$Ecotype[acc.info$Ecotype %in% c("Int", "VI")] <- "Other"
	acc.info$Ecotype[acc.info$Ecotype %in% c("Or-I", "Or-II", "Or-III")] <- "Wild"
  } else if (pop.list == 2) {
	acc.info$Ecotype[acc.info$Ecotype %in% c("Ind_Int", "IndI", "indica", "IndII")] <- "Ind"
	acc.info$Ecotype[acc.info$Ecotype == "Jap_Int"] <- "Jap"
	acc.info$Ecotype[acc.info$Ecotype %in% c("Int", "VI")] <- "Other"
	acc.info$Ecotype[acc.info$Ecotype %in% c("Or-I", "Or-II", "Or-III")] <- "Wild"
  } else if (pop.list == 1) {
	acc.info$Ecotype[acc.info$Ecotype %in% c("Ind_Int", "IndI", "indica", "IndII")] <- "Ind"
	acc.info$Ecotype[acc.info$Ecotype == "Jap_Int"] <- "Jap"
	acc.info$Ecotype[acc.info$Ecotype %in% c("Int", "VI")] <- "Other"
  }
  
  net.pie <- lapply(1:length(attr(h.sub, "index")), function(x){
    x.acc <- rownames(dat.mat)[attr(h.sub, "index")[[x]]]
    x.acc.info <- acc.info$Ecotype[acc.info$ID %in% x.acc]
    if (pop.list == 3) {
	return(c(length(which(x.acc.info=="Ind")), length(which(x.acc.info=="Jap")), 
	         length(which(x.acc.info=="Aus")), length(which(x.acc.info=="Wild")),
	         length(which(x.acc.info=="Other")))/length(x.acc.info))
    } else if (pop.list == 2) {
	return(c(length(which(x.acc.info=="Ind")), length(which(x.acc.info=="Jap")), 
		 length(which(x.acc.info=="TeJ")), length(which(x.acc.info=="TrJ")),
	         length(which(x.acc.info=="Aus")), length(which(x.acc.info=="Wild")),
	         length(which(x.acc.info=="Other")))/length(x.acc.info))
    } else if (pop.list == 1) {
	return(c(length(which(x.acc.info=="Ind")), length(which(x.acc.info=="Jap")), 
		 length(which(x.acc.info=="TeJ")), length(which(x.acc.info=="TrJ")),
	         length(which(x.acc.info=="Aus")), length(which(x.acc.info=="Or-I")),
		 length(which(x.acc.info=="Or-II")), length(which(x.acc.info=="Or-III")),
	         length(which(x.acc.info=="Other")))/length(x.acc.info))
    }
  })
  
  net.pie.df <- do.call(rbind, net.pie)
  if (pop.list == 3) {
	colnames(net.pie.df) <- c("Ind", "Jap", "Aus", "wild", "Other")
  } else if (pop.list == 2) {
	colnames(net.pie.df) <- c("Ind", "Jap", "TeJ", "TrJ", "Aus", "wild", "Other")
  } else if (pop.list == 1) {
	colnames(net.pie.df) <- c("Ind", "Jap", "TeJ", "TrJ", "Aus", "Or-I", "Or-II", "Or-III", "Other")
  }
  
  pegas:::plot.haploNet(net, size=attr(net, "freq"), 
       cex = 0.8, pie=net.pie.df, legend=FALSE, 
       threshold=0, ...)

  legend(x=legend.x, y=legend.y, legend=colnames(net.pie.df), fill=rainbow(ncol(net.pie.df)), border=NA, bty="n")
}

