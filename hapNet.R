

hapNet <- function(data=NULL, legend.x="topleft", legend.y=NULL, 
                   min.freq=10, max.freq=400, pop.list=1, snpSites = NULL, ...) {
  if (!is.null(snpSites) && length(snpSites)>=1) {
    data <- data[rownames(data) %in% snpSites, ]
  }
  dat.mat <- data
  dat.mat[is.na(dat.mat)] <- "-"
  dat.mat <- t(dat.mat)
  dat.mat.bin <- as.DNAbin(dat.mat)
  
  h <- haplotype(dat.mat.bin)
  h <- sort(h)
  h.sub <- subset(h, minfreq = min.freq, maxfreq = max.freq)
  dimnames(h.sub)[[1]] <- as.character(as.roman(1:nrow(h.sub)))
  net <- haploNet(h.sub)
  
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
  
  plot(net, size=attr(net, "freq"), 
       cex = 0.8, pie=net.pie.df, legend=FALSE, 
       threshold=0, ...)

  legend(x=legend.x, y=legend.y, legend=colnames(net.pie.df), fill=rainbow(ncol(net.pie.df)), border=NA, bty="n")
}

