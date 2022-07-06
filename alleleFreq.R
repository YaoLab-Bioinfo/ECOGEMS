
# A function to perform allele frequency analysis using SNP data in a specified genomic region.
# Change to the directory of ECOGEMS using the setwd function of R.
# Usage: type the next two lines in R Console without the leading #
# source("Global.R")
# allele.plot <- alleleFreq(snpSite = c("0602942293", "0138383182", "0329584501", "0316733111"), accGroup = c("Aus", "Indica", "TeJ", "TrJ", "Wild"), pieCols = c("cornflowerblue", "forestgreen"))
# Then the result plot would be displayed in a plotting device.
# For more info, please check the AlleleFreq menu of the ECOGEMS database.

alleleFreq <- function(snpSite = c("0602942293", "0138383182", "0329584501", "0316733111"),
                       accGroup = c("Aus", "Indica", "TeJ", "TrJ", "Wild"),
                       pieCols = c("cornflowerblue", "forestgreen") ) {
  myChr <- paste0("chr", substr(snpSite, 1, 2))
  myPos <- as.numeric(substr(snpSite, 3, 10))
  
  site.geno <- lapply(1:length(myChr), function(i) {
    return(fetchSnp(myChr[i], myPos[i], myPos[i]))
  })
  
  acc.group <- accGroup
  
  site.allele.freq <- lapply(site.geno, function(x){
    if (nrow(x[[2]]) == 0) {
      return(NULL)
    } else {
      maj.allele <- x[[2]][1, 1]
      min.allele <- x[[2]][1, 2]
      
      acc.grp.tab <- lapply(acc.group, function(i) {
        i.acc <- readLines(paste0("./data/", i, ".acc.txt"))
        x.i.allele.freq <- table(x[[1]][, colnames(x[[1]]) %in% i.acc])
        x.i.allele.freq <- x.i.allele.freq[c(maj.allele, min.allele)]
        if (is.na(x.i.allele.freq[2])) {
          x.i.allele.freq[2] <- 0
          names(x.i.allele.freq)[2] <- min.allele
        }
        
        return(x.i.allele.freq)
      })
      
      acc.grp.tab.df <- do.call(rbind, acc.grp.tab)
      rownames(acc.grp.tab.df) <- acc.group
      return(acc.grp.tab.df)
    }
  })
  
  snpSite <- snpSite[!sapply(site.allele.freq, is.null)]
  site.allele.freq[sapply(site.allele.freq, is.null)] <- NULL
  
  op <- par(mfrow=c(length(site.allele.freq), length(acc.group) + 1 ),
            oma = c(0, 0, 0, 0),
            mar = c(0, 1, 1, 0),
            mgp = c(0, 0, 0),
            xpd = NA)
  
  pie.cols <- pieCols
  
  for (i in 1:length(site.allele.freq)) {
    for (j in 1:nrow(site.allele.freq[[i]])) {
      if (i == 1) {
        if (j == 1) {
          pie(site.allele.freq[[i]][j, ], main=acc.group[j], col=pie.cols,
              ylab = snpSite[i], label=NA, radius=1)
        } else if (j == length(acc.group)) {
          pie(site.allele.freq[[i]][j, ], main=acc.group[j], col=pie.cols, label=NA, radius=1)
          plot(0, type = "n", axes = F, xlab="", ylab="")
          legend("center", legend = names(site.allele.freq[[i]][j, ]), 
                 fill = pie.cols, border = pie.cols, bty="n", cex = 2)
        } else {
          pie(site.allele.freq[[i]][j, ], main=acc.group[j], col=pie.cols, label=NA, radius=1)
        }
      } else {
        if (j == 1) {
          pie(site.allele.freq[[i]][j, ], col=pie.cols,
              ylab = snpSite[i], label=NA, radius=1)
        } else if (j == length(acc.group)) {
          pie(site.allele.freq[[i]][j, ], col=pie.cols, label=NA, radius=1)
          plot(0, type = "n", axes = F, xlab="", ylab="")
          legend("center", legend = names(site.allele.freq[[i]][j, ]), 
                 fill = pie.cols, border = pie.cols, bty="n", cex = 2)
        } else {
          pie(site.allele.freq[[i]][j, ], col=pie.cols, label=NA, radius=1)
        }
      }
    }
  }
  
  par(op)
  
  return(site.allele.freq)
}

