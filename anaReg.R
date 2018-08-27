

anaReg <- function(x=NULL) {
  x <- gsub("\\s+", "", x)
  if (grepl("chr", x)) {
    myChr <- gsub(":.+", "", x)
    myPos <- as.numeric(gsub("\\s","", strsplit(gsub(".+:", "", x),"-")[[1]]))
  } else {
    myChr <- paste0("chr", substr(x, 7,8))
    myPos <- c(gene.info$start[gene.info$id==x], gene.info$end[gene.info$id==x])
  }
  
  chr.size <- c(43268879, 35930381, 36406689, 35278225, 29894789, 31246789,
                29696629, 28439308, 23011239, 23134759, 28512666, 27497214)
  names(chr.size) <- paste0("chr", sprintf("%02d", 1:12))
  
  myPos[1] <- max(1, myPos[1])
  myPos[2] <- min(myPos[2], chr.size[myChr])
  
  return(list(chr=myChr, start=myPos[1], end=myPos[2]))
}

