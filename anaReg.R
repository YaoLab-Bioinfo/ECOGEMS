

anaReg <- function(x=NULL) {
  x <- gsub("\\s+", "", x)
  if (grepl("chr", x)) {
    myChr <- gsub(":.+", "", x)
    myPos <- as.numeric(gsub("\\s","", strsplit(gsub(".+:", "", x),"-")[[1]]))
  } else {
    myChr <- paste0("chr", substr(x, 7,8))
    myPos <- c(gene.info$start[gene.info$id==x], gene.info$end[gene.info$id==x])
  }
  
  return(list(chr=myChr, start=myPos[1], end=myPos[2]))
}

