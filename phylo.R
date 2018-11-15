

phylo <- function(chr="chr09", start=37800, end=46400, accession=NULL, mutType=NULL, snpSites = NULL) {
  start <- as.numeric(start)
  end <- as.numeric(end)
  reg.gr <- IRanges(start, end)
  snp.lst.chr <- snp.lst[snp.lst$chr==chr, ]
  snp.lst.gr <- IRanges(start=snp.lst.chr$start, end=snp.lst.chr$end)
  snp.fls <- snp.lst.chr$file[unique(queryHits(findOverlaps(snp.lst.gr, reg.gr)))]
  
  snp.data.lst <- lapply(snp.fls, function(x){
    load(x)
    return(snp.data.inter.Matrix)
  })
  snp.data <- do.call(rbind, snp.data.lst)
  snp.data <- snp.data[order(as.numeric(rownames(snp.data))), ]
  
  start <- as.numeric(paste0(substr(chr, 4, 5), sprintf("%08d", start)))
  end <- as.numeric(paste0(substr(chr, 4, 5), sprintf("%08d", end)))
  
  dat.res <- snp.data[as.numeric(rownames(snp.data))>=start & as.numeric(rownames(snp.data))<=end, , drop=FALSE]
  dat.res <- as.matrix(dat.res)
  
  accession <- gsub(",.+", "", accession)
  accession <- sapply(accession, function(x){
    if (x %in% c("Aus", "Indica", "IndicaI", "IndicaII", "Japonica", "TeJ", "TrJ", "Or-I", "Or-II", "Or-III")) {
      x.dat <- readLines(paste0("./data/", x, ".acc.txt"))
      return(x.dat)
    } else {
      return(x)
    }
  })
  accession <- unique(unlist(accession))
  
  if (!is.null(accession) && length(accession)>=2) {
    dat.res <- dat.res[, colnames(dat.res) %in% accession, drop=FALSE]
  }
  
  dat.res.row.c <- apply(dat.res, 1, function(x){
    length(unique(x[!is.na(x)]))
  })
  dat.res <- dat.res[dat.res.row.c>1, , drop=FALSE]
  
  if (!is.null(mutType) && length(mutType)>=1 && length(mutType)!=16) {
    eff.Rdata <- paste0("./data/", chr, ".snpeff.RData")
    load(eff.Rdata)
    snpeff.info <- snpeff[snpeff[, 1] %in% rownames(dat.res),]
    
    snpeff.info[,"eff"][grepl("IT", snpeff.info[,"eff"])] <- "Intergenic"
    snpeff.info[,"eff"][grepl("IR", snpeff.info[,"eff"])] <- "Intron"
    snpeff.info[,"eff"][grepl("IG", snpeff.info[,"eff"])] <- "Start_gained"
    snpeff.info[,"eff"][grepl("IL", snpeff.info[,"eff"])] <- "Start_lost"
    snpeff.info[,"eff"][grepl("SG", snpeff.info[,"eff"])] <- "Stop_gained"
    snpeff.info[,"eff"][grepl("SL", snpeff.info[,"eff"])] <- "Stop_lost"
    snpeff.info[,"eff"][grepl("Up", snpeff.info[,"eff"])] <- "Upstream"
    snpeff.info[,"eff"][grepl("Dn", snpeff.info[,"eff"])] <- "Downstream"
    snpeff.info[,"eff"][grepl("U3", snpeff.info[,"eff"])] <- "three_prime_UTR"
    snpeff.info[,"eff"][grepl("U5", snpeff.info[,"eff"])] <- "five_prime_UTR"
    snpeff.info[,"eff"][grepl("SSA", snpeff.info[,"eff"])] <- "Splice_site_acceptor"
    snpeff.info[,"eff"][grepl("SSD", snpeff.info[,"eff"])] <- "Splice_site_donor"
    snpeff.info[,"eff"][grepl("NSC", snpeff.info[,"eff"])] <- "Non_synonymous_coding"
    snpeff.info[,"eff"][grepl("NSS", snpeff.info[,"eff"])] <- "Non_synonymous_start"
    snpeff.info[,"eff"][grepl("SC", snpeff.info[,"eff"])] <- "Synonymous_coding"
    snpeff.info[,"eff"][grepl("SS", snpeff.info[,"eff"])] <- "Synonymous_stop"
    snpeff.info[,"eff"][grepl("IA", snpeff.info[,"eff"])] <- "Intergenic"
    
    snpeff.info <- snpeff.info[snpeff.info[, "eff"] %in% mutType, , drop=FALSE]
    
    dat.res <- dat.res[rownames(dat.res) %in% snpeff.info[, "id"], , drop=FALSE]
  }
  
  if (!is.null(snpSites) && length(snpSites)>=1) {
    dat.res <- dat.res[rownames(dat.res) %in% snpSites, , drop=FALSE]
  }
  
  #### calculate distance matrix
  dist.mat <- sapply(1:ncol(dat.res), function(x){
    colMeans(abs(dat.res - dat.res[,x]), na.rm=TRUE)
  })
  
  colnames(dist.mat) <- rownames(dist.mat)
  
  ### tree
  dist.mat <- as.dist(dist.mat)
  tre <- nj(dist.mat)
  
  p <- ggtree(tre, layout="circular", branch.length="none", size=0.01) + ggtitle("")
  p <- p + theme_void()
  p <- gheatmap(p, acc.tree, offset = 1, width=0.1, colnames = FALSE, color=NULL) +
    scale_fill_manual(breaks=c("Or-I", "Or-II", "Or-III", 
                               "Indica", "Japonica", "Other"), 
                      values=c("blue", "red", "black", 
                               "purple", "gold", "cyan"))
  figurecp <<- p
  treNwk <<- tre
  return(p)
}

