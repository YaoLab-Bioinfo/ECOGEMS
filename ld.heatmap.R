

ld.heatmap <- function(chr="chr09", start=37800, end=46400, snp.pos=c(1), 
                       gene=FALSE, ld.y=0.64, ld.w=0.80, flip=FALSE, accession=NULL, 
                       mutType=NULL, snpSites = NULL, ...){
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
  snp.allele.lst <- lapply(snp.fls, function(x){
    load(x)
    return(snp.data.allele)
  })
  snp.allele <- do.call(rbind, snp.allele.lst)
  snp.allele <- snp.allele[order(as.numeric(rownames(snp.allele))), ]
  
  start.c <- as.numeric(paste0(substr(chr, 4, 5), sprintf("%08d", start)))
  end.c <- as.numeric(paste0(substr(chr, 4, 5), sprintf("%08d", end)))
  
  dat <- snp.data[as.numeric(rownames(snp.data))>=start.c & as.numeric(rownames(snp.data))<=end.c, ]
  
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
    dat <- dat[, colnames(dat) %in% accession]
  }
  
  dat.row.c <- apply(dat, 1, function(x){
    length(unique(x[!is.na(x)]))
  })
  dat <- dat[dat.row.c>1, ]
  
  if (!is.null(mutType) && length(mutType)>=1) {
    eff.Rdata <- paste0("./data/", chr, ".snpeff.RData")
    load(eff.Rdata)
    snpeff.info <- snpeff[snpeff[, 1] %in% rownames(dat),]
    
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
    
    snpeff.info <- snpeff.info[snpeff.info[, "eff"] %in% mutType, ]
    
    dat <- dat[rownames(dat) %in% snpeff.info[, "id"], ]
  }
  
  if (!is.null(snpSites) && length(snpSites)>=1) {
    dat <- dat[rownames(dat) %in% snpSites, ]
  }
  
  snp.code <- snp.allele[rownames(dat), ]
  snp.code.pos <- as.numeric(substr(rownames(snp.code), 3, 10))
  
  dat <- as.matrix(dat)
  dat <- t(dat)
  dat[dat==1] <- 2
  
  dat.snp.mat <- as(dat, "SnpMatrix")
  
  if (gene) {
    ll <- LDheatmap(dat.snp.mat, snp.code.pos, 
                    flip=TRUE, title=NULL, ...)
    
    p1 <- geneStru(chr=chr, start=start, end=end)
    
    plot.new()
    llQplot2 <- LDheatmap.addGrob(ll, rectGrob(gp=gpar(col="white")), height=.3)
    pushViewport(viewport(x=0.483, y=ld.y, width=ld.w, height=.1))
    
    grid.draw(ggplotGrob(p1))
  } else {
    if (flip) {
      LDheatmap(dat.snp.mat, snp.code.pos, flip=TRUE, title=NULL, ...)
    } else {
      LDheatmap(dat.snp.mat, snp.code.pos, flip=FALSE, SNP.name = colnames(dat)[snp.pos], title=NULL, ...)
    }
  }
  
}

