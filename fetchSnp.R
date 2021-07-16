
# A function to extract SNP data in a specified genomic region.
# Change to the directory of ECOGEMS using the setwd function of R.
# Usage: type the next two lines in R Console without the leading #
# source("Global.R")
# snp.data <- fetchSnp(chr="chr07", start=29616705, end=29629223, accession=c("Aus", "TeJ"), mutType=NULL)
# The output snp.data is a list with two elements: the genotype matrix and the allele matrix

fetchSnp <- function(chr="chr07", start=29616705, end=29629223, accession=NULL, mutType=NULL){
  
  if (is.null(chr)) {
    return(NULL)
  } else {
    chr.size <- chrInfo$size[chrInfo$chr == chr]
    start <- max(0, start)
    end <- min(end, chr.size)
    
    start <- as.numeric(start)
    end <- as.numeric(end)
    reg.gr <- IRanges::IRanges(start, end)
    snp.lst.chr <- snp.lst[snp.lst$chr==chr, ]
    snp.lst.gr <- IRanges::IRanges(start=snp.lst.chr$start, end=snp.lst.chr$end)
    snp.fls <- snp.lst.chr$file[unique(S4Vectors::queryHits(IRanges::findOverlaps(snp.lst.gr, reg.gr)))]
    
    snp.fls.lst <- lapply(snp.fls, function(x){
      load(x)
      return(list(snp.data.inter.Matrix, snp.data.allele))
    })
    snp.fls.lst <- unlist(snp.fls.lst, recursive = FALSE)
    snp.data <- do.call(rbind, snp.fls.lst[seq(1, length(snp.fls.lst), by=2)])
    snp.data <- snp.data[order(as.numeric(rownames(snp.data))), ]
    snp.allele <- do.call(rbind, snp.fls.lst[seq(2, length(snp.fls.lst), by=2)])
    snp.allele <- snp.allele[order(as.numeric(rownames(snp.allele))), ]
    
    start <- as.numeric(paste0(substr(chr, 4, 5), sprintf("%08d", start)))
    end <- as.numeric(paste0(substr(chr, 4, 5), sprintf("%08d", end)))
    
    dat.res <- snp.data[as.numeric(rownames(snp.data))>=start & as.numeric(rownames(snp.data))<=end, , drop=FALSE]
    snp.code <- as.vector(t(snp.allele[rownames(dat.res), ]))
    allele.res <- snp.allele[rownames(dat.res), , drop=FALSE]
    
    dat.res.n <- seq_len(nrow(dat.res)) - 1
    dat.res.n <- rep(dat.res.n, each=ncol(dat.res))
    dat.res.n <- matrix(dat.res.n, ncol=ncol(dat.res), byrow = T)
    
    dat.res <- dat.res + 1 + dat.res.n * 2
    dat.res.mat <- matrix(snp.code[as.matrix(dat.res)], ncol=ncol(dat.res))
    rownames(dat.res.mat) <- rownames(dat.res)
    colnames(dat.res.mat) <- colnames(dat.res)
    
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
      dat.res.mat <- dat.res.mat[, colnames(dat.res.mat) %in% accession, drop=FALSE]
    }
    
    dat.res.mat.row.c <- apply(dat.res.mat, 1, function(x){
      length(unique(x[!is.na(x)]))
    })
    
    dat.res.mat <- dat.res.mat[dat.res.mat.row.c>1, , drop=FALSE]
    allele.res <- allele.res[dat.res.mat.row.c>1, , drop=FALSE]
    
    if (!is.null(mutType) && length(mutType)>=1 && length(mutType)!=16) {
      eff.Rdata <- paste0("./data/", chr, ".snpeff.RData")
      load(eff.Rdata)
      snpeff.info <- snpeff[snpeff[, 1] %in% rownames(dat.res.mat), , drop=FALSE]
      
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
      
      dat.res.mat <- dat.res.mat[rownames(dat.res.mat) %in% snpeff.info[, "id"], , drop=FALSE]
      allele.res <- allele.res[rownames(allele.res) %in% snpeff.info[, "id"], , drop=FALSE]
    }
    
    return(list(dat.res.mat, allele.res))
  }
  
}

