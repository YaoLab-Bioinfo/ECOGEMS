
snpInfo <- function(chr="chr07", start=29616705, end=29629223, accession=NULL, mutType=NULL) {
  snp.info <- fetchSnp(chr=chr, start=start, end=end, 
                       accession = accession)
  
  eff.Rdata <- paste0("./data/", chr, ".snpeff.RData")
  load(eff.Rdata)
  snpeff.info <- snpeff[snpeff[,1] %in% rownames(snp.info[[1]]), , drop=FALSE]
  
  snpeff.info[,"eff"] <- gsub("IT", "Intergenic", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("IR", "Intron", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("IG", "Start_gained", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("IL", "Start_lost", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("SG", "Stop_gained", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("SL", "Stop_lost", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("Up", "Upstream", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("Dn", "Downstream", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("U3", "three_prime_UTR", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("U5", "five_prime_UTR", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("SSA", "Splice_site_acceptor", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("SSD", "Splice_site_donor", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("NSC", "Non_synonymous_coding", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("NSS", "Non_synonymous_start", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("SC", "Synonymous_coding", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("SS", "Synonymous_stop", snpeff.info[,"eff"])
  snpeff.info[,"eff"] <- gsub("IA", "Intergenic", snpeff.info[,"eff"])
  
  colnames(snpeff.info) <- c("snpID", "reference", "alternative", "effect")
  
  if (!is.null(mutType) && length(mutType)>=1) {
    snpeff.info <- cbind(snpeff.info, eff="")
    snpeff.info[,"eff"][grepl("Intergenic", snpeff.info[,"effect"])] <- "Intergenic"
    snpeff.info[,"eff"][grepl("Intron", snpeff.info[,"effect"])] <- "Intron"
    snpeff.info[,"eff"][grepl("Start_gained", snpeff.info[,"effect"])] <- "Start_gained"
    snpeff.info[,"eff"][grepl("Start_lost", snpeff.info[,"effect"])] <- "Start_lost"
    snpeff.info[,"eff"][grepl("Stop_gained", snpeff.info[,"effect"])] <- "Stop_gained"
    snpeff.info[,"eff"][grepl("Stop_lost", snpeff.info[,"effect"])] <- "Stop_lost"
    snpeff.info[,"eff"][grepl("Upstream", snpeff.info[,"effect"])] <- "Upstream"
    snpeff.info[,"eff"][grepl("Downstream", snpeff.info[,"effect"])] <- "Downstream"
    snpeff.info[,"eff"][grepl("three_prime_UTR", snpeff.info[,"effect"])] <- "three_prime_UTR"
    snpeff.info[,"eff"][grepl("five_prime_UTR", snpeff.info[,"effect"])] <- "five_prime_UTR"
    snpeff.info[,"eff"][grepl("Splice_site_acceptor", snpeff.info[,"effect"])] <- "Splice_site_acceptor"
    snpeff.info[,"eff"][grepl("Splice_site_donor", snpeff.info[,"effect"])] <- "Splice_site_donor"
    snpeff.info[,"eff"][grepl("Non_synonymous_coding", snpeff.info[,"effect"])] <- "Non_synonymous_coding"
    snpeff.info[,"eff"][grepl("Non_synonymous_start", snpeff.info[,"effect"])] <- "Non_synonymous_start"
    snpeff.info[,"eff"][grepl("Synonymous_coding", snpeff.info[,"effect"])] <- "Synonymous_coding"
    snpeff.info[,"eff"][grepl("Synonymous_stop", snpeff.info[,"effect"])] <- "Synonymous_stop"
    snpeff.info[,"eff"][grepl("Intergenic", snpeff.info[,"effect"])] <- "Intergenic"
    
    snpeff.info <- snpeff.info[snpeff.info[, "eff"] %in% mutType, , drop=FALSE]
    snpeff.info <- snpeff.info[, c("snpID", "reference", "alternative", "effect"), drop=FALSE]
  }
  
  snp.allele <- as.data.frame(snp.info[[2]], stringsAsFactors=FALSE)
  snp.allele$snpID <- rownames(snp.allele)
  rownames(snp.allele) <- NULL
  
  dat.res <- merge(snp.allele, snpeff.info, by="snpID")
  return(list(snp.info, dat.res))
}

