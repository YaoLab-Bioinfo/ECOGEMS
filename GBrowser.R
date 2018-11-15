

GBrowser <- function(chr="chr07", start=29616705, end=29629223, accession=NULL, mutType=NULL) {
  start <- as.numeric(start)
  end <- as.numeric(end)
  
  set.seed(123)
  
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
  
  snp.reg <- fetchSnp(chr=chr, start=start, end=end, accession=accession, mutType=mutType)[[1]]
  
  eff.Rdata <- paste0("./data/", chr, ".snpeff.RData")
  load(eff.Rdata)
  
  snpeff.reg <- snpeff[snpeff[,1] %in% rownames(snp.reg), , drop=FALSE]
  snpeff.reg <- as.data.frame(snpeff.reg, stringsAsFactors = FALSE)
  snpeff.reg$chr <- paste0("chr", substr(snpeff.reg$id, 1, 2))
  snpeff.reg$pos <- as.numeric(substr(snpeff.reg$id, 3, 10))
  
  snpeff.reg$tag <- ""
  snpeff.reg$tag[grepl("IT", snpeff.reg$eff)] <- 1
  snpeff.reg$tag[grepl("IR", snpeff.reg$eff)] <- 4
  snpeff.reg$tag[grepl("SC", snpeff.reg$eff)] <- 7
  snpeff.reg$tag[grepl("SS", snpeff.reg$eff)] <- 8
  snpeff.reg$tag[grepl("SSA", snpeff.reg$eff)] <- 9
  snpeff.reg$tag[grepl("SSD", snpeff.reg$eff)] <- 10
  snpeff.reg$tag[grepl("NSC", snpeff.reg$eff)] <- 11
  snpeff.reg$tag[grepl("NSS", snpeff.reg$eff)] <- 12
  snpeff.reg$tag[grepl("IG", snpeff.reg$eff)] <- 13
  snpeff.reg$tag[grepl("IL", snpeff.reg$eff)] <- 14
  snpeff.reg$tag[grepl("SG", snpeff.reg$eff)] <- 15
  snpeff.reg$tag[grepl("SL", snpeff.reg$eff)] <- 16
  snpeff.reg$tag[grepl("Up", snpeff.reg$eff)] <- 2
  snpeff.reg$tag[grepl("Dn", snpeff.reg$eff)] <- 3
  snpeff.reg$tag[grepl("U3", snpeff.reg$eff)] <- 6
  snpeff.reg$tag[grepl("U5", snpeff.reg$eff)] <- 5
  
  snpeff.reg$eff <- gsub("IT", "Intergenic", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("IR", "Intron", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("SS", "Synonymous_stop", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("IG", "Start_gained", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("IL", "Start_lost", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("SG", "Stop_gained", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("SL", "Stop_lost", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("Up", "Upstream", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("Dn", "Downstream", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("U3", "three_prime_UTR", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("U5", "five_prime_UTR", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("SSA", "Splice_site_acceptor", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("SSD", "Splice_site_donor", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("NSC", "Non_synonymous_coding", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("NSS", "Non_synonymous_start", snpeff.reg$eff)
  snpeff.reg$eff <- gsub("SC", "Synonymous_coding", snpeff.reg$eff)
  
  
  snpeff.reg.1 <- snpeff.reg %>% group_by(id, chr, pos, ref, alt) %>% summarise(info=paste(eff, collapse="<br>"))
  snpeff.reg.2 <- snpeff.reg %>% group_by(id, chr, pos, ref, alt) %>% summarise(tag=max(tag))
  snpeff.reg.3 <- merge(snpeff.reg.1, snpeff.reg.2, by=c("id", "chr", "pos", "ref", "alt"))
  snpeff.reg.3$yr <- runif(nrow(snpeff.reg.3), min=-0.9, max=1.1)
  eff.tags <- c("Intergenic", "Intron", "Synonymous_coding","Synonymous_stop",
                "Splice_site_acceptor","Splice_site_donor",
                "Non_synonymous_coding","Non_synonymous_start", "Start_gained",
                "Start_lost","Stop_gained","Stop_lost","Upstream",
                "Downstream","three_prime_UTR","five_prime_UTR")
  names(eff.tags) <- c(1,4,7:16,2:3,6:5)
  snpeff.reg.3$tag <- eff.tags[snpeff.reg.3$tag]
  
  if (!is.null(mutType)) {
    snpeff.reg.3 <- snpeff.reg.3[snpeff.reg.3$tag %in% mutType, , drop=FALSE]
  }
  
  p1 <- ggplot(data=snpeff.reg.3) + geom_point(aes(x=pos, y=yr, color=tag, text=info, fill=tag), size=0.8, pch=25)
  
  
  gff$id <- gsub(":.+", "", gff$id)
  gff.mrna <- gff[gff$type == "mRNA", , drop=FALSE]
  gff.reg.mrna <- gff.mrna[gff.mrna$chr==chr & gff.mrna$start>=start & gff.mrna$end<=end, , drop=FALSE]
  gff.reg <- gff[gff$id %in% gff.reg.mrna$id, , drop=FALSE]
  
  gff.reg$anno <- paste(gff.reg$id, gff.reg$anno, sep=" <br> ")
  
  gff.reg.mrna.ir <- IRanges(gff.reg.mrna$start, gff.reg.mrna$end)
  gff.reg.mrna.op <- findOverlaps(gff.reg.mrna.ir, reduce(gff.reg.mrna.ir))
  gff.reg.mrna$grp <- subjectHits(gff.reg.mrna.op)
  
  gff.reg.mrna.1 <- gff.reg.mrna %>% group_by(grp) %>% mutate(y=row_number())
  
  gff.reg <- merge(gff.reg, gff.reg.mrna.1[, c("id", "y")], by="id")
  
  gff.reg$y <- gff.reg$y * 0.2 + 1
  
  plot.nm.lst <- lapply(unique(gff.reg$id), function(i){
    dat <- gff.reg[gff.reg$id == i, , drop=FALSE]
    i.strand <- dat$strand[1]
    
    if (i.strand == "-") {
      dat$y <- -dat$y
    }
    
    dat.nm <- dat[dat$type!="mRNA", , drop=FALSE]
    dat.nm <- dat.nm[-nrow(dat.nm), , drop=FALSE]
    
    if (nrow(dat.nm)>0) {
      dat.nm$ymin <- dat.nm$y+0.1
      dat.nm$ymax <- dat.nm$y+0.14
      dat.nm$ymin[dat.nm$type=="CDS"] <- dat.nm$ymin[dat.nm$type=="CDS"] - 0.02
      dat.nm$ymax[dat.nm$type=="CDS"] <- dat.nm$ymax[dat.nm$type=="CDS"] + 0.02
    }
    return(dat.nm)
  })
  plot.nm <- do.call(rbind, plot.nm.lst)
  if (!is.null(plot.nm) && nrow(plot.nm)>0) {
    p1 <- p1 + geom_rect(aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax), 
                         color="grey30", fill="grey30", data=plot.nm)
  }
  
  plot.mrna.lst <- lapply(unique(gff.reg$id), function(i){
    dat <- gff.reg[gff.reg$id == i, , drop=FALSE]
    i.strand <- dat$strand[1]
    
    if (i.strand == "-") {
      dat$y <- -dat$y
    }
    
    dat.mrna <- dat[dat$type=="mRNA", , drop=FALSE]
    return(dat.mrna)
  })
  plot.mrna <- do.call(rbind, plot.mrna.lst)
  if (!is.null(plot.mrna) && nrow(plot.mrna)>0) {
    p1 <- p1 + geom_rect(aes(xmin=start, xmax=end, ymin=y+0.118, ymax=y+0.122), 
                         color="grey30", fill="grey30", data=plot.mrna)
  }
  
  plot.tail.lst <- lapply(unique(gff.reg$id), function(i){
    dat <- gff.reg[gff.reg$id == i, , drop=FALSE]
    i.strand <- dat$strand[1]
    
    if (i.strand == "-") {
      dat$y <- -dat$y
    }
    
    dat.nm <- dat[dat$type!="mRNA", , drop=FALSE]
    
    i.anno <- dat$anno[1]
    i.id <- dat.nm$pare[1]
    
    tail.type <- dat.nm$type[nrow(dat.nm)]
    
    dat.tail <- data.frame(xx=rep(c(dat$start[nrow(dat)], 
                                    (dat$start[nrow(dat)] + dat$end[nrow(dat)])/2, dat$end[nrow(dat)]), each=2), 
                           stringsAsFactors = FALSE)
    if (i.strand == "-") {
      dat.tail$yy <- c(0.12, 0.12, 0.1, 0.14, 0.1, 0.14) + dat$y[1]
      dat.tail <- dat.tail[c(1,3,5,6,4,2), , drop=FALSE]
      dat.tail$pare <- i.id
      dat.tail$anno <- i.anno
      if (tail.type=="CDS") {
        dat.tail$yy[2:3] <- dat.tail$yy[2:3] - 0.02
        dat.tail$yy[4:5] <- dat.tail$yy[4:5] + 0.02
      }
    } else {
      dat.tail$yy <- c(0.1, 0.14, 0.1, 0.14, 0.12, 0.12) + dat$y[1]
      dat.tail <- dat.tail[c(1,3,5,6,4,2), , drop=FALSE]
      dat.tail$pare <- i.id
      dat.tail$anno <- i.anno
      if (tail.type=="CDS") {
        dat.tail$yy[1:2] <- dat.tail$yy[1:2] - 0.02
        dat.tail$yy[5:6] <- dat.tail$yy[5:6] + 0.02
      }
    }
    
    return(dat.tail)
  })
  plot.tail <- do.call(rbind, plot.tail.lst)
  if (!is.null(plot.tail) && nrow(plot.tail)>0) {
    p1 <- p1 + geom_polygon(aes(x=xx, y=yy, group=pare), color="grey30", fill="grey30", 
                            data=plot.tail)
  }
  
  
  p1 <- p1 + scale_y_continuous("", breaks=NULL)
  p1 <- p1 + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
    theme(panel.background=element_rect(fill="white",colour="white"))
  p1 <- p1 + xlab("Chromosome position")
  p1 <- p1 + guides(color=guide_legend(title=NULL))
  
  p1 <- p1 + guides(fill=FALSE)
  
  p3 <- p1 + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                   axis.line.y = element_blank())
  p3 <- ggplotly(p3, tooltip = c("pos", "info"))
  
  p3 <- p3 %>% layout(
    title = "",
    xaxis = list(
      rangeselector = list(),
      
      rangeslider = list(type = "category")
      ),
    
    yaxis = list(title = "")
  )
  
  return(list(p1, p3))
}

