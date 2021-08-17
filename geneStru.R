

geneStru <- function(chr="chr09", start=37800, end=46400){ 
  start <- as.numeric(start)
  end <- as.numeric(end)
  reg.gr <- IRanges::IRanges(start, end)
  snp.lst.chr <- snp.lst[snp.lst$chr==chr, ]
  snp.lst.gr <- IRanges::IRanges(start=snp.lst.chr$start, end=snp.lst.chr$end)
  snp.fls <- snp.lst.chr$file[unique(S4Vectors::queryHits(IRanges::findOverlaps(snp.lst.gr, reg.gr)))]
  
  snp.allele.lst <- lapply(snp.fls, function(x){
    load(x)
    return(snp.data.allele)
  })
  snp.allele <- do.call(rbind, snp.allele.lst)
  snp.allele <- snp.allele[order(as.numeric(rownames(snp.allele))), ]
  
  start.c <- as.numeric(paste0(substr(chr, 4, 5), sprintf("%08d", start)))
  end.c <- as.numeric(paste0(substr(chr, 4, 5), sprintf("%08d", end)))
  
  dat <- snp.allele[as.numeric(rownames(snp.allele))>=start.c & as.numeric(rownames(snp.allele))<=end.c, ]
  snp.code.pos <- as.numeric(substr(rownames(dat), 3, 10))
  
  
  gff$id <- gsub(":.+", "", gff$id)
  gff.mrna <- gff[gff$type == "mRNA", ]
  gff.reg.mrna <- gff.mrna[gff.mrna$chr==chr & gff.mrna$start>=start & gff.mrna$end<=end, ]
  gff.reg <- gff[gff$id %in% gff.reg.mrna$id, ]
  
  gff.reg$anno <- paste(gff.reg$id, gff.reg$anno, sep=" <br> ")
  
  gff.reg.mrna.ir <- IRanges::IRanges(gff.reg.mrna$start, gff.reg.mrna$end)
  gff.reg.mrna.op <- IRanges::findOverlaps(gff.reg.mrna.ir, IRanges::reduce(gff.reg.mrna.ir))
  gff.reg.mrna$grp <- S4Vectors::subjectHits(gff.reg.mrna.op)
  
  gff.reg.mrna.1 <- gff.reg.mrna %>% group_by(grp) %>% mutate(y=row_number())
  
  gff.reg <- merge(gff.reg, gff.reg.mrna.1[, c("id", "y")], by="id")
  
  gff.reg$y <- gff.reg$y * 0.2 + 1
  
  plot.mrna.lst <- lapply(unique(gff.reg$id), function(i){
    dat <- gff.reg[gff.reg$id == i, ]
    i.strand <- dat$strand[1]
    
    dat.mrna <- dat[dat$type=="mRNA", ]
    return(dat.mrna)
  })
  plot.mrna <- do.call(rbind, plot.mrna.lst)
  p1 <- ggplot2::ggplot(plot.mrna) + ggplot2::geom_rect(ggplot2::aes(xmin=start, xmax=end, ymin=y+0.118, ymax=y+0.122,
                                          text=anno), 
                                      color="grey30", fill="grey30")
  
  
  plot.nm.lst <- lapply(unique(gff.reg$id), function(i){
    dat <- gff.reg[gff.reg$id == i, ]
    i.strand <- dat$strand[1]
    
    dat.nm <- dat[dat$type!="mRNA", ]
    dat.nm <- dat.nm[-nrow(dat.nm), ]
    
    if (nrow(dat.nm)>0) {
      dat.nm$ymin <- dat.nm$y+0.1
      dat.nm$ymax <- dat.nm$y+0.14
      dat.nm$ymin[dat.nm$type=="CDS"] <- dat.nm$ymin[dat.nm$type=="CDS"] - 0.02
      dat.nm$ymax[dat.nm$type=="CDS"] <- dat.nm$ymax[dat.nm$type=="CDS"] + 0.02
    }
    return(dat.nm)
  })
  plot.nm <- do.call(rbind, plot.nm.lst)
  if (nrow(plot.nm)>0) {
    p1 <- p1 + ggplot2::geom_rect(ggplot2::aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, text=anno), 
                         color="grey30", fill="grey30", data=plot.nm)
  }
  
  
  plot.tail.lst <- lapply(unique(gff.reg$id), function(i){
    dat <- gff.reg[gff.reg$id == i, ]
    i.strand <- dat$strand[1]
    
    dat.nm <- dat[dat$type!="mRNA", ]
    
    i.anno <- dat$anno[1]
    i.id <- dat.nm$pare[1]
    
    tail.type <- dat.nm$type[nrow(dat.nm)]
    
    dat.tail <- data.frame(xx=rep(c(dat$start[nrow(dat)], 
                                    (dat$start[nrow(dat)] + dat$end[nrow(dat)])/2, dat$end[nrow(dat)]), each=2), 
                           stringsAsFactors = FALSE)
    if (i.strand == "-") {
      dat.tail$yy <- c(0.12, 0.12, 0.1, 0.14, 0.1, 0.14) + dat$y[1]
      dat.tail <- dat.tail[c(1,3,5,6,4,2), ]
      dat.tail$pare <- i.id
      dat.tail$anno <- i.anno
      if (tail.type=="CDS") {
        dat.tail$yy[2:3] <- dat.tail$yy[2:3] - 0.02
        dat.tail$yy[4:5] <- dat.tail$yy[4:5] + 0.02
      }
    } else {
      dat.tail$yy <- c(0.1, 0.14, 0.1, 0.14, 0.12, 0.12) + dat$y[1]
      dat.tail <- dat.tail[c(1,3,5,6,4,2), ]
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
  p1 <- p1 + ggplot2::geom_polygon(ggplot2::aes(x=xx, y=yy, group=pare), color="grey30", fill="grey30", 
                          data=plot.tail)
  
  snp.pos.df <- data.frame(x=snp.code.pos, ymin=1.23, ymax=1.25, stringsAsFactors = FALSE)
  p1 <- p1 + ggplot2::geom_linerange(ggplot2::aes(x=x, ymin=ymin, ymax=ymax), data=snp.pos.df)
  p1 <- p1 + ggplot2::geom_segment(ggplot2::aes(x=min(snp.code.pos), xend=max(snp.code.pos), y=1.25, yend=1.25))
  
  p1 <- p1 + ggplot2::scale_y_continuous("", breaks=NULL)
  p1 <- p1 + ggplot2::theme(panel.grid.major=ggplot2::element_blank(),panel.grid.minor=ggplot2::element_blank()) + 
    ggplot2::theme(panel.background=ggplot2::element_rect(fill="white",colour="white"))
  p1 <- p1 + ggplot2::scale_x_continuous("", breaks=NULL)
  
  return(p1)
}


