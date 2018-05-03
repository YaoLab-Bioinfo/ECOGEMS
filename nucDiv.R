

nucDiv <- function(chr="chr07", nuc.start=2800000, nuc.end=2900000, step=10, 
                   groups=c("Wild", "Cultivar"), numerator="Wild", denominator="Cultivar",
                   mutType = NULL, snpSites = NULL) {
  data <- fetchSnp(chr=chr, start=nuc.start, end=nuc.end, mutType = mutType)[[1]]
  
  if (!is.null(snpSites) && length(snpSites)>=1) {
    data <- data[rownames(data) %in% snpSites, ]
  }
  
  dat.mat <- t(data)
  dat.mat[is.na(dat.mat)] <- "-"
  dat.bin <- as.DNAbin(dat.mat)
  
  div.group <- lapply(unique(c(groups, numerator, denominator)), function(x){
    x.accession <- readLines(paste0("./data/", x, ".acc.txt"))
    dat <- dat.bin[rownames(dat.bin) %in% x.accession, ]
    
    nuc.div <- lapply(seq(1, ncol(dat.bin), by=step), function(i){
      dat.i <- dat[, i:min(i+step-1, ncol(dat.bin))]
      
      # if(!is.matrix(dat.i)) {return(NULL)}
      div <- nuc.div(dat.i, pairwise.deletion = TRUE)
      
      return(div)
    })
    
    nuc.div.df <- do.call(rbind, nuc.div)
    nuc.div.df <- data.frame(nuc.div.df, stringsAsFactors = FALSE)
  })
  
  div.group.df <- do.call(cbind, div.group)
  names(div.group.df) <- unique(c(groups, numerator, denominator))
  dat.pos <- as.numeric(substr(colnames(dat.mat), 3, 10))
  nuc.pos <- dat.pos[seq(1, ncol(dat.mat), by=step)][1:nrow(div.group.df)]
  div.group.df$pos <- nuc.pos
  
  div.group.df.1 <- div.group.df[, c("pos", groups)]
  div.group.df.2 <- div.group.df[, c("pos", numerator, denominator)]
  
  div.group.df.1.long <- gather(div.group.df.1, group, diversity, -pos)
  div.group.df.2.long <- gather(div.group.df.2, group, diversity, -pos)
  
  nuc.chr <- substr(chr, 4, 5)
  nuc.gene.info <- gene.info
  nuc.gene.info$chr <- substr(nuc.gene.info$id, 7, 8)
  nuc.gene.info <- nuc.gene.info[nuc.gene.info$chr==nuc.chr & 
                                 nuc.gene.info$start>=as.numeric(nuc.start) &
                                 nuc.gene.info$end<=as.numeric(nuc.end), ]
  
  p1 <- ggplot(div.group.df.1.long) + geom_line(aes(x=pos, y=diversity, color=group))
  p1 <- p1 + xlab("") + ylab("Nucleotide diversity")
  p1 <- p1 + theme_classic() + ylim(-0.14, NA)
  p1 <- p1 + theme(legend.title = element_blank())
  p1 <- p1 + theme(legend.position="top")
  
  if (nrow(nuc.gene.info)>=1) {
    p1 <- p1 + geom_rect(aes(xmin=start, xmax=end, ymin=-0.05, ymax=-0.07), color="grey40", data=nuc.gene.info)
    p1 <- p1 + geom_text(aes(x=(start+end)/2, y=-0.12, label=id), angle=50, size=2.5, data=nuc.gene.info)
  }
  
  p1 <- p1 + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                   axis.line.x = element_blank())
  
  div.group.df.2$value <- div.group.df.2[,numerator]/div.group.df.2[,denominator]
  
  p2 <- ggplot(div.group.df.2) + geom_line(aes(x=pos, y=value))
  p2 <- p2 + xlab("genomic position") + ylab(paste0(numerator, "/", denominator))
  p2 <- p2 + theme_classic()
  
  gp1 <- ggplotGrob(p1)
  gp2 <- ggplotGrob(p2)
  maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
  gp1$widths[2:5] <- as.list(maxWidth)
  gp2$widths[2:5] <- as.list(maxWidth)
  grid.draw(
    grid.arrange(gp1, gp2, ncol=1, heights=c(2.3, 1))
  )
}

