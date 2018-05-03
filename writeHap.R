
 
writeHap <- function(haplotype=NULL) {
  hap.con <- haplotype
  
  nex.con <- c("#NEXUS", "", "BEGIN TAXA;", paste0("DIMENSIONS NTAX=", length(hap.con[[1]]), ";"), "", "TAXLABELS",
               paste0("Hap_", names(hap.con[[1]])), ";", "END;", "", "", "[Hap#  Freq. Sequences]")
  
  for (i in 1:length(hap.con[[2]])) {
    con <- paste(paste0("[Hap_", names(hap.con[[2]])[i], ":"), length(hap.con[[2]][[i]]), paste(hap.con[[2]][[i]], collapse = " ") )
    con <- paste0(con, "]")
    nex.con <- c(nex.con, con)
  }
  
  nex.con <- c(nex.con, "", "", "BEGIN CHARACTERS;", paste0("DIMENSIONS NCHAR=", nchar(hap.con[[1]])[1], ";"),
               "FORMAT DATATYPE=DNA  MISSING=? GAP=- MATCHCHAR=.;", "MATRIX", "")
  
  for (i in 1:length(hap.con[[1]])) {
    con <- paste(paste0("Hap_", names(hap.con[[1]])[i]), hap.con[[1]][[i]] )
    nex.con <- c(nex.con, con)
  }
  
  nex.con <- c(nex.con, ";", "END;", "", "BEGIN TRAITS;")
  nex.con <- c(nex.con, paste0("  Dimensions NTRAITS=", ncol(hap.con[[3]]), ";"))
  nex.con <- c(nex.con, "  Format labels=yes missing=? separator=Comma;")
  nex.con <- c(nex.con, paste0("  TraitLabels ", paste(colnames(hap.con[[3]]), collapse = " "), ";"))
  nex.con <- c(nex.con, "  Matrix")
  for (i in 1:nrow(hap.con[[3]])) {
    nex.con <- c(nex.con, paste(paste0("Hap_", names(hap.con[[1]])[i]), paste(hap.con[[3]][i, ], collapse = ",")))
  }
  nex.con <- c(nex.con, ";", "", "END;")
  
  return(nex.con)
}

