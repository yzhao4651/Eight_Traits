###function for finding the genes anotation 
Findneargenes <- function(mySNPs, search_radius){
  saved_genes <- list()
  length(saved_genes) <- nrow(mySNPs)
  names(saved_genes) <- as.character(mySNPs$SNP)
  chromnames <- paste("Chr", formatC(mySNPs$Chromosome, width = 2, flag = '0'), sep = "")
  for(i in 1:nrow(mySNPs)){
    gr <- GRanges(chromnames[i],
                  IRanges(mySNPs$Position[i] - search_radius,
                          mySNPs$Position[i] + search_radius))
    mygenes <- transcriptsByOverlaps(txdb, gr)
    saved_genes[[i]] <- mygenes$tx_name
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        Gene = unlist(saved_genes))
  return(gene_df)
}