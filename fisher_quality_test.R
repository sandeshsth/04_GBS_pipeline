#!/usr/bin/env Rscript

ARGV = commandArgs(trailingOnly=TRUE)

#options(scipen=999) # to skip e notation
options(scipen=0) # to skip e notation

df <- read.table(file=ARGV[1], header = TRUE, stringsAsFactors = FALSE, sep = "\t")

all_P <- vector(mode = "numeric", length = 0)
quality <- vector(mode = "numeric", length = 0)


for(i in 1:nrow(df)){
  rw <- df[i,]
  het <- as.integer(rw[5])
  ref <- as.integer(rw[6])
  alt <- as.integer(rw[7])
  nothing <- as.integer(rw[8])
  
  #print(nothing)
  
  matrix_table <- matrix(c(het, ref, alt, nothing), nrow = 2)
  #print(matrix_table)
  #break
  fisher_P <- fisher.test(matrix_table, alternative = "two.sided")$p.value
  
  # if the p-value is 0, -log10 will be Inf so changed to 1e-100
    if(fisher_P == 0){
        fisher_P = 1e-100
        #print(fisher_P)
    }
  
  all_P[i] <- fisher_P
  
  # quality as log p-value
  quality[i] <- -log10(fisher_P)
 
}

new_df <- cbind(df, "F_test_p" = all_P, "QUALITYSCORE" = quality)
another <- new_df[, c("CHROM", "POS", "QUALITYSCORE")] 

# write.table(new_df, "with_P_quality_snps.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(another, file=paste0("fisher_quality_", ARGV[1]), sep = "\t", quote = FALSE, row.names = FALSE)

