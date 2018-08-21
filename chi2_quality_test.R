options(scipen=0) # to skip e notation

ARGV = commandArgs(trailingOnly=TRUE)

df <- read.table(file=ARGV[1], header = TRUE, stringsAsFactors = FALSE, sep = "\t")
all_chi2 <- vector(mode = "numeric", length = 0)
all_chiP <- vector(mode = "numeric", length = 0)
all_log_chiP <- vector(mode = "numeric", length = 0)
quality <- vector(mode = "numeric", length = 0)

for(i in 1:nrow(df)){
  rw <- df[i,]
  het <- as.integer(rw[5])
  ref <- as.integer(rw[6])
  alt <- as.integer(rw[7])
  total <- het+ref+alt
  
  observed  = c(ref, het, alt)
  
  allele1 = (ref+(het/2))/total
  allele2 = (alt+(het/2))/total
  
  # 94% inbred: 
  exp_ref = allele1^2 + 2*allele1*allele2*0.47
  exp_het = 2*allele1*allele2*0.06
  exp_alt = allele2^2 + 2*allele1*allele2*0.47
  
  expected.prop  = c(exp_ref, exp_het, exp_alt)
  expected.count = sum(observed)*expected.prop
  chi2 = sum((observed- expected.count)^2/ expected.count)
  all_chi2[i] <- chi2
 
# critical value, 2df: For alpha,0.01 =9.21  For alpha, 0.001= 13.816
  
  if(is.na(chi2)){
    quality[i] <- 0
  }
  else if(chi2 < 9.21){
    quality[i] <- 100
  }
  else{
    quality[i] <- 0
 }
}
new_df <- cbind(df, "Chi2" = all_chi2, "QUALITYSCORE" = quality)
another <- new_df[, c("CHROM", "POS", "QUALITYSCORE")] 

#table(quality)
#write.table(new_df, "with_chi2_snps.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(another, file=paste0("chi2_quality_", ARGV[1]), sep = "\t", quote = FALSE, row.names = FALSE)


