# setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
# load("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/data_Saved_R_Objects/miRNA_target/putative_binary_interaction.rda")

setwd("/home/MARQNET/0099dod/encRNA/combine_target_databases")
load("putative_binary_interaction.rda")



putative_encRNA = matrix(ncol = 3)
colnames(putative_encRNA) = c("lncRNA","mRNA","miRNA")

ptm = proc.time()
i = 1;
sapply(common_miRNAs, function(miRNA){
  #print(i); i <<- i + 1
  mRNAs = df1[which(df1$miRNA == miRNA),c("putative_mRNA")]; 
  lncRNAs = df2[which(df2$miRNA == miRNA),c("putative_lncRNAs")]; 
  encRNA = expand.grid(x = lncRNAs, y = mRNAs)
  colnames(encRNA) = c("lncRNA","mRNA")
  encRNA$miRNA = rep(miRNA)
  putative_encRNA <<- rbind(putative_encRNA,encRNA)
  #print(dim(putative_encRNA))
})
putative_encRNA = putative_encRNA[-1,]
ptm = proc.time() - ptm; ptm
save(putative_encRNA, file = "putative_encRNA.rda")

common_miRNAs


## -------------- parallelized version ----------------------------------------------------
common_miRNAs = sort(intersect(starbase_mRNA_miRNA_pairs$miRNA,miRcode_lncRNA_pairs$miRNA))
length(common_miRNAs)

df1 = starbase_mRNA_miRNA_pairs[which(starbase_mRNA_miRNA_pairs$miRNA %in% common_miRNAs),]
df2 = miRcode_lncRNA_pairs[which(miRcode_lncRNA_pairs$miRNA %in% common_miRNAs),]

putative_encRNA = matrix(ncol = 3)
colnames(putative_encRNA) = c("lncRNA","mRNA","miRNA")

ptm = proc.time()
require(parallel)
num_cores = 3
cl <- makeCluster(num_cores)
output = parLapply(cl, common_miRNAs[1:10], function(miRNA,df1,df2){
  mRNAs = df1[which(df1$miRNA == miRNA),c("putative_mRNA")]; 
  lncRNAs = df2[which(df2$miRNA == miRNA),c("putative_lncRNAs")];
  encRNA = expand.grid(x = lncRNAs, y = mRNAs)
  print(nrow(encRNA))
  # colnames(encRNA) = c("lncRNA","mRNA")
  # encRNA$miRNA = rep(miRNA)
  return(encRNA)
},df1,df2)
putative_encRNA = do.call(rbind,output)
stopCluster(cl)
ptm = proc.time() - ptm; ptm


