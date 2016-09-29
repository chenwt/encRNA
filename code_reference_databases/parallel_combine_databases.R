# setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
# load("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/data_Saved_R_Objects/miRNA_target/putative_binary_interaction.rda")
setwd("/home/MARQNET/0099dod/encRNA/combine_target_databases")
load("putative_binary_interaction.rda")
require(foreach);require(doParallel);require(doSNOW)

common_miRNAs = sort(intersect(starbase_mRNA_miRNA_pairs$miRNA,miRcode_lncRNA_pairs$miRNA))
length(common_miRNAs)

df1 = starbase_mRNA_miRNA_pairs[which(starbase_mRNA_miRNA_pairs$miRNA %in% common_miRNAs),]
df2 = miRcode_lncRNA_pairs[which(miRcode_lncRNA_pairs$miRNA %in% common_miRNAs),]

putative_encRNA = matrix(ncol = 3)
colnames(putative_encRNA) = c("lncRNA","mRNA","miRNA")
## ------------- try with parApply --------------------------------------------------------

# ptm = proc.time()
# require(parallel)
# num_cores = 3
# cl <- makeCluster(num_cores)
# output = parLapply(cl, common_miRNAs[1:3], function(miRNA,df1,df2){
#   mRNAs = df1[which(df1$miRNA == miRNA),c("putative_mRNA")]; 
#   lncRNAs = df2[which(df2$miRNA == miRNA),c("putative_lncRNAs")];
#   encRNA = expand.grid(x = lncRNAs, y = mRNAs)
#   print(nrow(encRNA))
#   colnames(encRNA) = c("lncRNA","mRNA")
#   encRNA$miRNA = rep(miRNA)
#   return(encRNA)
# },df1,df2)
# putative_encRNA = do.call(rbind,output)
# save(putative_encRNA, file = "parallel_putative_encRNA.rda")
# stopCluster(cl)
# ptm = proc.time() - ptm; ptm

## ------------- try with foreach --------------------------------------------------------
gc()
num_cluster = detectCores()
ptm = proc.time()
cl <- snow::makeCluster(num_cluster, outfile = "")
doSNOW::registerDoSNOW(cl)
putative_encRNA_batch_121_176 = foreach(i = 121:176, .combine = rbind) %dopar%{
  print(paste("Iteration number: ", i, sep = ""))
  mRNAs = df1[which(df1$miRNA == common_miRNAs[i]),c("putative_mRNA")]; 
  lncRNAs = df2[which(df2$miRNA == common_miRNAs[i]),c("putative_lncRNAs")];
  encRNA = expand.grid(x = lncRNAs, y = mRNAs)
  #print(nrow(encRNA))
  colnames(encRNA) = c("lncRNA","mRNA")
  encRNA$miRNA = rep(common_miRNAs[i])
  return(encRNA)
}
snow::stopCluster(cl)
ptm = proc.time() - ptm
putative_encRNA_batch_121_176 = data.table(putative_encRNA_batch_121_176)
print(paste("Total running time: ", ptm[3], sep = ""))
save(putative_encRNA_batch_121_176, file = "data_Saved_R_Objects/miRNA_target/putative_encRNA_batch_121_176.rda")
# sprintf("Total running time with %g cores: %g", num_cluster ,ptm[3])


## ------------- using data.table to combine resulst -------------------------------------
load("data_Saved_R_Objects/miRNA_target/putative_encRNA_batch_121_176.rda")
load("data_Saved_R_Objects/miRNA_target/putative_encRNA_batch_1_80.rda")
load("data_Saved_R_Objects/miRNA_target/putative_encRNA_batch_81_120.rda")

putative_encRNA_batch_1_120 = data.table::rbindlist(
  list(putative_encRNA_batch_1_80,putative_encRNA_batch_81_120))
putative_encRNA = data.table::rbindlist(
  list(putative_encRNA_batch_1_120,putative_encRNA_batch_121_176))

## ------------- using data.table to combine resulst -------------------------------------

load("data_Saved_R_Objects/brca_df.rda")
nrow(putative_encRNA) # [1] 395,384,991

brca_putative_encRNA = putative_encRNA[which(putative_encRNA$lncRNA %in% rownames(brca_lncRNA_df))]
nrow(brca_putative_encRNA) # 108,973,139
brca_putative_encRNA = brca_putative_encRNA[which(brca_putative_encRNA$mRNA %in% rownames(brca_mRNA_df))]
nrow(brca_putative_encRNA) # 104,588,830
brca_putative_encRNA = brca_putative_encRNA[which(brca_putative_encRNA$miRNA %in% rownames(brca_miRNA_df))]
nrow(brca_putative_encRNA) # 60,862,939

length(which(brca_putative_encRNA$lncRNA %in% rownames(brca_lncRNA_df)))
length(which(brca_putative_encRNA$mRNA %in% rownames(brca_mRNA_df)))
length(which(brca_putative_encRNA$miRNA %in% rownames(brca_miRNA_df)))

save(brca_putative_encRNA, file = "data_Saved_R_Objects/miRNA_target/brca_putative_encRNA.rda")
