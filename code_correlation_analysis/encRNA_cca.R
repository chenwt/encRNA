setwd("/media/ducdo/UUI1/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("get_normal_tumor_expression_matrices.R")
# load("data_Saved_R_Objects/brca_df.rda")

# setwd("/home/MARQNET/0099dod/encRNA")
# load("brca_df.rda")

# require(CCA)
# 
# normal_indices = 1:79
# tumor_indices = 90:ncol(brca_lncRNA_df)
# 
# mRNA_normal = brca_mRNA_df[,normal_indices];
# miRNA_normal = brca_miRNA_df[,normal_indices];
# lncRNA_normal = brca_lncRNA_df[,normal_indices];
# 
# mRNA_tumor = brca_mRNA_df[,tumor_indices];
# miRNA_tumor = brca_miRNA_df[,tumor_indices];
# lncRNA_tumor = brca_lncRNA_df[,tumor_indices]
# 
# dim(mRNA_normal) # 17,613 * 79
# dim(lncRNA_normal) # 4,828 * 79

# estim.regul(X = t(mRNA_normal[1:100,]), 
#             Y = t(lncRNA_normal[1:100,]), 
#             grid1 = c(0.01,0.5),
#             grid2 = c(0.1,0.2,0.3),
#             plt = TRUE)

# mRNA_lncRNA_normal_cc = rcc(X = t(mRNA_normal), 
#                             Y = t(lncRNA_normal), 
#                             lambda1 = 0.01, 
#                             lambda2 = 0.1)
# save(mRNA_lncRNA_normal_cc, file = "mRNA_lncRNA_normal_cc.rda")
# load("code_correlation_analysis/mRNA_lncRNA_normal_cc.rda"); gc()
# mRNA_lncRNA_normal_cc$cor[1:10]


######## -------------- Compute Sparse CCA -----------------------------------------------------------

### from all normal and all tumor
require(PMA)
load("data_Saved_R_Objects/cca/mRNA_lncRNA_normal_cca.rda")
load("data_Saved_R_Objects/cca/mRNA_lncRNA_tumor_cca.rda")

mRNA_lncRNA_normal_scca = get_scca_result(x = t(mRNA_normal), z = t(lncRNA_normal), nperms = 100)
mRNA_lncRNA_normal_perm_out = mRNA_lncRNA_normal_scca$perm
normal_out = mRNA_lncRNA_normal_scca$out

mRNA_lncRNA_tumor_scca = get_scca_result(x = t(mRNA_tumor), z = t(lncRNA_tumor), nperms = 100)
mRNA_lncRNA_tumor_perm_out = mRNA_lncRNA_tumor_scca$perm
tumor_out = mRNA_lncRNA_tumor_scca$out

### from matched normal_tumor
load("data_Saved_R_Objects/brca_expression_matched.rda")
source("get_normal_tumor_expression_matrices.R")
mRNA_lncRNA_normal_matched_scca = get_scca_result(x = t(mRNA_normal_matched), 
                                                  z = t(lncRNA_normal_matched), 
                                                  nperms = 100)
save(mRNA_lncRNA_normal_matched_scca, 
     file = "data_Saved_R_Objects/cca/mRNA_lncRNA_normal_matched_scca.rda")

######## -------------- Analyze Sparse CCA -----------------------------------------------------------
setwd("/media/ducdo/UUI1/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("get_normal_tumor_expression_matrices.R")
load("data_Saved_R_Objects/cca/mRNA_lncRNA_normal_cca.rda")
load("data_Saved_R_Objects/cca/mRNA_lncRNA_tumor_cca.rda")

#### how many mRNAs and lncRNAs are kept 
# normal
length(which(normal_out$u != 0)) # mRNAs: 9,251 out of 17,613
length(which(normal_out$v != 0)) # lncRNAs: 2720 out of 4,828  
# tumor
length(which(tumor_out$u != 0)) # mRNAs: 723 out of 17,613
length(which(tumor_out$v != 0)) # lncRNAs: 195 out of 4,828  


######## -------------- Sparse CCA function -----------------------------------------------------------

get_scca_result = function(x, z, nperms = 100){
  require(PMA)
  perm_out <- CCA.permute(x = x,
                          z = z,
                          typex = "standard",
                          typez = "standard",
                          nperms = nperms)
  out <- CCA(x = x,
             z = z, 
             typex = "standard", 
             typez = "standard",
             penaltyx = perm_out$bestpenaltyx,
             penaltyz = perm_out$bestpenaltyz,
             v = perm_out$v.init,
             K = 1)
  l = list(perm = perm_out, out = out)
  return(l)
}




