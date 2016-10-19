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

require(PMA)
load("data_Saved_R_Objects/cca/mRNA_lncRNA_normal_cca.rda")
load("data_Saved_R_Objects/cca/mRNA_lncRNA_tumor_cca.rda")
gc()
mRNA_lncRNA_normal_perm_out <- CCA.permute(x = t(mRNA_normal),
                        z = t(lncRNA_normal),
                        typex="standard",
                        typez="standard",
                        nperms=100)
normal_out <- CCA(x = t(mRNA_normal),
                  z = t(lncRNA_normal), 
                  typex = "standard", 
                  typez = "standard",
                  penaltyx = mRNA_lncRNA_normal_perm_out$bestpenaltyx,
                  penaltyz = mRNA_lncRNA_normal_perm_out$bestpenaltyz,
                  v = mRNA_lncRNA_normal_perm_out$v.init,
                  K = 1) 
save(mRNA_lncRNA_normal_perm_out, normal_out, file = "data_Saved_R_Objects/cca/mRNA_lncRNA_normal_cca.rda")
gc()
mRNA_lncRNA_tumor_perm_out <- CCA.permute(x = t(mRNA_tumor),
                                           z = t(lncRNA_tumor),
                                           typex="standard",
                                           typez="standard",
                                           nperms=100)
tumor_out <- CCA(x = t(mRNA_normal),
                  z = t(lncRNA_normal), 
                  typex = "standard", 
                  typez = "standard",
                  penaltyx = mRNA_lncRNA_tumor_perm_out$bestpenaltyx,
                  penaltyz = mRNA_lncRNA_tumor_perm_out$bestpenaltyz,
                  v = mRNA_lncRNA_tumor_perm_out$v.init,
                  K = 1) 
save(mRNA_lncRNA_tumor_perm_out,tumor_out, file = "data_Saved_R_Objects/cca/mRNA_lncRNA_tumor_cca.rda")

######## -------------- Analyze Sparse CCA -----------------------------------------------------------
setwd("/media/ducdo/UUI1/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("get_normal_tumor_expression_matrices.R")
load("data_Saved_R_Objects/cca/mRNA_lncRNA_normal_cca.rda")
load("data_Saved_R_Objects/cca/mRNA_lncRNA_tumor_cca.rda")








