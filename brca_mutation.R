# -------------- Basic setup ---------------------------------------------------------------------
setwd("/media/ducdo/UUI1/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
load("brca_mutation.rda")
load("Saved_R_Objects/brca_common_objects3.rda")
require(TCGA2STAT)

# brca_mutation <- getTCGA(disease="BRCA", 
#                          data.type="Mutation", 
#                          type="somatic", 
#                          clinical = TRUE)
# 
# save(brca_mutation, file = "brca_mutation.rda")

ncol(brca_lncRNA_common); colnames(brca_lncRNA_common)[1:10]
ncol(brca_mRNA_common); colnames(brca_mRNA_common)[1:10]
ncol(brca_miRNA_common); colnames(brca_miRNA_common)[1:10]
# ----------------- Explore ----------------------------------------------------------------------
dim(brca_mutation$dat) # [1] 16806   977
head(brca_mutation$dat)[1:5,1:5]
