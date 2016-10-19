setwd("/media/ducdo/UUI1/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
load("data_Saved_R_Objects/brca_df.rda")
require(limma)

# create status object 
status = c(rep("Normal", 79), rep("Tumor", 457))
# create model matrix
model_matrix = model.matrix(~ as.factor(status))

## ------------- DEA on lncRNA ---------------------------------------------------------------------------
dim(brca_lncRNA_df) # [1] 4828  536

brca_lncRNA_fit_limma = lmFit(brca_lncRNA_df, brca_lncRNA_model_matrix)
brca_lncRNA_ebayes_limma = eBayes(brca_lncRNA_fit_limma)
brca_lncRNA_limma_all= topTable(brca_lncRNA_ebayes_limma, number = dim(brca_lncRNA_df)[1])

sum(brca_lncRNA_limma_all$adj.P.Val < 0.05) # 3148
head(brca_lncRNA_limma_all,3)

# check how many genes are up-regulated, how are down-regulated
brca_lncRNA_limma_all_subset = subset(brca_lncRNA_limma_all, brca_lncRNA_limma_all$adj.P.Val < 0.05)
dim(brca_lncRNA_limma_all_subset) # [1] 3148    6
sum(brca_lncRNA_limma_all_subset$logFC < 0) # 2023 lncRNA genes are down-regulated
sum(brca_lncRNA_limma_all_subset$logFC > 0) # 1125 lncRNA genes are up-regulated 

de_brca_lncRNA_df = brca_lncRNA_df[rownames(brca_lncRNA_limma_all_subset),]


## ------------- DEA on mRNA ---------------------------------------------------------------------------
dim(brca_mRNA_df) # [1] 17613 536

brca_mRNA_fit_limma = lmFit(brca_mRNA_df, model_matrix)
brca_mRNA_ebayes_limma = eBayes(brca_mRNA_fit_limma)
brca_mRNA_limma_all= topTable(brca_mRNA_ebayes_limma, number = dim(brca_mRNA_df)[1])

sum(brca_mRNA_limma_all$adj.P.Val < 0.05) # 13220
head(brca_mRNA_limma_all,3)

# check how many genes are up-regulated, how are down-regulated
brca_mRNA_limma_all_subset = subset(brca_mRNA_limma_all, brca_mRNA_limma_all$adj.P.Val < 0.05)
dim(brca_mRNA_limma_all_subset) # [1] 13220    6
sum(brca_mRNA_limma_all_subset$logFC < 0) # 6181 lncRNA genes are down-regulated
sum(brca_mRNA_limma_all_subset$logFC > 0) # 7039 lncRNA genes are up-regulated 

de_brca_mRNA_df = brca_mRNA_df[rownames(brca_mRNA_limma_all_subset),]

## ------------- DEA on miRNA --------------------------------------------------------#-------------------

dim(brca_miRNA_df) # [1] 343 536

brca_miRNA_fit_limma = lmFit(brca_miRNA_df, model_matrix)
brca_miRNA_ebayes_limma = eBayes(brca_miRNA_fit_limma)
brca_miRNA_limma_all= topTable(brca_miRNA_ebayes_limma, number = dim(brca_miRNA_df)[1])

sum(brca_miRNA_limma_all$adj.P.Val < 0.05) # 262
head(brca_miRNA_limma_all,3)

# check how many genes are up-regulated, how are down-regulated
brca_miRNA_limma_all_subset = subset(brca_miRNA_limma_all, brca_miRNA_limma_all$adj.P.Val < 0.05)
dim(brca_miRNA_limma_all_subset) # [1] 262    6
sum(brca_miRNA_limma_all_subset$logFC < 0) # 121 lncRNA genes are down-regulated
sum(brca_miRNA_limma_all_subset$logFC > 0) # 141 lncRNA genes are up-regulated 

de_brca_miRNA_df = brca_miRNA_df[rownames(brca_miRNA_limma_all_subset),]

## saving
save(brca_lncRNA_limma_all, brca_mRNA_limma_all_subset, brca_miRNA_limma_all_subset,
     de_brca_lncRNA_df, de_brca_mRNA_df, de_brca_miRNA_df,
     file = "Saved_R_Objects/differential_analysis/de_lncRNA_mRNA_miRNA.rda")

##############

load("data_Saved_R_Objects/differential_analysis/de_lncRNA_mRNA_miRNA.rda")
head(brca_miRNA_limma_all_subset)

"hsa-mir-22" %in% rownames(brca_miRNA_limma_all_subset) # true

brca_miRNA_limma_all_subset["hsa-mir-22",]











s