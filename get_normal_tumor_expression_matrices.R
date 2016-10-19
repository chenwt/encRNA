setwd("/media/ducdo/UUI1/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
load("data_Saved_R_Objects/brca_df.rda")

normal_indices = 1:79
tumor_indices = 90:ncol(brca_lncRNA_df)

mRNA_normal = brca_mRNA_df[,normal_indices];
miRNA_normal = brca_miRNA_df[,normal_indices];
lncRNA_normal = brca_lncRNA_df[,normal_indices];

mRNA_tumor = brca_mRNA_df[,tumor_indices];
miRNA_tumor = brca_miRNA_df[,tumor_indices];
lncRNA_tumor = brca_lncRNA_df[,tumor_indices]
