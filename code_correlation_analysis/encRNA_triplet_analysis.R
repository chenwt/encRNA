setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")

load("data_Saved_R_Objects/corr_matrices/normal_encRNA_850356.rda")
load("data_Saved_R_Objects/corr_matrices/tumor_encRNA_850356.rda")
load("data_Saved_R_Objects/miRNA_target/brca_putative_encRNA.rda"); gc()

########################################################################################
####       Statistics on putative triplet                       ########################
########################################################################################

dim(brca_putative_encRNA) # 60,862,939        3
length(unique(brca_putative_encRNA$mRNA)); gc()
length(unique(brca_putative_encRNA$lncRNA)); gc()
length(unique(brca_putative_encRNA$miRNA)); gc()

# ---check triplets statisfying encRNA hypothesis and included in ------------------------
# ----putative binding information--------------------------------------------------------

# Recall that the normal_encRNA and tumor_enRNA are obtained from the top 1% sensitivity of 
# top lncRNA-mRNA pairs with all 343 miRNAs

#### question: how many triplets in normal_encRNA satistifying encRNA hypothesis?
#### answer: 347,072 for normal samples; 844,475 for tumor samples

# normal samples
length(which(normal_encRNA$lncRNA_miRNA_corr < 0 & normal_encRNA$mRNA_miRNA_corr < 0)) # 347,072
indices = which(normal_encRNA$lncRNA_miRNA_corr < 0 & normal_encRNA$mRNA_miRNA_corr < 0, arr.ind = TRUE)
normal_encRNA_goodCorr = normal_encRNA[indices,]; nrow(normal_encRNA_goodCorr) # 347,072

# length(which(normal_encRNA$lncRNA_miRNA_corr < -0.5 & normal_encRNA$mRNA_miRNA_corr < -0.5)) # 347,072

length(unique(normal_encRNA_goodCorr$mRNA)) 
length(unique(normal_encRNA_goodCorr$miRNA)) 
length(unique(normal_encRNA_goodCorr$lncRNA)) 
length(unique(normal_encRNA_goodCorr$encRNA_pair))

# tumor samples
length(which(tumor_encRNA$lncRNA_miRNA_corr < 0 & tumor_encRNA$mRNA_miRNA_corr < 0)) 
indices = which(tumor_encRNA$lncRNA_miRNA_corr < 0 & tumor_encRNA$mRNA_miRNA_corr < 0, arr.ind = TRUE)
tumor_encRNA_goodCorr = tumor_encRNA[indices,]; nrow(tumor_encRNA_goodCorr) 

length(which(tumor_encRNA$lncRNA_miRNA_corr < -0.5 & tumor_encRNA$mRNA_miRNA_corr < -0.5)) 
length(unique(tumor_encRNA_goodCorr$mRNA)) 
length(unique(tumor_encRNA_goodCorr$miRNA)) 
length(unique(tumor_encRNA_goodCorr$lncRNA)) 
length(unique(tumor_encRNA_goodCorr$encRNA_pair))

# check overlap
length(intersect(unique(normal_encRNA_goodCorr$encRNA_triple), unique(tumor_encRNA_goodCorr$encRNA_triple)))
length(intersect(unique(normal_encRNA_goodCorr$mRNA), unique(tumor_encRNA_goodCorr$mRNA)))
length(intersect(unique(normal_encRNA_goodCorr$miRNA), unique(tumor_encRNA_goodCorr$miRNA)))
length(intersect(unique(normal_encRNA_goodCorr$lncRNA), unique(tumor_encRNA_goodCorr$lncRNA)))
length(intersect(unique(normal_encRNA_goodCorr$encRNA_pair), unique(tumor_encRNA_goodCorr$encRNA_pair)))

#################################################################################################
#### question: out of all triplets, how many included in the putative binding information? ######
#### ans: 21,553 for normal samples; 11,516 for tumor                                      ######  
#################################################################################################

time = proc.time()
normal_encRNA_sensitivity_bound = get_putative_encRNA_interaction(normal_encRNA)
time = proc.time() - time; print(time) # 40s 

time = proc.time()
tumor_encRNA_sensitivity_bound = get_putative_encRNA_interaction(tumor_encRNA)
time = proc.time() - time; print(time) # 120s

save(normal_encRNA_sensitivity_bound, tumor_encRNA_sensitivity_bound, 
     file = "data_Saved_R_Objects/miRNA_target/predicted_normal_tumor_sensitivity_binding.rda")

dim(normal_encRNA_sensitivity_bound) # 21,553   12
length(unique(normal_encRNA_sensitivity_bound$mRNA)) 
length(unique(normal_encRNA_sensitivity_bound$miRNA)) 
length(unique(normal_encRNA_sensitivity_bound$lncRNA)) 
length(unique(normal_encRNA_sensitivity_bound$encRNA_pair))

dim(tumor_encRNA_sensitivity_bound) # [1] 11516    13
length(unique(tumor_encRNA_sensitivity_bound$mRNA)) 
length(unique(tumor_encRNA_sensitivity_bound$miRNA))
length(unique(tumor_encRNA_sensitivity_bound$lncRNA)) 
length(unique(tumor_encRNA_sensitivity_bound$encRNA_pair))


# check overlap
length(intersect(unique(normal_encRNA_sensitivity_bound$encRNA_triple), unique(tumor_encRNA_sensitivity_bound$encRNA_triple)))

length(intersect(unique(normal_encRNA_sensitivity_bound$mRNA), unique(tumor_encRNA_sensitivity_bound$mRNA)))
length(intersect(unique(normal_encRNA_sensitivity_bound$miRNA), unique(tumor_encRNA_sensitivity_bound$miRNA)))
length(intersect(unique(normal_encRNA_sensitivity_bound$lncRNA), unique(tumor_encRNA_sensitivity_bound$lncRNA)))
length(intersect(unique(normal_encRNA_sensitivity_bound$encRNA_pair), unique(tumor_encRNA_sensitivity_bound$encRNA_pair)))


#### follow-up question: out of those encRNA whose interactions included in the database, how many satisfying the encRNA hypothesis?
# normal samples
indices = which(normal_encRNA_sensitivity_bound$lncRNA_miRNA_corr < 0 & normal_encRNA_sensitivity_bound$mRNA_miRNA_corr < 0, arr.ind = TRUE)
normal_encRNA_sensitivity_bound_goodCoor = normal_encRNA_sensitivity_bound[indices,]; 
nrow(normal_encRNA_sensitivity_bound_goodCoor)
length(unique(normal_encRNA_sensitivity_bound_goodCoor$mRNA)) 
length(unique(normal_encRNA_sensitivity_bound_goodCoor$miRNA)) 
length(unique(normal_encRNA_sensitivity_bound_goodCoor$lncRNA)) 
length(unique(normal_encRNA_sensitivity_bound_goodCoor$encRNA_pair))
# tumor samples
indices = which(tumor_encRNA_sensitivity_bound$lncRNA_miRNA_corr < 0 & tumor_encRNA_sensitivity_bound$mRNA_miRNA_corr < 0, arr.ind = TRUE)
tumor_encRNA_sensitivity_bound_goodCoor = tumor_encRNA_sensitivity_bound[indices,]; 
nrow(tumor_encRNA_sensitivity_bound_goodCoor)
length(unique(tumor_encRNA_sensitivity_bound_goodCoor$mRNA)) 
length(unique(tumor_encRNA_sensitivity_bound_goodCoor$miRNA)) 
length(unique(tumor_encRNA_sensitivity_bound_goodCoor$lncRNA)) 
length(unique(tumor_encRNA_sensitivity_bound_goodCoor$encRNA_pair))

# check overlap
length(intersect(unique(normal_encRNA_sensitivity_bound_goodCoor$encRNA_triple), unique(tumor_encRNA_sensitivity_bound_goodCoor$encRNA_triple)))
length(intersect(unique(normal_encRNA_sensitivity_bound_goodCoor$mRNA), unique(tumor_encRNA_sensitivity_bound_goodCoor$mRNA)))
length(intersect(unique(normal_encRNA_sensitivity_bound_goodCoor$miRNA), unique(tumor_encRNA_sensitivity_bound_goodCoor$miRNA)))
length(intersect(unique(normal_encRNA_sensitivity_bound_goodCoor$lncRNA), unique(tumor_encRNA_sensitivity_bound_goodCoor$lncRNA)))
length(intersect(unique(normal_encRNA_sensitivity_bound_goodCoor$encRNA_pair), unique(tumor_encRNA_sensitivity_bound_goodCoor$encRNA_pair)))

save(normal_encRNA_sensitivity_bound_goodCoor, tumor_encRNA_sensitivity_bound_goodCoor, 
     file = "data_Saved_R_Objects/miRNA_target/predicted_normal_tumor_goodCorr.rda")

########################################################################################
####       Select all sensitivity > 0                           ########################
########################################################################################

setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")
load("data_Saved_R_Objects/corr_matrices/normal_sensitivity_full.rda"); gc()
