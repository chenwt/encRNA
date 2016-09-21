setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")

load("data_Saved_R_Objects/brca_df.rda")
# load("Saved_R_Objects/corr_matrices/corr_matrices.rda")
# load("Saved_R_Objects/corr_matrices/normal_tumor_corr_matrices.rda")
load("data_Saved_R_Objects/corr_matrices/all_corr_matrices_with_names.rda")
# load("data_Saved_R_Objects/corr_matrices/tumor_normal_lncRNA_mRNA_pair9999.rda")
# load("data_Saved_R_Objects/corr_matrices/sensitivity_matricies_9999.rda")
# load("Saved_R_Objects/corr_matrices/sensitivity_matrix.rda")
load("data_Saved_R_Objects/corr_matrices/normal_encRNA_850356.rda")
load("data_Saved_R_Objects/corr_matrices/tumor_encRNA_850356.rda")
load("data_Saved_R_Objects/corr_matrices/normal_sensitivity_full.rda")
load("data_Saved_R_Objects/corr_matrices/tumor_sensivitivity_full.rda")
load("data_Saved_R_Objects/corr_matrices/normal_tumor_lncRNA_mRNA_pairs.rda")
require(ppcor); require(rlist);
require(foreach); require(doParallel);
require(ggplot2)

# -------------------- Dimension check --------------------------------------------------
dim(lncRNA_mRNA_corr_matrix) # [1]  4828 17613
dim(miRNA_lncRNA_corr_matrix) #[1]  343 4828
dim(miRNA_mRNA_corr_matrix) #[1]   343 17613

dim(normal_sensitivity_full) # 850356    343


# ------------------- Visualization ----------------------------------------------------

# plot density of all correaltion between lncRNA and mRNA
plot(density(tumor_lncRNA_mRNA_corr_matrix), 
     col = "red", lwd = 3, ylim = c(0,4), 
     main = "mRNA - lncRNA correlation (normal and tumor)")
lines(density(normal_lncRNA_mRNA_corr_matrix), col = "green", lwd = 3)
legend(x = 'topright', pch = 15,
       col = c("green","red"),
       legend = c("normal", "tumor"))



# ------------------- Get lncRNA mRNA pair ----------------------------------------------
# DO NOT NEED TO RUN THIS CHUNK AGAIN
#select top pairs lncRNA-mRNA from normal dataset
quantile(normal_lncRNA_mRNA_corr_matrix, c(0.99)) # 0.710914
gc()
normal_lncRNA_mRNA_pairs = get_lncRNA_mRNA_pairs(normal_lncRNA_mRNA_corr_matrix,
                                                     quantile(normal_lncRNA_mRNA_corr_matrix, c(0.99)))

# select top pairs lncRNA-mRNA from tumor dataset
quantile(tumor_lncRNA_mRNA_corr_matrix, c(0.99)) # 0.3476711; 
gc()
tumor_lncRNA_mRNA_pairs = get_lncRNA_mRNA_pairs(tumor_lncRNA_mRNA_corr_matrix,
                                                quantile(tumor_lncRNA_mRNA_corr_matrix, c(0.99)))
dim(tumor_lncRNA_mRNA_pairs)
save(normal_lncRNA_mRNA_pairs, tumor_lncRNA_mRNA_pairs, file = "data_Saved_R_Objects/corr_matrices/normal_tumor_lncRNA_mRNA_pairs.rda")

# density of top 1% mRNA - lncRNA correlation
plot(density(tumor_lncRNA_mRNA_pairs$corr), 
     col = "red", lwd = 3, ylim = c(0,16), 
     xlab = "Correlations from top 850,356 lncRNA-mRNA pairs",
     main = "Top 1% mRNA - lncRNA correlation (normal and tumor)")
lines(density(normal_lncRNA_mRNA_pairs$corr), col = "green", lwd = 3)

legend(x = 'topright', pch = 15,
       col = c("green","red"),
       legend = c("normal (cut-off = 0.71)", "tumor (cut-off = 0.35)"))


# -------- check overlap between normal and tumor lncRNA - mRNA pair ---------------------

# question: how many highly correlated lncRNA-mRNA pairs are in both normal and tumor datasets?
normal_lncRNA_mRNA_pairs_vector = paste(normal_lncRNA_mRNA_pairs$lncRNA, normal_lncRNA_mRNA_pairs$mRNA, 
                                        sep = "-")
tumor_lncRNA_mRNA_pairs_vector = paste(tumor_lncRNA_mRNA_pairs$lncRNA, tumor_lncRNA_mRNA_pairs$mRNA, 
                                        sep = "-")
length(normal_lncRNA_mRNA_pairs_vector) # 850,356
length(intersect(normal_lncRNA_mRNA_pairs_vector,tumor_lncRNA_mRNA_pairs_vector)) # 92,725
# which is 10.9 % of 850356 pairs


# --------------- build sensitivity matrix from normal cells------------------------------

# get only subset expression data from normal cells
# normal_indices = 1:79
# mRNA_normal = brca_mRNA_df[,normal_indices];
# miRNA_normal = brca_miRNA_df[,normal_indices];
# lncRNA_normal = brca_lncRNA_df[,normal_indices]
# 
# # get sensitivity matrix for normal samples
# normal_pairs = normal_lncRNA_mRNA_pairs
# 
# num_cluster = 8
# cl <- makePSOCKcluster(num_cluster)
# registerDoParallel(cl)
# 
# ptm_normal = proc.time()
# normal_sensitivity_full = foreach(i = 1:nrow(normal_pairs), .combine = rbind) %dopar% {
#   require(ppcor)
#   the_pair = normal_pairs[i,]
#   mRNA_vector = mRNA_normal[the_pair$mRNA_index,]
#   lncRNA_vector = lncRNA_normal[the_pair$lncRNA_index,]
#   sensitivity_corr_vector = c()
#   for (j in 1:nrow(miRNA_normal)){
#     partial_corr = pcor.test(mRNA_vector, lncRNA_vector, miRNA_normal[j,])$estimate
#     sensitivity_corr_vector = append(sensitivity_corr_vector, the_pair$corr - partial_corr)
#   }
#   return(sensitivity_corr_vector)
# }
# ptm_normal = proc.time() - ptm_normal;
# ptm_normal
# stopCluster(cl)
# 
# sprintf("Running time with %g cores: %g", num_cluster ,ptm_normal[3])
# dim(normal_sensitivity_full)
# save(normal_sensitivity_full, file = "normal_sensivitivity_testFull.rda")


# dim(test); View(test)

# dim(normal_sensivitivity_test)
# 
# normal_lncRNA_mRNA_names = paste(normal_pairs$lncRNA, normal_pairs$mRNA, sep = "-")
# rownames(normal_sensivitivity_test) = normal_lncRNA_mRNA_names
# colnames(normal_sensivitivity_test) = rownames(brca_miRNA_df)
#


# # --------------- build sensitivity matrix from cancer cells------------------------------
# # get only subset expression data from tumor cells
# tumor_indices = 80:ncol(brca_miRNA_df)
# mRNA_tumor = brca_mRNA_df[,tumor_indices];
# miRNA_tumor = brca_miRNA_df[,tumor_indices];
# lncRNA_tumor = brca_lncRNA_df[,tumor_indices]
# 
# # get sensitivity matrix for tumor samples
# tumor_pairs = tumor_lncRNA_mRNA_pairs
# 
# num_cluster = 8
# cl <- makePSOCKcluster(num_cluster)
# registerDoParallel(cl)
# 
# ptm_tumor = proc.time()
# tumor_sensitivity_full = foreach(i = 1:nrow(tumor_pairs), .combine = rbind) %dopar% {
#   require(ppcor)
#   the_pair = tumor_pairs[i,]
#   mRNA_vector = mRNA_tumor[the_pair$mRNA_index,]
#   lncRNA_vector = lncRNA_tumor[the_pair$lncRNA_index,]
#   sensitivity_corr_vector = c()
#   for (j in 1:nrow(miRNA_tumor)){
#     partial_corr = pcor.test(mRNA_vector, lncRNA_vector, miRNA_tumor[j,])$estimate
#     sensitivity_corr_vector = append(sensitivity_corr_vector, the_pair$corr - partial_corr)
#   }
#   return(sensitivity_corr_vector)
# }
# ptm_tumor = proc.time() - ptm_tumor;
# ptm_tumor
# stopCluster(cl)
# 
# sprintf("Running time with %g cores: %g", num_cluster ,ptm_tumor[3])
# dim(tumor_sensitivity_full)
# 
# tumor_lncRNA_mRNA_names = paste(tumor_pairs$lncRNA, tumor_pairs$mRNA, sep = "-")
# rownames(tumor_sensitivity_full) = tumor_lncRNA_mRNA_names
# colnames(tumor_sensitivity_full) = rownames(brca_miRNA_df)
# 
# save(tumor_sensitivity_full, file = "tumor_sensivitivity_full.rda")


# ---check triplets statisfying encRNA hypothesis and included in ------------------------
# ----putative binding information--------------------------------------------------------

## NORMAL SAMPLES
# now, get the summary matrix from the sensitivity matrix, with threshold 0.99 for sensitivity matrix
# DO NOT RUN THIS CHUNK AGAIN
dim(normal_sensitivity_full) # 850356    343  --> 291,672,108 total triplets (lncRNA - mRNA - miRNA)
normal_99quantile = quantile(as.vector(normal_sensitivity_full),c(0.99)); gc() # 0.2454465
ptm = proc.time()
normal_encRNA = get_encRNA(matrix = normal_sensitivity_full,
                           lncRNA_mRNA_corr = normal_lncRNA_mRNA_corr_matrix,
                           miRNA_mRNA_corr = normal_miRNA_mRNA_corr_matrix,
                           miRNA_lncRNA_corr = normal_miRNA_lncRNA_corr_matrix,
                           threshold = normal_99quantile) 
ptm = proc.time() - ptm; ptm #  1103.398 seconds ~ 18 mins

# save(normal_encRNA, file = "data_Saved_R_Objects/corr_matrices/normal_encRNA_850356.rda")
dim(normal_encRNA) # [1] 2916723     12 --> total 2,916,723 triplets
length(unique(normal_encRNA$mRNA)) # 6819
length(unique(normal_encRNA$miRNA)) # 54
length(unique(normal_encRNA$lncRNA)) # 947
length(unique(normal_encRNA$encRNA_pair)) # 552576
length(unique(normal_encRNA$mRNA_miRNA_pair)) # 58721
length(unique(normal_encRNA$lncRNA_miRNA_pair)) # 7683
# FROM HERE CAN RUN AGAIN

## TUMOR SAMPLES
dim(tumor_sensitivity_full) # 850356    343 --> good
tumor_99quantile = quantile(as.vector(tumor_sensitivity_full),c(0.99)); gc() # 0.08703023 --> quite low
ptm = proc.time()
tumor_encRNA = get_encRNA(matrix = tumor_sensitivity_full,
                          lncRNA_mRNA_corr = tumor_lncRNA_mRNA_corr_matrix,
                          miRNA_mRNA_corr = tumor_miRNA_mRNA_corr_matrix,
                          miRNA_lncRNA_corr = tumor_miRNA_lncRNA_corr_matrix,
                          threshold = tumor_99quantile) # corresponding to 99 percentile
gc()
ptm = proc.time() - ptm; ptm #  770.898 ~ 12.84 mins

save(tumor_encRNA, file = "data_Saved_R_Objects/corr_matrices/tumor_encRNA_850356.rda")


length(unique(tumor_encRNA$mRNA)) 
length(unique(tumor_encRNA$lncRNA)) 
length(unique(tumor_encRNA$miRNA)) 
length(unique(tumor_encRNA$encRNA_pair)) 
length(unique(tumor_encRNA$mRNA_miRNA_pair)) 
length(unique(tumor_encRNA$lncRNA_miRNA_pair)) 


#### question: test the overlap between encRNA of tumor and normal 
length(intersect(unique(normal_encRNA$mRNA), unique(tumor_encRNA$mRNA)))
length(intersect(unique(normal_encRNA$lncRNA), unique(tumor_encRNA$lncRNA)))
length(intersect(unique(normal_encRNA$miRNA), unique(tumor_encRNA$miRNA)))
length(intersect(unique(normal_encRNA$encRNA_pair), unique(tumor_encRNA$encRNA_pair)))
length(intersect(unique(normal_encRNA$mRNA_miRNA_pair), unique(tumor_encRNA$mRNA_miRNA_pair)))
length(intersect(unique(normal_encRNA$lncRNA_miRNA_pair), unique(tumor_encRNA$lncRNA_miRNA_pair)))
length(intersect(unique(normal_encRNA$encRNA_triple), unique(tumor_encRNA$encRNA_triple)))

## -------- Visualizeation the distribution of mRNA-miRNA's and lncRNA-miRNA's correlation in encRNA -----------
# aim: can be used to filter encRNA

# plot the sensitivity matrix
# plot density of all sensitivity values (291,672,108 data points)
plot(density(as.vector(normal_sensitivity_full)), 
     col = "green", lwd = 3, ylim = c(0,160), 
     main = "Sensitivity correlations from 850,356 encRNA pairs vs 343 miRNAs")
lines(density(as.vector(tumor_sensitivity_full)), col = "red", lwd = 3)
legend(x = 'topright', pch = 15,
       col = c("green","red"),
       legend = c("normal", "tumor"))

# for the top 1% sensitivity cut-off, plot all the remaining sensitivity
plot(density(normal_encRNA$sensitivity), 
     main = "Sensivity correlation of top 1% of 850,356 encRNA pairs vs 343 miRNAs",
     xlab = "Sensivity correlation of 2,916,723 encRNAs-miRNA pairs",
     xlim = c(0,1),
     ylim = c(0,31),
     col = "green", lwd = 3)
lines(density(tumor_encRNA$sensitivity), 
      col = "red", lwd = 3)
legend(x = 'topright', pch = 15,
       col = c("red","green"),
       legend = c("tumor", "normal"))


# plot lncRNA-miRNA and mRNA-miRNA correlation correlation of normal_encRNA and tumor_encRNA
par(mfrow = c(2,1))
plot(density(normal_encRNA$lncRNA_miRNA_corr), 
     main = "Normal samples - encRNA pairs correlation (S-based top 1% encRNA triplets)",
     xlab = "Correlation of 2,916,723 RNA-miRNA pairs",
     ylim = c(0,7),
     col = "orange", lwd = 3)
lines(density(normal_encRNA$mRNA_miRNA_corr), 
      col = "blue", lwd = 3)
legend(x = 'topleft', pch = 15,
       col = c("orange","blue"),
       legend = c("lncRNA-miRNA correlation", "mRNA-miRNA correlation"))

plot(density(tumor_encRNA$lncRNA_miRNA_corr), 
     main = "Tumor samples - encRNA pairs correlation (S-based top 1% encRNA triplets)",
     xlab = "Correlation of 2,916,723 RNA-miRNA pairs",
     ylim = c(0,7),
     col = "orange", lwd = 3)
lines(density(tumor_encRNA$mRNA_miRNA_corr), 
      col = "blue", lwd = 3)
legend(x = 'topleft', pch = 15,
       col = c("orange","blue"),
       legend = c("lncRNA-miRNA correlation", "mRNA-miRNA correlation"))

plot(density(tumor_encRNA$sensitivity), 
     col = "red", lwd = 3, ylim = c(0,16), 
     main = "Top 99 percentile mRNA - lncRNA correlation (normal and tumor)")
lines(density(normal_lncRNA_mRNA_pairs$corr), col = "green", lwd = 3)


# plot mRNA-lncRNA correlation correlation of normal_encRNA and tumor_encRNA
plot(density(normal_encRNA$lncRNA_mRNA_corr), 
     main = "mRNA-lncRNA correlation (S-based top 1% encRNA triplets)",
     xlab = "mRNA-lncRNA correlation of 2,916,723 encRNAs-miRNA pairs",
     xlim = c(0,1),
     ylim = c(0,31),
     col = "green", lwd = 3)
lines(density(tumor_encRNA$lncRNA_mRNA_corr), 
      col = "red", lwd = 3)
legend(x = 'topright', pch = 15,
       col = c("red","green"),
       legend = c("tumor", "normal"))


summary(normal_encRNA$mRNA_miRNA_corr)
summary(normal_encRNA$lncRNA_miRNA_corr)

#### question: how many triplets in normal_encRNA satistifying encRNA hypothesis?
#### answer: 347,072 for normal samples; 844,475 for tumor samples

# normal samples
length(which(normal_encRNA$lncRNA_miRNA_corr < 0 & normal_encRNA$mRNA_miRNA_corr < 0)) # 347,072
indices = which(normal_encRNA$lncRNA_miRNA_corr < 0 & normal_encRNA$mRNA_miRNA_corr < 0, arr.ind = TRUE)
normal_encRNA_goodCorr = normal_encRNA[indices,]; nrow(normal_encRNA_goodCorr) # 347,072

length(which(normal_encRNA$lncRNA_miRNA_corr < -0.5 & normal_encRNA$mRNA_miRNA_corr < -0.5)) # 347,072


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


#### question: out of all triplets, how many included in the putative binding information? 
#### ans: 21,553 for normal samples; 11,516 for tumor
# normal samples
ptm = proc.time()
normal_encRNA_sensitivity_bound = get_matched_enRNA_sensitivity_with_putative_binding(normal_encRNA)
ptm = proc.time() - ptm; ptm # 28 seconds
dim(normal_encRNA_sensitivity_bound) # 21,553   12
length(unique(normal_encRNA_sensitivity_bound$mRNA)) # 3563
length(unique(normal_encRNA_sensitivity_bound$miRNA)) # 17
length(unique(normal_encRNA_sensitivity_bound$lncRNA)) # 142
length(unique(normal_encRNA_sensitivity_bound$encRNA_pair))
#tumor samples
ptm = proc.time()
tumor_encRNA_sensitivity_bound = get_matched_enRNA_sensitivity_with_putative_binding(tumor_encRNA)
ptm = proc.time() - ptm; ptm # 18 seconds
dim(tumor_encRNA_sensitivity_bound) # [1] 11516    13
length(unique(tumor_encRNA_sensitivity_bound$mRNA)) # 3220
length(unique(tumor_encRNA_sensitivity_bound$miRNA))# 51
length(unique(tumor_encRNA_sensitivity_bound$lncRNA)) # 312
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
normal_encRNA_sensitivity_bound_goodCoor = normal_encRNA_sensitivity_bound[indices,]; nrow(normal_encRNA_sensitivity_bound_goodCoor) 
nrow(normal_encRNA_sensitivity_bound_goodCoor)
length(unique(normal_encRNA_sensitivity_bound_goodCoor$mRNA)) 
length(unique(normal_encRNA_sensitivity_bound_goodCoor$miRNA)) 
length(unique(normal_encRNA_sensitivity_bound_goodCoor$lncRNA)) 
length(unique(normal_encRNA_sensitivity_bound_goodCoor$encRNA_pair))
# tumor samples
indices = which(tumor_encRNA_sensitivity_bound$lncRNA_miRNA_corr < 0 & tumor_encRNA_sensitivity_bound$mRNA_miRNA_corr < 0, arr.ind = TRUE)
tumor_encRNA_sensitivity_bound_goodCoor = tumor_encRNA_sensitivity_bound[indices,]; nrow(tumor_encRNA_sensitivity_bound_goodCoor) 
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



# #### follow-up question: out of those encRNA triplets satisfying encRNA hypothesis, how many included in putattive bindining information?
# #### answer: 1955 triplets for normal pairs, 
# # normal pairs
# ptm = proc.time()
# normal_encRNA_sensitivity_goodCorr_bound = get_matched_enRNA_sensitivity_with_putative_binding(normal_encRNA_goodCorr)
# ptm = proc.time() - ptm; ptm # 11 seconds 
# dim(normal_encRNA_sensitivity_goodCorr_bound) # 1955   12 --> 1955 triplets included in miRcode and starBase
# length(unique(normal_encRNA_sensitivity_goodCorr_bound$mRNA)) 
# length(unique(normal_encRNA_sensitivity_goodCorr_bound$miRNA)) 
# length(unique(normal_encRNA_sensitivity_goodCorr_bound$lncRNA)) 
# length(unique(normal_encRNA_sensitivity_goodCorr_bound$encRNA_pair))
# 
# # cancer pairs
# ptm = proc.time()
# tumor_encRNA_sensitivity_goodCorr_bound = get_matched_enRNA_sensitivity_with_putative_binding(tumor_encRNA_goodCorr)
# ptm = proc.time() - ptm; ptm # 11 seconds 
# dim(tumor_encRNA_sensitivity_goodCorr_bound) 
# length(unique(tumor_encRNA_sensitivity_goodCorr_bound$mRNA)) 
# length(unique(tumor_encRNA_sensitivity_goodCorr_bound$miRNA)) 
# length(unique(tumor_encRNA_sensitivity_goodCorr_bound$lncRNA)) 
# length(unique(tumor_encRNA_sensitivity_goodCorr_bound$encRNA_pair))
# 
# # check overlap
# length(intersect(unique(normal_encRNA_sensitivity_goodCorr_bound$mRNA), unique(tumor_encRNA_sensitivity_goodCorr_bound$mRNA)))
# length(intersect(unique(normal_encRNA_sensitivity_goodCorr_bound$miRNA), unique(tumor_encRNA_sensitivity_goodCorr_bound$miRNA)))
# length(intersect(unique(normal_encRNA_sensitivity_goodCorr_bound$lncRNA), unique(tumor_encRNA_sensitivity_goodCorr_bound$lncRNA)))
# length(intersect(unique(normal_encRNA_sensitivity_goodCorr_bound$encRNA_pair), unique(tumor_encRNA_sensitivity_goodCorr_bound$encRNA_pair)))
# 
# 
# length(intersect(unique(normal_encRNA_sensitivity_goodCorr_bound$mRNA), unique(tumor_encRNA_sensitivity_goodCorr_bound$mRNA)))


# RECAP: the 1,955 triples satisfy 3 conditions:
# 1. top 99 percentile sensitivity matrix
# 2. following encRNA hypothesis (negative correlation lncRNA-miRNA, mRNA-miRNA, positive correalation lncRNA-mRNA)
# 3. the RNA-miRNA interaction is included in the putative binding databases



# Question: in the 1,955 triples, is there any mRNA-lnRNA included in the normal-tumor overlapped 92,725 
# mRNA-lncRNA pairs (in the beginning of the analysis)?
overlapped_mRNA_lncRNA_two_conditions = intersect(normal_lncRNA_mRNA_pairs_vector,tumor_lncRNA_mRNA_pairs_vector)
length(intersect(normal_encRNA_sensitivity_goodCorr_bound$encRNA_pair,overlapped_mRNA_lncRNA_two_conditions)) # 58 

# -------------- building mRNA-mRNA and lncRNA-lncRNA pairs ---------------------------------------------
load("data_Saved_R_Objects/corr_matrices/normal_tumor_lncRNA_mRNA_pairs.rda")
# idea: based on the normal_lncRNA_mRNA_pairs and tumor_lncRNA_mRNA_pairs --> build 4 correlation matricies of 
# mRNA-mRNA and lncRNA-lncRNA in tumor and normal cases
normal_lncRNA_mRNA_pairs
tumor_lncRNA_mRNA_pairs

normal_mRNA = unique(normal_lncRNA_mRNA_pairs$mRNA); length(normal_mRNA) # 12182
normal_lncRNA = unique(normal_lncRNA_mRNA_pairs$lncRNA); length(normal_lncRNA) # 2379
tumor_mRNA = unique(tumor_lncRNA_mRNA_pairs$mRNA); length(tumor_mRNA) # 16320
tumor_lncRNA = unique(tumor_lncRNA_mRNA_pairs$lncRNA); length(tumor_lncRNA) # 4575

normal_indices = 1:79
mRNA_normal = brca_mRNA_df[normal_mRNA,normal_indices]; dim(mRNA_normal) # [1] 12182    79
lncRNA_normal = brca_lncRNA_df[normal_lncRNA,normal_indices]; dim(lncRNA_normal) # 2379   79

ptm = proc.time()
normal_mRNA_mRNA_corr_matrix = cor(t(mRNA_normal)); gc()
normal_mRNA_mRNA_corr_matrix[lower.tri(normal_mRNA_mRNA_corr_matrix,diag = TRUE)] <- NA; gc()
normal_lncRNA_lncRNA_corr_matrix = cor(t(lncRNA_normal)); gc()
normal_lncRNA_lncRNA_corr_matrix[lower.tri(normal_lncRNA_lncRNA_corr_matrix,diag = TRUE)] <- NA; gc()
ptm = proc.time() - ptm; ptm # 


originalCorrMatrix = normal_mRNA_mRNA_corr_matrix
corrThreshold = 0

get_mRNA_mRNA_pairs = function(originalCorrMatrix, corrThreshold){
  correlation_pairs = which(!is.na(originalCorrMatrix), arr.ind = TRUE); gc()
  mRNA1 = rownames(originalCorrMatrix)[correlation_pairs[,1]]; gc()
  mRNA2 = colnames(originalCorrMatrix)[correlation_pairs[,2]]; gc()
  dataframe = as.data.frame(cbind(correlation_pairs, mRNA1, mRNA2))
  colnames(dataframe) = c("mRNA_index1", "mRNA_index2","mRNA1", "mRNA2")
  dataframe$mRNA_index1 = as.numeric(as.character(dataframe$mRNA_index1))
  dataframe$mRNA_index2 = as.numeric(as.character(dataframe$mRNA_index2))
  # correlation_vector = normal_lncRNA_mRNA_corr_matrix[dataframe$lncRNA_index,dataframe$mRNA_index]
  # dataframe = cbind(dataframe, originalCorrMatrix[which(originalCorrMatrix > corrThreshold)])
  # colnames(dataframe) = c("mRNA_index1", "mRNA_index2","mRNA1", "mRNA2","corr")
  # dataframe = dataframe[order(-dataframe$corr),]
  # return(dataframe)
}

