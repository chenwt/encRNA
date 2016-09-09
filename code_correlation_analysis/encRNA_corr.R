setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
load("Saved_R_Objects/brca_df.rda")
# load("Saved_R_Objects/corr_matrices/corr_matrices.rda")
# load("Saved_R_Objects/corr_matrices/normal_tumor_corr_matrices.rda")
load("Saved_R_Objects/corr_matrices/all_corr_matrices_with_names.rda")
load("Saved_R_Objects/corr_matrices/sensitivity_matrix.rda")

# test: results must be different - checked!
lncRNA_mRNA_corr_matrix[1,1]
tumor_lncRNA_mRNA_corr_matrix[1,1]
normal_miRNA_lncRNA_corr_matrix[1,1]

dim(lncRNA_mRNA_corr_matrix) # [1]  4828 17613
dim(miRNA_lncRNA_corr_matrix) #[1]  343 4828
dim(miRNA_mRNA_corr_matrix) #[1]   343 17613

rownames(lncRNA_mRNA_corr_matrix) = rownames(brca_lncRNA_df)
colnames(lncRNA_mRNA_corr_matrix) = rownames(brca_mRNA_df)

rownames(miRNA_mRNA_corr_matrix) = rownames(brca_miRNA_df)
colnames(miRNA_mRNA_corr_matrix) = rownames(brca_mRNA_df)

rownames(miRNA_lncRNA_corr_matrix) = rownames(brca_miRNA_df)
colnames(miRNA_lncRNA_corr_matrix) = rownames(brca_lncRNA_df)


rownames(normal_lncRNA_mRNA_corr_matrix) = rownames(brca_lncRNA_df)
colnames(normal_lncRNA_mRNA_corr_matrix) = rownames(brca_mRNA_df)
rownames(tumor_lncRNA_mRNA_corr_matrix) = rownames(brca_lncRNA_df)
colnames(tumor_lncRNA_mRNA_corr_matrix) = rownames(brca_mRNA_df)

rownames(normal_miRNA_mRNA_corr_matrix) = rownames(brca_miRNA_df)
colnames(normal_miRNA_mRNA_corr_matrix) = rownames(brca_mRNA_df)
rownames(tumor_miRNA_mRNA_corr_matrix) = rownames(brca_miRNA_df)
colnames(tumor_miRNA_mRNA_corr_matrix) = rownames(brca_mRNA_df)

rownames(normal_miRNA_lncRNA_corr_matrix) = rownames(brca_miRNA_df)
colnames(normal_miRNA_lncRNA_corr_matrix) = rownames(brca_lncRNA_df)
rownames(tumor_miRNA_lncRNA_corr_matrix) = rownames(brca_miRNA_df)
colnames(tumor_miRNA_lncRNA_corr_matrix) = rownames(brca_lncRNA_df)

save(lncRNA_mRNA_corr_matrix, miRNA_mRNA_corr_matrix, miRNA_lncRNA_corr_matrix,
     normal_lncRNA_mRNA_corr_matrix, tumor_lncRNA_mRNA_corr_matrix,
     normal_miRNA_mRNA_corr_matrix, tumor_miRNA_mRNA_corr_matrix,
     normal_miRNA_lncRNA_corr_matrix, tumor_miRNA_lncRNA_corr_matrix,
     file = "Saved_R_Objects/corr_matrices/all_corr_matrices_with_names.rda")

## ------- both normal and tumor -------------------------------------------------------

## lncRNA_mRNA 
# convert to vector
lncRNA_mRNA_corr_vector = as.vector(lncRNA_mRNA_corr_matrix)
summary(lncRNA_mRNA_corr_vector)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.80060 -0.06800  0.01551  0.02119  0.10570  0.99930 
quantile(lncRNA_mRNA_corr_vector, c(0.3,0.5,0.7,0.9,0.99))
# 30%         50%         70%         90%         99% 
# -0.04921599  0.01550982  0.08454999  0.19998947  0.40084424 

# how many mRNA-lncRNA pairs whose correlation is higher than: 0.6, 0.7, 0.8, 0.9
nrow(which(lncRNA_mRNA_corr_matrix > 0.6, arr.ind = TRUE)) # 76061 pairs
nrow(which(lncRNA_mRNA_corr_matrix > 0.7, arr.ind = TRUE)) # 19589 pairs
nrow(which(lncRNA_mRNA_corr_matrix > 0.8, arr.ind = TRUE)) # 3617 pairs
nrow(which(lncRNA_mRNA_corr_matrix > 0.9, arr.ind = TRUE)) # 502 pairs

## miRNA_mRNA

miRNA_mRNA_corr_vector = as.vector(miRNA_mRNA_corr_matrix)
quantile(miRNA_mRNA_corr_vector, c(0.3,0.5,0.7,0.9,0.99))
# 30%          50%          70%          90%          99% 
# -0.071368945 -0.001705558  0.069046615  0.189977536  0.410233097 
nrow(which(miRNA_mRNA_corr_matrix > 0.6, arr.ind = TRUE)) # 4394 pairs
nrow(which(miRNA_mRNA_corr_matrix > 0.7, arr.ind = TRUE)) # 502 pairs
nrow(which(miRNA_mRNA_corr_matrix > 0.8, arr.ind = TRUE)) # 64 pairs
nrow(which(miRNA_mRNA_corr_matrix > 0.9, arr.ind = TRUE)) # 1 pair

which(miRNA_mRNA_corr_matrix > 0.9, arr.ind = TRUE)

### testing 

# get the column (correspond to mRNA) satisfy the condition:
mRNA_indices_0.9 = unname(unique(which(lncRNA_mRNA_corr_matrix > 0.9, arr.ind = TRUE)[,2]))
mRNA_indices_0.8 = unname(unique(which(lncRNA_mRNA_corr_matrix > 0.8, arr.ind = TRUE)[,2]))
mRNA_indices_0.7 = unname(unique(which(lncRNA_mRNA_corr_matrix > 0.7, arr.ind = TRUE)[,2]))
mRNA_indices_0.6 = unname(unique(which(lncRNA_mRNA_corr_matrix > 0.6, arr.ind = TRUE)[,2]))
mRNA_indices_0.2 = unname(unique(which(lncRNA_mRNA_corr_matrix > 0.2, arr.ind = TRUE)[,2]))


normal_mRNA_indices_0.9 = unname(unique(which(normal_lncRNA_mRNA_corr_matrix > 0.9, arr.ind = TRUE)[,2]))
normal_mRNA_indices_0.8 = unname(unique(which(normal_lncRNA_mRNA_corr_matrix > 0.8, arr.ind = TRUE)[,2]))
normal_mRNA_indices_0.7 = unname(unique(which(normal_lncRNA_mRNA_corr_matrix > 0.7, arr.ind = TRUE)[,2]))
normal_mRNA_indices_0.6 = unname(unique(which(normal_lncRNA_mRNA_corr_matrix > 0.6, arr.ind = TRUE)[,2]))
normal_mRNA_indices_0.2 = unname(unique(which(normal_lncRNA_mRNA_corr_matrix > 0.2, arr.ind = TRUE)[,2]))


tumor_mRNA_indices_0.9 = unname(unique(which(tumor_lncRNA_mRNA_corr_matrix > 0.9, arr.ind = TRUE)[,2]))
tumor_mRNA_indices_0.8 = unname(unique(which(tumor_lncRNA_mRNA_corr_matrix > 0.8, arr.ind = TRUE)[,2]))
tumor_mRNA_indices_0.7 = unname(unique(which(tumor_lncRNA_mRNA_corr_matrix > 0.7, arr.ind = TRUE)[,2]))
tumor_mRNA_indices_0.6 = unname(unique(which(tumor_lncRNA_mRNA_corr_matrix > 0.6, arr.ind = TRUE)[,2]))
tumor_mRNA_indices_0.2 = unname(unique(which(tumor_lncRNA_mRNA_corr_matrix > 0.2, arr.ind = TRUE)[,2]))


nrow(which(normal_lncRNA_mRNA_corr_matrix > 0.6, arr.ind = TRUE)) # 3221969 pairs
nrow(which(normal_lncRNA_mRNA_corr_matrix > 0.7, arr.ind = TRUE)) # 996604 pairs
nrow(which(normal_lncRNA_mRNA_corr_matrix > 0.8, arr.ind = TRUE)) # 156770 pairs
nrow(which(normal_lncRNA_mRNA_corr_matrix > 0.9, arr.ind = TRUE)) # 502 pairs

nrow(which(tumor_lncRNA_mRNA_corr_matrix > 0.6, arr.ind = TRUE)) # 28094 pairs
nrow(which(tumor_lncRNA_mRNA_corr_matrix > 0.7, arr.ind = TRUE)) # 7857 pairs
nrow(which(tumor_lncRNA_mRNA_corr_matrix > 0.8, arr.ind = TRUE)) # 1735 pairs
nrow(which(tumor_lncRNA_mRNA_corr_matrix > 0.9, arr.ind = TRUE)) # 406 pairs


# all miRNA-mRNA correlation
plot(density(miRNA_mRNA_corr_matrix), 
     main = "miRNA-mRNA correlation", 
     ylim = c(0,4),
     lwd = 3, cex.main = 0.9, col = "blue")
lines(density(normal_miRNA_mRNA_corr_matrix),  col = "black", lwd = 3)
lines(density(tumor_miRNA_mRNA_corr_matrix),  col = "red", lwd = 3)
legend(x = 'topright', pch = 15,
       col = c("blue", "black", "red"),
       legend = c("both normal and tumor", "normal only", "tumor only"))

# normal cells, contrainst on all_mRNA-lncRNA correaltion
subset_normal_miRNA_mRNA_all_corr_matrix_0.9 = normal_miRNA_mRNA_corr_matrix[,mRNA_indices_0.9]
subset_normal_miRNA_mRNA_all_corr_matrix_0.8 = normal_miRNA_mRNA_corr_matrix[,mRNA_indices_0.8]
subset_normal_miRNA_mRNA_all_corr_matrix_0.7 = normal_miRNA_mRNA_corr_matrix[,mRNA_indices_0.7]
subset_normal_miRNA_mRNA_all_corr_matrix_0.6 = normal_miRNA_mRNA_corr_matrix[,mRNA_indices_0.6]


plot(density(normal_miRNA_mRNA_corr_matrix), 
     main = "normal_miRNA-mRNA correlation (mRNAs selected based on all_lncRNA_mRNA correlation", 
     lwd = 3, cex.main = 0.9)
#lines(density(subset_normal_miRNA_mRNA_corr_matrix_0.2), col = "gray", lwd = 3)
lines(density(subset_normal_miRNA_mRNA_all_corr_matrix_0.6), col = "green", lwd = 3)
lines(density(subset_normal_miRNA_mRNA_all_corr_matrix_0.7), col = "orange", lwd = 3)
lines(density(subset_normal_miRNA_mRNA_all_corr_matrix_0.8), col = "blue", lwd = 3)
lines(density(subset_normal_miRNA_mRNA_all_corr_matrix_0.9),  col = "red", lwd = 3)

legend(x = 'topright', pch = 15,
       col = c("black","green", "orange", "blue","red"),
       legend = c("no constraint","r>0.6", "r>0.7", "r>0.8", "r>0.9"))


# normal cells, contrainst on normal_mRNA-lncRNA correaltion
subset_normal_miRNA_mRNA_normal_corr_matrix_0.9 = normal_miRNA_mRNA_corr_matrix[,normal_mRNA_indices_0.9]
subset_normal_miRNA_mRNA_normal_corr_matrix_0.8 = normal_miRNA_mRNA_corr_matrix[,normal_mRNA_indices_0.8]
subset_normal_miRNA_mRNA_normal_corr_matrix_0.7 = normal_miRNA_mRNA_corr_matrix[,normal_mRNA_indices_0.7]
subset_normal_miRNA_mRNA_normal_corr_matrix_0.6 = normal_miRNA_mRNA_corr_matrix[,normal_mRNA_indices_0.6]


plot(density(normal_miRNA_mRNA_corr_matrix), 
     main = "normal_miRNA-mRNA correlation (mRNAs selected based on normal_lncRNA_mRNA correlation)", 
     lwd = 3, cex.main = 0.9)
#lines(density(subset_normal_miRNA_mRNA_corr_matrix_0.2), col = "gray", lwd = 3)
lines(density(subset_normal_miRNA_mRNA_normal_corr_matrix_0.6), col = "green", lwd = 3)
lines(density(subset_normal_miRNA_mRNA_normal_corr_matrix_0.7), col = "orange", lwd = 3)
lines(density(subset_normal_miRNA_mRNA_normal_corr_matrix_0.8), col = "blue", lwd = 3)
lines(density(subset_normal_miRNA_mRNA_normal_corr_matrix_0.9),  col = "red", lwd = 3)

legend(x = 'topright', pch = 15,
       col = c("black","green", "orange", "blue","red"),
       legend = c("no constraint","r>0.6", "r>0.7", "r>0.8", "r>0.9"))

# tumor cells, contrainst on mRNA-lncRNA correaltion
subset_tumor_miRNA_mRNA_normal_corr_matrix_0.9 = tumor_miRNA_mRNA_corr_matrix[,normal_mRNA_indices_0.9]
subset_tumor_miRNA_mRNA_normal_corr_matrix_0.8 = tumor_miRNA_mRNA_corr_matrix[,normal_mRNA_indices_0.8]
subset_tumor_miRNA_mRNA_normal_corr_matrix_0.7 = tumor_miRNA_mRNA_corr_matrix[,normal_mRNA_indices_0.7]
subset_tumor_miRNA_mRNA_normal_corr_matrix_0.6 = tumor_miRNA_mRNA_corr_matrix[,normal_mRNA_indices_0.6]


plot(density(tumor_miRNA_mRNA_corr_matrix), 
     main = "tumor_miRNA-mRNA correlation (mRNAs selected based on normal_lncRNA_mRNA correlation)", 
     lwd = 3, cex.main = 0.9)
#lines(density(subset_normal_miRNA_mRNA_corr_matrix_0.2), col = "gray", lwd = 3)
lines(density(subset_tumor_miRNA_mRNA_normal_corr_matrix_0.6), col = "green", lwd = 3)
lines(density(subset_tumor_miRNA_mRNA_normal_corr_matrix_0.7), col = "orange", lwd = 3)
lines(density(subset_tumor_miRNA_mRNA_normal_corr_matrix_0.8), col = "blue", lwd = 3)
lines(density(subset_tumor_miRNA_mRNA_normal_corr_matrix_0.9),  col = "red", lwd = 3)

legend(x = 'topright', pch = 15,
       col = c("black","green", "orange", "blue","red"),
       legend = c("no constraint","r>0.6", "r>0.7", "r>0.8", "r>0.9"))


# --------------- Paci et all ----------------------------------------------------------

# ------------- helper functions ------------------------------------------------------

get_lncRNA_mRNA_pairs = function(originalCorrMatrix, corrThreshold){
  correlation_pairs = which(originalCorrMatrix > corr_threshold, arr.ind = TRUE)
  lncRNA = rownames(originalCorrMatrix)[correlation_pairs[,1]]
  mRNA = colnames(originalCorrMatrix)[correlation_pairs[,2]]
  dataframe = as.data.frame(cbind(correlation_pairs, lncRNA, mRNA))
  colnames(dataframe) = c("lncRNA_index", "mRNA_index","lncRNA", "mRNA")
  dataframe$lncRNA_index = as.numeric(as.character(dataframe$lncRNA_index))
  dataframe$mRNA_index = as.numeric(as.character(dataframe$mRNA_index))
  # correlation_vector = normal_lncRNA_mRNA_corr_matrix[dataframe$lncRNA_index,dataframe$mRNA_index]
  dataframe = cbind(dataframe, originalCorrMatrix[which(originalCorrMatrix > corr_threshold)])
  colnames(dataframe) = c("lncRNA_index", "mRNA_index","lncRNA", "mRNA", "corr")
  dataframe = dataframe[order(-dataframe$corr),]
  return(dataframe)
}


# --------------- build sensitivity matrix from normal cells------------------------------
require(ppcor)

# select top pairs lncRNA-mRNA from normal dataset
corr_threshold = 0.9
lncRNA_mRNA_pairs = get_lncRNA_mRNA_pairs(normal_lncRNA_mRNA_corr_matrix, corr_threshold)

# get only subset expression data from normal cells 
normal_indices = 1:79
mRNA_normal = brca_mRNA_df[,normal_indices]; 
miRNA_normal = brca_miRNA_df[,normal_indices];
lncRNA_normal = brca_lncRNA_df[,normal_indices]

lncRNA_mRNA_pairs[1:3,]

# lncRNA_index mRNA_index            lncRNA      mRNA      corr
# 555           907       1686 ENSG00000229645.4 C14orf139 0.9997353
# 2831         3657       9931 ENSG00000265142.2      MYH4 0.9992861
# 3935          162      14252 ENSG00000203875.6     SNHG5 0.9992630

pairs = lncRNA_mRNA_pairs[,]

sensitivity_matrix = matrix(ncol = 343)

ptm = proc.time()
for (i in 1:nrow(pairs)){
  the_pair = pairs[i,]
  mRNA_vector = mRNA_normal[the_pair$mRNA_index,]
  lncRNA_vector = lncRNA_normal[the_pair$lncRNA_index,]
  sensitivity_corr_vector = c()
  for (j in 1:nrow(miRNA_normal)){
    partial_corr = pcor.test(mRNA_vector, lncRNA_vector, miRNA_normal[j,])$estimate
    sensitivity_corr_vector = append(sensitivity_corr_vector, the_pair$corr - partial_corr)
  }
  sensitivity_matrix <<- rbind(sensitivity_matrix, sensitivity_corr_vector)
}
ptm = proc.time() - ptm

sensitivity_matrix = sensitivity_matrix[-1,]
dim(sensitivity_matrix) # 4801 343 

sensitivity_matrix_0.9 = sensitivity_matrix
lncRNA_mRNA_pairs_0.9 = lncRNA_mRNA_pairs

lncRNA_mRNA_names = paste(lncRNA_mRNA_pairs_0.9$lncRNA, lncRNA_mRNA_pairs_0.9$mRNA, sep = "-")

rownames(sensitivity_matrix_0.9) = lncRNA_mRNA_names 
colnames(sensitivity_matrix_0.9) = rownames(brca_miRNA_df)
save(lncRNA_mRNA_pairs, sensitivity_matrix_0.9, file = "Saved_R_Objects/corr_matrices/sensitivity_matrix_0.9.rda")

require(gplots)

# --------------- build sensitivity matrix from cancer cells------------------------------
corr_threshold = 0.8
lncRNA_mRNA_pairs = get_lncRNA_mRNA_pairs(tumor_lncRNA_mRNA_corr_matrix, corr_threshold)








