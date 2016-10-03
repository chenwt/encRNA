setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")
load("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/data_Saved_R_Objects/brca_df.rda")

setwd("/home/MARQNET/0099dod/encRNA/biweight_correlation")
load("/home/MARQNET/0099dod/encRNA/brca_df.rda")
# load("data_Saved_R_Objects/corr_matrices/all_corr_matrices_with_names.rda")
# load("data_Saved_R_Objects/corr_matrices/normal_tumor_lncRNA_mRNA_pairs.rda")
# load("data_Saved_R_Objects/corr_matrices/normal_sensitivity_full.rda")

### ----------- SET UP DATA FRAME OF GENE EXPRESSION -------------------------------------
# get only subset expression data from normal cells
normal_indices = 1:79
mRNA_normal = brca_mRNA_df[,normal_indices];
miRNA_normal = brca_miRNA_df[,normal_indices];
lncRNA_normal = brca_lncRNA_df[,normal_indices]

tumor_indices = 80:ncol(brca_mRNA_df)
mRNA_tumor = brca_mRNA_df[,tumor_indices];
miRNA_tumor = brca_miRNA_df[,tumor_indices];
lncRNA_tumor = brca_lncRNA_df[,tumor_indices]

## ------------ TEST BI-WEIGHT CORRELATION ------------------------------------------------
require(WGCNA); require(doParallel); require(foreach); require(doSNOW)

matrix1 = lncRNA_tumor
matrix2 = mRNA_tumor[1:10,]

num_cluster = parallel::detectCores()
#cl <- snow::makeCluster(num_cluster, outfile = "")
cl <- snow::makeCluster(num_cluster)
doSNOW::registerDoSNOW(cl)
running_time  = proc.time()
result = foreach(i = 1:nrow(matrix1), .combine = rbind) %dopar% {
  #print(paste("completion:", round(i*100/4828,3), "%", sep = " "))
  require(WGCNA)
  matrix1_vector = matrix1[i,]
  corr = sapply(1:nrow(matrix2), function(j){
    WGCNA::bicor(matrix1_vector, matrix2[j,])
  })
  return(corr)
}
snow::stopCluster(cl)
running_time  = proc.time() - running_time; 
sprintf("Running time with %g cores: %g", num_cluster ,running_time [3])
# dim(result)
# nrow(miRNA_normal); nrow(lncRNA_normal)
rownames(result) = rownames(matrix1); colnames(result) = rownames(matrix2)

tumor_bicor_lncRNA_mRNA = result
save(tumor_bicor_lncRNA_mRNA, file = "biweight_correlation/tumor_bicor_lncRNA_mRNA.rda")


## ------------ TEST MUTUAL INFORMATION ------------------------------------------------
require(infotheo)

# discretize data frame 
mRNA_normal = infotheo::discretize(mRNA_normal)
miRNA_normal = infotheo::discretize(miRNA_normal)
lncRNA_normal = infotheo::discretize(lncRNA_normal)

mRNA_tumor = infotheo::discretize(mRNA_tumor)
miRNA_tumor = infotheo::discretize(miRNA_tumor)
lncRNA_tumor = infotheo::discretize(lncRNA_tumor)

mutualInfoVector = c()
for (i in 1:nrow(mRNA_normal)){
  print(i)
  I = infotheo::mutinformation(miRNA_normal[2,], mRNA_normal[i,])
  mutualInfoVector = append(mutualInfoVector, I)
}


#---------------------------------------------------------------------------------------





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

# require(ppcor); require(foreach); require(doParallel);
# 
# num_cluster = 4
# cl <- makePSOCKcluster(num_cluster)
# registerDoParallel(cl)
# # matrix = matrix(NA, ncol = 343)
# running_time  = proc.time()
# #iter = nrow(normal_pairs)
# iter = 1000
# normal_sensitivity = foreach(i = 1:iter, .combine = rbind) %dopar% {
#   the_pair = normal_pairs[i,]
#   mRNA_name = as.character(the_pair$mRNA);
#   lncRNA_name = as.character(the_pair$lncRNA);
#   sensitivity_corr_vector = c()
#   for (j in 1:nrow(miRNA_normal)){
#     miRNA_name = as.character(rownames(miRNA_normal)[j]);
#     # par_cor_list = parCorList(mRNA_name,lncRNA_name,miRNA_name, status = "normal");
#     #sensitivity = parCorList(mRNA_name,lncRNA_name,miRNA_name, status = "normal");
#     
#     cor_mRNA_lncRNA = normal_lncRNA_mRNA_corr_matrix[lncRNA_name,mRNA_name]; 
#     cor_mRNA_miRNA = normal_miRNA_mRNA_corr_matrix[miRNA_name,mRNA_name]; 
#     cor_lncRNA_miRNA = normal_miRNA_lncRNA_corr_matrix[miRNA_name,lncRNA_name]; 
#     #print(cor_mRNA_lncRNA); print(cor_mRNA_miRNA);   print(cor_lncRNA_miRNA)
#     
#     pc_mRNA_lncRNA_miRNA = (cor_mRNA_lncRNA - cor_mRNA_miRNA * cor_lncRNA_miRNA) / 
#       sqrt((1-cor_mRNA_miRNA^2)*(1-cor_lncRNA_miRNA^2))
#     
#     pc_lncRNA_miRNA_mRNA = (cor_lncRNA_miRNA - cor_mRNA_lncRNA * cor_mRNA_miRNA) / 
#       sqrt((1-cor_mRNA_lncRNA^2)*(1-cor_mRNA_miRNA^2))
#     
#     pc_mRNA_miRNA_lncRNA = (cor_mRNA_miRNA - cor_mRNA_lncRNA * cor_lncRNA_miRNA) / 
#       sqrt((1-cor_mRNA_lncRNA^2)*(1-cor_mRNA_lncRNA^2))
#     
#     epsilon = 1/3*(pc_mRNA_lncRNA_miRNA/cor_mRNA_lncRNA 
#                    + pc_mRNA_miRNA_lncRNA/cor_mRNA_miRNA 
#                    + pc_lncRNA_miRNA_mRNA/cor_lncRNA_miRNA)
#     
#     # work from here
#     condition_satisfied = (abs(cor_mRNA_lncRNA) < abs(epsilon*cor_mRNA_miRNA)) & (abs(cor_mRNA_lncRNA) < abs(epsilon*cor_lncRNA_miRNA))
#     
#     sensitivity = cor_mRNA_lncRNA - pc_mRNA_lncRNA_miRNA
#     
#     sensitivity_corr_vector = append(sensitivity_corr_vector, sensitivity)
#   }
#   return(sensitivity_corr_vector)
# }
# stopCluster(cl)
# 
# running_time  = proc.time() - running_time ;
# sprintf("Running time with %g cores: %g", num_cluster ,running_time [3])
# 
# 
# 
# # function to compute partial correlation partialCor(x,y|z) 
# parCorList = function(mRNA_name,lncRNA_name,miRNA_name, status = "normal"){
#   #load("data_Saved_R_Objects/corr_matrices/all_corr_matrices_with_names.rda")
#   cor_mRNA_lncRNA = cor_mRNA_miRNA = cor_lncRNA_miRNA = 0
#   if (status == "normal"){
#     cor_mRNA_lncRNA = normal_lncRNA_mRNA_corr_matrix[lncRNA_name,mRNA_name]; 
#     cor_mRNA_miRNA = normal_miRNA_mRNA_corr_matrix[miRNA_name,mRNA_name]; 
#     cor_lncRNA_miRNA = normal_miRNA_lncRNA_corr_matrix[miRNA_name,lncRNA_name]; 
#   }else if (status == "tumor"){
#     cor_mRNA_lncRNA = tumor_lncRNA_mRNA_corr_matrix(lncRNA_name,mRNA_name); 
#     cor_mRNA_miRNA = tumor_miRNA_mRNA_corr_matrix(miRNA_name,mRNA_name); 
#     cor_lncRNA_miRNA = tumor_miRNA_lncRNA_corr_matrix(miRNA_name,lncRNA_name);  
#   }
#   
#   #print(cor_mRNA_lncRNA); print(cor_mRNA_miRNA);   print(cor_lncRNA_miRNA)
#   
#   pc_mRNA_lncRNA_miRNA = (cor_mRNA_lncRNA - cor_mRNA_miRNA * cor_lncRNA_miRNA) / 
#     sqrt((1-cor_mRNA_miRNA^2)*(1-cor_lncRNA_miRNA^2))
#   
#   pc_lncRNA_miRNA_mRNA = (cor_lncRNA_miRNA - cor_mRNA_lncRNA * cor_mRNA_miRNA) / 
#     sqrt((1-cor_mRNA_lncRNA^2)*(1-cor_mRNA_miRNA^2))
#   
#   pc_mRNA_miRNA_lncRNA = (cor_mRNA_miRNA - cor_mRNA_lncRNA * cor_lncRNA_miRNA) / 
#     sqrt((1-cor_mRNA_lncRNA^2)*(1-cor_mRNA_lncRNA^2))
#   
#   epsilon = 1/3*(pc_mRNA_lncRNA_miRNA/cor_mRNA_lncRNA 
#                  + pc_mRNA_miRNA_lncRNA/cor_mRNA_miRNA 
#                  + pc_lncRNA_miRNA_mRNA/cor_lncRNA_miRNA)
#   
#   # work from here
#   condition_satisfied = (abs(cor_mRNA_lncRNA) < abs(epsilon*cor_mRNA_miRNA)) & (abs(cor_mRNA_lncRNA) < abs(epsilon*cor_lncRNA_miRNA))
#   
#   sensitivity = cor_mRNA_lncRNA - pc_mRNA_lncRNA_miRNA
#   
#   return(sensitivity)
  
  # return(list(
  #   condition_satisfied = condition_satisfied,
  #   sensitivity = sensitivity,
  #   epsilon = epsilon,
  #   pc_mRNA_lncRNA_miRNA = pc_mRNA_lncRNA_miRNA,
  #   pc_lncRNA_miRNA_mRNA = pc_lncRNA_miRNA_mRNA,
  #   pc_mRNA_miRNA_lncRNA = pc_mRNA_miRNA_lncRNA))
#}


