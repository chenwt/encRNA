load("data_Saved_R_Objects/brca_df.rda")
load("data_Saved_R_Objects/corr_matrices/all_corr_matrices_with_names.rda")
load("data_Saved_R_Objects/corr_matrices/normal_tumor_lncRNA_mRNA_pairs.rda")


# get only subset expression data from normal cells
# normal_indices = 1:79
mRNA_normal = brca_mRNA_df[,normal_indices];
miRNA_normal = brca_miRNA_df[,normal_indices];
lncRNA_normal = brca_lncRNA_df[,normal_indices]

# get sensitivity matrix for normal samples
normal_pairs = normal_lncRNA_mRNA_pairs

require(ppcor)

ptm_normal = proc.time()
for(i in 1:1){
  the_pair = normal_pairs[i,]
  mRNA_vector = mRNA_normal[the_pair$mRNA_index,]
  lncRNA_vector = lncRNA_normal[the_pair$lncRNA_index,]
  mRNA_name = as.character(the_pair$mRNA);
  lncRNA_name = as.character(the_pair$lncRNA);
  for (j in 1:nrow(miRNA_normal)){
    miRNA_vector = miRNA_normal[j,]
    miRNA_name = as.character(rownames(miRNA_normal)[j]);
    
    par_cor_list = parCorList(mRNA_name,lncRNA_name,miRNA_name, status = "normal");
  }
}
ptm_normal = proc.time() - ptm_normal;



# function to compute partial correlation partialCor(x,y|z) 
parCorList = function(mRNA_name,lncRNA_name,miRNA_name, status = "normal"){
  #load("data_Saved_R_Objects/corr_matrices/all_corr_matrices_with_names.rda")
  cor_mRNA_lncRNA = cor_mRNA_miRNA = cor_lncRNA_miRNA = 0
  if (status == "normal"){
    cor_mRNA_lncRNA = normal_lncRNA_mRNA_corr_matrix[lncRNA_name,mRNA_name]; 
    cor_mRNA_miRNA = normal_miRNA_mRNA_corr_matrix[miRNA_name,mRNA_name]; 
    cor_lncRNA_miRNA = normal_miRNA_lncRNA_corr_matrix[miRNA_name,lncRNA_name]; 
  }else if (status == "tumor"){
    cor_mRNA_lncRNA = tumor_lncRNA_mRNA_corr_matrix(lncRNA_name,mRNA_name); 
    cor_mRNA_miRNA = tumor_miRNA_mRNA_corr_matrix(miRNA_name,mRNA_name); 
    cor_lncRNA_miRNA = tumor_miRNA_lncRNA_corr_matrix(miRNA_name,lncRNA_name);  
  }
  
  pc_mRNA_lncRNA_miRNA = (cor_mRNA_lncRNA - cor_mRNA_miRNA * cor_lncRNA_miRNA) / 
    sqrt((1-cor_mRNA_miRNA^2)*(1-cor_lncRNA_miRNA^2))
  
  pc_mRNA_miRNA_lncRNA = (cor_mRNA_miRNA - cor_mRNA_lncRNA * cor_lncRNA_miRNA) / 
    sqrt((1-cor_mRNA_lncRNA^2)*(1-cor_mRNA_lncRNA^2))
  
  pc_lncRNA_miRNA_mRNA = (cor_lncRNA_miRNA - cor_mRNA_lncRNA * cor_mRNA_miRNA) / 
    sqrt((1-cor_mRNA_lncRNA^2)*(1-cor_mRNA_miRNA^2))
  
  epsilon = 1/3*(pc_mRNA_lncRNA_miRNA/cor_mRNA_lncRNA 
                 + pc_mRNA_miRNA_lncRNA/cor_mRNA_miRNA 
                 + pc_lncRNA_miRNA_mRNA/cor_lncRNA_miRNA)
  
  # work from here
  (abs(mRNA_lncRNA_cor) < abs(epsilon*mRNA_miRNA_cor) & abs(mRNA_lncRNA_cor) < abs(epsilon*lncRNA_miRNA_cor))
  
  return(list(
    epsilon = epsilon,
    pc_mRNA_lncRNA_miRNA = pc_mRNA_lncRNA_miRNA,
    pc_lncRNA_miRNA_mRNA = pc_lncRNA_miRNA_mRNA,
    pc_mRNA_miRNA_lncRNA = pc_mRNA_miRNA_lncRNA))
  
  
  #value = (cor_xy - cor_xz * cor_yz) / sqrt((1-cor_xz^2)*(1-cor_yz^2))
  #return(value)
}








num_cluster = 8
cl <- makePSOCKcluster(num_cluster)
registerDoParallel(cl)

ptm_normal = proc.time()
normal_sensitivity_full = foreach(i = 1:nrow(normal_pairs), .combine = rbind) %dopar% {
  require(ppcor)
  the_pair = normal_pairs[i,]
  mRNA_vector = mRNA_normal[the_pair$mRNA_index,]
  lncRNA_vector = lncRNA_normal[the_pair$lncRNA_index,]
  sensitivity_corr_vector = c()
  for (j in 1:nrow(miRNA_normal)){
    partial_corr = pcor.test(mRNA_vector, lncRNA_vector, miRNA_normal[j,])$estimate
    sensitivity_corr_vector = append(sensitivity_corr_vector, the_pair$corr - partial_corr)
  }
  return(sensitivity_corr_vector)
}
ptm_normal = proc.time() - ptm_normal;
ptm_normal
stopCluster(cl)

sprintf("Running time with %g cores: %g", num_cluster ,ptm_normal[3])
dim(normal_sensitivity_full)
save(normal_sensitivity_full, file = "normal_sensivitivity_testFull.rda")
