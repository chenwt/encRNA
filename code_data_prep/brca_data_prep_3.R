setwd("/media/ducdo/UUI1/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
# load("brca_common_objects2.rda")
load("data_Saved_R_ObjectsSaved_R_Objects/brca_common_objects3.rda")
load("data_Saved_R_Objects/brca_df.rda")

colnames(brca_mRNA_common)[1:10]
colnames(brca_miRNA_common)[1:10]
colnames(brca_lncRNA_common)[1:10]




## --------------- BASIC SETUP ---------------------------------------

# remove other variables if needed  
rm(list = ls()); gc()

# set directory
setwd("/media/ducdo/UUI1/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
load("brca2_TCGA2STAT.rda")
require(TCGA2STAT)

## --------------- SAMPLE MATCHING ---------------------------------------

### mRNA ########
dim(brca_mRNA$dat) #[1] 20502   878

brca_mRNA_bytype = SampleSplit(brca_mRNA$dat)
brca_mRNA_normal_samples = colnames(brca_mRNA_bytype$normal); 
length(brca_mRNA_normal_samples) # 100
brca_mRNA_tumor_samples = colnames(brca_mRNA_bytype$primary.tumor) # 1093
length(brca_mRNA_tumor_samples) # 775


brca_mRNA_normal_samples[1] # [1] "TCGA-A7-A0D9-11A-53R-A089-07-2"

brca_mRNA_normal_samples_v2 =  strsplit(brca_mRNA_normal_samples, split = "[-]")
brca_mRNA_normal_samples_v2 = unlist(lapply(brca_mRNA_normal_samples_v2, function(name) paste("BRCA","Normal",name[1],name[2],name[3],sep=".")))
brca_mRNA_tumor_samples_v2 =  strsplit(brca_mRNA_tumor_samples, split = "[-]")
brca_mRNA_tumor_samples_v2 = unlist(lapply(brca_mRNA_tumor_samples_v2, function(name) paste("BRCA","Tumor",name[1],name[2],name[3],sep=".")))

### find overlapped samples with miRNA and lncNRA #######
brca_miRNA <- getTCGA(disease="BRCA", data.type="miRNASeq", type = "rpmmm", clinical = T)
save(brca_mRNA, brca_miRNA, file = "brca2_TCGA2STAT.rda")

### get miRNA normal and tumor samples, then manipulate the sameple names ######
brca.miRNA.bytype = SampleSplit(brca_miRNA$dat)
brca_miRNA_normal_samples = colnames(brca.miRNA.bytype$normal) 
length(brca_miRNA_normal_samples) # 87
brca_miRNA_tumor_samples = colnames(brca.miRNA.bytype$primary.tumor) 
length(brca_miRNA_tumor_samples) # 755

brca_miRNA_normal_samples[1] # TCGA-A7-A13E-11A-61R-A12O-13"
brca_miRNA_normal_samples_v2 =  strsplit(brca_miRNA_normal_samples, split = "[-]")
brca_miRNA_normal_samples_v2 = unlist(lapply(brca_miRNA_normal_samples_v2, function(name) paste("BRCA","Normal",name[1],name[2],name[3],sep=".")))
brca_miRNA_tumor_samples_v2 =  strsplit(brca_miRNA_tumor_samples, split = "[-]")
brca_miRNA_tumor_samples_v2 = unlist(lapply(brca_miRNA_tumor_samples_v2, function(name) paste("BRCA","Tumor",name[1],name[2],name[3],sep=".")))

#### load lncRNA data ########
brca_lncRNA_df = read.csv(file = "data_expression_files/lncRNA_Tanric_062616.csv", header = TRUE)
rownames(brca_lncRNA_df) = brca_lncRNA_df[,1]
brca_lncRNA_df = brca_lncRNA_df[,-1]
dim(brca_lncRNA_df) # [1] 12727   942

brca_lncRNA_normal_samples = colnames(brca_lncRNA_df)[grep("Normal",colnames(brca_lncRNA_df))] # 105
brca_lncRNA_tumor_samples = colnames(brca_lncRNA_df)[grep("Tumor",colnames(brca_lncRNA_df))] # 837; note that all normal samples have paired tumor samples

length(brca_lncRNA_normal_samples) # 105
length(brca_lncRNA_tumor_samples) # 837

# find the overlapped names
length(intersect(intersect(brca_mRNA_normal_samples_v2,brca_miRNA_normal_samples_v2), 
                 brca_lncRNA_normal_samples)) # 79
length(intersect(intersect(brca_mRNA_tumor_samples_v2,brca_miRNA_tumor_samples_v2),
                 brca_lncRNA_tumor_samples)) # 457

# thus, accross all mRNA, miRNA, lnRNA, there are common 79 normal and common 457 samples
brca_normal_common_samples = intersect(intersect(brca_mRNA_normal_samples_v2,brca_miRNA_normal_samples_v2), brca_lncRNA_normal_samples)
brca_tumor_common_samples = intersect(intersect(brca_mRNA_tumor_samples_v2,brca_miRNA_tumor_samples_v2), brca_lncRNA_tumor_samples)
length(brca_normal_common_samples); length(brca_tumor_common_samples)

#### subset lncRNA
brca_lncRNA_common = brca_lncRNA_df[,which(colnames(brca_lncRNA_df) %in% c(brca_normal_common_samples,brca_tumor_common_samples))]
dim(brca_lncRNA_common) # [1] 12727   536

#### reformat mRNA and miRNA
brca_normal_common_samples_v3 = strsplit(brca_normal_common_samples, split = "[.]")
brca_tumor_common_samples_v3 = strsplit(brca_tumor_common_samples, split = "[.]")

brca_normal_common_samples_v3 = unlist(lapply(
  brca_normal_common_samples_v3, function(name) paste(name[3],name[4],name[5],"11", sep="-")))

brca_tumor_common_samples_v3 = unlist(lapply(
  brca_tumor_common_samples_v3, function(name) paste(name[3],name[4],name[5],"01", sep="-")))

#### subset mRNA
brca_mRNA_common = brca_mRNA$dat[,c(brca_mRNA_normal_samples, brca_mRNA_tumor_samples)]
# brca_mRNA_common = brca_rnaV2$dat[,(which( substr(colnames(brca_mRNA$dat),1,15) %in% c(brca_normal_common_samples_v3, brca_tumor_common_samples_v3) ))]
colnames(brca_mRNA_common) = c(brca_mRNA_normal_samples_v2, brca_mRNA_tumor_samples_v2)
brca_mRNA_common = brca_mRNA_common[,colnames(brca_lncRNA_common)]


#### subset miRNA
brca_miRNA_common = brca_miRNA$dat[,c(brca_miRNA_normal_samples, brca_miRNA_tumor_samples)]
colnames(brca_miRNA_common) = c(brca_miRNA_normal_samples_v2, brca_miRNA_tumor_samples_v2)
brca_miRNA_common = brca_miRNA_common[,colnames(brca_lncRNA_common)]

### check matched normal-tumor

length(intersect(intersect(colnames(brca_mRNA_common),colnames(brca_lncRNA_common)), colnames(brca_miRNA_common)))
normal_samples = colnames(brca_mRNA_common)[1:79]; normal_samples = substr(normal_samples, 13, 25)
tumor_samples = colnames(brca_mRNA_common)[80:536]; tumor_samples = substr(tumor_samples, 12, 24)
length(intersect(normal_samples, tumor_samples))

load("data_Saved_R_Objects/brca_df.rda")
length(intersect(intersect(colnames(brca_mRNA_df),colnames(brca_lncRNA_df)), colnames(brca_miRNA_df)))
normal_samples = colnames(brca_mRNA_df)[1:79]; normal_samples = substr(normal_samples, 13, 25)
tumor_samples = colnames(brca_mRNA_df)[80:536]; tumor_samples = substr(tumor_samples, 12, 24)
length(intersect(normal_samples, tumor_samples))

common_name_core = intersect(normal_samples, tumor_samples) 
normal_samples = paste("BRCA.Normal.", common_name_core, sep = "")
tumor_samples = paste("BRCA.Tumor.", common_name_core, sep = "")
sample_names = c(normal_samples, tumor_samples)

brca_mRNA_matched = brca_mRNA_df[,sample_names]; dim(brca_mRNA_matched)
brca_miRNA_matched = brca_miRNA_df[,sample_names]; dim(brca_miRNA_matched)
brca_lncRNA_matched = brca_lncRNA_df[,sample_names]; dim(brca_lncRNA_matched)

save(brca_mRNA_matched,brca_miRNA_matched,brca_lncRNA_matched,file = "data_Saved_R_Objects/brca_expression_matched.rda")

# check dimension matching
dim(brca_mRNA_common) # 20502   536
dim(brca_miRNA_common) #  1046  536
dim(brca_lncRNA_common) #  12727   536

save(brca_mRNA_common,brca_miRNA_common,brca_lncRNA_common, file = "Saved_R_Objects/brca_common_objects3.rda")

par(mfrow = c(3,1))
plot(density(brca_mRNA_common[,1]),lwd=3,ylim=c(0,1), main = "brca_mRNA_common")
for(i in 2:ncol(brca_mRNA_common)){lines(density(brca_mRNA_common[,i]),lwd=3)}
plot(density(brca_miRNA_common[,1]),lwd=3,ylim=c(0,5), main = "brca_miRNA_common")
for(i in 2:ncol(brca_miRNA_common)){lines(density(brca_miRNA_common[,i]),lwd=3)}
plot(density(brca_lncRNA_common[,1]),lwd=3,ylim=c(0,20), main = "brca_lncRNA_common")
for(i in 2:ncol(brca_lncRNA_common)){lines(density(brca_lncRNA_common[,i]),lwd=3)}

## --------------- GENE CLEANING (Remove low count gene) ---------------------------------------

## ------ lncRNA ---------
# First, remove all genes whose variance = 0 accross al samples
# check genes with variance = 0
dim(brca_lncRNA_common) # 12727   536
brca_lncRNA_var = apply(brca_lncRNA_common, 1, var)
length(which(brca_lncRNA_var  == 0)) # 41

# those genes also have expression = 0 accross all samples
brca_lncRNA_df = brca_lncRNA_common[-which(brca_lncRNA_var  == 0),]
dim(brca_lncRNA_df) # [1] 12686   536

#### second, look at all cross all normal samples, select a gene if more than 75% of sample
# having expression more than 0. Among the selected genes, look at all tumor samples, 
# select those that have expression in more than 75% cases

selected_threshold = 0.75
# split into normal and tumor subset
brca_lncRNA_df_normal = brca_lncRNA_df[,1:79]
brca_lncRNA_df_tumor = brca_lncRNA_df[,80:ncol(brca_lncRNA_df)]

brca_total_normal = ncol(brca_lncRNA_df_normal) # 79
brca_total_tumor= ncol(brca_lncRNA_df_tumor) #457

brca_lncRNA_is_kept = c(); 

brca_lncRNA_is_kept = sapply(1:nrow(brca_lncRNA_df), function(i) 
{
  sum(brca_lncRNA_df_normal[i,] > 0) > selected_threshold*brca_total_normal
}
)

sum(brca_lncRNA_is_kept) # [1]  5768

# create vector to retain names of kept genes
brca_lncRNA_kept_genes = rownames(brca_lncRNA_df)[brca_lncRNA_is_kept]

# work with tumor subset
brca_lncRNA_df_tumor = brca_lncRNA_df_tumor[brca_lncRNA_is_kept, ]
brca_lncRNA_is_kept = sapply(1:nrow(brca_lncRNA_df_tumor), function(i) 
{
  sum(brca_lncRNA_df_tumor[i,] > 0) > selected_threshold*brca_total_tumor
}
)
sum(brca_lncRNA_is_kept) # 4828

brca_lncRNA_kept_genes = brca_lncRNA_kept_genes[brca_lncRNA_is_kept]
brca_lncRNA_df = subset(brca_lncRNA_df, rownames(brca_lncRNA_df) %in% brca_lncRNA_kept_genes)
brca_lncRNA_df = log2(brca_lncRNA_df + 1)
dim(brca_lncRNA_df) # [1] 4828  536

## ------ mRNA ---------

# First, remove all genes whose variance = 0 accross al samples
# check genes with variance = 0
dim(brca_mRNA_common) # 20502   536
brca_mRNA_var = apply(brca_mRNA_common, 1, var)
length(which(brca_mRNA_var  == 0)) # 17

# those genes also have expression = 0 accross all samples
brca_mRNA_df = brca_mRNA_common[-which(brca_mRNA_var  == 0),]
dim(brca_mRNA_df) # [1] 20485   536

#### second, look at all cross all normal samples, select a gene if more than 75% of sample
# having expression more than 0. Among the selected genes, look at all tumor samples, 
# select those that have expression in more than 75% cases

selected_threshold = 0.75
# work with normal subset
brca_mRNA_df_normal = brca_mRNA_df[,1:79]
brca_mRNA_df_tumor = brca_mRNA_df[,80:ncol(brca_mRNA_df)]

brca_total_normal = ncol(brca_mRNA_df_normal) # 79
brca_total_tumor= ncol(brca_mRNA_df_tumor) #457

brca_mRNA_is_kept = c(); 

brca_mRNA_is_kept = sapply(1:nrow(brca_mRNA_df), function(i) 
{
  sum(brca_mRNA_df_normal[i,] > 0) > selected_threshold*brca_total_normal
}
)

sum(brca_mRNA_is_kept) # [1]  18013

# create vector to retain names of kept genes
brca_mRNA_kept_genes = rownames(brca_mRNA_df)[brca_mRNA_is_kept]

# work with tumor subset
brca_mRNA_df_tumor = brca_mRNA_df_tumor[brca_mRNA_is_kept, ]
brca_mRNA_is_kept = sapply(1:nrow(brca_mRNA_df_tumor), function(i) 
{
  sum(brca_mRNA_df_tumor[i,] > 0) > selected_threshold*brca_total_tumor
}
)
sum(brca_mRNA_is_kept) # 17613

brca_mRNA_kept_genes = brca_mRNA_kept_genes[brca_mRNA_is_kept]
brca_mRNA_df = subset(brca_mRNA_df, rownames(brca_mRNA_df) %in% brca_mRNA_kept_genes)
brca_mRNA_df = log2(brca_mRNA_df + 1)
dim(brca_mRNA_df) # [1] 17613  536

## ------ miRNA ---------

# First, remove all genes whose variance = 0 accross al samples
# check genes with variance = 0
dim(brca_miRNA_common) # [1] 1046  536
brca_miRNA_var = apply(brca_miRNA_common, 1, var)
length(which(brca_miRNA_var  == 0)) # 169

# those genes also have expression = 0 accross all samples
brca_miRNA_df = brca_miRNA_common[-which(brca_miRNA_var  == 0),]
dim(brca_miRNA_df) # [1] 877   536

#### second, look at all cross all normal samples, select a gene if more than 75% of sample
# having expression more than 0. Among the selected genes, look at all tumor samples, 
# select those that have expression in more than 75% cases

selected_threshold = 0.75
# work with normal subset
brca_miRNA_df_normal = brca_miRNA_df[,1:79]
brca_miRNA_df_tumor = brca_miRNA_df[,80:ncol(brca_miRNA_df)]

brca_total_normal = ncol(brca_miRNA_df_normal) # 79
brca_total_tumor= ncol(brca_miRNA_df_tumor) #457

brca_miRNA_is_kept = c(); 

brca_miRNA_is_kept = sapply(1:nrow(brca_miRNA_df), function(i) 
{
  sum(brca_miRNA_df_normal[i,] > 0) > selected_threshold*brca_total_normal
}
)

sum(brca_miRNA_is_kept) # [1]  360

# create vector to retain names of kept genes
brca_miRNA_kept_genes = rownames(brca_miRNA_df)[brca_miRNA_is_kept]

# work with tumor subset
brca_miRNA_df_tumor = brca_miRNA_df_tumor[brca_miRNA_is_kept, ]
brca_miRNA_is_kept = sapply(1:nrow(brca_miRNA_df_tumor), function(i) 
{
  sum(brca_miRNA_df_tumor[i,] > 0) > selected_threshold*brca_total_tumor
}
)
sum(brca_miRNA_is_kept) # 343

brca_miRNA_kept_genes = brca_miRNA_kept_genes[brca_miRNA_is_kept]
brca_miRNA_df = subset(brca_miRNA_df, rownames(brca_miRNA_df) %in% brca_miRNA_kept_genes)
brca_miRNA_df = log2(brca_miRNA_df + 1)
dim(brca_miRNA_df) # [1] 343  536

save(brca_lncRNA_df, brca_mRNA_df, brca_miRNA_df, file = "Saved_R_Objects/brca_df3.rda")

par(mfrow = c(3,1))
plot(density(brca_mRNA_df[,1]),lwd=3, main = "brca_mRNA_filtered_log2")
for(i in 2:ncol(brca_mRNA_df)){lines(density(brca_mRNA_df[,i]),lwd=3)}
plot(density(brca_miRNA_df[,1]),lwd=3, main = "brca_miRNAfiltered_log2")
for(i in 2:ncol(brca_miRNA_df)){lines(density(brca_miRNA_df[,i]),lwd=3)}
plot(density(brca_lncRNA_df[,1]),lwd=3, main = "brca_lncRNA_common_filtered_log2")
for(i in 2:ncol(brca_lncRNA_df)){lines(density(brca_lncRNA_df[,i]),lwd=3)}

par(mfrow = c(3,1))
hist(brca_mRNA_df)
hist(brca_miRNA_df)
hist(as.matrix(brca_lncRNA_df))

# ----------------- Correlation -----------------------------------------------
load("Saved_R_Objects/corr_matrices/corr_matrices.rda")
load("Saved_R_Objects/corr_matrices/normal_tumor_corr_matrices")

# test_mRNA = brca_mRNA_df
# test_miRNA = brca_miRNA_df
brca_lncRNA_df = as.matrix(brca_lncRNA_df)

dim(test_mRNA) #[1] 17613   536
dim(test_miRNA) #[1] 343 536
dim(test_lncRNA) #[1] 4828  536

getCorrMatrix = function(matrix1, matrix2){
  corr_matrix = c()
  for (i in 1:nrow(matrix1)){
    row_matrix1 = matrix1[i,]
    corr = sapply(1:nrow(matrix2), function(j) cor(row_matrix1, matrix2[j,]))
    corr_matrix = rbind(corr_matrix, corr)
  }
  rownames(corr_matrix) = NULL
  return(corr_matrix)
}


##### both normal and tumor ###############################

miRNA_mRNA_corr_matrix = c()
miRNA_lncRNA_corr_matrix = c()
lncRNA_mRNA_corr_matrix = c()

#### miRNA-mRNA correlation
# for (i in 1:nrow(test_miRNA)){
#   miRNA = test_miRNA[i,]
#   corr = sapply(1:nrow(test_mRNA), function(j) cor(miRNA, test_mRNA[j,]))
#   miRNA_mRNA_corr_matrix = rbind(miRNA_mRNA_corr_matrix, corr)
# }
# 
# rownames(miRNA_mRNA_corr_matrix) = NULL

miRNA_mRNA_corr_matrix = getCorrMatrix(brca_miRNA_df, brca_mRNA_df);

hist(miRNA_mRNA_corr_matrix)
summary(miRNA_mRNA_corr_matrix)

#### miRNA-lncRNA correlation
# for (i in 1:nrow(test_miRNA)){
#   miRNA = test_miRNA[i,]
#   corr = sapply(1:nrow(test_lncRNA), function(j) cor(miRNA, test_lncRNA[j,]))
#   miRNA_lncRNA_corr_matrix = rbind(miRNA_lncRNA_corr_matrix, corr)
# }
# 
# rownames(miRNA_lncRNA_corr_matrix) = NULL

miRNA_lncRNA_corr_matrix = getCorrMatrix(brca_miRNA_df, as.matrix(brca_mRNA_df));
hist(miRNA_lncRNA_corr_matrix)
summary(as.vector(miRNA_lncRNA_corr_matrix))

save(miRNA_mRNA_corr_matrix, miRNA_lncRNA_corr_matrix, file = "Saved_R_Objects/corr_matrices/corr_matrices.rda")

#### lncRNA-mRNA correlation
lncRNA_mRNA_corr_matrix = getCorrMatrix(brca_lncRNA_df, brca_mRNA_df)

# load("Saved_R_Objects/corr_matrices/corr_matrices.rda")
save(miRNA_mRNA_corr_matrix, miRNA_lncRNA_corr_matrix, lncRNA_mRNA_corr_matrix, file = "Saved_R_Objects/corr_matrices/corr_matrices.rda")

##### normal only ###############################


mRNA_normal = brca_mRNA_df[,1:79]; mRNA_tumor =  brca_mRNA_df[,80:ncol(brca_mRNA_df)];
miRNA_normal = brca_miRNA_df[,1:79]; miRNA_tumor =  brca_miRNA_df[,80:ncol(brca_miRNA_df)];
lncRNA_normal = as.matrix(brca_lncRNA_df[,1:79]); 
lncRNA_tumor = as.matrix(brca_lncRNA_df[,80:ncol(brca_lncRNA_df)]);

normal_lncRNA_mRNA_corr_matrix = getCorrMatrix(lncRNA_normal, mRNA_normal)
tumor_lncRNA_mRNA_corr_matrix = getCorrMatrix(lncRNA_tumor, mRNA_tumor)

normal_miRNA_mRNA_corr_matrix = getCorrMatrix(miRNA_normal, mRNA_normal)
tumor_miRNA_mRNA_corr_matrix = getCorrMatrix(miRNA_tumor, mRNA_tumor)

normal_miRNA_lncRNA_corr_matrix = getCorrMatrix(miRNA_normal, lncRNA_normal)
tumor_miRNA_lncRNA_corr_matrix = getCorrMatrix(miRNA_tumor, lncRNA_tumor)

save(normal_lncRNA_mRNA_corr_matrix,
     tumor_lncRNA_mRNA_corr_matrix,
     normal_miRNA_mRNA_corr_matrix, 
     tumor_miRNA_mRNA_corr_matrix, 
     normal_miRNA_lncRNA_corr_matrix,
     tumor_miRNA_lncRNA_corr_matrix, 
     file = "Saved_R_Objects/corr_matrices/normal_tumor_corr_matrices.rda")






