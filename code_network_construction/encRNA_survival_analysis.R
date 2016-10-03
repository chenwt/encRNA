#########################################################
## ------ perform survival analysis----------------------
setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_network_construction/network_helper_functions.R")
load("brca2_TCGA2STAT.rda")
load("data_Saved_R_Objects/brca_df.rda")
load("data_Saved_R_Objects/miRNA_target/predicted_normal_tumor_goodCorr.rda")
require(survival)

clinical = brca_miRNA$clinical

## select sample names also showed up in filtered normal and tumor sample 
colnames(brca_miRNA_df)

normal_samples = colnames(brca_miRNA_df)[1:79]; length(normal_samples) # 79
tumor_samples = colnames(brca_miRNA_df)[80:536]; length(tumor_samples) # 457

normal_samples =  strsplit(normal_samples, split = "[.]")
normal_samples = unlist(lapply(normal_samples, function(name){
  paste(name[3],name[4],name[5],sep="-")
}))

tumor_samples =  strsplit(tumor_samples, split = "[.]")
tumor_samples = unlist(lapply(tumor_samples, function(name){
  paste(name[3],name[4],name[5],sep="-")
}))

length(intersect(normal_samples, tumor_samples)) # 76 out of 79 normal samples having matched tumor samples
length(intersect(rownames(clinical), tumor_samples)) # 45 -->
# thus all tumor samples are included in the clinical dataset

# subset clinical data
clinical = clinical[tumor_samples,c("vitalstatus","daystodeath","daystolastfollowup")]
clinical = as.data.frame(clinical)
nrow(clinical) # 457
length(which(is.na(clinical$daystodeath)))
length(which(is.na(clinical$daystolastfollowup)))
# there is one sample does not have neither daystodeath or daystolastfollowup
which(is.na(clinical$daystodeath) & is.na(clinical$daystolastfollowup))
# TCGA-E9-A245 
# 413 
clinical = clinical[-which(is.na(clinical$daystodeath) & is.na(clinical$daystolastfollowup)),]
nrow(clinical) # 456

# create a new column which aggregate daystodeath and daystolastfollowup
clinical$time = ifelse(is.na(clinical$daystodeath), clinical$daystolastfollowup, clinical$daystodeath)
sum(is.na(clinical$time)) # 0 --> good!
clinical$SurvObj <- with(clinical, Surv(time, vitalstatus == 1))

# TO DO: adding expression data as column into the clinical data_frame. 
# epxression data are selected from the modules 

# EXTRA: use Cox-Lasso as an idenpedent steps to select distinct genes in term of survival
