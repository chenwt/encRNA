#########################################################
## ------ perform survival analysis----------------------
setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_network_construction/network_helper_functions.R")
load("brca2_TCGA2STAT.rda")
load("data_Saved_R_Objects/brca_df.rda")
load("data_Saved_R_Objects/miRNA_target/predicted_normal_tumor_goodCorr.rda")
require(survival)

## ------- sample name manipulation ----------------------------------------------
clinical_df = brca_miRNA$clinical

## select sample names also showed up in filtered normal and tumor sample 
colnames(brca_miRNA_df)[1:10]

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
length(intersect(rownames(clinical_df), tumor_samples)) # 457 -->
# thus all tumor samples are included in the clinical_df dataset

## ------- subset clinical_df data ----------------------------------------------
clinical_df = clinical_df[tumor_samples,c("vitalstatus","daystodeath","daystolastfollowup", "pathologyTstage")]
clinical_df = as.data.frame(clinical_df)
clinical_df$vitalstatus = as.factor(as.numeric(clinical_df$vitalstatus) - 1)
clinical_df$daystodeath = as.numeric(as.character(clinical_df$daystodeath))
clinical_df$daystolastfollowup = as.numeric(as.character(clinical_df$daystolastfollowup))
clinical_df$pathologyTstage = as.integer(substr(clinical_df$pathologyTstage, start = 2, stop = 2))

temp =  strsplit(rownames(clinical_df), split = "[-]")
rownames(clinical_df) = unlist(lapply(temp, function(name){
  paste("BRCA","Tumor",name[1],name[2],name[3],sep=".")
}))

nrow(clinical_df) # 457
length(which(is.na(clinical_df$daystodeath)))
length(which(is.na(clinical_df$daystolastfollowup)))
# there is one sample does not have neither daystodeath or daystolastfollowup
which(is.na(clinical_df$daystodeath) & is.na(clinical_df$daystolastfollowup))
# TCGA-E9-A245 
# 413 
# remove those samples
clinical_df = clinical_df[-which(is.na(clinical_df$daystodeath) & is.na(clinical_df$daystolastfollowup)),]
nrow(clinical_df) # 456

# divide patients into 2 groups based on pathological stage
clinical_df$risk_level = ifelse(clinical_df$pathologyTstage < 3, "low-risk", "high-risk")

## -------  create a new column which aggregate daystodeath and daystolastfollowup -------------------
clinical_df$time = ifelse(is.na(clinical_df$daystodeath), 
                          as.numeric(as.character(clinical_df$daystolastfollowup)), 
                          as.numeric(as.character(clinical_df$daystodeath)))

clinical_df$time = ifelse(is.na(clinical_df$daystodeath), 
                          clinical_df$daystolastfollowup, 
                          clinical_df$daystodeath)
sum(is.na(clinical_df$time)) # 0 --> good!

## ------- start Suvival anlysis --------------------------------------------------

# create Survival object
require(survival)

# plot the Kaplan-Meier curve between two high-risk and low-risk group
clinical_df$SurvObj <- with(clinical_df, Surv(time, vitalstatus == 1))
fit <- survfit(clinical_df$SurvObj ~ clinical_df$risk_level)
plot(fit, mark.time = T,  col = c("red","blue"))
# perform log-rank test
survdiff(formula = clinical_df$SurvObj ~ clinical_df$risk_level)

coxph(formula = clinical_df$SurvObj ~ factor(clinical_df$risk_level))


# get expression data
sample_names = rownames(clinical_df)


### ... run network_sensitivity_based to get hub genes
summary(normal_lncRNA_degree)
selected_lncRNAs = names(which(tumor_lncRNA_degree >= quantile(tumor_lncRNA_degree, c(0.80))))
# cox propotional hazard model based on hub lncRNA and mRNA
lncRNA_expression = brca_lncRNA_df[selected_lncRNAs,sample_names]
# add expression data into the clinical matrix
clinical_df = cbind(clinical_df, t(lncRNA_expression))

hazard_ratio = pvalue_hazard_ratio = c()
for (i in 1:length(selected_lncRNAs)){
  phm_fit =  coxph(formula = as.formula(paste("SurvObj ~ risk_level + ", selected_lncRNAs[i])), data = clinical_df)
  hazard_ratio = append(hazard_ratio, summary(phm_fit)$coefficient[2,2])
  pvalue_hazard_ratio = append(pvalue_hazard_ratio, summary(phm_fit)$coefficient[2,5])
}

hazard_ratio[which(pvalue_hazard_ratio < 0.05)]

# TO DO: adding expression data as column into the clinical_df data_frame. 
# epxression data are selected from the modules 

# EXTRA: use Cox-Lasso as an idenpedent steps to select distinct genes in term of survival
