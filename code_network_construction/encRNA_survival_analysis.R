#########################################################
## ------ perform survival analysis----------------------
setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_network_construction/network_helper_functions.R")
load("brca2_TCGA2STAT.rda")
load("data_Saved_R_Objects/brca_df.rda")
load("data_Saved_R_Objects/miRNA_target/predicted_normal_tumor_goodCorr.rda")
require(survival); require(ggplot2)

## ----------- helper_function --------------------------------------------------------------------------

get_survival_stat = function(clinical_expression_df, RNA){
  hazard_ratio = pvalue_hazard_ratio = c()
  for (i in 1:length(RNA)){
    phm_fit =  coxph(formula = as.formula(paste("SurvObj ~ ", RNA[i])), data = clinical_expression_df)
    hazard_ratio = append(hazard_ratio, summary(phm_fit)$coefficient[1,2])
    pvalue_hazard_ratio = append(pvalue_hazard_ratio, summary(phm_fit)$coefficient[1,5])
  }
  names(hazard_ratio) = names(pvalue_hazard_ratio) = RNA
  return(list(hazard_ratio = hazard_ratio, pvalue_hazard_ratio = pvalue_hazard_ratio))
}


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

# get sample names s
sample_names = rownames(clinical_df)

## ------- start Suvival anlysis --------------------------------------------------

# create Survival object
require(survival)

# plot the Kaplan-Meier curve between two high-risk and low-risk group
clinical_df$SurvObj <- with(clinical_df, Surv(time, vitalstatus == 1))
save(clinical_df, file = "data_Saved_R_Objects/clinical_df.rda")

fit <- survfit(clinical_df$SurvObj ~ clinical_df$risk_level)
plot(fit, mark.time = T,  col = c("red","blue"),
     xlab = "Days",ylab = "Survival Proability",
     main = "Kaplan-Meier plot between low-risk vs high-risk group")
legend(x = "topright", c("low-risk","high-risk"), lty=c(1,1), col = c("blue","red"))
#text(x = 3.5,y = 100, labels = "p log_rank = 0.0375")
# perform log-rank test
survdiff(formula = clinical_df$SurvObj ~ clinical_df$risk_level)

coxph(formula = clinical_df$SurvObj ~ factor(clinical_df$risk_level))

### ... run network_sensitivity_based to get hub genes
# summary(normal)
# selected_lncRNAs = names(which(tumor_lncRNA_degree >= quantile(tumor_lncRNA_degree, c(0.80))))

# select all mRNAs and lncRNA from normal and tumor encRNA triplet
normal_mRNAs = unique(normal_encRNA_sensitivity_bound_goodCoor$mRNA)
normal_lncRNAs = unique(normal_encRNA_sensitivity_bound_goodCoor$lncRNA)
tumor_mRNAs = unique(tumor_encRNA_sensitivity_bound_goodCoor$mRNA); tumor_mRNAs = tumor_mRNAs[-241]
tumor_lncRNAs = unique(tumor_encRNA_sensitivity_bound_goodCoor$lncRNA)

### cox propotional hazard model based on those lncRNA and mRNA

# add expression data into the clinical matrix
load("data_Saved_R_Objects/clinical_df.rda")
normal_mRNAs_clinical_expression_df = cbind(clinical_df, t(brca_mRNA_df[normal_mRNAs,sample_names]))
normal_mRNAs_cox = get_survival_stat(clinical_expression_df = normal_mRNAs_clinical_expression_df,
                         RNA = normal_mRNAs)
length(which(normal_mRNAs_cox$pvalue_hazard_ratio < 0.05)) # 22

load("data_Saved_R_Objects/clinical_df.rda")
normal_lncRNAs_clinical_expression_df = cbind(clinical_df, t(brca_lncRNA_df[normal_lncRNAs,sample_names]))
normal_lncRNAs_cox = get_survival_stat(clinical_expression_df = normal_lncRNAs_clinical_expression_df,
                                     RNA = normal_lncRNAs)
length(which(normal_lncRNAs_cox$pvalue_hazard_ratio < 0.05)) # 3

load("data_Saved_R_Objects/clinical_df.rda")
tumor_mRNAs_clinical_expression_df = cbind(clinical_df, t(brca_mRNA_df[tumor_mRNAs,sample_names]))
tumor_mRNAs_cox = get_survival_stat(clinical_expression_df = tumor_mRNAs_clinical_expression_df,
                                     RNA = tumor_mRNAs)
length(which(tumor_mRNAs_cox$pvalue_hazard_ratio < 0.05)) # 44

load("data_Saved_R_Objects/clinical_df.rda")
tumor_lncRNAs_clinical_expression_df = cbind(clinical_df, t(brca_lncRNA_df[tumor_lncRNAs,sample_names]))
tumor_lncRNAs_cox = get_survival_stat(clinical_expression_df = tumor_lncRNAs_clinical_expression_df,
                                       RNA = tumor_lncRNAs)
length(which(tumor_lncRNAs_cox$pvalue_hazard_ratio < 0.05)) # 26

# create violin plot
mRNA_normal_hazard_ratio <- data.frame(group = "mRNA_normal_hazard_ratio", value = 
                  normal_mRNAs_cox$hazard_ratio[which(normal_mRNAs_cox$pvalue_hazard_ratio < 0.05)])
lncRNA_normal_hazard_ratio <- data.frame(group = "lncRNA_normal_hazard_ratio", value = 
                  tumor_mRNAs_cox$hazard_ratio[which(tumor_mRNAs_cox$pvalue_hazard_ratio < 0.05)])
mRNA_tumor_hazard_ratio <- data.frame(group = "mRNA_tumor_hazard_ratio", value = 
                                        tumor_mRNAs_cox$hazard_ratio[which(tumor_mRNAs_cox$pvalue_hazard_ratio < 0.05)])
lncRNA_tumor_hazard_ratio <- data.frame(group = "lncRNA_tumor_hazard_ratio", value = 
                                        tumor_lncRNAs_cox$hazard_ratio[which(tumor_lncRNAs_cox$pvalue_hazard_ratio < 0.05)])

plot.data = rbind(mRNA_normal_hazard_ratio, lncRNA_normal_hazard_ratio, 
                  mRNA_tumor_hazard_ratio, lncRNA_tumor_hazard_ratio)


ggplot(plot.data, aes(x=group, y=value, fill=group)) +  geom_violin()
ggplot(plot.data, aes(x=group, y=value, fill=group)) +  geom_boxplot()




