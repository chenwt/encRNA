# -------------- Basic setup ---------------------------------------------------------------------
setwd("/media/ducdo/UUI1/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("code_data_prep/get_normal_tumor_expression_matrices.R")
source("code_data_prep/helper_functions.R")
load("data_Saved_R_Objects/mutation_copyNumber/brca_mutation.rda")
load("data_Saved_R_Objects/mutation_copyNumber/brca_cnasnp.rda")
load("data_Saved_R_Objects/brca2_TCGA2STAT.rda")
load("data_Saved_R_Objects/brca_df.rda")
load("data_Saved_R_Objects/annotation/mRNA_lncRNA_miRNA_annotation.rda")
require(TCGA2STAT)

### ----- DOWNLOAD DATA --------------------------------------------------
# brca_mutation <- getTCGA(disease="BRCA", 
#                          data.type="Mutation", 
#                          type="somatic", 
#                          clinical = TRUE)
# 
# save(brca_mutation, file = "brca_mutation.rda")

# brca_cnasnp <- getTCGA(disease="BRCA", data.type="CNA_SNP")
# save(brca_cnasnp, file = "data_Saved_R_Objects/mutation_copyNumber/brca_cnasnp.rda")

#### ----- BASIC CHECK --------------------------------------------------
brca_mutation_df = brca_mutation$dat; 
dim(brca_mutation_df) # 16806   977
View(brca_mutation_df[1:100,])

nrow(brca_mRNA_df); rownames(brca_mRNA_df)[1:10] # 17,613
ncol(brca_mRNA_df); colnames(brca_mRNA_df)[1:10] # 16,806

nrow(brca_mutation_df); rownames(brca_mutation_df)[1:10] # 16,806
ncol(brca_mutation_df); colnames(brca_mutation_df)[1:10] # 16,806

# how many gene names are overlapped from mRNA and mutation data? 
length(intersect(rownames(brca_mRNA_df), rownames(brca_mutation_df))) # 13,595

# how many gene names in the mutation data are also annoated in ensembl?
length(intersect(rownames(brca_mutation_df), brca_mRNA_name_annotation$hgnc_symbol)) # 13,018

sample_mutation_split = SampleSplit(brca_mutation$dat)
dim(sample_mutation_split$primary.tumor)

# ----------------- Explore Mutation------------------------------------------------------
dim(brca_mutation$dat) # [1] 16806   977
head(brca_mutation$dat)[1:5,1:5]
brca_mutation_df = brca_mutation$dat; 
dim(brca_mutation_df) # 16806   977

# compare with the original mRNA matrix
require(TCGA2STAT)
brca_mRNA_bytype = SampleSplit(brca_mRNA$dat)
brca_mRNA_normal_samples = colnames(brca_mRNA_bytype$normal)
brca_mRNA_normal_samples = substr(brca_mRNA_normal_samples, start = 1, stop = 12)
brca_mRNA_tumor_samples = colnames(brca_mRNA_bytype$primary.tumor); 
brca_mRNA_tumor_samples = substr(brca_mRNA_tumor_samples, start = 1, stop = 12)

length(intersect(brca_mRNA_normal_samples, colnames(brca_mutation_df))) # 99
length(intersect(brca_mRNA_tumor_samples, colnames(brca_mutation_df))) # 754

s1 = intersect(brca_mRNA_normal_samples, colnames(brca_mutation_df))
s2 = intersect(brca_mRNA_tumor_samples, colnames(brca_mutation_df))
length(intersect(s1,s2)) # 96

# compare with the sample-filtered mRNA matrix
source("code_data_prep/helper_functions.R")
brca_mRNA_normal_samples = changeSampleNameVersion(sampleName = colnames(mRNA_normal), original =  ".")
brca_mRNA_tumor_samples = changeSampleNameVersion(sampleName = colnames(mRNA_tumor), original = ".")
s1 = intersect(brca_mRNA_normal_samples, colnames(brca_mutation_df)); length(s1) # 78
s2 = intersect(brca_mRNA_tumor_samples, colnames(brca_mutation_df)); length(s2) # 443
length(intersect(s1,s2)) # 75
  
# ----------------- Explore Copy Number---------------------------------------------------
source("code_data_prep/helper_functions.R")
brca_mRNA_normal_samples = changeSampleNameVersion(sampleName = colnames(mRNA_normal), original =  ".")
brca_mRNA_tumor_samples = changeSampleNameVersion(sampleName = colnames(mRNA_tumor), original = ".")

require(TCGA2STAT)
View(brca_cnasnp$dat[2000:3000,])
brca_cnasnp_sample_split = SampleSplit(brca_cnasnp$dat)

brca_cnasnp_normal = brca_cnasnp_sample_split$normal; dim(brca_cnasnp_normal) # 22618  1113
brca_cnasnp_tumor = brca_cnasnp_sample_split$primary.tumor; dim(brca_cnasnp_tumor) # 22618  1089

brca_cnasnp_normal_samples = changeSampleNameVersion(sampleNames = colnames(brca_cnasnp_normal), 
                                                     original = "-",
                                                     type = "normal")
brca_cnasnp_tumor_samples = changeSampleNameVersion(sampleNames = colnames(brca_cnasnp_tumor), 
                                                     original = "-",
                                                     type = "tumor")

s1 = intersect(colnames(mRNA_normal), brca_cnasnp_normal_samples); length(s1) # 76
s2 = intersect(colnames(mRNA_tumor), brca_cnasnp_tumor_samples); length(s2) # 446

length(intersect(changeSampleNameVersion(s1,original = "."),changeSampleNameVersion(s2,original = "."))) # 72

# ----------------- Test Code ---------------------------------------------------
setwd("~/Desktop")
require(TCGA2STAT)
brca_cn  = read.table(file = "brca_scna_hg19.txt", header = T)
brca_cn_samples = as.character(unique(brca_cn$Sample)); length(brca_cn_samples)

brca_cnasnp_sample_split = SampleSplit(brca_cnasnp$dat)
brca_cnasnp_normal = brca_cnasnp_sample_split$normal; dim(brca_cnasnp_normal) # 22618  1113
brca_cnasnp_tumor = brca_cnasnp_sample_split$primary.tumor; dim(brca_cnasnp_tumor) # 22618  1089

t = append(colnames(brca_cnasnp_normal),colnames(brca_cnasnp_tumor))
intersect(t, brca_cn_samples)
intersect(colnames(brca_cnasnp_normal), brca_cn_samples)
intersect(colnames(brca_cnasnp_tumor), brca_cn_samples)

t1 = brca_cnasnp_tumor[,"TCGA-Z7-A8R6-01A-11D-A41E-01"]
t2 = brca_cn[which(brca_cn$Sample %in% "TCGA-Z7-A8R6-01A-11D-A41E-01"),]
length(t1)
dim(t2)
View(t1)
length(unique(t1))
length(unique(t2$Segment_Mean))
t2[which(t2$Segment_Mean == -0.1365),]
View(brca_cnasnp_tumor[which(t1 == -0.1365),])




