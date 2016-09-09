## --------------- BASIC SETUP -----------------------------------------------------------

# remove other variables if needed  
rm(list = ls()); gc()

# set directory
setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")

# load libraries
require(TCGA2STAT); require(limma)

# load de_brca_mRNA_starbase_df, de_brca_lncRNA_mircode_df
load("brca_df2.rda")

load("de_brca_miRcode_starbase.rda") 
dim(de_brca_lncRNA_mircode_df) # [1] 2100  536
dim(de_brca_mRNA_starbase_df) # [1] [1] 9389  536

load("de_brca_enRNAs_candidate_list.rda")

# Saved objects in this scripts
save(de_brca_mRNA_starbase_common_miRNA_df, de_brca_lncRNA_mircode_common_miRNA_df,
     file = "de_brca_lncRNA_mRNA_common_miRNA_dfs.rda" )
save(de_brca_lncRNA_cor_matrix, de_brca_mRNA_cor_matrix, 
     file = "de_brca_lncRNA_mRNA_cor_matrices.rda")
save(brca_miRNA_df,de_brca_mRNA_starbase_common_miRNA_df,
     de_brca_lncRNA_mircode_common_miRNA_df,
     file = "de_brca_mRNA_lncRNA_common_miRNA_Aug0816.rda")

## --------------- SUBSET ----------------------------------------------------------------

# length(de_brca_enRNAs_candidate_mRNAs_list) # 27 --> 27 common miRNA
de_brca_enRNAs_candidate_mRNAs_df = reshape2::melt(de_brca_enRNAs_candidate_mRNAs_list)
de_brca_enRNAs_candidate_mRNAs_name = unique(de_brca_enRNAs_candidate_mRNAs_df$value)
de_brca_enRNAs_candidate_lncRNAs_df = reshape2::melt(de_brca_enRNAs_candidate_lncRNAs_list)
de_brca_enRNAs_candidate_lncRNAs_name = unique(de_brca_enRNAs_candidate_lncRNAs_df$value)
# subset de_brca_lncRNA_mircode_df and de_brca_mRNA_starbase_df
de_brca_mRNA_starbase_common_miRNA_df = subset(de_brca_mRNA_starbase_df, 
                                                rownames(de_brca_mRNA_starbase_df) %in% 
                                                 de_brca_enRNAs_candidate_mRNAs_name)
de_brca_lncRNA_mircode_common_miRNA_df = subset(de_brca_lncRNA_mircode_df, 
                                                rownames(de_brca_lncRNA_mircode_df) %in% 
                                                  de_brca_enRNAs_candidate_lncRNAs_name)

save(de_brca_mRNA_starbase_common_miRNA_df, de_brca_lncRNA_mircode_common_miRNA_df, 
     file = "de_brca_lncRNA_mRNA_common_miRNA_dfs.rda" )


## --------------- VISUALIZATION----------------------------------------------------------

heatmap(as.matrix(de_brca_mRNA_starbase_common_miRNA_df)[1:6,1:4])

require(gplots)
heatmap.2(as.matrix(de_brca_mRNA_starbase_common_miRNA_df)[1:6,1:4], col=topo.colors(75), scale="none",
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)


d <- density(de_brca_mRNA_starbase_common_miRNA_df[,1])
plot(d,main="Density plot of differentially expressed genes",ylim=c(0,.45),xlab="Gene expression", ylab="Density")
for (s in 2:ncol(de_brca_mRNA_starbase_common_miRNA_df)){
  d = density(de_brca_mRNA_starbase_common_miRNA_df[,s])
  lines(d)
}

## --------------- CORRELATION ----------------------------------------------------------

dim(de_brca_mRNA_starbase_common_miRNA_df) # [1] 6451  536
dim(de_brca_lncRNA_mircode_common_miRNA_df) # [1] 1838  536

# find pair 

de_brca_lncRNA_cor_matrix = cor(t(de_brca_lncRNA_mircode_common_miRNA_df))
de_brca_mRNA_cor_matrix = cor(t(de_brca_mRNA_starbase_common_miRNA_df))
save(de_brca_lncRNA_cor_matrix, de_brca_mRNA_cor_matrix, 
     file = "de_brca_lncRNA_mRNA_cor_matrices.rda")




#################### reorder samples across 3 tables to have the same sample ordering

colnames(de_brca_mRNA_starbase_common_miRNA_df)[1:10]
colnames(de_brca_lncRNA_mircode_common_miRNA_df)[1:10]
colnames(brca_miRNA_df)[1:10]

brca.miRNA.bytype = SampleSplit(brca_miRNA$dat)
brca_miRNA_normal_samples = colnames(brca.miRNA.bytype$normal) 
brca_miRNA_tumor_samples = colnames(brca.miRNA.bytype$primary.tumor) 

idx = which(colnames(brca_miRNA_df) %in% brca_miRNA_normal_samples)
temp = colnames(brca_miRNA_df)[idx]
temp = strsplit(temp, split = "[-]")
temp = unlist(lapply(temp, function(name) paste("BRCA","Normal",name[1],name[2],name[3],sep=".")))
colnames(brca_miRNA_df)[idx] = temp
idx = which(colnames(brca_miRNA_df) %in% brca_miRNA_tumor_samples)
temp = colnames(brca_miRNA_df)[idx]
temp = strsplit(temp, split = "[-]")
temp = unlist(lapply(temp, function(name) paste("BRCA","Tumor",name[1],name[2],name[3],sep=".")))
colnames(brca_miRNA_df)[idx] = temp

brca_mRNA_bytype = SampleSplit(brca_mRNA$dat)
brca_mRNA_normal_samples = substr(colnames(brca_mRNA_bytype$normal),1,28); 
length(brca_mRNA_normal_samples) # 100
brca_mRNA_tumor_samples = substr(colnames(brca_mRNA_bytype$primary.tumor),1,28) # 1093
length(brca_mRNA_tumor_samples) # 775

idx = which(colnames(de_brca_mRNA_starbase_common_miRNA_df) %in% brca_mRNA_normal_samples)
temp = colnames(de_brca_mRNA_starbase_common_miRNA_df)[idx]
temp = strsplit(temp, split = "[-]")
temp = unlist(lapply(temp, function(name) paste("BRCA","Normal",name[1],name[2],name[3],sep=".")))
colnames(de_brca_mRNA_starbase_common_miRNA_df)[idx] = temp
idx = which(colnames(de_brca_mRNA_starbase_common_miRNA_df) %in% brca_mRNA_tumor_samples)
temp = colnames(de_brca_mRNA_starbase_common_miRNA_df)[idx]
temp = strsplit(temp, split = "[-]")
temp = unlist(lapply(temp, function(name) paste("BRCA","Tumor",name[1],name[2],name[3],sep=".")))
colnames(de_brca_mRNA_starbase_common_miRNA_df)[idx] = temp

brca_miRNA_df = brca_miRNA_df[,colnames(de_brca_lncRNA_mircode_common_miRNA_df)]
dim(brca_miRNA_df)
brca_common_miRNA = brca_miRNA_df[rownames(brca_miRNA_df) %in% names(de_brca_enRNAs_candidate_mRNAs_list),]
de_brca_mRNA_starbase_common_miRNA_df = 
  de_brca_mRNA_starbase_common_miRNA_df[,colnames(de_brca_lncRNA_mircode_common_miRNA_df)]
dim(de_brca_mRNA_starbase_common_miRNA_df)

save(brca_miRNA_df,de_brca_mRNA_starbase_common_miRNA_df,
     de_brca_lncRNA_mircode_common_miRNA_df,
     file = "de_brca_mRNA_lncRNA_common_miRNA_Aug0816.rda")

#################### testing ###########################

names(de_brca_enRNAs_candidate_mRNAs_list)
which(rownames(brca_miRNA_df) %in% names(de_brca_enRNAs_candidate_mRNAs_list))
# [1]  20  92  99 132 216 241

names(de_brca_enRNAs_candidate_mRNAs_list)[3] # "hsa-mir-761"
idx = which(names(de_brca_enRNAs_candidate_lncRNAs_list) %in% "hsa-mir-761")

lncRNAs = de_brca_enRNAs_candidate_lncRNAs_list[[idx]]
mRNAs = de_brca_enRNAs_candidate_mRNAs_list[[3]]

lncRNAs[1] # "ENSG00000117242.7"
mRNAs[1] # PPME1

load("de_brca_mRNA_lncRNA_common_miRNA_Aug0816.rda")
dim(brca_miRNA_df) # [1] 343 536

mRNA_PPME1 = de_brca_mRNA_starbase_common_miRNA_df[rownames(de_brca_mRNA_starbase_common_miRNA_df) ==
                                                     "PPME1",]
lncRNA_ENSG00000117242.7 = de_brca_lncRNA_mircode_common_miRNA_df[
  rownames(de_brca_lncRNA_mircode_common_miRNA_df) == "ENSG00000117242.7",]

miRNA_hsa_mir_761 = brca_miRNA_df[rownames(brca_miRNA_df) == "hsa-mir-761",]




