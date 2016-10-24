setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
load("brca2_TCGA2STAT.rda")
load("data_Saved_R_Objects/brca_df.rda")

# get clinical data from downloaded TCGA objects
clinical = brca_miRNA$clinical

# get samples names from processed data
sample_names = colnames(brca_miRNA_df)
# split string by [.]
sample_names = strsplit(sample_names, split = "[.]")
sample_names = unlist(lapply(sample_names, function(name) 
  paste(name[3],name[4],name[5], sep="-")))
# check for duplication
sum(duplicated(sample_names)) 
sum(sample_names %in% rownames(clinical))

# get clinical data
clinical_data = clinical[match(sample_names, rownames(clinical)),]
# temp =  strsplit(rownames(clinical_data), split = "[-]")
# rownames(clinical_data)[1:79] = unlist(lapply(temp[1:79], function(name) paste("BRCA","Normal",name[1],name[2],name[3],sep=".")))
# rownames(clinical_data)[80:536] = unlist(lapply(temp[80:536], function(name) paste("BRCA","Tumor",name[1],name[2],name[3],sep=".")))

clinical_data = clinical_data[which(!duplicated(rownames(clinical_data))),]
dim(clinical_data)
