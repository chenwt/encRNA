setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")

load("data_Saved_R_Objects/brca_df.rda")
# load("Saved_R_Objects/corr_matrices/corr_matrices.rda")
# load("Saved_R_Objects/corr_matrices/normal_tumor_corr_matrices.rda")
load("data_Saved_R_Objects/corr_matrices/all_corr_matrices_with_names.rda")
load("data_Saved_R_Objects/corr_matrices/tumor_normal_lncRNA_mRNA_pair9999.rda")
load("data_Saved_R_Objects/corr_matrices/sensitivity_matricies_9999.rda")
# load("Saved_R_Objects/corr_matrices/sensitivity_matrix.rda")

require(ppcor); require(rlist);
require(foreach); require(doParallel);

# -------------------- Dimension check --------------------------------------------------
dim(lncRNA_mRNA_corr_matrix) # [1]  4828 17613
dim(miRNA_lncRNA_corr_matrix) #[1]  343 4828
dim(miRNA_mRNA_corr_matrix) #[1]   343 17613

# ------------------- Visualization ----------------------------------------------------
plot(density(tumor_lncRNA_mRNA_corr_matrix), 
     col = "red", lwd = 3, ylim = c(0,4), 
     main = "mRNA - lncRNA correlation (normal and tumor)")
lines(density(normal_lncRNA_mRNA_corr_matrix), col = "green", lwd = 3)
legend(x = 'topright', pch = 15,
       col = c("green","red"),
       legend = c("normal", "tumor"))
# ------------------- Get lncRNA mRNA pair ----------------------------------------------

# # get 99% quantile from normal and cancer correlation matrices
# normal_quantile99 = quantile(normal_lncRNA_mRNA_corr_matrix, c(0.99))
# # 99% 
# # 0.710914 
# tumor_quantile99 = quantile(tumor_lncRNA_mRNA_corr_matrix, c(0.99))
# # 99% 
# # 0.3476711 
# 
# # select top pairs lncRNA-mRNA from normal dataset
# normal_lncRNA_mRNA_pairs99 = get_lncRNA_mRNA_pairs(normal_lncRNA_mRNA_corr_matrix, 
#                                                  normal_quantile99)
# dim(normal_lncRNA_mRNA_pairs99)
# # [1] 850356      5 --> too many
# 
# tumor_lncRNA_mRNA_pairs99 = get_lncRNA_mRNA_pairs(tumor_lncRNA_mRNA_corr_matrix, 
#                                                   tumor_quantile99)
# dim(tumor_lncRNA_mRNA_pairs99)
# # [1] 850356      5 --> similarly too many, thus apply a more stringent 

# thus, make the threshold for normal sample more stringent --> 99.99%
normal_quantile9999 = quantile(normal_lncRNA_mRNA_corr_matrix, c(0.9999))
# 99.99% 
# 0.8866987 
# select top pairs lncRNA-mRNA from normal dataset
normal_lncRNA_mRNA_pairs9999 = get_lncRNA_mRNA_pairs(normal_lncRNA_mRNA_corr_matrix, 
                                                 normal_quantile9999)
dim(normal_lncRNA_mRNA_pairs9999)
# [1] 8504    5

tumor_quantile9999 = quantile(tumor_lncRNA_mRNA_corr_matrix, c(0.9999))
# 99.99% 
# 0.694144 
tumor_lncRNA_mRNA_pairs9999 = get_lncRNA_mRNA_pairs(tumor_lncRNA_mRNA_corr_matrix, 
                                                  tumor_quantile9999)
dim(tumor_lncRNA_mRNA_pairs9999) 
# [1] 8504  5

save(tumor_lncRNA_mRNA_pairs9999, normal_lncRNA_mRNA_pairs9999, 
     file = "Saved_R_Objects/corr_matrices/tumor_normal_lncRNA_mRNA_pair9999.rda")

# --------------- build sensitivity matrix from normal cells------------------------------

# get only subset expression data from normal cells 
normal_indices = 1:79
mRNA_normal = brca_mRNA_df[,normal_indices]; 
miRNA_normal = brca_miRNA_df[,normal_indices];
lncRNA_normal = brca_lncRNA_df[,normal_indices]

# get sensitivity matrix for normal samples 
pairs = normal_lncRNA_mRNA_pairs9999
normal_sensitivity_matrix9999 = matrix(ncol = 343)

# expect to run in 1.5 hour
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
  normal_sensitivity_matrix9999 <<- rbind(normal_sensitivity_matrix9999, sensitivity_corr_vector)
}
ptm = proc.time() - ptm

dim(normal_sensitivity_matrix9999)
normal_sensitivity_matrix9999 = normal_sensitivity_matrix9999[-1,]

normal_lncRNA_mRNA_names = paste(normal_lncRNA_mRNA_pairs9999$lncRNA, normal_lncRNA_mRNA_pairs9999$mRNA, sep = "-")
rownames(normal_sensitivity_matrix9999) = normal_lncRNA_mRNA_names 
colnames(normal_sensitivity_matrix9999) = rownames(brca_miRNA_df)

# --------------- build sensitivity matrix from cancer cells------------------------------
# get only subset expression data from cancer cells 
tumor_indices = 80:ncol(brca_mRNA_df)
mRNA_tumor = brca_mRNA_df[,tumor_indices]; 
miRNA_tumor = brca_miRNA_df[,tumor_indices];
lncRNA_tumor = brca_lncRNA_df[,tumor_indices]

# get sensitivity matrix for normal samples 
pairs = tumor_lncRNA_mRNA_pairs9999
tumor_sensitivity_matrix9999 = matrix(ncol = 343)

# expect to run in 1.5 hour
ptm = proc.time()
for (i in 1:nrow(pairs)){
  the_pair = pairs[i,]
  mRNA_vector = mRNA_tumor[the_pair$mRNA_index,]
  lncRNA_vector = lncRNA_tumor[the_pair$lncRNA_index,]
  sensitivity_corr_vector = c()
  for (j in 1:nrow(miRNA_tumor)){
    partial_corr = pcor.test(mRNA_vector, lncRNA_vector, miRNA_tumor[j,])$estimate
    sensitivity_corr_vector = append(sensitivity_corr_vector, the_pair$corr - partial_corr)
  }
  tumor_sensitivity_matrix9999 <<- rbind(tumor_sensitivity_matrix9999, sensitivity_corr_vector)
}
ptm = proc.time() - ptm

dim(tumor_sensitivity_matrix9999) # [1] 8505  343
tumor_sensitivity_matrix9999 = tumor_sensitivity_matrix9999[-1,]

tumor_lncRNA_mRNA_names = paste(tumor_lncRNA_mRNA_pairs9999$lncRNA, tumor_lncRNA_mRNA_pairs9999$mRNA, sep = "-")
rownames(tumor_sensitivity_matrix9999) = tumor_lncRNA_mRNA_names 
colnames(tumor_sensitivity_matrix9999) = rownames(brca_miRNA_df)

save(normal_sensitivity_matrix9999, tumor_sensitivity_matrix9999, 
     file = "Saved_R_Objects/corr_matrices/sensitivity_matricies_9999.rda")


# ------------ visualazation of senstivity matrix ------------------------------------

plot(density(tumor_sensitivity_matrix9999), 
     lwd = 3, col = "red",
     main = "sensitivity matrix")
lines(density(normal_sensitivity_matrix9999),
      lwd = 3, col = "green")
legend(x = 'topright', pch = 15,
       col = c("green","red"),
       legend = c("normal", "tumor"))

require(gplots)

png(file = "Saved_R_Objects/corr_matrices/heatmap_normal.png")
heatmap.2(normal_sensitivity_matrix9999, 
          col=redgreen(75), dendrogram = "none", trace = 'none',
          labRow = F, colRow = F)
dev.off()  ## 10Kb file in current dir

heatmap.2(tumor_sensitivity_matrix9999, 
          col=redgreen(75), dendrogram = "none", trace = 'none')

# ------------ extract triple from sensitivity matrix ------------------------------------

## NORMAL SENSITIVITY

head(normal_lncRNA_mRNA_pairs9999)
# lncRNA_index mRNA_index            lncRNA      mRNA      corr
# 893           907       1686 ENSG00000229645.4 C14orf139 0.9997353
# 4898         3657       9931 ENSG00000265142.2      MYH4 0.9992861
# 6870          162      14252 ENSG00000203875.6     SNHG5 0.9992630
quantile(normal_sensitivity_matrix9999,c(0.99))
# 99% 
# 0.1272952
normal_encRNA = get_encRNA(matrix = normal_sensitivity_matrix9999,
                           miRNA_mRNA_corr = normal_miRNA_mRNA_corr_matrix,
                           miRNA_lncRNA_corr = normal_miRNA_lncRNA_corr_matrix,
                           threshold = quantile(normal_sensitivity_matrix9999,c(0.99))) 
length(which(normal_encRNA$lncRNA_miRNA_corr < 0 & normal_encRNA$mRNA_miRNA_corr < 0)) # 1793
length(which(normal_encRNA$lncRNA_miRNA_corr < 0 & normal_encRNA$mRNA_miRNA_corr < 0)) / nrow(normal_encRNA) * 100
# 6.146937

dim(normal_encRNA) # [1] 29169     8
head(normal_encRNA,3)
#                   lncRNA mRNA        miRNA      sensitivity    encRNA_pair                      encRNA_triple encRNA_pair_index miRNA_index
# 6573  ENSG00000203499.6  DSP hsa-mir-200a 0.6741329  ENSG00000203499.6-DSP  ENSG00000203499.6-DSP-hsa-mir-1-2              5800         118
# 13048 ENSG00000254148.3 CDH1 hsa-mir-200c 0.6631619 ENSG00000254148.3-CDH1 ENSG00000254148.3-CDH1-hsa-mir-1-2              4420         120
# 6088  ENSG00000203499.6 CDH1 hsa-mir-200a 0.6547107 ENSG00000203499.6-CDH1 ENSG00000203499.6-CDH1-hsa-mir-1-2              4943         118
length(unique(normal_encRNA$mRNA)) # 1456
length(unique(normal_encRNA$lncRNA)) # 200
length(unique(normal_encRNA$miRNA)) # 33
length(unique(normal_encRNA$encRNA_pair)) # 5010

## TUMOR SENSITIVITY 
quantile(tumor_sensitivity_matrix9999, c(0.99))
# 99% 
# 0.07240933 
tumor_encRNA = get_encRNA(matrix = tumor_sensitivity_matrix9999,
                           miRNA_mRNA_corr = tumor_miRNA_mRNA_corr_matrix,
                           miRNA_lncRNA_corr = tumor_miRNA_lncRNA_corr_matrix,
                           threshold = quantile(tumor_sensitivity_matrix9999,c(0.99))) 

length(which(tumor_encRNA$lncRNA_miRNA_corr < 0 & tumor_encRNA$mRNA_miRNA_corr < 0)) # 593
length(which(tumor_encRNA$lncRNA_miRNA_corr < 0 & tumor_encRNA$mRNA_miRNA_corr < 0)) / nrow(tumor_encRNA) * 100
# 2.033

length(unique(tumor_encRNA$mRNA)) # 593
length(unique(tumor_encRNA$lncRNA)) # 124
length(unique(tumor_encRNA$miRNA)) # 44
length(unique(tumor_encRNA$encRNA_pair)) # 5249

## INTERSECTION
intersect(unique(tumor_encRNA$miRNA), unique(normal_encRNA$miRNA))
# [1] "hsa-mir-196a-1" "hsa-mir-135b"   "hsa-mir-146a"   "hsa-mir-139"   
# [5] "hsa-mir-378"    "hsa-mir-22"     "hsa-mir-452"    "hsa-mir-200c"  
# [9] "hsa-mir-146b"   "hsa-mir-141"    "hsa-mir-375"    "hsa-mir-224"   
intersect(unique(tumor_encRNA$encRNA_pair), unique(normal_encRNA$encRNA_pair))
# [1] "ENSG00000231367.1-ROPN1"    "ENSG00000249042.1-MLPH"    
# [3] "ENSG00000260807.2-ROPN1"    "ENSG00000228639.2-GABRP"   
# [5] "ENSG00000254148.3-ROPN1B"   "ENSG00000245812.2-AQP7"    
# [7] "ENSG00000243350.1-GATA3"    "ENSG00000232638.1-GATA3"   
# [9] "ENSG00000254148.3-GABRP"    "ENSG00000229645.4-CHRDL1"  
# [11] "ENSG00000245812.2-KCNIP2"   "ENSG00000245812.2-NPR1"    
# [13] "ENSG00000228639.2-ROPN1B"   "ENSG00000245812.2-CHRDL1"  
# [15] "ENSG00000245812.2-ITIH5"    "ENSG00000249042.1-CCDC125" 
# [17] "ENSG00000245812.2-AOC3"     "ENSG00000245812.2-TMEM132C"
# [19] "ENSG00000245812.2-LOC90586" "ENSG00000228639.2-MIA"     
# [21] "ENSG00000254148.3-MIA"      "ENSG00000245812.2-CAV1"    
# [23] "ENSG00000236333.3-ADIPOQ"   "ENSG00000236333.3-LOC90586"
# [25] "ENSG00000228639.2-SOX10"    "ENSG00000227591.1-GPD1"    
# [27] "ENSG00000254148.3-SOX10"    "ENSG00000227591.1-RDH5"    
# [29] "ENSG00000236333.3-CIDEA"    "ENSG00000236333.3-AQP7" 
intersect(unique(tumor_encRNA$encRNA_triple), unique(normal_encRNA$encRNA_triple))
character(0)

save(normal_encRNA, tumor_encRNA, file = "Saved_R_Objects/corr_matrices/normal_tumor_encRNA_triplets.rda" )

# normal_encRNA = get_encRNA(matrix = normal_sensitivity_matrix9999,
#                            threshold = 0.3) 
# dim(normal_encRNA) # [1] 3783    8

# ------------ check with putative target binding database -----------------------------

load("mircode_objects.rda") # mircode; mircode_lncRNA -- from 
load("starbase_mRNA_miRNA_interactions.rda") # from starBase_Human_Pan-Cancer_mirna_mRNA_interaction.csv

##### NORMAL SAMPLES ######

## putative lncRNA - miRNA interaction from miRcode 
length(unique(mircode_lncRNA$gene_id)) # [1] 10349
length(intersect(unique(mircode_lncRNA$gene_id), unique(normal_encRNA$lncRNA))) # 81
length(intersect(unique(mircode_lncRNA$gene_id), unique(tumor_encRNA$lncRNA))) # 98

# question: for the given lncRNAs found in the triple RNAs, look for what miRNAs binding to it from the reference datatabse (miRcode)
normal_overlapped_lncRNAs = intersect(unique(mircode_lncRNA$gene_id), unique(normal_encRNA$lncRNA))
df = mircode_lncRNA[which(mircode_lncRNA$gene_id %in% normal_overlapped_lncRNAs), 
                    c("gene_id","microrna")]
putative_miRNA_lncRNAs = get_putative_lncRNA_miRNA(df)
normal_encRNA_subset = normal_encRNA[which(normal_encRNA$lncRNA %in% normal_overlapped_lncRNAs),]
# now, find overlap between sequence-based predicted miRNAs with sensitivity-based miRNAs
length(intersect(unique(putative_miRNA_lncRNAs$lncRNA_miRNA_pair),unique(normal_encRNA_subset$lncRNA_miRNA_pair))) # 106

## putative mRNA - miRNA interaction from starbase
colnames(starbase_mrna_mirna_interaction)[1:2] = c("miRNA", "mRNA")
starbase_mrna_mirna_interaction$mRNA_miRNA_pair = paste(starbase_mrna_mirna_interaction$mRNA,
                                                        starbase_mrna_mirna_interaction$miRNA,
                                                        sep = "-")
tumor_overlapped_mRNAs = intersect(unique(starbase_mrna_mirna_interaction$mRNA), unique(normal_encRNA$mRNA))
length(tumor_overlapped_mRNAs) # 1098
length(intersect(unique(starbase_mrna_mirna_interaction$mRNA_miRNA_pair),unique(normal_encRNA_subset$mRNA_miRNA_pair))) # 106

processed_starbase_mrna_mirna_interaction = 


##### TUMOR SAMPLES ######



# ------------ check correlation signs (miRNA-mRNA and miRNA-lncRNA) ------------------

View(normal_encRNA)

normal_miRNA_lncRNA_corr_matrix[which(rownames(normal_miRNA_lncRNA_corr_matrix) == "hsa-mir-200c"),
                              which(colnames(normal_miRNA_lncRNA_corr_matrix) == "ENSG00000261183.1")]
normal_miRNA_mRNA_corr_matrix[which(rownames(normal_miRNA_mRNA_corr_matrix) == "hsa-mir-200c"),
                              which(colnames(normal_miRNA_mRNA_corr_matrix) == "SERPINB5")]

test = normal_encRNA[1:10, c("mRNA", "miRNA", "lncRNA")]
test$k = normal_miRNA_lncRNA_corr_matrix[as.character(test$miRNA),as.character(test$lncRNA)]

for (i in 1:nrow(normal_encRNA)){ 
  print(i)
  print(normal_miRNA_lncRNA_corr_matrix[as.character(normal_encRNA$miRNA[i]),
                                        as.character(normal_encRNA$lncRNA[i]) ]); 
}
# ------------- helper functions ------------------------------------------------------

get_lncRNA_mRNA_pairs = function(originalCorrMatrix, corrThreshold){
  correlation_pairs = which(originalCorrMatrix > corrThreshold, arr.ind = TRUE)
  lncRNA = rownames(originalCorrMatrix)[correlation_pairs[,1]]
  mRNA = colnames(originalCorrMatrix)[correlation_pairs[,2]]
  dataframe = as.data.frame(cbind(correlation_pairs, lncRNA, mRNA))
  colnames(dataframe) = c("lncRNA_index", "mRNA_index","lncRNA", "mRNA")
  dataframe$lncRNA_index = as.numeric(as.character(dataframe$lncRNA_index))
  dataframe$mRNA_index = as.numeric(as.character(dataframe$mRNA_index))
  # correlation_vector = normal_lncRNA_mRNA_corr_matrix[dataframe$lncRNA_index,dataframe$mRNA_index]
  dataframe = cbind(dataframe, originalCorrMatrix[which(originalCorrMatrix > corrThreshold)])
  colnames(dataframe) = c("lncRNA_index", "mRNA_index","lncRNA", "mRNA", "corr")
  dataframe = dataframe[order(-dataframe$corr),]
  return(dataframe)
}

matrix = normal_sensitivity_matrix9999;
miRNA_mRNA_corr = normal_miRNA_mRNA_corr_matrix
miRNA_lncRNA_corr = normal_miRNA_lncRNA_corr_matrix
threshold = quantile(normal_sensitivity_matrix9999,c(0.99))

get_encRNA = function(matrix, miRNA_mRNA_corr, miRNA_lncRNA_corr, threshold){
  triple = which(matrix > threshold, arr.ind = TRUE)
  encRNAs = rownames(matrix)[triple[,1]]
  miRNA = colnames(matrix)[triple[,2]]
  dataframe = as.data.frame(cbind(triple, encRNAs, miRNA))
  colnames(dataframe) = c("encRNA_pair_index", "miRNA_index","encRNA_pair", "miRNA")
  dataframe$encRNA_pair_index = as.numeric(as.character(dataframe$encRNA_pair_index))
  dataframe$miRNA_index = as.numeric(as.character(dataframe$miRNA_index))
  dataframe = cbind(dataframe, matrix[which(matrix > threshold)])
  colnames(dataframe) = c("encRNA_pair_index", "miRNA_index","encRNA_pair", "miRNA", "sensitivity")
  dataframe = dataframe[order(-dataframe$sensitivity),]
  # dataframe$lncRNA = substr(dataframe$encRNA_pair, start = 1, stop = 17)
  # dataframe$mRNA = substr(dataframe$encRNA_pair, start = 19, stop = nchar(as.character(dataframe$encRNA_pair)))
  
  pos = regexpr("-",as.character(dataframe$encRNA_pair))
  dataframe$lncRNA = substr(as.character(dataframe$encRNA_pair), start = 1, stop = pos - 1)
  dataframe$mRNA = substr(as.character(dataframe$encRNA_pair), start = pos + 1, stop = nchar(as.character(dataframe$encRNA_pair))) 
  
  dataframe$lncRNA_miRNA_corr = rep(NA, nrow(dataframe))
  dataframe$mRNA_miRNA_corr = rep(NA, nrow(dataframe))
  
  for (i in 1:nrow(dataframe)){
    dataframe$lncRNA_miRNA_corr[i] = miRNA_lncRNA_corr[as.character(dataframe$miRNA[i]),as.character(dataframe$lncRNA[i])]
    dataframe$mRNA_miRNA_corr[i] = miRNA_mRNA_corr[as.character(dataframe$miRNA[i]),as.character(dataframe$mRNA[i])]
  }
  
  dataframe$encRNA_triple = paste(dataframe$encRNA_pair, miRNA, sep = "-")
  dataframe$lncRNA_miRNA_pair = paste(dataframe$lncRNA, miRNA, sep = "-")
  dataframe$mRNA_miRNA_pair = paste(dataframe$mRNA, miRNA, sep = "-")
  dataframe = dataframe[,c("lncRNA", "mRNA", "miRNA", "sensitivity", 
                           "lncRNA_miRNA_corr", "mRNA_miRNA_corr",
                           "encRNA_pair", "encRNA_triple", 
                           "lncRNA_miRNA_pair" , "mRNA_miRNA_pair", 
                           "encRNA_pair_index", "miRNA_index")]
  return(dataframe)
}



# for each lncRNA ensemble ids, find potential miRNA interaction
getMiRNAs = function(miRNA_family = NULL){
  if (is.null(miRNA_family)){
    print("must provide miRNA_family")
    break;
  }
  part1 = substr(miRNA_family, start = 1, stop = 3)
  part1 = paste("hsa",tolower(part1),sep = "-")
  
  # substring anything after miR till the end, divided by "/"
  part2 = substr(miRNA_family,start = 5, stop = nchar(miRNA_family))
  part2 = unlist(strsplit(x = part2, split = "/"))
  
  # foreach element, remove 3p and 5p parts
  part2 = gsub("-3p","",part2)
  part2 = gsub("-5p","",part2)
  
  # return individual mircRNA, 
  # example: 106abc will be disconstructed into 106a, 106b, 106c
  part2 =  sapply(part2, function(element){
    if (grepl("\\D",element)){
      digit_part = gsub(pattern = "\\D", replacement = "", x = element)
      character_parts = gsub(pattern = "\\d", replacement = "", x = element)
      character_parts = unlist(strsplit(x = character_parts,split = ""))
      returned_value = paste(digit_part, character_parts,sep = "")
    }else{
      element
    }
  })
  part2 = unname(unlist(part2))
  return(paste(part1,part2,sep="-"))
}

get_putative_lncRNA_miRNA = function(dataframe){
  l = list()
  apply(dataframe, 1, function(r){
    k = getMiRNAs(r[2])
    l <<- list.append(l, k)
  })
  names(l) = dataframe$gene_id
  df = reshape2::melt(l)
  colnames(df) = c("miRNA", "putative_lncRNAs")
  df$lncRNA_miRNA_pair = paste(df$putative_lncRNAs, df$miRNA, sep = "-")
  return(df)
}





