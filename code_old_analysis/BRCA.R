setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")

library(TCGAbiolinks)
library(TCGA2STAT)


## ---------------- Utility function ---------------------------

numOverlappedSamples <- function(cancerType = "brca"){
  query <- TCGAquery(tumor = cancerType, level = 3)
  matSamples <- TCGAquery_integrate(query)
  methylation450_index <- grep(x = colnames(matSamples), pattern = "HumanMethylation450")
  methylation27_index <- grep(x = colnames(matSamples), pattern = "HumanMethylation27")
  rnaSeqV2_index <- grep(x = colnames(matSamples), pattern = "IlluminaHiSeq_RNASeqV2")
  miRNASeq_index <- grep(x = colnames(matSamples), pattern = "IlluminaHiSeq_miRNASeq")
  return(matSamples[c(rnaSeqV2_index,methylation450_index,methylation27_index,miRNASeq_index),
                    c(rnaSeqV2_index,methylation450_index,methylation27_index,miRNASeq_index)])
}

reformatTCGASampleName = function(givenVector, startIndex = 6, stopIndex = 12, 
                                  oldPattern = "-", newPattern = "."){
  givenVector = substr(givenVector, start = startIndex, stop = stopIndex); 
  givenVector = gsub(oldPattern,newPattern,givenVector)
  return(givenVector)
}


####################### Called from TCGA2Biolinks ######################
brca_samples = numOverlappedSamples(cancerType = "brca"); brca_samples;
#                                       IlluminaHiSeq_RNASeqV2 HumanMethylation450 HumanMethylation27 IlluminaHiSeq_miRNASeq
# IlluminaHiSeq_RNASeqV2                   1218                 790                314                    758
# HumanMethylation450                       790                 930                  9                    629
# HumanMethylation27                        314                   9                349                    134
# IlluminaHiSeq_miRNASeq                    758                 629                134                    869
#----------- BRCA DATA --------------------
brca_mRNA_datQuery <- TCGAquery(tumor = "BRCA",
                                platform = "IlluminaHiSeq_RNASeqV2",
                                level = "3")
brca_mRNA_lsSample <- TCGAquery_samplesfilter(query = brca_mRNA_datQuery);
brca_mRNA_dataSmTP <- TCGAquery_SampleTypes(barcode = brca_mRNA_lsSample$IlluminaHiSeq_RNASeqV2, typesample = "TP")
brca_mRNA_dataSmNT <- TCGAquery_SampleTypes(barcode = brca_mRNA_lsSample$IlluminaHiSeq_RNASeqV2, typesample ="NT")
length(brca_mRNA_dataSmTP) # 1097
length(brca_mRNA_dataSmNT) # 114

brca_miRNA_datQuery <- TCGAquery(tumor = "BRCA",
                                 platform = "IlluminaHiSeq_miRNASeq",
                                 level = "3")
brca_miRNA_lsSample <- TCGAquery_samplesfilter(query = brca_miRNA_datQuery);
brca_miRNA_dataSmTP <- TCGAquery_SampleTypes(barcode = brca_miRNA_lsSample$IlluminaHiSeq_miRNASeq, typesample = "TP")
brca_miRNA_dataSmNT <- TCGAquery_SampleTypes(barcode = brca_miRNA_lsSample$IlluminaHiSeq_miRNASeq, typesample ="NT")
length(brca_miRNA_dataSmTP) # 775
length(brca_miRNA_dataSmNT) # 87


brca_methyl_datQuery <- TCGAquery(tumor = "BRCA",
                                  platform = "HumanMethylation27",
                                  level = "3")
brca_methyl_lsSample <- TCGAquery_samplesfilter(query = brca_methyl_datQuery);
brca_methyl_dataSmTP <- TCGAquery_SampleTypes(barcode = brca_methyl_lsSample$HumanMethylation27, typesample = "TP")
brca_methyl_dataSmNT <- TCGAquery_SampleTypes(barcode = brca_methyl_lsSample$HumanMethylation27, typesample ="NT")
length(brca_methyl_dataSmTP) # 316
length(brca_methyl_dataSmNT) # 27

####################### Called from TCGA2STAT ##############
########
brca_rnaV2 <- getTCGA(disease="BRCA", data.type="RNASeq2", clinical = TRUE)
dim(brca_rnaV2$dat) # [1] 20501   1212
brca_miRNA <- getTCGA(disease="BRCA", data.type="miRNASeq", type="count")
dim(brca_miRNA$dat) # [1] 1046   849
brca_methyl <- getTCGA(disease="BRCA", data.type="Methylation", type="27K")
dim(brca_methyl$dat) # [1] 27578   343

brca_rnaV2_sName <- reformatTCGASampleName(colnames(brca_rnaV2$dat), startIndex = 1, stopIndex = 12,
                                           oldPattern = "-", newPattern = "-")
brca_miRNA_sName <- reformatTCGASampleName(colnames(brca_miRNA$dat), startIndex = 1, stopIndex = 12,
                                           oldPattern = "-", newPattern = "-")
brca_methyl_sName <- reformatTCGASampleName(colnames(brca_methyl$dat), startIndex = 1, stopIndex = 12,
                                            oldPattern = "-", newPattern = "-")

length(intersect(brca_rnaV2_sName, brca_miRNA_sName)) # 754
length(intersect(brca_rnaV2_sName, brca_methyl_sName)) # 312
length(intersect(brca_miRNA_sName, brca_methyl_sName)) # 132

save(brca_rnaV2_sName, brca_miRNA_sName,brca_methyl_sName, file = "brca_TCGA2STAT_shortname.rda" )

# Split the OMICs data by sample types
brca.rnaseq2.bytype <- SampleSplit(brca_rnaV2$dat)
brca.miRNA.bytype <- SampleSplit(brca_miRNA$dat)
brca.methyl.bytype <- SampleSplit(brca_methyl$dat)

brca.rnaV2_tumor_sName = reformatTCGASampleName( colnames(brca.rnaseq2.bytype$primary.tumor), startIndex = 1, stopIndex = 12,
                                                 oldPattern = "-", newPattern = "-")
length(brca.rnaV2_tumor_sName) # 1093
brca.methyl_tumor_sName = reformatTCGASampleName( colnames(brca.methyl.bytype$primary.tumor), startIndex = 1, stopIndex = 12,
                                                  oldPattern = "-", newPattern = "-")  
length(brca.methyl_tumor_sName) # 314
brca.miRNA_tumor_sName = reformatTCGASampleName( colnames(brca.miRNA.bytype$primary.tumor), startIndex = 1, stopIndex = 12,
                                                 oldPattern = "-", newPattern = "-") 
length(brca.miRNA_tumor_sName) # 755

brca.rnaV2_normal_sName = reformatTCGASampleName( colnames(brca.rnaseq2.bytype$normal), startIndex = 1, stopIndex = 12,
                                                  oldPattern = "-", newPattern = "-")
length(brca.rnaV2_normal_sName) # 112
brca.methyl_normal_sName = reformatTCGASampleName( colnames(brca.methyl.bytype$normal), startIndex = 1, stopIndex = 12,
                                                   oldPattern = "-", newPattern = "-")  
length(brca.methyl_normal_sName) # 27
brca.miRNA_normal_sName = reformatTCGASampleName( colnames(brca.miRNA.bytype$normal), startIndex = 1, stopIndex = 12,
                                                  oldPattern = "-", newPattern = "-") 
length(brca.miRNA_normal_sName) # 87

length(intersect(brca.rnaV2_tumor_sName, brca.miRNA_tumor_sName)) # 753
length(intersect(brca.methyl_tumor_sName, brca.methyl_tumor_sName)) # 314
length(intersect(brca.miRNA_tumor_sName, brca.methyl_tumor_sName)) # 132

# matched paired-sampl
length(intersect(brca.rnaV2_tumor_sName, brca.rnaV2_normal_sName)) # 112
length(intersect(brca.methyl_tumor_sName, brca.methyl_normal_sName)) # 27 --> good!
length(intersect(brca.miRNA_tumor_sName, brca.miRNA_normal_sName)) # 86 --> good!

# tumor samples having all mRNA, miRNA, and methylation
brca.rnaV2_miRNA_tumor_sName = intersect(brca.rnaV2_tumor_sName, brca.miRNA_tumor_sName); 
length(brca.rnaV2_miRNA_tumor_sName) # 753
length(intersect(brca.rnaV2_miRNA_tumor_sName,brca.methyl_tumor_sName)) # 132 --> very good!
