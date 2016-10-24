#### --------- LOAD SAVED OBJECTS ----------------------------------------------
setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
#load("data_Starbase_miRcode_files/starBase_Human_Pan-Cancer_miRNA-LncRNA_Interactions2016-07-20_18-03.csv")
  
#### --------- IMPORT DATA FROM STARBASE V2.0 ----------------------------------

# ## lncRNA-miRNA interaction
# lncrna_mirna_interaction = read.csv( file = "starBase_Human_Pan-Cancer_miRNA-LncRNA_Interactions2016-07-20_18-03.csv", 
#                                      sep = ",", 
#                                      header = TRUE)
# head(lncrna_mirna_interaction,2)
# # name mirAccession      geneName targetSites bioComplex clipReadNum cancerNum
# # 1 hsa-miR-200b-3p MIMAT0000318   CTA-204B4.6           1          1           1        -1
# # 2 hsa-miR-200b-3p MIMAT0000318 CTD-2630F21.1           1          8         190        -1
# 
# ## mRNA-miRNA interaction
# mrna_mirna_interaction = read.csv(file = "starBase_Human_Pan-Cancer_mirna_mRNA_interaction.csv", 
#                                     sep = ",", 
#                                     header = TRUE)
# head(mrna_mirna_interaction,2)
# # name                  geneName targetScanSites picTarSites RNA22Sites PITASites miRandaSites CancerNum   X X.1 X.2   X.3 X.4
# # 1 hsa-miR-200b-3p      DEK             0[0          0]        0[0        0]          0[0        0] 0[0  0] 1[5 1384]   3
# # 2 hsa-miR-200b-3p     PELO             0[0          0]        0[0        0]          0[0        0] 0[0  0] 1[1    6]  10
# 
# # because of tab seperation problem, we need to do some cleanup
# temp_targetScanSites = mrna_mirna_interaction$targetScanSites
# temp_picTarSites = mrna_mirna_interaction$picTarSites
# temp_RNA22Sites = mrna_mirna_interaction$RNA22Sites
# temp_PITASites = mrna_mirna_interaction$PITASites
# temp_miRandaSites = mrna_mirna_interaction$miRandaSites
# temp_CancerNum = mrna_mirna_interaction$CancerNum
# temp_X = mrna_mirna_interaction$X
# temp_X1 = mrna_mirna_interaction$X.1
# temp_X2 = mrna_mirna_interaction$X.2
# temp_X3 = mrna_mirna_interaction$X.3
# temp_X4 = mrna_mirna_interaction$X.4
# 
# mrna_mirna_interaction$targetScanSites = paste(temp_targetScanSites,",",temp_picTarSites, sep = "")
# mrna_mirna_interaction$picTarSites = paste(temp_RNA22Sites, ",", temp_PITASites, sep = "")
# mrna_mirna_interaction$RNA22Sites = paste(temp_miRandaSites, ",", temp_CancerNum, sep = "")
# mrna_mirna_interaction$PITASites = paste(temp_X, ",", temp_X1, sep="")
# mrna_mirna_interaction$miRandaSites = paste(temp_X2, ",", temp_X3, sep= "")
# mrna_mirna_interaction$CancerNum = temp_X4
# 
# mrna_mirna_interaction =mrna_mirna_interaction[,-c(9,10,11,12,13)] 
# head(mrna_mirna_interaction,2)
# 
# save(lncrna_mirna_interaction, mrna_mirna_interaction, file = "starbase_lncRNA_objects.rda")

#### --------- COMPARE WITH miRcode ----------------------------------

test_miRNA = "miR-199"

starbase_indices = grep(pattern = test_miRNA, x = lncrna_mirna_interaction$name)
length(starbase_indices) #75
starbase_mRNA_target = lncrna_mirna_interaction[starbase_indices,]$geneName

mircode_indices = grep(pattern = test_miRNA, x = mircode_lncRNA$microrna)
length(mircode_indices) # [1] 2314
mircode_mRNA_target = mircode_lncRNA[mircode_indices,]$gene_symbol
  
length(intersect(starbase_mRNA_target,mircode_mRNA_target) )
# [1] "XIST"     "SNHG1"    "GAS5"     "FGD5-AS1"

# overall intersection of total lncRNA from two datasets
starbase_lncRNA = unique(lncrna_mirna_interaction$geneName) # 1127
mircode_lncRNA = unique(mircode_lncRNA$gene_symbol) # 10349

length(intersect(starbase_lncRNA,mircode_lncRNA)) # [1] 89

# look at mRNA-mirna database, extract all the genes
starbase_mRNA = unique(mrna_mirna_interaction$geneName) #13802
length(intersect(starbase_mRNA,starbase_lncRNA)) # 50
length(intersect(starbase_mRNA,mircode_lncRNA)) # [1] 19








  
  
  