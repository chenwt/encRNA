## --------------- BASIC SETUP ---------------------------------------

# remove other variables if needed
rm(list = ls()); gc()

# set directory
setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")

# load original brca_lncRNA if wanted
load("original_brca_lncRNA.rda"); # load original_brca_lncRNA: 12727 * 943

# Load processed brca_lncRNA data (12246 rows) and also load the original ensemblID
load("brca_lncRNA.rda"); # load processed brca_lncRNA: 12246 * 942
load("original_ensemblID.rda") # load original_ensemblID: 12727
# Load annotation data
load("lncRNA_extAnnot.rda") # load lncRNA_extAnnot: 12246 * 7
# Load differentially expression lncRNA 
load("limma_all.rda") # load limma_all: 12246 * 7
# Load mircode object 
load("mircode.rda") # load: mircode, mircode_lncRNA, de_lncRNA_miRNA
# mircode: 1329071 * 12
# mircode_lncRNA: 138883 * 12 --> only includes overlapping lncRNA and intergenic lncRNA
# de_lncRNA_miRNA: differentially expressed lncRNA which also includes in the miRcode database --> maybe we use try better database 
# which have more interactions

## --------------- Load libraries ---------------------------------------
library(TCGA2STAT) 
library(dplyr)
library(limma)

## --------------- DNA methylation, mRNA, miRNA analysis --------------------------
brca.mRNA_methyl <- OMICSBind(dat1 = brca_rnaV2$dat, dat2 = brca_methyl$dat)

dim(brca_rnaV2$dat) # 20501  1212
dim(brca_rnaV2$clinical) # 1097   18

brca_mRNA_PT = brca.rnaseq2.bytype$primary.tumor; dim(brca_mRNA_PT) # [1] 20501  1093
brca_mRNA_NT = brca.rnaseq2.bytype$normal; dim(brca_mRNA_NT) # [1] 20501   112

head(colnames(brca_mRNA_PT));
head(brca.rnaV2_tumor_sName); 
length(brca.rnaV2_normal_sName) 
head(rownames(brca_rnaV2$clinical))

brca.rnaV2_tumor_sName = reformatTCGASampleName( colnames(brca.rnaseq2.bytype$primary.tumor), startIndex = 1, stopIndex = 12,
                                                 oldPattern = "-", newPattern = "-")

# check methylation

# brca.rnaseq2.bytype <- SampleSplit(brca_rnaV2$dat)
# brca.miRNA.bytype <- SampleSplit(brca_miRNA$dat)
# brca.methyl.bytype <- SampleSplit(brca_methyl$dat)

dim(brca_methyl$dat) # [1] 27578   343
length(brca.methyl_tumor_sName) # 314
length(brca.methyl_normal_sName)# 27
head(brca.methyl_normal_sName)
# [1] "TCGA-BH-A0B7" "TCGA-BH-A0BL" "TCGA-BH-A0BO" "TCGA-BH-A0BQ" "TCGA-BH-A0BW" "TCGA-BH-A0DE"
head(brca.methyl_tumor_sName)
# [1] "TCGA-A1-A0SD" "TCGA-A2-A04N" "TCGA-A2-A04P" "TCGA-A2-A04Q" "TCGA-A2-A04T" "TCGA-A2-A04U"
length(intersect(brca.methyl_tumor_sName,brca.methyl_normal_sName)) # 27
head(colnames(brca_methyl$dat))

brca.tumor_index = grep("-01A-", x = colnames(brca_methyl$dat))
length(brca.tumor_index) # 310
brca.normal_index = grep("-11A-", x = colnames(brca_methyl$dat))
length(brca.normal_index) # 27

## --------------- Methylation Analysis --------------------------

# Density plot of normal vs tumor indices
i=brca.tumor_index[1]
plot(density(na.omit(brca_methyl$dat[,i]),from=0,to=1),
     main="",ylim=c(0,15),type="n")

for(i in brca.tumor_index){
  lines(density(na.omit(brca_methyl$dat[,i]),from=0,to=1),col=2)
}

for(i in brca.normal_index){
  lines(density(na.omit(brca_methyl$dat[,i]),from=0,to=1),col=1)
}

# remove sample labeled as recurrent tumor
brca_methyl_data = as.data.frame(brca_methyl$dat[,-grep("-01B-", x = colnames(brca_methyl$dat))])
dim(brca_methyl_data) # [1] 27578   339

# create a temp pdata so that we can perform t-test to find methylated region
pdata_brca = data.frame(row.names = colnames(brca_methyl_data))
pdata_brca$Status = ifelse(grepl("-01A", x = rownames(pdata_brca)), "Cancer", "Normal")
pdata_brca$Status = as.factor(pdata_brca$Status)

library(limma)
X<-model.matrix(~pdata_brca$Status )
fit<-lmFit(brca_methyl_data,X)
eb <- ebayes(fit)

library(rafalib)

splot(fit$coef[,2],-log10(eb$p.value[,2]),xlab="Effect size",ylab="-log10 p-value")

View(brca_methyl$cpgs)

# To do
# Big task 1: create GRange object from brca_methyl$cpgs, need to remove NA data if needed
# Big task 2: perform regional methylated 

## --------------- Read data from TANRIC lncRNA (BRCA) --------------------------

brca_lncRNA = read.csv(file = "lncRNA_Tanric_062616.csv", header = TRUE)
head(brca_lncRNA)[1:3]
dim(brca_lncRNA) # 12727 *  943 samples

## --------------- get annotation from bioMart -----------------
##### annotation with biomaRT
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

listMarts()
####### select a dataset
# connect to specific BioMart database, this must be a valid name given by listMarts 
ensembl = useMart("ensembl")
# BiomMarts databasese can contain several datasets, for Ensembl every species is a 
# different datasset
listDatasets(ensembl)

# use a dataset
ensembl = useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
# ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

###### Build a query
# attributes defines the value that we are interested in retrieving 
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

lncRNA_ensemblID = brca_lncRNA[,1]
length(lncRNA_ensemblID) # [1] 12727 : total number of lncRNA in the dataset
# remove decimal part of the ensembl Id
lncRNA_ensemblID = substr(lncRNA_ensemblID, start = 1, stop = 15); 
# lncRNA_basicAnnot = getBM(attributes = c("chromosome_name","hgnc_symbol","ensembl_gene_id", "clone_based_vega_transcript_name"),
#              filters = "ensembl_gene_id",
#              values = lncRNA_ensemblID, 
#              mart = ensembl)
# lncRNA_basicAnnot$hgnc_symbol[(lncRNA_basicAnnot$hgnc_symbol == "")] = NA
# dim(lncRNA_basicAnnot) # [1] 12246     3 ==> thus, there is 12727 - 12246 = 481 lncRNA that 
# # ensemble database appears to not be able to find annotation
# sum(is.na(lncRNA_basicAnnot$hgnc_symbol)) # [1] 9800: ==>? out of 12246 annotated lncRNA, 9800 of them
# # does not have annotated hgnc_symbol

original_brca_lncRNA = brca_lncRNA;


lncRNA_extAnnot = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","external_gene_name",
                                       "chromosome_name","start_position","end_position","strand"),
                          filters = "ensembl_gene_id",
                          values = lncRNA_ensemblID, 
                          mart = ensembl)
dim(lncRNA_extAnnot); # [1] 12246     7
head(lncRNA_extAnnot)
# ensembl_gene_id hgnc_symbol external_gene_name chromosome_name start_position end_position strand
# 1 ENSG00000005206      SPPL2B             SPPL2B              19        2328615      2354806      1
# 2 ENSG00000083622                     AC000111.6               7      117604791    117647415     -1
# 3 ENSG00000088970         KIZ                KIZ              20       21125983     21246622      1
# 4 ENSG00000099869     IGF2-AS            IGF2-AS              11        2140501      2148666      1
# 5 ENSG00000100181      TPTEP1             TPTEP1              22       16601887     16698742      1
# 6 ENSG00000104691       UBXN8              UBXN8               8       30732247     30767006      1

length(lncRNA_extAnnot$external_gene_name) # 12246
length(unique(lncRNA_extAnnot$external_gene_name)) # 12239
duplicated_cases = which(duplicated(lncRNA_extAnnot$external_gene_name))
lncRNA_extAnnot$external_gene_name[duplicated_cases]
# [1] "PGM5-AS1"      "LINC01347"     "PABPC1L2B-AS1" "LINC01422"     "PRICKLE2-AS1"  "LINC01481"     "PROX1-AS1" 
d1 = which(lncRNA_extAnnot$external_gene_name == "PGM5-AS1") # [1]  958 1171
d2 = which(lncRNA_extAnnot$external_gene_name == "LINC01347") # [1]  426 2922
d3 = which(lncRNA_extAnnot$external_gene_name == "PABPC1L2B-AS1") #[1] 1478 3069
d4 = which(lncRNA_extAnnot$external_gene_name == "LINC01422") # [1]  595 4057
d5 = which(lncRNA_extAnnot$external_gene_name == "PRICKLE2-AS1") # [1] 5124 5171
d6 = which(lncRNA_extAnnot$external_gene_name == "LINC01481") # [1] 8203 8253
d7 = which(lncRNA_extAnnot$external_gene_name == "PROX1-AS1") # [1] [1]  2604 11761

d1 = lncRNA_extAnnot$ensembl_gene_id[d1] # [1] "ENSG00000224958" "ENSG00000225655"
d2 = lncRNA_extAnnot$ensembl_gene_id[d2] # [1] "ENSG00000214837" "ENSG00000231512"
d3 = lncRNA_extAnnot$ensembl_gene_id[d3] # [1] "ENSG00000226725" "ENSG00000231963"
d4 = lncRNA_extAnnot$ensembl_gene_id[d4] # [1] "ENSG00000223704" "ENSG00000235271"
d5 = lncRNA_extAnnot$ensembl_gene_id[d5] # [1] "ENSG00000241111" "ENSG00000241572"
d6 = lncRNA_extAnnot$ensembl_gene_id[d6] # [1] "ENSG00000257613" "ENSG00000257815"
d7 = lncRNA_extAnnot$ensembl_gene_id[d7] # [1] "ENSG00000230461" "ENSG00000272167"

duplicated_id = c(d1,d2,d3,d4,d5,d6,d7)

## --------------- reassign the name for variable brca_lncRNA -----------------
# reassign the name for variable brca_lncRNA
dim(lncRNA_extAnnot) # [1] 12246     7
dim(brca_lncRNA) # [1] 12727   943

# check what id is gin lncRNA_extAnnot but not in brca_lncRNA
length(lncRNA_extAnnot$ensembl_gene_id) # 12246
length(lncRNA_ensemblID) # 12727
length(intersect(lncRNA_extAnnot$ensembl_gene_id, lncRNA_ensemblID)) # 12246
brca_lncRNA$Gene_ID = lncRNA_ensemblID; length(brca_lncRNA$Gene_ID) # 12727
selected_indices = which(brca_lncRNA$Gene_ID %in% lncRNA_extAnnot$ensembl_gene_id)

# subsetting brca_lncRNA to contain only gene having symbol annotation --> POTENTIAL PITFALL here, 
# WORTH REVISITING IT IN THE FUTURE!!!

brca_lncRNA = brca_lncRNA[selected_indices,];
dim(brca_lncRNA) # [1] 12246   943
# [1] 12246

# re-name rows
rownames(brca_lncRNA) = brca_lncRNA$Gene_ID

# remove gene_id column
brca_lncRNA = brca_lncRNA[,-1]
brca_lncRNA[1:4,1:4]

# save processed brca_lncRNA
#save(brca_lncRNA, file = "brca_lncRNA.rda")

## --------------- DEA on 12246 lncRNAs -----------------
reduced_names = substr(colnames(brca_lncRNA), start = 13, stop = 24)
status = substr(colnames(brca_lncRNA), start = 13, stop = 24)
normal_col_indices = grep(colnames(brca_lncRNA), pattern = ".Normal.")
tumor_col_indices = grep(colnames(brca_lncRNA), pattern = ".Tumor.")

status = seq(0, length = ncol(brca_lncRNA))
for (i in 1:ncol(brca_lncRNA)){
  status[i] = unlist(strsplit(colnames(brca_lncRNA)[i], split = "[.]"))[2]
}

require(limma)
mod = model.matrix(~ as.factor(status))
fit_limma = lmFit(brca_lncRNA, mod)
ebayes_limma = eBayes(fit_limma)
limma_all = topTable(ebayes_limma, number = dim(brca_lncRNA)[1])

# note that limma_all is all sorted by pvalue with the first row as the smallest  

hist(limma_all$t, col = 4)
names(limma_all)
dim(limma_all) # [1] 12246     6
head(limma_all)
#                 logFC       AveExpr         t       P.Value     adj.P.Val        B
# ENSG00000241684 -0.48491118 0.10079111 -40.49784 2.345821e-208 2.095620e-204 466.2363
# ENSG00000234456 -2.79028164 1.11936252 -40.47226 3.422538e-208 2.095620e-204 465.8586
# ENSG00000264016 -0.59261984 0.21133252 -39.41198 2.291504e-201 9.353919e-198 450.1461
# ENSG00000254862 -0.24386254 0.05623317 -39.13744 1.364670e-199 4.177937e-196 446.0603
# ENSG00000256508 -0.08697416 0.01736069 -37.27244 1.856850e-187 4.547798e-184 418.1297
# ENSG00000262097 -0.44456487 0.11464853 -36.55601 9.172257e-183 1.872058e-179 407.3254

fp_bh_limma = p.adjust(limma_all$P.Value, method = "BH")
length(fp_bh_limma) # 12426
hist(fp_bh_limma, col = 4)
quantile(fp_bh_limma)

sum(fp_bh_limma < 0.05) # [1] 5046
# actually fp_bh_limma is exactly simmilar to limma_all$adj.P.Val: 
# sum(fp_bh_limma - limma_all$adj.P.Val) # return 0

original_de_brca_lncRNA = rownames(limma_all)[1:5046]
#save(limma_all, de_brca_lncRNA, file = "limma_all.rda")

# TO-DO: plot taking account effect size
intersect(original_de_brca_lncRNA, de_brca_lncRNA) 
which(original_de_brca_lncRNA %in% intersect(duplicated_id, de_brca_lncRNA) )
# [1] "ENSG00000224958" "ENSG00000225655" "ENSG00000223704" "ENSG00000235271" "ENSG00000241572" "ENSG00000257613" "ENSG00000257815"
# [8] "ENSG00000230461" "ENSG00000272167"

# remove those cases
de_brca_lncRNA= original_de_brca_lncRNA[-which(original_de_brca_lncRNA %in% intersect(duplicated_id, de_brca_lncRNA))]
length(de_brca_lncRNA) # 1] 5037

# convert to gene symbol
de_brca_lncRNA = lncRNA_extAnnot$external_gene_name[lncRNA_extAnnot$ensembl_gene_id %in% de_brca_lncRNA]

#save(limma_all, de_brca_lncRNA, file = "limma_all.rda")

## --------------- get interaction from miRcode ---------------------------------

mircode = read.table(file = "mircode_highconsfamilies.txt", header = TRUE, sep = "\t")
# save(mircode, file = "mircode.rda")
head(mircode)

length(unique(mircode$gene_symbol)) # 47458
unique(mircode$gene_class)
length(unique(mircode$gene_class))
table(mircode$gene_class)[1]
# coding 
# 1018974 
table(mircode$gene_class)[2]
# lncRNA, intergenic 
# 74646 
table(mircode$gene_class)[3]
# lncRNA, overlapping 
# 64237 
table(mircode$gene_class)[4]
# other 
# 40266 
table(mircode$gene_class)[5]
# pseudogene 
# 130948 

mircode_lncRNA_row_indices = which(mircode$gene_class %in% c("lncRNA, overlapping", "lncRNA, intergenic"))
mircode_lncRNA = mircode[mircode_lncRNA_row_indices,]

head(mircode_lncRNA,1)
length(unique(mircode_lncRNA$gene_symbol)) # 10349

intersect(mircode_lncRNA$gene_symbol, de_brca_lncRNA)

# we have 5037 differentially expressed lncRNA, we have 10349 unique lncRNA from mrcode
length(intersect(unique(mircode_lncRNA$gene_symbol), de_brca_lncRNA)) # only 293
de_lncRNA_miRNA = intersect(mircode_lncRNA$gene_symbol, de_brca_lncRNA)
save(mircode, mircode_lncRNA, de_lncRNA_miRNA, file = "mircode.rda")

# try intersect with the whole unique gene symbol in mrcode to see if we have more intersected elements thn 293
length(intersect(unique(mircode$gene_symbol), de_brca_lncRNA)) # 407, just little better than 293

#should try out with ensemble instead of symbolized de_brca_lncRNA
length(intersect(unique(mircode$gene_id), unique(original_ensemblID))) # 8741 --> look promising!




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