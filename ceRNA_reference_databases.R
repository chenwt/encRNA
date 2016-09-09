setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")

## ----- miRSponge ----------------------------------------

mirSponge_prediction_miRNA_target = read.table(file = "ceRNA_TargetDatabases/miRSponge/Prediction_of_miRNA_targets_in_miRSponge.txt", 
                                               sep = "\t",  header = FALSE, fill = TRUE)
head(mirSponge_prediction_miRNA_target, 4)
# V1    V2                          V3                         V4                          V5    V6    V7  V8     V9
# 1 let-7a  VCAN  3 uugauaugUUGGAUGAUGGAGu 5             ||:| |||||:||   5 atttgctgAATCAACTACTTCc 3  5424  5445 142 -14.10
# 2 let-7a IGF1R 3 uuGAUAUGUUGGAU-GAUGGAGu 5      | | || |::|| |||||||  5 ccCCAAAC-ATTTATCTACCTCa 3  7749  7770 153 -17.10
# 3 let-7a IGF1R 3 uugauaUGUUGGAU-GAUGGAGu 5          | |:::|: |||||||  5 tttgccAGAGTTTGTCTACCTCt 3 11790 11812 147 -14.54
# 4 let-7a IGF1R  3 uugauauguuggaUGAUGGAGu 5                  :|||||||   5 acggtgatcgaggGCTACCTCc 3  1200  1221 141 -14.40

mirSponge_experimentally_validated_ceRNA = 
  read.csv(file = "ceRNA_TargetDatabases/miRSponge/Experimentally_validated_miRNA_sponges_and_ceRNAs_in_miRSponge.csv",
           sep = ",", header = TRUE, fill = TRUE) 
head(mirSponge_experimentally_validated_ceRNA)
View(mirSponge_experimentally_validated_ceRNA)

mirSponge_experimentally_validated_miRNA_target = read.table(
  file = "ceRNA_TargetDatabases/miRSponge/Experimentally_validated_miRNA_targets_in_miRSponge.txt", 
  sep = "\t",  header = TRUE, fill = TRUE)
head(mirSponge_experimentally_validated_ceRNA)
View(mirSponge_experimentally_validated_ceRNA)

## ----- lnCeDb ----------------------------------------

# miRNA-lncRNA,mRNA -- to be checked!
lnCeDb_mir_lnc_parclip72 = read.table(file = "ceRNA_TargetDatabases/lnCeDb/mir_lnc_parclip72.txt", 
                                      sep = "\t",  header = FALSE)
head(lnCeDb_mir_lnc_parclip72,3)
# V1                V2 V3                 V4
# 1 hsa-miR-34b-5p ENST00000505030.1  2 7-mer-A1,7-mer-A1,
# 2 hsa-miR-34b-5p ENST00000500197.2  2 7-mer-A1,7-mer-A1,
# 3 hsa-miR-34b-5p ENST00000504246.1  2 7-mer-A1,7-mer-A1,

# miRNA-lncRNA,mRNA -- to be checked!
lnCeDb_mir_lnc_targets_all = read.table(file = "ceRNA_TargetDatabases/lnCeDb/mir_lnc_targets_all.txt", 
                                        sep = "\t",  header = FALSE)

head(lnCeDb_mir_lnc_targets_all,2)
# V1                V2              V3                           V4
# 1 hsa-miR-17-5p ENST00000395900.1  4 7-mer-m8,7-mer-m8,7-mer-m8,7-mer-m8,
# 2 hsa-miR-17-5p ENST00000421848.1  2                   7-mer-m8,7-mer-m8,

lnCeDb_miR_mRNA_lnc_ceRNA = scan(file = "ceRNA_TargetDatabases/lnCeDb/miR_mRNA_lnc_ceRNA.txt", 
                                 what = "", sep = "\n")

lnCeDb_miR_mRNA_lnc_ceRNA = strsplit(lnCeDb_miR_mRNA_lnc_ceRNA, ";")
miRNA_names = c()

# the returned lnCeDb_miR_mRNA_lnc_ceRNA is a list, whose names of list's element are the miRNA names
lnCeDb_miR_mRNA_lnc_ceRNA = lapply(lnCeDb_miR_mRNA_lnc_ceRNA, function(miRNA_group){
  miRNA_group_vector = unlist(strsplit(miRNA_group, "\t"))
  miRNA_names <<- append(miRNA_names, miRNA_group_vector[1]) # <<- for calling global namess
  miRNA_group_vector[-1]
  # names(miRNA_group) = miRNA_group[1]
  # miRNA_group = miRNA_group[-1] 
})
names(lnCeDb_miR_mRNA_lnc_ceRNA) = miRNA_names

## ----- LncACTdb ----------------------------------------

lncACTdb_cancer_associated = read.table(file = "ceRNA_TargetDatabases/LncACTdb/Cancer_associated_lncACTs.txt", 
                                        sep = "\t",  header = TRUE)

dim(lncACTdb_cancer_associated) # [1] 878   5
head(lncACTdb_cancer_associated,2)
#     lncRNAEnsgID   lncRNAName       miRNA Genename Disease
# 1 ENSG00000124915 DKFZP434K028 hsa-miR-328     CD44    kirp
# 2 ENSG00000124915 DKFZP434K028 hsa-miR-145  C11orf9    kirp

# miRNA - mRNA
lncACTdb_experimental_validated_miRNA_targets = read.table(file = "ceRNA_TargetDatabases/LncACTdb/Experimental_validated_miRNA_targets.txt", 
                                                           sep = "\t",  header = FALSE, fill = TRUE)
dim(lncACTdb_experimental_validated_miRNA_targets) # [1] 43497     3

head(lncACTdb_experimental_validated_miRNA_targets)
# V1              V2    V3
# 1    hsa-let-7 ENSG00000149948 HMGA2
# 2    hsa-let-7 ENSG00000009307  NRAS

# main table
lncACTdb_functionally_activated_lncACTS = read.table(file = "ceRNA_TargetDatabases/LncACTdb/Functionally_activated_lncACTS.txt", 
                                                     sep = "\t",  header = TRUE)
dim(lncACTdb_functionally_activated_lncACTS) # 5119 
head(lncACTdb_functionally_activated_lncACTS)
# lncRNA.Ensembl.ID    lncRNA.Name        miRNA Gene.Name Gene.Ensembl.ID
# 1   ENSG00000100181         TPTEP1 hsa-miR-133b    BCL2L2 ENSG00000129473
# 2   ENSG00000100181         TPTEP1 hsa-miR-133b     PTBP2 ENSG00000117569

lncACTdb_lncRNA_datasource = read.table(file = "ceRNA_TargetDatabases/LncACTdb/Long_noncoding_RNAs_datasource.gtf.txt", 
                                        sep = ";",  header = TRUE)
head(lncACTdb_lncRNA_datasource)

save(mirSponge_prediction_miRNA_target, mirSponge_experimentally_validated_ceRNA, mirSponge_experimentally_validated_miRNA_target,
     lnCeDb_mir_lnc_parclip72, lnCeDb_mir_lnc_targets_all, lnCeDb_miR_mRNA_lnc_ceRNA,
     lncACTdb_cancer_associated, lncACTdb_experimental_validated_miRNA_targets, lncACTdb_functionally_activated_lncACTS,
     file = "ceRNA_TargetDatabases/ceRNA_databases.rda")

## ----- miRCode ----------------------------------------
load("mircode_objects.rda")
head(mircode)
#     gene_id             gene_symbol gene_class       microrna      seed_pos seed_type repeat. total_cons_. primates_cons_. mammals_cons_. vertebrates_cons_.    tr_region
# 1 ENSG00000000003.9      TSPAN6     coding miR-132/212/212-3p chrX:99884666     8-mer       0           13              67              0                  0       3pUTR,
# 2 ENSG00000000003.9      TSPAN6     coding         miR-133abc chrX:99884907  7-mer-A1       0           11              33              9                  0       3pUTR,

## ------------- Original data ----------------------------------------------------------------

load("brca_common_objects2.rda") # brca_lncRNA_common, brca_miRNA_common, brca_mRNA_common (before gene filtering, after sample filtering)
nrow(brca_mRNA_common) # 20501
rownames(brca_mRNA_common)[1:10] # [1] "A1BG"   "A1CF"   "A2BP1"  "A2LD1"  "A2ML1"  "A2M"    "A4GALT" "A4GNT"  "AAA1"   "AAAS"  
nrow(brca_miRNA_common) # 1046
rownames(brca_miRNA_common)[10:15] # [1] "hsa-let-7g"    "hsa-let-7i"    "hsa-mir-1-1"   "hsa-mir-1-2"   "hsa-mir-100"   "hsa-mir-101-1"
nrow(brca_lncRNA_common) # 12727
rownames(brca_lncRNA_common)[1:5] # 1] "ENSG00000005206.12" "ENSG00000031544.10" "ENSG00000083622.8"  "ENSG00000088970.11" "ENSG00000099869.6" 

draw.single.venn(area = 22, category = "Test")

load("brca_df2.rda") # after low-count filter and log transform, load: brca_lncRNA_df, brca_miRNA_df, brca_mRNA_df
nrow(brca_mRNA_df) # 16669
nrow(brca_miRNA_df) # 343
nrow(brca_lncRNA_df) # 4828



## ------------- Analysis ----------------------------------------------------------------

########  miRSponge ######################################################################
library(VennDiagram)
library(biomaRt)
load("ceRNA_TargetDatabases/ceRNA_databases.rda")

##### mirSponge_prediction_miRNA_target
head(mirSponge_prediction_miRNA_target)
dim(mirSponge_prediction_miRNA_target) #[1] 220603      9
# V1: miRNA, V2: gene
length(unique(mirSponge_prediction_miRNA_target$V1)) # 687 miRNAs
length(unique(mirSponge_prediction_miRNA_target$V2)) # 173 mRNA or lncRNAs

plot.new(); frame()

# interaction of mRNA from brca_mRNA_common with mRNAs included in mirSponge_prediction_miRNA_target
draw.pairwise.venn(area1 = nrow(brca_mRNA_common), 
                   area2 = length(unique(mirSponge_prediction_miRNA_target$V2)),
                   length(intersect(rownames(brca_mRNA_common),unique(mirSponge_prediction_miRNA_target$V2))),
                   category= c("brca_mRNA", "genes or lncRNAs from mirSponge_prediction_miRNA_target"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0))

# use biomart to annotate lncRNA in brca_lncRNA_common
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl = useMart("ensembl")
ensembl = useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

brca_lncRNA_common_annotated = substr(rownames(brca_lncRNA_common), start = 1, stop = 15)
brca_lncRNA_common_annotated = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","external_gene_name",
                                                    "chromosome_name","start_position","end_position","strand"),
                                     filters = "ensembl_gene_id",
                                     values = brca_lncRNA_common_annotated, 
                                     mart = ensembl)
length(unique(brca_lncRNA_common_annotated$external_gene_name)) # [1] 12148 < 12727
brca_lncRNA_common_annotated_unique = unique(brca_lncRNA_common_annotated$external_gene_name)

# intersection of lncRNA from annotated brca_lncRNA_common with lncRNAs included in mirSponge_prediction_miRNA_target
draw.pairwise.venn(area1 = length(brca_lncRNA_common_annotated_unique), 
                   area2 = length(unique(mirSponge_prediction_miRNA_target$V2)),
                   length(intersect(brca_lncRNA_common_annotated_unique,unique(mirSponge_prediction_miRNA_target$V2))),
                   category= c("brca_lncRNA", "lncRNA from mirSponge_prediction_miRNA_target"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0))

# intersection of miRNA from brca_miRNA_common with miRNA from mirSponge_prediction_miRNA_target
mirSponge_prediction_miRNA_target_miRNA = tolower(paste("hsa-",mirSponge_prediction_miRNA_target$V1,
                                                        sep = ""))
draw.pairwise.venn(area1 = nrow(brca_miRNA_common), 
                   area2 = length(unique(mirSponge_prediction_miRNA_target_miRNA)),
                   length(intersect(rownames(brca_miRNA_common),unique(mirSponge_prediction_miRNA_target_miRNA))),
                   category= c("brca_miRNA", "miRNAfrom mirSponge_prediction_miRNA_target"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0))

######## mirSponge_experimentally_validated_miRNA_target
head(mirSponge_experimentally_validated_miRNA_target)

# V1: miRNA.name, V3: Target.Name
length(unique(mirSponge_experimentally_validated_miRNA_target$miRNA.Name)) # 490
length(unique(mirSponge_experimentally_validated_miRNA_target$Target.Name)) # 1695

# interaction of mRNA from brca_mRNA_common with mRNAs included in mirSponge_experimentally_validated_miRNA_target
draw.pairwise.venn(area1 = nrow(brca_mRNA_common), 
                   area2 = length(unique(mirSponge_experimentally_validated_miRNA_target$Target.Name)),
                   length(intersect(rownames(brca_mRNA_common),unique(mirSponge_experimentally_validated_miRNA_target$Target.Name))),
                   category= c("brca_mRNA", "genes or lncRNAs from mirSponge_experimentally_validated_miRNA_target"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0)) # 1411
# intersection of lncRNA from annotated brca_lncRNA_common with lncRNAs included in mirSponge_prediction_miRNA_target
draw.pairwise.venn(area1 = length(brca_lncRNA_common_annotated_unique), 
                   area2 = length(unique(mirSponge_experimentally_validated_miRNA_target$Target.Name)),
                   length(intersect(brca_lncRNA_common_annotated_unique,unique(mirSponge_experimentally_validated_miRNA_target$Target.Name))),
                   category= c("brca_lncRNA", "lncRNA from mirSponge_experimentally_validated_miRNA_target"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0)) # 2
# intersection of miRNA from brca_miRNA_common with miRNA from mirSponge_prediction_miRNA_target
mirSponge_experimentally_validated_miRNA_target_miRNA = tolower(paste("hsa-",mirSponge_experimentally_validated_miRNA_target$miRNA.Name,
                                                                      sep = ""))
draw.pairwise.venn(area1 = nrow(brca_miRNA_common), 
                   area2 = length(unique(mirSponge_experimentally_validated_miRNA_target_miRNA)),
                   length(intersect(rownames(brca_miRNA_common),unique(mirSponge_experimentally_validated_miRNA_target_miRNA))),
                   category= c("brca_miRNA", "miRNAfrom mirSponge_prediction_miRNA_target"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0)) # 151

######## mirSponge_experimentally_validated_ceRNA (To do)


########  lnCeDb #########################################################################

####### lnCeDb_mir_lnc_parclip72 #########
head(lnCeDb_mir_lnc_parclip72) 
# V1: miRNA -- Example: hsa-miR-34b-5p
# V2: lncRNA -- Example: ENST00000505030.1


# check intersection lncRNA from brca_lncRNA_common with lncRNA from lnCeDb_mir_lnc_parclip72 (to do)
lnCeDb_mir_lnc_parclip72$V2[1:4]
rownames(brca_lncRNA_common)[1:3]

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl = useMart("ensembl")
ensembl = useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

nchar(as.character(lnCeDb_mir_lnc_parclip72$V2[1]))
lnCeDb_mir_lnc_parclip72_lncRNAs =substr(lnCeDb_mir_lnc_parclip72$V2, start = 1, stop = 15)
test = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","external_gene_name",
                            "chromosome_name","start_position","end_position","strand"),
             filters = "ensembl_gene_id",
             values = unique(lnCeDb_mir_lnc_parclip72_lncRNAs), 
             mart = ensembl)



# check intersection miRNA from brca_miRNA_common with miRNA from lnCeDb_mir_lnc_parclip72 
lnCeDb_mir_lnc_parclip72_miRNA = tolower(lnCeDb_mir_lnc_parclip72$V1)
length(unique(tolower(lnCeDb_mir_lnc_parclip72$V1))) # 41 distince miRNAs from lnCeDb_mir_lnc_parclip72
draw.pairwise.venn(area1 = nrow(brca_miRNA_common), 
                   area2 = length(unique(lnCeDb_mir_lnc_parclip72_miRNA)),
                   length(intersect(rownames(brca_miRNA_common),unique(lnCeDb_mir_lnc_parclip72_miRNA))),
                   category= c("brca_miRNA", "miRNAfrom lnCeDb_mir_lnc_parclip72"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0)) # 24


####### lnCeDb_mir_lnc_targets_all #########

head(lnCeDb_mir_lnc_targets_all)
# V1: miRNAs, example hsa-miR-17-5p ; 
# V2: lncRNA transcripts, example: ENST00000395900.1

# check intersection miRNA from brca_miRNA_common with miRNA from lnCeDb_mir_lnc_targets_all 
lnCeDb_mir_lnc_targets_all_miRNAs = tolower(lnCeDb_mir_lnc_targets_all$V1)
length(unique(lnCeDb_mir_lnc_targets_all_miRNAs)) # 2042
draw.pairwise.venn(area1 = nrow(brca_miRNA_common), 
                   area2 = length(unique(lnCeDb_mir_lnc_targets_all_miRNAs)),
                   length(intersect(rownames(brca_miRNA_common),unique(lnCeDb_mir_lnc_targets_all_miRNAs))),
                   category= c("brca_miRNA", "miRNAfrom lnCeDb_mir_lnc_targets_all_miRNAs"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0)) # 486 


# check intersection lncRNA or mRNA from brca_lncRNA_common with lncRNA or from lnCeDb_mir_lnc_targets_all (TO DO)
length(unique(lnCeDb_mir_lnc_targets_all$V2)) # 22615

####### lnCeDb_miR_mRNA_lnc_ceRNA ######### (TO DO)
head(lnCeDb_miR_mRNA_lnc_ceRNA)

########  LncACTdb #########################################################################

#### lncACTdb_experimental_validated_miRNA_targets ###########
head(lncACTdb_experimental_validated_miRNA_targets)
# V1: miRNAs, example: hsa-miR-215
# V2: miRNA_target(mRNA or lncRNA?); example: ENSG00000180008
# note: may need to remove some nonsense element

# check intersection miRNA from brca_miRNA_common with miRNA from lncACTdb_experimental_validated_miRNA_targets 
lncACTdb_experimental_validated_miRNA_targets_miRNA = tolower(lncACTdb_experimental_validated_miRNA_targets$V1)
length(unique(lncACTdb_experimental_validated_miRNA_targets_miRNA)) # 697
draw.pairwise.venn(area1 = nrow(brca_miRNA_common), 
                   area2 = length(unique(lncACTdb_experimental_validated_miRNA_targets_miRNA)),
                   length(intersect(rownames(brca_miRNA_common),unique(lncACTdb_experimental_validated_miRNA_targets_miRNA))),
                   category= c("brca_miRNA", "miRNAfrom lncACTdb_experimental_validated_miRNA_targets_miRNA"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0)) # 142

# check intersection lncRNA or mRNA from brca_lncRNA_common with lncRNA or from lncACTdb_experimental_validated_miRNA_targets (TO DO)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl = useMart("ensembl")
ensembl = useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

lncACTdb_experimental_validated_miRNA_targets_encRNAs = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","external_gene_name",
                                                                             "chromosome_name","start_position","end_position","strand"),
                                                              filters = "ensembl_gene_id",
                                                              values = unique(lncACTdb_experimental_validated_miRNA_targets$V2), 
                                                              mart = ensembl)
length(unique(lncACTdb_experimental_validated_miRNA_targets_encRNAs$external_gene_name)) # 12133

draw.pairwise.venn(area1 = nrow(brca_mRNA_common), 
                   area2 = length(unique(lncACTdb_experimental_validated_miRNA_targets_encRNAs$external_gene_name)),
                   length(intersect(rownames(brca_mRNA_common),unique(lncACTdb_experimental_validated_miRNA_targets_encRNAs$external_gene_name))),
                   category= c("brca_mRNA", "mRNA and lncRNA from lncACTdb_experimental_validated_miRNA_targets_encRNAs"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0)) # 11026

# use biomart to annotate lncRNA in brca_lncRNA_common
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl = useMart("ensembl")
ensembl = useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

brca_lncRNA_common_annotated = substr(rownames(brca_lncRNA_common), start = 1, stop = 15)
brca_lncRNA_common_annotated = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","external_gene_name",
                                                    "chromosome_name","start_position","end_position","strand"),
                                     filters = "ensembl_gene_id",
                                     values = brca_lncRNA_common_annotated, 
                                     mart = ensembl)
length(unique(brca_lncRNA_common_annotated$external_gene_name)) # [1] 12148 < 12727
brca_lncRNA_common_annotated_unique = unique(brca_lncRNA_common_annotated$external_gene_name)

draw.pairwise.venn(area1 = length(brca_lncRNA_common_annotated_unique), 
                   area2 = length(unique(lncACTdb_experimental_validated_miRNA_targets_encRNAs$external_gene_name)),
                   length(intersect(brca_lncRNA_common_annotated_unique,unique(lncACTdb_experimental_validated_miRNA_targets_encRNAs$external_gene_name))),
                   category= c("brca_lncRNA_common_annotated_unique", "mRNA and lncRNA from lncACTdb_experimental_validated_miRNA_targets_encRNAs"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0)) # 83? too small?

#### lncACTdb_functionally_activated_lncACTS ########### --> Good! 
head(lncACTdb_functionally_activated_lncACTS)
# lncRNA.Name: lncRNA.Name, example: RP11-467L20.10
# miRNA, example: hsa-miR-145
# Gene.Name: mRNA, example: C11orf9

# miRNA intersection
draw.pairwise.venn(area1 =nrow(brca_miRNA_common), 
                   area2 = length(unique(lncACTdb_functionally_activated_lncACTS$miRNA)),
                   length(intersect(rownames(brca_miRNA_common),unique(tolower(lncACTdb_functionally_activated_lncACTS$miRNA)))),
                   category= c("brca_miRNA_common", "lncACTdb_functionally_activated_lncACTS$miRNA"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0)) # 129

# mRNA intersection
draw.pairwise.venn(area1 =nrow(brca_mRNA_common), 
                   area2 = length(unique(lncACTdb_functionally_activated_lncACTS$Gene.Name)),
                   length(intersect(rownames(brca_mRNA_common),unique(lncACTdb_functionally_activated_lncACTS$Gene.Name))),
                   category= c("brca_mRNA_common", "lncACTdb_functionally_activated_lncACTS$Gene.Name"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0)) # 1280

# lncRNA intersection
draw.pairwise.venn(area1 = length(brca_lncRNA_common_annotated_unique), 
                   area2 = length(unique(lncACTdb_functionally_activated_lncACTS$lncRNA.Name)),
                   length(intersect(brca_lncRNA_common_annotated_unique,unique(lncACTdb_functionally_activated_lncACTS$lncRNA.Name))),
                   category= c("brca_lncRNA_common_annotated_unique", "lncACTdb_functionally_activated_lncACTS$lncRNA.Name"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0)) # 168

####  lncACTdb_cancer_associated ###########
head(lncACTdb_cancer_associated)
# lncRNAEnsgID   lncRNAName       miRNA Genename Disease

# miRNA intersection
draw.pairwise.venn(area1 =nrow(brca_miRNA_common), 
                   area2 = length(unique(lncACTdb_cancer_associated$miRNA)),
                   length(intersect(rownames(brca_miRNA_common),unique(tolower(lncACTdb_cancer_associated$miRNA)))),
                   category= c("brca_miRNA_common", "lncACTdb_cancer_associated$miRNA"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0)) # 83

# mRNA intersection
draw.pairwise.venn(area1 =nrow(brca_mRNA_common), 
                   area2 = length(unique(lncACTdb_cancer_associated$Genename)),
                   length(intersect(rownames(brca_mRNA_common),unique(lncACTdb_cancer_associated$Genename))),
                   category= c("brca_mRNA_common", "lncACTdb_cancer_associated$Genename"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0)) # 268

# lncRNA intersection
draw.pairwise.venn(area1 = length(brca_lncRNA_common_annotated_unique), 
                   area2 = length(unique(lncACTdb_cancer_associated$lncRNAName)),
                   length(intersect(brca_lncRNA_common_annotated_unique,unique(lncACTdb_cancer_associated$lncRNAName))),
                   category= c("brca_lncRNA_common_annotated_unique", "lncACTdb_cancer_associated$lncRNAName"),
                   fill =  c("light blue", "pink"),
                   cat.pos = c(0,0)) # 151





























