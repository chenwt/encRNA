## --------------- BASIC SETUP ---------------------------------------

# remove other variables if needed  
rm(list = ls()); gc()

# set directory
setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")

### Required inputs

load("mircode_objects.rda") 
# mircode: 1329071 * 12; but length(unique(mircode$gene_id)) is 50665
# mircode_lncRNA: 138883 * 12 --> only includes overlapping lncRNA and intergenic lncRNA;
#                 but length(unique(mircode_lncRNA$gene_id)) is 10349

load("brca_lncRNA_df.rda")
# [1] brca_lncRNA_df

load("de_brca_lncRNA_objects.rda")
# [1] de_brca_lncRNA_df
# [2] de_brca_lncRNA_mircode_df

load("brca_common_objects.rda")
# [1] "brca_lncRNA_common" 
# [2] "brca_miRNA_common"   
# [3] "brca_mRNA_common"

### Output:
load("de_brca_lncRNA_annotation.rda")

## --------------- LOAD LIBRARIES ----------------------------------------------
library(limma)
library(biomaRt)

## --------------- BASIC CHECK -------------------------------------------------

head(mircode)
dim(mircode) # [1] 1329071      12
dim(mircode_lncRNA) # [1] 138883     12
length(unique(mircode_lncRNA)) # 10349
length(rownames(de_brca_lncRNA_mircode_df) %in% unique(mircode_lncRNA$gene_id)) # 1093
# save(mircode, mircode_lncRNA, file = "mircode_objects.rda")

## --------------- SUBSETTING miRcode ------------------------------------------

mircode_de_lncRNA = mircode_lncRNA[which(mircode_lncRNA$gene_id %in% 
                                           rownames(de_brca_lncRNA_mircode_df)),]

dim(mircode_de_lncRNA) #[1] 31952    12
head(mircode_de_lncRNA)

# --------------- BIOMART ANNOTATION   -----------------------------------------

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl = useMart("ensembl")
ensembl = useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

# filters = listFilters(ensembl)
# attributes = listAttributes(ensembl)

# remove decimal part of the ensembl Id
de_brca_lncRNA_mircode_samples = substr(rownames(de_brca_lncRNA_mircode_df), start = 1, stop = 15)
length(de_brca_lncRNA_mircode_samples) # [1] 1903
length(unique(de_brca_lncRNA_mircode_samples)) # [1] 1903

de_brca_lncRNA_annotation = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","external_gene_name",
                                       "chromosome_name","start_position","end_position","strand"),
                        filters = "ensembl_gene_id",
                        values = de_brca_lncRNA_mircode_samples, 
                        mart = ensembl)
dim(de_brca_lncRNA_annotation) # [1] 1837    7 --> so, (1903 - 1837) = 66 ensemblID is not annotated
save(de_brca_lncRNA_annotation, file = "de_brca_lncRNA_annotation.rda")



















