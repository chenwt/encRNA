# set directory
setwd("/media/ducdo/UUI1/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
require(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl = useMart("ensembl")
ensembl = useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

load("data_Saved_R_Objects/brca_df.rda")
brca_lncRNA_names = substr(rownames(brca_lncRNA_df), start = 1, stop = 15) # 4828 and all uniques --> good!
brca_mRNA_names = rownames(brca_mRNA_df) # 17,613
brca_miRNA_names = rownames(brca_miRNA_df) # 17,613

##### --------------- lncRNA -----------------------------------------------------------------
rownames(brca_lncRNA_df)[1:3] # "ENSG00000005206.12" "ENSG00000088970.11" "ENSG00000100181.17"

brca_lncRNA_name_annotation = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","external_gene_name","gene_biotype",
                                                 "chromosome_name","start_position","end_position","strand"),
                                  filters = "ensembl_gene_id",
                                  values = brca_lncRNA_names, 
                                  mart = ensembl)
dim(brca_lncRNA_name_annotation) # 4528    7

##### --------------- mRNA -----------------------------------------------------------------
brca_mRNA_name_annotation = getBM(attributes = c("hgnc_symbol","ensembl_gene_id","external_gene_name","gene_biotype",
                                              "chromosome_name","start_position","end_position","strand"),
                               filters = "hgnc_symbol",
                               values = brca_mRNA_names, 
                               mart = ensembl)
dim(brca_mRNA_name_annotation) # 17196     7
length(intersect(brca_mRNA_names, brca_mRNA_name_annotation$hgnc_symbol)) # 15507
setdiff(brca_mRNA_names, brca_mRNA_name_annotation$hgnc_symbol)
setdiff(brca_mRNA_name_annotation$hgnc_symbol, brca_mRNA_names)

##### --------------- miRNA -----------------------------------------------------------------

library(biomaRt)
library(GenomicRanges)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
miRNA = getBM(attributes=c("mirbase_id","ensembl_gene_id","chromosome_name","start_position", "end_position","gene_biotype","external_gene_name"),
              values = brca_miRNA_names,
              mart=ensembl) 
length(miRNA$mirbase_id)
miRNA = miRNA[which(miRNA$mirbase_id != ""),]
length(intersect(miRNA$mirbase_id,brca_miRNA_names)) # 320
setdiff(brca_miRNA_names,miRNA$mirbase_id)
