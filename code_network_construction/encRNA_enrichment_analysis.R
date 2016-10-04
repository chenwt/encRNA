setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_network_construction/network_helper_functions.R")

require(igraph); require(RCy3); require(plyr); require(dplyr); require(data.table); require(rlist)

load("brca2_TCGA2STAT.rda")
load("data_Saved_R_Objects/miRNA_target/predicted_normal_tumor_goodCorr.rda")
load("data_Saved_R_Objects/brca_df.rda")

load("data_Saved_R_Objects/corr_matrices/normal_encRNA_850356.rda")
load("data_Saved_R_Objects/corr_matrices/tumor_encRNA_850356.rda")

##########################################################################################
########            GENES OR LNCRNAS OF SIMILAR MIRNAS                            ########
##########################################################################################

normal_unique_miRNAs = as.character(unique(normal_encRNA_sensitivity_bound_goodCoor$miRNA))
tumor_unique_miRNAs = as.character(unique(tumor_encRNA_sensitivity_bound_goodCoor$miRNA))

length(normal_unique_miRNAs) # 7
length(tumor_unique_miRNAs) # 21

#### obtain list of mRNA and lncRNA sharing similar miRNAs
normal_RNA_list = get_RNAs_sharing_common_miRNAs(normal_encRNA_sensitivity_bound_goodCoor,normal_unique_miRNAs)
normal_mRNA_list = normal_RNA_list[['mRNA']]; normal_lncRNA_list = normal_RNA_list[['lncRNA']]
length(normal_mRNA_list) #7
tumor_RNA_list = get_RNAs_sharing_common_miRNAs(tumor_encRNA_sensitivity_bound_goodCoor,tumor_unique_miRNAs)
tumor_mRNA_list = tumor_RNA_list[['mRNA']]; tumor_lncRNA_list = tumor_RNA_list[['lncRNA']]
length(tumor_mRNA_list) # 21

#### GOTerm enrichment analysis on those normal_mRNA_list
# obtain all ensemble gene symbols
all_ensembl_gene_symbols = get_all_human_ensemble()

# filter the RNA list to include only those having annotated in ensemble
normal_mRNA_list = get_mRNAs_in_Ensembl(geneList = normal_mRNA_list,
                                        ensemble = all_ensembl_gene_symbols,
                                        minimumSize = 4)
length(normal_mRNA_list) # 7
tumor_mRNA_list = get_mRNAs_in_Ensembl(geneList = tumor_mRNA_list,
                                       ensemble = all_ensembl_gene_symbols,
                                       minimumSize = 4)
length(tumor_mRNA_list) # 8 (decreased from 21)

#### obtain result from topGo
normal_mRNA_topGo_common_miRNAs = list()
for (i in 1:length(normal_mRNA_list)){
  result = get_TopGo_result(all_ensembl_gene_symbols = all_ensembl_gene_symbols,
                            gene_list = normal_mRNA_list[[i]],
                            topNodes = 20)
  normal_mRNA_topGo_common_miRNAs = list.append(normal_mRNA_topGo_common_miRNAs, result)
}

# plot Go Term
# http://stackoverflow.com/questions/6901405/r-ggplot-ordering-bars-in-barplot-like-plot
require(ggplot2)
result$Term = reorder(result$Term, -log10(as.numeric(result$classicFisher)))
ggplot2::ggplot(data = result,
       mapping = aes(x = result$Term, y = -log10(as.numeric(result$classicFisher)), fill = result$Group)) + geom_bar(stat = "identity") + coord_flip()
  

##########################################################################################
########            Module Detection and Analysis                                 ########
##########################################################################################

normal_encRNA_graph = build_iGraph(df = normal_encRNA_sensitivity_bound_goodCoor)
tumor_encRNA_graph = build_iGraph(df = tumor_encRNA_sensitivity_bound_goodCoor)

## ----- basic network statistic ------------------------- 

# vertice and edge count
vcount(normal_encRNA_graph); ecount(normal_encRNA_graph) # 359, 1907
V(normal_encRNA_graph)[which(V(normal_encRNA_graph)$type == "mRNA")] # 301
V(normal_encRNA_graph)[which(V(normal_encRNA_graph)$type == "lncRNA")] # 58

vcount(tumor_encRNA_graph); ecount(tumor_encRNA_graph) # [1] 1302, 1218
V(tumor_encRNA_graph)[which(V(tumor_encRNA_graph)$type == "mRNA")] # 956
V(tumor_encRNA_graph)[which(V(tumor_encRNA_graph)$type == "lncRNA")] # 346

# select hub genes
normal_mRNA_degree = igraph::degree(normal_encRNA_graph)[which(V(normal_encRNA_graph)$type == "mRNA")]
normal_lncRNA_degree = igraph::degree(normal_encRNA_graph)[which(V(normal_encRNA_graph)$type == "lncRNA")]
tumor_mRNA_degree = igraph::degree(tumor_encRNA_graph)[which(V(tumor_encRNA_graph)$type == "mRNA")]
tumor_lncRNA_degree = igraph::degree(tumor_encRNA_graph)[which(V(tumor_encRNA_graph)$type == "lncRNA")]

par(mfrow = c(2,2))
hist(normal_mRNA_degree, breaks = length(unique(normal_mRNA_degree)))
hist(normal_lncRNA_degree, breaks = length(unique(normal_lncRNA_degree)))
hist(tumor_mRNA_degree, breaks = length(unique(tumor_mRNA_degree)))
hist(tumor_lncRNA_degree, breaks = length(unique(tumor_lncRNA_degree)))
par(mfrow = c(1,1))

# graph density
graph.density(normal_encRNA_graph) # [1] 0.02967585
graph.density(tumor_encRNA_graph) # [1] 0.0014381

# clustering coefficent

## ----- obtain list of mRNAs ------------------------- 
module_normal_encRNA_graph = cluster_louvain(normal_encRNA_graph,weights = E(normal_encRNA_graph)$weight)
module_tumor_encRNA_graph = cluster_louvain(tumor_encRNA_graph,weights = E(tumor_encRNA_graph)$weight)

normal_encRNA_graph_membership = getAllSubgraphsFromMembership(igraphObject = normal_encRNA_graph,
                                                               membership = module_normal_encRNA_graph$membership)
tumor_encRNA_graph_membership = getAllSubgraphsFromMembership(igraphObject = tumor_encRNA_graph,
                                                               membership = module_tumor_encRNA_graph$membership)

#### get the list of mRNAs
normal_encRNA_mRNA_list = lapply(normal_encRNA_graph_membership, function(graph){
  names(V(graph)[which(V(graph)$type == "mRNA")])
})

tumor_encRNA_mRNA_list = lapply(tumor_encRNA_graph_membership, function(graph){
  names(V(graph)[which(V(graph)$type == "mRNA")])
})

## ------ perform enrichment analysis on modules-------------------
## (must perform above steps)
setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
load("brca2_TCGA2STAT.rda")
load("data_Saved_R_Objects/miRNA_target/predicted_normal_tumor_goodCorr.rda")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_network_construction/network_helper_functions.R")
require(topGO); require(biomaRt)
# general idea: for each mRNA list, perform topGo enrichment

# get all possible ids
all_ensembl_gene_symbols = get_all_human_ensemble()
# get the genes from the module which are also includedi the ensemble gene symbols
# this step would need some refinement
normal_encRNA_mRNA_list = get_mRNAs_in_Ensembl(geneList = normal_encRNA_mRNA_list,
                                                ensemble = all_ensembl_gene_symbols, 
                                               minimumSize = 4)
tumor_encRNA_mRNA_list = get_mRNAs_in_Ensembl(geneList = tumor_encRNA_mRNA_list,
                                                ensemble = all_ensembl_gene_symbols, 
                                              minimumSize = 4)

tumor_encRNA_mRNA_list = tumor_encRNA_mRNA_list[-which(length_gene_list <= 1)]

# perform GO enrichment analysis for each element in the list
normal_encRNA_module_enrichment = list()
tumor_encRNA_module_enrichment = list()

require(rlist); require(org.Hs.eg.db)

lapply(normal_encRNA_mRNA_list, function(mRNA_group){
  allRes = get_TopGo_result(all_ensembl_gene_symbols = all_ensembl_gene_symbols, 
                            gene_list = mRNA_group, 
                            topNodes = 20)
  normal_encRNA_module_enrichment <<- list.append(normal_encRNA_module_enrichment,allRes)
})

lapply(tumor_encRNA_mRNA_list, function(mRNA_group){
  allRes = get_TopGo_result(all_ensembl_gene_symbols = all_ensembl_gene_symbols, 
                            gene_list = mRNA_group, 
                            topNodes = 20)
  tumor_encRNA_module_enrichment <<- list.append(tumor_encRNA_module_enrichment,allRes)
})

save(normal_encRNA_module_enrichment, tumor_encRNA_module_enrichment, file = "data_Saved_R_Objects/louvain_enrichment.rda")

