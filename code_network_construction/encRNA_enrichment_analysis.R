setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_network_construction/network_helper_functions.R")

require(igraph); require(RCy3); require(plyr); require(dplyr); require(data.table); require(rlist)
require(ggplot2); require(gridExtra); require(GO.db)

load("brca2_TCGA2STAT.rda")
load("data_Saved_R_Objects/miRNA_target/predicted_normal_tumor_goodCorr.rda")
load("data_Saved_R_Objects/brca_df.rda")

load("data_Saved_R_Objects/corr_matrices/normal_encRNA_850356.rda")
load("data_Saved_R_Objects/corr_matrices/tumor_encRNA_850356.rda")

# ------------ helper function ----------------------------------------------------------

plot_ea = function(data, title){
  require(ggplot2)
  data$Term = reorder(data$Term, -log10(as.numeric(data$classicFisher)))
  ggplot2::ggplot(data = data,
                  mapping = aes(x = Term, y = -log10(as.numeric(classicFisher)), fill = Group)) + geom_bar(stat = "identity") + xlab("GO Term")+  ylab("-log10(p-value)") + coord_flip()+ ggtitle(title)
}

#### -------------- obtain all ensemble gene symbols --------------------------------
all_ensembl_gene_symbols = get_all_human_ensemble()


##########################################################################################
########            GENES OR LNCRNAS OF SIMILAR MIRNAS                            ########
##########################################################################################

normal_miRNAs = as.character(unique(normal_encRNA_sensitivity_bound_goodCoor$miRNA))
tumor_miRNAs = as.character(unique(tumor_encRNA_sensitivity_bound_goodCoor$miRNA))
normal_mRNAs = unique(normal_encRNA_sensitivity_bound_goodCoor$mRNA)
normal_lncRNAs = unique(normal_encRNA_sensitivity_bound_goodCoor$lncRNA)
tumor_mRNAs = unique(tumor_encRNA_sensitivity_bound_goodCoor$mRNA); # tumor_mRNAs = tumor_mRNAs[-241]
tumor_lncRNAs = unique(tumor_encRNA_sensitivity_bound_goodCoor$lncRNA)

length(normal_miRNAs) # 7
length(tumor_miRNAs) # 21

write(tumor_miRNAs, file = "data_Saved_R_Objects/encRNA_network/tumor_miRNA.txt")

##########################################################################################
#### Enrichment analysis for all mRNAs
##########################################################################################

par(mfrow = c(1,2))
normal_mRNA_ea = get_TopGo_result(all_ensembl_gene_symbols = all_ensembl_gene_symbols,
                          gene_list = normal_mRNAs,
                          topNodes = 100)
tumor_mRNA_ea = get_TopGo_result(all_ensembl_gene_symbols = all_ensembl_gene_symbols,
                                  gene_list = tumor_mRNAs,
                                  topNodes = 100)
p1 = plot_ea(normal_mRNA_ea, "mRNA normal")
p2 = plot_ea(tumor_mRNA_ea, "mRNA tumor")
grid.arrange(p1, p2, ncol=2)

##########################################################################################
#### GOTerm enrichment analysis on all mRNA sharing similar miRNA
##########################################################################################

#### obtain list of mRNA and lncRNA sharing similar miRNAs
normal_RNA_list = get_RNAs_sharing_common_miRNAs(normal_encRNA_sensitivity_bound_goodCoor,normal_miRNAs)
normal_mRNA_list = normal_RNA_list[['mRNA']]; normal_lncRNA_list = normal_RNA_list[['lncRNA']]
length(normal_mRNA_list) #7
tumor_RNA_list = get_RNAs_sharing_common_miRNAs(tumor_encRNA_sensitivity_bound_goodCoor,tumor_miRNAs)
tumor_mRNA_list = tumor_RNA_list[['mRNA']]; tumor_lncRNA_list = tumor_RNA_list[['lncRNA']]
length(tumor_mRNA_list) # 21

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
names(normal_mRNA_topGo_common_miRNAs) = names(normal_mRNA_list)

go_term = c()
sapply(normal_mRNA_topGo_common_miRNAs, function(df){
  go_term <<- append(go_term, df$GO.ID)
})

Term(go_term[duplicated(go_term)])

tumor_mRNA_topGo_common_miRNAs = list()
for (i in 1:length(normal_mRNA_list)){
  result = get_TopGo_result(all_ensembl_gene_symbols = all_ensembl_gene_symbols,
                            gene_list = tumor_mRNA_list[[i]],
                            topNodes = 20)
  tumor_mRNA_topGo_common_miRNAs = list.append(tumor_mRNA_topGo_common_miRNAs, result)
}
names(tumor_mRNA_topGo_common_miRNAs) = names(tumor_mRNA_list)

##########################################################################################
#### GOTerm enrichment analysis on all on mRNA group based on connection with lncRNA
##########################################################################################

normal_mRNA_connected_lncRNA_list = get_mRNAs_connected_to_lncRNA(normal_encRNA_graph,names(normal_lncRNA_degree))
tumor_mRNA_connected_lncRNA_list = get_mRNAs_connected_to_lncRNA(tumor_encRNA_graph,names(tumor_lncRNA_degree))

normal_encRNA_mRNA_list = get_mRNAs_in_Ensembl(geneList = normal_mRNA_connected_lncRNA_list,
                                               ensemble = all_ensembl_gene_symbols, 
                                               minimumSize = 10)
tumor_encRNA_mRNA_list = get_mRNAs_in_Ensembl(geneList = normal_mRNA_connected_lncRNA_list,
                                               ensemble = all_ensembl_gene_symbols, 
                                               minimumSize = 10)


##########################################################################################
#### GOTerm enrichment analysis on all on miRNAs
##########################################################################################
normal_miRNAs = gsub(x = normal_miRNAs, pattern = "mir", replacement = "miR")
tumor_miRNAs = gsub(x = tumor_miRNAs, pattern = "mir", replacement = "miR")
write(normal_miRNAs, file = "data_Saved_R_Objects/encRNA_network/normal_miRNAs.txt")
write(tumor_miRNAs, file = "data_Saved_R_Objects/encRNA_network/tumor_miRNAs.txt")

##########################################################################################
########            Module Detection and Analysis                                 ########
##########################################################################################

normal_encRNA_graph = build_iGraph(df = normal_encRNA_sensitivity_bound_goodCoor)
tumor_encRNA_graph = build_iGraph(df = tumor_encRNA_sensitivity_bound_goodCoor)

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

