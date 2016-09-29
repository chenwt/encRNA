
setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_network_construction/network_helper_functions.R")

require(data.table); require(igraph)

load("brca2_TCGA2STAT.rda")
load("data_Saved_R_Objects/miRNA_target/predicted_normal_tumor_goodCorr.rda")
load("data_Saved_R_Objects/brca_df.rda")

load("data_Saved_R_Objects/corr_matrices/normal_encRNA_850356.rda")
load("data_Saved_R_Objects/corr_matrices/tumor_encRNA_850356.rda")

########################################################################################3
# test code: construct graph of mRNA-miRNA from encRNA_df
# https://www.r-bloggers.com/network-visualization-part-5-cytoscape-an-update-rcy3/

require(igraph); require(RCy3); require(plyr); require(dplyr); require(data.table)

# for testing purpose, subsetting it to only first 10 rows
normal_encRNA_sensitivity_bound = get_matched_enRNA_sensitivity_with_putative_binding(normal_encRNA)
tumor_encRNA_sensitivity_bound = get_matched_enRNA_sensitivity_with_putative_binding(tumor_encRNA)
dim(normal_encRNA_sensitivity_bound) # 21553    12


##########################################################################################
########            Module Detection and Analysis                                 ########
##########################################################################################

normal_encRNA_graph = build_iGraph(df = normal_encRNA_sensitivity_bound_goodCoor)
tumor_encRNA_graph = build_iGraph(df = tumor_encRNA_sensitivity_bound_goodCoor)

graph.density(normal_encRNA_graph)
graph.density(tumor_encRNA_graph)

## ----- obtain list of mRNAs ------------------------- 
module_normal_encRNA_graph = cluster_louvain(normal_encRNA_graph,weights = E(normal_encRNA_graph)$weight)
module_tumor_encRNA_graph = cluster_louvain(tumor_encRNA_graph,weights = E(tumor_encRNA_graph)$weight)

normal_encRNA_graph_membership = getAllSubgraphsFromMembership(igraphObject = normal_encRNA_graph,
                                                               membership = module_normal_encRNA_graph$membership)
test = normal_encRNA_graph_membership[[1]]
vcount(test)
V(test)$type

#### get the list of mRNAs
normal_encRNA_mRNA_list = lapply(normal_encRNA_graph_membership, function(graph){
  names(V(graph)[which(V(graph)$type == "mRNA")])
})

#########################################################
## ------ perform enrichment analysis -------------------
load("brca2_TCGA2STAT.rda")
load("data_Saved_R_Objects/miRNA_target/predicted_normal_tumor_goodCorr.rda")

require(topGO); require(biomaRt)
# general idea: for each mRNA list, perform topGo enrichment

# get all possible ids
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 85)
listDatasets(ensembl)
mart = useDataset("hsapiens_gene_ensembl", mart= ensembl)
all_ensembl_gene_symbols = getBM(attributes = c("hgnc_symbol","external_gene_name"),
                                 values = "*", mart = mart)
length(all_ensembl_gene_symbols$hgnc_symbol) # 35,524

# get the genes from the module which are also includedi the ensemble gene symbols
# this step would need some refinement

normal_encRNA_mRNA_list = lapply(normal_encRNA_mRNA_list, function(gene_list){
  gene_list = gene_list[which(gene_list %in% all_ensembl_gene_symbols$hgnc_symbol)]; 
})

# perform GO enrichment analysis for each element in the list
normal_encRNA_module_enrichment = list()
require(rlist); require(org.Hs.eg.db)

lapply(normal_encRNA_mRNA_list[1:2], function(mRNA_group){
  universe = factor(as.integer(all_ensembl_gene_symbols$hgnc_symbol %in% mRNA_group))
  names(universe) = all_ensembl_gene_symbols$hgnc_symbol
  
  GOdata = new('topGOdata',
               ontology="MF",
               allGenes=universe,
               annot=annFUN.org,
               mapping="org.Hs.eg.db", 
               ID="symbol")
  
  #fischer
  resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  #top 20 significant nodes determined by KS method, can compare to classic 
  #fischer and weight
  allRes <- GenTable(GOdata, classicFisher = resultFis, 
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 20)
  normal_encRNA_module_enrichment <<- list.append(normal_encRNA_module_enrichment,allRes)
})

#########################################################
## ------ perform survival analysis----------------------
setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_network_construction/network_helper_functions.R")
load("brca2_TCGA2STAT.rda")
load("data_Saved_R_Objects/brca_df.rda")
load("data_Saved_R_Objects/miRNA_target/predicted_normal_tumor_goodCorr.rda")
require(survival)

clinical = brca_miRNA$clinical

## select sample names also showed up in filtered normal and tumor sample 
colnames(brca_miRNA_df)

normal_samples = colnames(brca_miRNA_df)[1:79]; length(normal_samples) # 79
tumor_samples = colnames(brca_miRNA_df)[80:536]; length(tumor_samples) # 457

normal_samples =  strsplit(normal_samples, split = "[.]")
normal_samples = unlist(lapply(normal_samples, function(name){
  paste(name[3],name[4],name[5],sep="-")
}))

tumor_samples =  strsplit(tumor_samples, split = "[.]")
tumor_samples = unlist(lapply(tumor_samples, function(name){
  paste(name[3],name[4],name[5],sep="-")
}))

length(intersect(normal_samples, tumor_samples)) # 76 out of 79 normal samples having matched tumor samples
length(intersect(rownames(clinical), tumor_samples)) # 45 -->
# thus all tumor samples are included in the clinical dataset

# subset clinical data
clinical = clinical[tumor_samples,c("vitalstatus","daystodeath","daystolastfollowup")]
clinical = as.data.frame(clinical)
nrow(clinical) # 457
length(which(is.na(clinical$daystodeath)))
length(which(is.na(clinical$daystolastfollowup)))
# there is one sample does not have neither daystodeath or daystolastfollowup
which(is.na(clinical$daystodeath) & is.na(clinical$daystolastfollowup))
# TCGA-E9-A245 
# 413 
clinical = clinical[-which(is.na(clinical$daystodeath) & is.na(clinical$daystolastfollowup)),]
nrow(clinical) # 456

# create a new column which aggregate daystodeath and daystolastfollowup
clinical$time = ifelse(is.na(clinical$daystodeath), clinical$daystolastfollowup, clinical$daystodeath)
sum(is.na(clinical$time)) # 0 --> good!
clinical$SurvObj <- with(clinical, Surv(time, vitalstatus == 1))

# TO DO: adding expression data as column into the clinical data_frame. 
# epxression data are selected from the modules 

# EXTRA: use Cox-Lasso as an idenpedent steps to select distinct genes in term of survival


########################################################################################
###               Export to Cytoscape                                         ##########
########################################################################################
# list multiple edges
E(graph_encRNA)[which_multiple(graph_encRNA,  eids = E(graph_encRNA))]

graphNel_obj <- igraph::as_graphnel(graph_encRNA); gc()
#save(graphNel_obj, file = "graphNel_obj.rda")
# check to see 
graph::nodeData(graphNel_obj, igraph::V(graph_encRNA)$name, 'color')
graph::nodeData(graphNel_obj, igraph::V(graph_encRNA)$name, 'type')
graph::edgeData(graphNel_obj, 
                as.character(encRNA_df$lncRNA), 
                as.character(encRNA_df$mRNA),
                'weight')

# init attributes
graphNel_obj <- RCy3::initNodeAttribute(graphNel_obj, 'color', 'char', 'a') 
graphNel_obj <- RCy3::initNodeAttribute(graphNel_obj, 'type', 'char', 'a') 
graphNel_obj <- RCy3::initEdgeAttribute (graphNel_obj, "weight", 'numeric', 0)

# Next, we will create a new graph window in cytoscape
# cytoscape_window <- RCy3::CytoscapeWindow("encRNA_normal", graph = graphNel_obj, overwriteWindow = TRUE)
cytoscape_window <- RCy3::CytoscapeWindow("encRNA_tumor", graph = graphNel_obj, overwriteWindow = TRUE)

# We can display graph, with defaults color/size scheme
RCy3::displayGraph(cytoscape_window)

# set node and edge attributes
RCy3::setNodeAttributesDirect(cytoscape_window, 'type', 'char', 
                              igraph::V(graph_encRNA)$name, 
                              igraph::V(graph_encRNA)$type)
RCy3::setNodeAttributesDirect(cytoscape_window, 'color', 'char', 
                              igraph::V(graph_encRNA)$name, 
                              igraph::V(graph_encRNA)$color)


edge_names = names(RCy3::cy2.edge.names (cytoscape_window@graph))

# test$pair = as.factor(test$pair)
# test$pair = gdata::reorder.factor(test$pair, new.order=edge_names)
# test %>% dplyr::arrange(pair)
encRNA_df = encRNA_df[match(edge_names, encRNA_df$pair),]

RCy3::setEdgeAttributesDirect(obj = cytoscape_window, 
                              attribute.name = 'weight', 
                              attribute.type = 'numeric', 
                              edge.names = as.character (RCy3::cy2.edge.names (cytoscape_window@graph)), 
                              values = encRNA_df$sensitivity)


# If you also want to choose a layout from R, a list  of available layouts can be accessed as follow:
cy <- RCy3::CytoscapeConnection()
hlp <-RCy3::getLayoutNames(cy)
# hlp
# [1] "attribute-circle"      "stacked-node-layout"   "degree-circle"         "circular"              "attributes-layout"     "kamada-kawai"         
# [7] "force-directed"        "grid"                  "hierarchical"          "fruchterman-rheingold" "isom"      
# We'll select the "fruchterman-rheingold" layout. This layout is the layout number 10 
# To see properties for the given layout, use:
# RCy3::getLayoutPropertyNames(cy, hlp[10])
# We can choose any property we want and provide them as a list
#RCy3::setLayoutProperties (cytoscape_window, hlp[10], list (gravity_multiplier = 'similarity', nIterations = 1))
RCy3::setLayoutProperties (cytoscape_window, hlp[10], list (gravity_multiplier = 'similarity', nIterations = 1))
RCy3::layoutNetwork(cytoscape_window, hlp[10])

RCy3::setDefaultNodeSize(cytoscape_window, 50)
RCy3::setDefaultNodeFontSize(cytoscape_window, 10)


######################################################################################
