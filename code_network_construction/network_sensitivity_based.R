
setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")

load("data_Saved_R_Objects/brca_df.rda")
load("data_Saved_R_Objects/corr_matrices/normal_encRNA_850356.rda")

ptm = proc.time()
test =  sapply(1:10000, function(index){
  normal_lncRNA_mRNA_corr_matrix[as.character(normal_encRNA$lncRNA[index]),as.character(normal_encRNA$mRNA[index])]
})
ptm = proc.time() - ptm

########################################################################################3
# test code: construct graph of mRNA-miRNA from normal_encRNA_df
# https://www.r-bloggers.com/network-visualization-part-5-cytoscape-an-update-rcy3/

require(igraph); require(RCy3); require(plyr); require(dplyr); require(data.table)

# for testing purpose, subsetting it to only first 10 rows
#normal_encRNA_sensitivity_bound = get_matched_enRNA_sensitivity_with_putative_binding(normal_encRNA)
dim(normal_encRNA_sensitivity_bound) # 21553    12
normal_encRNA_df = normal_encRNA_sensitivity_bound[,c("mRNA","lncRNA","sensitivity")]

normal_encRNA_df$pair = paste(normal_encRNA_df$mRNA, normal_encRNA_df$lncRNA, sep = "~")

# there are duplicated pairs in normal_encRNA_df. In such a case, compute the average 
# sensitivity of rows from a similar pair
keys = colnames(normal_encRNA_df)[!grepl('sensitivity', colnames(normal_encRNA_df))]
X = as.data.table(normal_encRNA_df)
normal_encRNA_df =X[,list(sensitivity = mean(sensitivity)),keys] # this code is ought be changed based on how sensitivity will be computed 
normal_encRNA_df = as.data.frame(normal_encRNA_df)

# start constructing the igraph object
node_type_df = data.frame(node = c(unique(normal_encRNA_df$mRNA), unique(normal_encRNA_df$lncRNA)))
node_type_df$type = c(rep("mRNA", length(unique(normal_encRNA_df$mRNA))),
                      rep("lncRNA", length(unique(normal_encRNA_df$lncRNA)))) 
graph_normal_encRNA = graph.data.frame(d = normal_encRNA_df[,c("mRNA","lncRNA")], 
                                       vertices = node_type_df,
                                       directed = F)
V(graph_normal_encRNA)$color <- ifelse(V(graph_normal_encRNA)$type == "mRNA", "red", "blue")
E(graph_normal_encRNA)$weight = normal_encRNA_df[,"sensitivity"]

plot.igraph(graph_normal_encRNA)
ecount(graph_normal_encRNA) # 19402
vcount(graph_normal_encRNA) # 3705

# list multiple edges
E(graph_normal_encRNA)[which_multiple(graph_normal_encRNA, 
                                                        eids = E(graph_normal_encRNA))]


graphNel_obj <- igraph::as_graphnel(graph_normal_encRNA); gc()
#save(graphNel_obj, file = "graphNel_obj.rda")
# check to see 
graph::nodeData(graphNel_obj, igraph::V(graph_normal_encRNA)$name, 'color')
graph::nodeData(graphNel_obj, igraph::V(graph_normal_encRNA)$name, 'type')
graph::edgeData(graphNel_obj, 
                as.character(normal_encRNA_df$lncRNA), 
                as.character(normal_encRNA_df$mRNA),
                'weight')

# init attributes
graphNel_obj <- RCy3::initNodeAttribute(graphNel_obj, 'color', 'char', 'a') 
graphNel_obj <- RCy3::initNodeAttribute(graphNel_obj, 'type', 'char', 'a') 
graphNel_obj <- RCy3::initEdgeAttribute (graphNel_obj, "weight", 'numeric', 0)

# Next, we will create a new graph window in cytoscape
cytoscape_window <- RCy3::CytoscapeWindow("encRNA_normal", graph = graphNel_obj, overwriteWindow = TRUE)
# We can display graph, with defaults color/size scheme
RCy3::displayGraph(cytoscape_window)

# set node and edge attributes
RCy3::setNodeAttributesDirect(cytoscape_window, 'type', 'char', 
                              igraph::V(graph_normal_encRNA)$name, 
                              igraph::V(graph_normal_encRNA)$type)
RCy3::setNodeAttributesDirect(cytoscape_window, 'color', 'char', 
                              igraph::V(graph_normal_encRNA)$name, 
                              igraph::V(graph_normal_encRNA)$color)


edge_names = names(RCy3::cy2.edge.names (cytoscape_window@graph))

# test$pair = as.factor(test$pair)
# test$pair = gdata::reorder.factor(test$pair, new.order=edge_names)
# test %>% dplyr::arrange(pair)
normal_encRNA_df = normal_encRNA_df[match(edge_names, normal_encRNA_df$pair),]

RCy3::setEdgeAttributesDirect(obj = cytoscape_window, 
                              attribute.name = 'weight', 
                              attribute.type = 'numeric', 
                              edge.names = as.character (RCy3::cy2.edge.names (cytoscape_window@graph)), 
                              values = normal_encRNA_df$sensitivity)


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

RCy3::setDefaultNodeSize(cytoscape_window, 100)
RCy3::setDefaultNodeFontSize(cytoscape_window, 10)
RCy3::s(cytoscape_window, 10)


######################################################################################
