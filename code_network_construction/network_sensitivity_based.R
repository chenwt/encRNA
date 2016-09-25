
setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")

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

# encRNA_df = normal_encRNA_sensitivity_bound[,c("mRNA","lncRNA","sensitivity")]
encRNA_df = tumor_encRNA_sensitivity_bound[,c("mRNA","lncRNA","sensitivity")]

encRNA_df$pair = paste(encRNA_df$mRNA, encRNA_df$lncRNA, sep = "~")

# there are duplicated pairs in encRNA_df. In such a case, compute the average 
# sensitivity of rows from a similar pair
keys = colnames(encRNA_df)[!grepl('sensitivity', colnames(encRNA_df))]
X = as.data.table(encRNA_df)
encRNA_df =X[,list(sensitivity = mean(sensitivity)),keys] # this code is ought be changed based on how sensitivity will be computed 
encRNA_df = as.data.frame(encRNA_df)

# start constructing the igraph object
node_type_df = data.frame(node = c(unique(encRNA_df$mRNA), unique(encRNA_df$lncRNA)))
node_type_df$type = c(rep("mRNA", length(unique(encRNA_df$mRNA))),
                      rep("lncRNA", length(unique(encRNA_df$lncRNA)))) 
graph_encRNA = graph.data.frame(d = encRNA_df[,c("mRNA","lncRNA")], 
                                       vertices = node_type_df,
                                       directed = F)
V(graph_encRNA)$color <- ifelse(V(graph_encRNA)$type == "mRNA", "red", "blue")
E(graph_encRNA)$weight = encRNA_df[,"sensitivity"]

plot.igraph(graph_encRNA)
ecount(graph_encRNA) # 19402
vcount(graph_encRNA) # 3705

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
