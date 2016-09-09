## --------------- LOAD LIBRARIES -------------------------------------------------
library('igraph')
library('gplots')
library('ggplot2')
library('knitr')
library('limma')
library('reshape2')
library('RColorBrewer')
library('WGCNA')
library('gdata')

## --------------- LOAD OBJECTS ----------------------------------------------------------

# remove other variables if needed  
rm(list = ls()); gc()

# set directory
setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")

# load("de_brca_lncRNA_objects.rda")
load("de_brca_miRNA_to_encRNAs.rda") # de_brca_miRNA_to_mRNA_list, de_brca_miRNA_to_lncRNA_list
load("de_brca_enRNAs_candidate_list.rda") # de_brca_enRNAs_candidate_mRNAs_list, de_brca_enRNAs_candidate_lncRNAs_list -- common miRNAs

## --------------- IGHRAPH FROM DATA FRAME------------------------------------------------ 
load("de_brca_enRNAs_candidate_list.rda")

length(de_brca_enRNAs_candidate_mRNAs_list) # 27
total_lncRNAs = sapply(de_brca_enRNAs_candidate_mRNAs_list, function(miRNA_target) length(miRNA_target))
sum(total_lncRNAs) # 23553

de_brca_enRNAs_candidate_mRNAs_df = reshape2::melt(de_brca_enRNAs_candidate_mRNAs_list)
dim(de_brca_enRNAs_candidate_mRNAs_df) # [1] 23553     2 -- correct!
head(de_brca_enRNAs_candidate_mRNAs_df,3)
# value          L1
# 1   FBN2 hsa-mir-429
# 2 RNF160 hsa-mir-429
# 3 SEC24A hsa-mir-429
length(unique(de_brca_enRNAs_candidate_mRNAs_df$value)) # 6451
length(unique(de_brca_enRNAs_candidate_mRNAs_df$L1)) # 27
# thus, we expect the derived graph would have 6451 + 27 = 6478 nodes and 23553 edges 

graph_de_brca_miRNA_to_mRNA = graph.data.frame(d = de_brca_enRNAs_candidate_mRNAs_df, directed = F)
vcount(graph_de_brca_miRNA_to_mRNA) # 6478 good!
ecount(graph_de_brca_miRNA_to_mRNA) # 23553 good!

de_brca_enRNAs_candidate_lncRNAs_df = reshape2::melt(de_brca_enRNAs_candidate_lncRNAs_list)
dim(de_brca_enRNAs_candidate_lncRNAs_df) # [1] 9215      2
length(unique(de_brca_enRNAs_candidate_lncRNAs_df$value)) # 1838
length(unique(de_brca_enRNAs_candidate_lncRNAs_df$L1)) # 27
# thus, we expect the derived graph would have 1838 + 27 = 1865 nodes and 9215 edges 

graph_de_brca_miRNA_to_lncRNA = graph.data.frame(d = de_brca_enRNAs_candidate_lncRNAs_df, directed = F)
vcount(graph_de_brca_miRNA_to_lncRNA) # 1865 good!
ecount(graph_de_brca_miRNA_to_lncRNA) # 9215 good!

graph_de_brca_encRNA = graph_de_brca_miRNA_to_mRNA + graph_de_brca_miRNA_to_lncRNA
vcount(graph_de_brca_encRNA) # 8316 == 6478 + 1838 --> good
ecount(graph_de_brca_encRNA) # 32768 = 23553 + 9215 --> good!








## --------------- SIMMILARITY MEASURES ----------------------------------------

# Similarity measure which combines elements from Pearson correlation and
# Euclidean distance.
cordist <- function(dat) {
  cor_matrix  <- cor(t(dat))
  dist_matrix <- as.matrix(dist(dat, diag=TRUE, upper=TRUE))
  dist_matrix <- log1p(dist_matrix)
  dist_matrix <- 1 - (dist_matrix / max(dist_matrix))
  sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix)/ 2)
}

# extract simmilarity matrix
de_brca_lncRNA_sim_matrix = cordist(de_brca_lncRNA_mircode_df)

# Let's see what our similarity matrix looks like at this point. 
# Because the heatmap.2 function (which includes a biclustering step) can be
# pretty slow, we will use a sub-sample of our data --
# for visualization purposes this is fine.

heatmap_indices = sample(nrow(de_brca_lncRNA_sim_matrix), 500)

heatmap.2(t(de_brca_lncRNA_sim_matrix[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          labRow=NA, labCol=NA, 
          trace='none', dendrogram='row',
          xlab='Gene', ylab='Gene',
          main='Similarity matrix',
          density.info='none', revC=TRUE)

## Construct adjacency matrix
de_brca_lncRNA_adj_matrix = adjacency.fromSimilarity(de_brca_lncRNA_sim_matrix, power=3, type='signed')

# Convert to matrix
de_brca_lncRNA_gene_ids = rownames(de_brca_lncRNA_adj_matrix)
de_brca_lncRNA_adj_matrix = matrix(de_brca_lncRNA_adj_matrix, 
                                    nrow=nrow(de_brca_lncRNA_adj_matrix))
rownames(de_brca_lncRNA_adj_matrix) = de_brca_lncRNA_gene_ids
colnames(de_brca_lncRNA_adj_matrix) = de_brca_lncRNA_gene_ids

heatmap.2(t(de_brca_lncRNA_adj_matrix[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          labRow=NA, labCol=NA, 
          trace='none', dendrogram='row',
          xlab='Gene', ylab='Gene',
          main='Adjacency matrix',
          density.info='none', revC=TRUE)

## --------------- COEXPRESSION ANALYSIS ----------------------------------------

# Cluster gene expression profiles; the flashClust function from
# the authors of WGCNA is another options for larger datasets.
# For input, we use the reciprocal of the adjacency matrix; hierarchical
# clustering works by comparing the _distance_ between objects instead of the
# _similarity_.
de_brca_lncRNA_gene_tree <- hclust(as.dist(1 - de_brca_lncRNA_adj_matrix), 
                                   method="average")

# we will use the cuttreeDynamicTree method to break apart the hc dendrogram
# into separate modules
de_brca_lncRNA_module_labels <- cutreeDynamicTree(dendro=de_brca_lncRNA_gene_tree, 
                                                  minModuleSize=15,
                                                  deepSplit=TRUE)

# assign a color to each module for easier visualization and referencing
de_brca_lncRNA_module_colors <- labels2colors(de_brca_lncRNA_module_labels)


# create geneinfo 
de_brca_lncRNA_geneInfo = data.frame(gene_ids = de_brca_lncRNA_gene_ids, gene_module = de_brca_lncRNA_module_labels, 
                                     gene_rgb = col2hex(de_brca_lncRNA_module_colors))

summary(upperTriangle(de_brca_lncRNA_adj_matrix, diag=FALSE, byrow=FALSE))
de_brca_lncRNA_graph <- export_network_to_graphml(de_brca_lncRNA_adj_matrix, filename='~/de_brca_lncRNA_network.graphml',
                               threshold=0.9, nodeAttrDataFrame=de_brca_lncRNA_geneInfo)

## --------------- EXPORT NETWORK FUNCTIONS ------------------------------------

export_network_to_graphml <- function (adj_mat, filename=NULL, weighted=TRUE,
                                       threshold=0.5, max_edge_ratio=3,
                                       nodeAttr=NULL, nodeAttrDataFrame=NULL,
                                       edgeAttributes=NULL, verbose=FALSE) {
  library('igraph')
  
  # Determine filename to use
  if (is.null(filename)) {
    filename='network.graphml'
  }
  
  # TODO 2015/04/09
  # Add option to rescale correlations for each module before applying
  # threshold (this is simpler than the previous approach of trying to
  # determine a different threshold for each module)
  #
  # Still, modules with very low correlations should be given somewhat
  # less priority than those with very high correlations.
  
  #module_colors <- unique(nodeAttrDataFrame$color)
  #module_genes <- which(nodeAttrDataFrame$color == color)
  #module_adjmat <- adj_mat[module_genes,]
  #num_genes <- length(module_genes)
  
  # Adjust threshold if needed to limit remaining edges
  max_edges <- max_edge_ratio * nrow(adj_mat)
  
  edge_to_total_ratio <- max_edges / length(adj_mat)
  edge_limit_cutoff <- as.numeric(quantile(abs(adj_mat), 1 - edge_to_total_ratio))
  
  # Also choose a minimum threshold to make sure that at least some edges
  # are left
  min_threshold <- as.numeric(quantile(abs(adj_mat), 0.9999))
  
  threshold <- min(min_threshold, max(threshold, edge_limit_cutoff))
  
  # Remove edges with weights lower than the cutoff
  adj_mat[abs(adj_mat) < threshold] <- 0
  
  # Drop any genes with no edges (TODO: Make optional)
  orphaned <- (colSums(adj_mat) == 0)
  adj_mat <- adj_mat[!orphaned, !orphaned]
  
  # Also remove annotation entries
  if (!is.null(nodeAttr)) {
    nodeAttr <- nodeAttr[!orphaned]
  }
  if (!is.null(nodeAttrDataFrame)) {
    nodeAttrDataFrame <- nodeAttrDataFrame[!orphaned,]
  }
  
  # Keep track of non-positive edges and rescale to range 0,1
  is_zero     <- adj_mat == 0
  is_negative <- adj_mat < 0
  
  adj_mat <- (abs(adj_mat) - threshold) / (max(adj_mat) - threshold)
  adj_mat[is_zero] <- 0
  adj_mat[is_negative] <- -adj_mat[is_negative]
  
  if (verbose) {
    message(sprintf("Outputting matrix with %d nodes and %d edges", 
                    nrow(adj_mat), sum(adj_mat > 0)))
  }
  
  # Create a new graph and add vertices
  # Weighted graph
  if (weighted) {
    g <- graph.adjacency(adj_mat, mode='undirected', weighted=TRUE, diag=FALSE)
  } else {
    adj_mat[adj_mat != 0] <- 1
    g <- graph.adjacency(adj_mat, mode='undirected', diag=FALSE)
  }
  
  # Add single node annotation from vector
  if (!is.null(nodeAttr)) {
    g <- set.vertex.attribute(g, "attr", value=nodeAttr)
  }
  
  # Add node one or more node annotations from a data frame
  if (!is.null(nodeAttrDataFrame)) {
    for (colname in colnames(nodeAttrDataFrame)) {
      g <- set.vertex.attribute(g, colname, value=nodeAttrDataFrame[,colname])
    }
  }
  
  edge_correlation_negative <- c()
  
  # neg_correlations[edge_list]
  edge_list <- get.edgelist(g)
  
  for (i in 1:nrow(edge_list)) {
    from <- edge_list[i, 1]    
    to   <- edge_list[i, 2]    
  }
  
  # Save graph to a file
  write.graph(g, filename, format='graphml')
  
  # return igraph
  return(g)
}












