## ----- helper functions -------------------------
build_iGraph = function(df){
  require(data.table)
  encRNA_df = df[,c("mRNA","lncRNA","sensitivity")]
  nrow(encRNA_df)
  encRNA_df$pair = paste(encRNA_df$mRNA, encRNA_df$lncRNA, sep = "~")
  
  # there are duplicated pairs in encRNA_df. In such a case, compute the average 
  # sensitivity of rows from a similar pair
  keys = colnames(encRNA_df)[!grepl('sensitivity', colnames(encRNA_df))]
  X = as.data.table(encRNA_df)
  encRNA_df =X[,list(sensitivity = mean(sensitivity)),keys] # this code is ought be changed based on how sensitivity will be computed 
  encRNA_df = as.data.frame(encRNA_df)
  nrow(encRNA_df)
  
  # start constructing the igraph object
  node_type_df = data.frame(node = c(unique(encRNA_df$mRNA), unique(encRNA_df$lncRNA)))
  node_type_df$type = c(rep("mRNA", length(unique(encRNA_df$mRNA))),
                        rep("lncRNA", length(unique(encRNA_df$lncRNA)))) 
  graph_encRNA = graph.data.frame(d = encRNA_df[,c("mRNA","lncRNA")], 
                                  vertices = node_type_df,
                                  directed = F)
  V(graph_encRNA)$color <- ifelse(V(graph_encRNA)$type == "mRNA", "red", "blue")
  E(graph_encRNA)$weight = encRNA_df[,"sensitivity"]
  return(graph_encRNA)
}

getAllSubgraphsFromMembership = function(igraphObject = NULL, membership = NULL){
  list_of_all_subgraphs = list()
  for (i in sort(unique(membership))){
    list_of_all_subgraphs[[i]] = induced.subgraph(igraphObject,
                                                  which(membership==i))
  }
  return(list_of_all_subgraphs)
}

get_all_human_ensemble = function(){
  require(biomaRt)
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 85)
  listDatasets(ensembl)
  mart = useDataset("hsapiens_gene_ensembl", mart= ensembl)
  all_ensembl_gene_symbols = getBM(attributes = c("hgnc_symbol"),
                                   values = "*", mart = mart)
  length(all_ensembl_gene_symbols$hgnc_symbol) # 56820
  return(all_ensembl_gene_symbols)
}

get_mRNAs_in_Ensembl = function(geneList, ensemble, minimumSize){
  encRNA_mRNA_list = lapply(geneList, function(gene_list){
    gene_list = gene_list[which(gene_list %in% ensemble$hgnc_symbol)]; 
  })
  encRNA_mRNA_size = unlist(lapply(encRNA_mRNA_list, function(gene_list){
    length(gene_list)
  }))
  if(length(which(encRNA_mRNA_size < minimumSize)) != 0){
    encRNA_mRNA_list = encRNA_mRNA_list[-which(encRNA_mRNA_size < minimumSize)]
  }
  return(encRNA_mRNA_list)
}

get_RNAs_sharing_common_miRNAs = function(encRNA_df, miRNAs){
  require(rlist)
  mRNA_list = lncRNA_list = list()
  for (i in 1:length(miRNAs)){
    mRNAs = encRNA_df$mRNA[which(encRNA_df$miRNA == miRNAs[i])]
    mRNA_list = rlist::list.append(mRNA_list,unique(mRNAs))
    lncRNAs = encRNA_df$lncRNA[which(encRNA_df$miRNA == miRNAs[i])]
    lncRNA_list = rlist::list.append(lncRNA_list,unique(lncRNAs))
  }
  names(mRNA_list) = names(lncRNA_list) = miRNAs
  list = list(mRNA_list, lncRNA_list)
  names(list) = c("mRNA", "lncRNA")
  return(list)
}

get_TopGo_result = function(all_ensembl_gene_symbols, gene_list, topNodes){
  require(topGO)
  universe = factor(as.integer(all_ensembl_gene_symbols$hgnc_symbol %in% gene_list))
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
  allRes <- GenTable(GOdata, 
                     classicFisher = resultFis, 
                     orderBy = "classicFisher", 
                     ranksOf = "classicFisher", 
                     topNodes = topNodes)
  allRes$Group = AnnotationDbi::Ontology(allRes$GO.ID)
  return(allRes)
}

get_node_betweeness = function(igraph_object, RNA_type){
  betweenness = igraph::betweenness(igraph_object)
  return(sort(betweenness[which(V(igraph_object)$type == RNA_type)],decreasing = T))
}

get_node_degree = function(igraph_object, RNA_type){
  degree = igraph::degree(igraph_object)[which(V(igraph_object)$type == RNA_type)]
  return(sort(degree,decreasing = T))
}

get_mRNAs_connected_to_lncRNA = function(igraph_object, lncRNA_vector){
  require(rlist)
  result = list()
  for(i in 1:length(lncRNA_vector)){
    intra_edges = E(igraph_object) [ from(lncRNA_vector[i])]
    mRNAs = unique(names(V(igraph_object)[get.edges(igraph_object,intra_edges)[,1]]))
    result = list.append(result, mRNAs)
  }
  names(result) = lncRNA_vector
  return(result)
}

