setwd("/media/ducdo/UUI1/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
source("/media/ducdo/UUI1/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_correlation_analysis/helper_functions.R")
source("/media/ducdo/UUI1/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616/code_network_construction/network_helper_functions.R")

require(igraph); require(RCy3); require(plyr); require(dplyr); require(data.table); require(rlist); require(ggplot2)

load("data_Saved_R_Objects/encRNA_network/encRNA_network_statistic.rda")
load("brca2_TCGA2STAT.rda")
load("data_Saved_R_Objects/miRNA_target/predicted_normal_tumor_goodCorr.rda")
load("data_Saved_R_Objects/brca_df.rda")

load("data_Saved_R_Objects/corr_matrices/normal_encRNA_850356.rda")
load("data_Saved_R_Objects/corr_matrices/tumor_encRNA_850356.rda")



slices = sort(table(df1$miRNA[which(df1$miRNA > 0)]), decreasing = T)
lbls <- names(sort(table(df1$miRNA[which(df1$miRNA > 0)]), decreasing = T))

slices = sort(table(df2$miRNA[which(df2$miRNA > 0)]), decreasing = T)
lbls <- names(sort(table(df2$miRNA[which(df2$miRNA > 0)]), decreasing = T))
pie(slices, labels = lbls, main="miRNAs regulating normal ceRNA", col=brewer.pal(length(lbls),"Set3"))
pie(slices, labels = lbls[1:7], main="miRNAs regulating normal ceRNA", col=brewer.pal(length(lbls),"Set3"))

##########################################################################################
#########                       Network statistics                               ######### 
##########################################################################################

normal_encRNA_graph = build_iGraph(df = normal_encRNA_sensitivity_bound_goodCoor)
vcount(normal_encRNA_graph); ecount(normal_encRNA_graph)
tumor_encRNA_graph = build_iGraph(df = tumor_encRNA_sensitivity_bound_goodCoor)
vcount(tumor_encRNA_graph); ecount(tumor_encRNA_graph)
## ----- basic network statistic ------------------------- 

load("data_Saved_R_Objects/encRNA_network/encRNA_network_statistic.rda")

# vertice and edge count
vcount(normal_encRNA_graph); ecount(normal_encRNA_graph) # 359, 1907
V(normal_encRNA_graph)[which(V(normal_encRNA_graph)$type == "mRNA")] # 301
V(normal_encRNA_graph)[which(V(normal_encRNA_graph)$type == "lncRNA")] # 58

vcount(tumor_encRNA_graph); ecount(tumor_encRNA_graph) # [1] 1302, 1218
V(tumor_encRNA_graph)[which(V(tumor_encRNA_graph)$type == "mRNA")] # 956
V(tumor_encRNA_graph)[which(V(tumor_encRNA_graph)$type == "lncRNA")] # 346

# select hub genes
normal_mRNA_degree = get_node_degree(normal_encRNA_graph, "mRNA")
normal_lncRNA_degree = get_node_degree(normal_encRNA_graph, "lncRNA")
tumor_mRNA_degree = get_node_degree(tumor_encRNA_graph, "mRNA")
tumor_lncRNA_degree = get_node_degree(tumor_encRNA_graph, "lncRNA")

normal_mRNA = as.data.frame(normal_mRNA_degree)
ggplot(normal_mRNA, aes(normal_mRNA$normal_mRNA_degree)) + 
    geom_histogram(fill = "darkgreen", bins = length(unique(normal_mRNA$normal_mRNA_degree)), binwidth = 0.5) +
    xlab("Degree") + ylab("count")



par(mfrow = c(2,2))
hist(normal_mRNA_degree, breaks = length(unique(normal_mRNA_degree)))
hist(normal_lncRNA_degree, breaks = length(unique(normal_lncRNA_degree)))
hist(tumor_mRNA_degree, breaks = length(unique(tumor_mRNA_degree)))
hist(tumor_lncRNA_degree, breaks = length(unique(tumor_lncRNA_degree)))
par(mfrow = c(1,1))

# graph density
graph.density(normal_encRNA_graph) # [1] 0.02967585
graph.density(tumor_encRNA_graph) # [1] 0.0014381

# node betweenness
normal_mRNA_betweeness = get_node_betweeness(normal_encRNA_graph, "mRNA")
normal_lncRNA_betweeness = get_node_betweeness(normal_encRNA_graph, "lncRNA")
tumor_mRNA_betweeness = get_node_betweeness(tumor_encRNA_graph, "mRNA")
tumor_lncRNA_betweeness = get_node_betweeness(tumor_encRNA_graph, "lncRNA")

par(mfrow = c(2,2))
hist(normal_mRNA_betweeness, breaks = length(unique(normal_mRNA_betweeness)))
hist(normal_lncRNA_betweeness, breaks = length(unique(normal_lncRNA_betweeness)))
hist(tumor_mRNA_betweeness, breaks = length(unique(tumor_mRNA_betweeness)))
hist(tumor_lncRNA_betweeness, breaks = length(unique(tumor_lncRNA_betweeness)))
par(mfrow = c(1,1))


# save(normal_mRNA_degree, normal_lncRNA_degree,tumor_mRNA_degree,tumor_lncRNA_degree, 
#      normal_mRNA_betweeness, normal_lncRNA_betweeness, tumor_mRNA_betweeness, tumor_lncRNA_betweeness,
#      file = "data_Saved_R_Objects/encRNA_network/encRNA_network_statistic.rda")

intra_edges = E(normal_encRNA_graph) [ from("ENSG00000229852.2") ]
V(normal_encRNA_graph)[get.edges(normal_encRNA_graph,intra_edges)[,1]]

################# get differential mRNAS and lncRNA ####################################
setwd("/media/ducdo/UUI1/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
load("data_Saved_R_Objects/differential_analysis/de_lncRNA_mRNA_miRNA.rda")

dim(brca_lncRNA_df) # [1] 4828    6
dim(brca_lncRNA_limma_all) # [1] 4828    6
sum(brca_lncRNA_limma_all$adj.P.Val < 0.05) # 3148
brca_lncRNA_limma_all_subset = subset(brca_lncRNA_limma_all, brca_lncRNA_limma_all$adj.P.Val < 0.05)
dim(brca_lncRNA_limma_all_subset) # [1] 3148    6
sum(brca_lncRNA_limma_all_subset$logFC < 0) # 2023 lncRNA genes are down-regulated
sum(brca_lncRNA_limma_all_subset$logFC > 0) # 1125 lncRNA genes are up-regulated 

# check overlap between differentially regulated lncRNAs with high degree lncRNAs

# check overlap between differentially regulated lncRNAs with high betweeness lncRNAs

# check overlap between differentially regulated lncRNAs with high hazard lncRNAs

# check if high degree lncRNA relates with its gene expression trend compared with the whole population
hist(normal_lncRNA_degree)
# box plot of expression of all lncRNA in the normal_degree
expression_all_lncRNA_network = brca_lncRNA_df[names(normal_lncRNA_degree),]

avg_expression_all_lncRNA = apply(brca_lncRNA_df, 1, mean)
avg_expression_all_lncRNA_network = apply(expression_all_lncRNA_network, 1, mean)
avg_expression_all_de_lncRNA = apply(brca_lncRNA_df[rownames(brca_lncRNA_limma_all_subset),], 1, mean)
avg_expression_all_up_lncRNA = apply(brca_lncRNA_df[rownames(brca_lncRNA_limma_all_subset)[which(brca_lncRNA_limma_all_subset$logFC > 0)],],1,mean)

d1 = data.frame(group = "all lncRNA", value = unname(avg_expression_all_lncRNA))
d2 = data.frame(group = "lncRNAs in network", value = unname(avg_expression_all_lncRNA_network))
d3 = data.frame(group = "d.e. lncRNAs ", value = unname(avg_expression_all_de_lncRNA))
d4 = data.frame(group = "up lncRNAs ", value = unname(avg_expression_all_up_lncRNA))
d = rbind(d1,d2,d3,d4)

require(ggplot2)
ggplot(d, aes(x = group, y = value, fill = group)) + geom_boxplot()

##########################################################################################
#########                       Check prevalent miRNAs                           ######### 
##########################################################################################

normal_encRNA_graph = build_iGraph(df = normal_encRNA_sensitivity_bound_goodCoor)
vcount(normal_encRNA_graph); ecount(normal_encRNA_graph)
tumor_encRNA_graph = build_iGraph(df = tumor_encRNA_sensitivity_bound_goodCoor)
vcount(tumor_encRNA_graph); ecount(tumor_encRNA_graph)

df1 = normal_encRNA_sensitivity_bound_goodCoor; df1$miRNA = as.character(df1$miRNA)
length(unique(df1$miRNA))
table(df1$miRNA[which(df1$miRNA > 0)])

df2 = tumor_encRNA_sensitivity_bound_goodCoor; df2$miRNA = as.character(df2$miRNA)
length(unique(df2$miRNA))
table(df2$miRNA[which(df2$miRNA > 0)])

intersect(names(sort(table(df1$miRNA[which(df1$miRNA > 0)]))),names(sort(table(df2$miRNA[which(df2$miRNA > 0)]))))

