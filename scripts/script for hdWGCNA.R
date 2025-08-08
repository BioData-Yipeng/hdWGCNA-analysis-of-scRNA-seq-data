
# https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html
# This script is to run WGCNA with single cell data, codes were modified from above link.
#prepare the environment
getwd()
setwd("F:/data science/R/6_GSE 133433 hdWGCNA")
library(dplyr)
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)

#install some libraries if needed
install.packages(c('UCell', 'GeneOverlap'))
install.packages("devtools")
devtools::install_github("smorabit/hdWGCNA")
remotes::install_github("carmonalab/UCell")
BiocManager::install("GeneOverlap")

# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# set random seed for reproducibility
set.seed(12345)
# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

#the showcase dataset is from GSE133433
#run standard analysis of scRNA seq data until the step of UMAP
wtfemale.data <- Read10X(data.dir ="data/wtfemale/" )
fivexfemale.data <- Read10X(data.dir ="data/5xfemale/" )
wtfemale <- CreateSeuratObject(counts = wtfemale.data,project = "wtfemale")
fivexfemale <- CreateSeuratObject(counts = fivexfemale.data,project = "5xfemale")
dim(wtfemale)
dim(fivexfemale)

merged_seurat <- merge(wtfemale, y=list(fivexfemale),
                     add.cell.ids=c("wtfemale", "5xfemale"),
                     project = "merged")
table(merged_seurat$orig.ident)

merged_seurat$condition <- ifelse(merged_seurat$orig.ident %in% c("wtfemale"), "wt", "5x") #add metadata to distinguish groups
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
variable_features <- VariableFeatures(merged_seurat)
#check variable features if needed
#top10 <- head(variable_features,10)
#p1 <- VariableFeaturePlot(merged_seurat)
#p2 <- LabelPoints(p1,top10,repel = T)
#p2
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat) #pca is run on the 2000 most variable features
ElbowPlot(merged_seurat,ndims=50) #check variance of each pc

merged_seurat <- FindNeighbors(merged_seurat,
                             dims = 1:15,
                             k.param = 30,
                             reduction = "pca")
merged_seurat <- FindClusters(merged_seurat,
                            resolution = 0.4,
                            algorithm = 4,
                            random.seed = 123,
                            cluster.name = "MergedCluster")
merged_seurat <- RunUMAP(merged_seurat,
                       dims=1:15,
                       reduction = "pca",
                       reduction.name = "MergedUMAP")

levels(Idents(merged_seurat)) #check how many clusters were returned

DimPlot(merged_seurat,
        reduction = "MergedUMAP",
        group.by = "condition") #this is used to visualize batch effects

Integrated_seurat <- IntegrateLayers(object = merged_seurat,
                                   method = HarmonyIntegration,
                                   orig.reduction = "pca",
                                   new.reduction = "harmony",
                                   verbose=TRUE) #harmony method of integration

Integrated_seurat <- FindNeighbors(Integrated_seurat,
                                 reduction = "harmony",
                                 dims = 1:15,
                                 k.param = 30) #dims, k.param, algorithm all need to adjust to produce reasonable cluster

Integrated_seurat <- FindClusters(Integrated_seurat,
                                resolution=0.1,
                                algorithm = 4,
                                random.seed = 123,
                                cluster.name = "harmonycluster",
                                verbose=TRUE)
levels(Idents(Integrated_seurat))
Integrated_seurat <- RunUMAP(Integrated_seurat,
                           reduction = "harmony",
                           dims = 1:15,
                           reduction.name = "harmonyUMAP")

#visualization clusters in wt vs 5X
DimPlot(Integrated_seurat,
        reduction = "harmonyUMAP",
        group.by = "harmonycluster",
        split.by ="condition",
        label = TRUE) + ggtitle("clusters wt vs 5X")

table(Integrated_seurat$orig.ident)

#visualization and compare effect of integration
beforeintegrate <- DimPlot(merged_seurat,
        reduction = "MergedUMAP",
        group.by = "condition") + ggtitle("before integration")
afterintegrate <- DimPlot(Integrated_seurat,
        reduction = "harmonyUMAP",
        group.by = "condition") + ggtitle("after integration")
beforeintegrate+afterintegrate #after integration, cells from the same cluster gather together.

#now start WGCNA analysis steps
Integrated_seurat <- JoinLayers(Integrated_seurat)
wgcna_obj <- SetupForWGCNA(Integrated_seurat,
                           gene_select = "fraction",
                           fraction = 0.2, #gene with express in 50% cells were selected
                           wgcna_name = "microglia")

# construct metacells  in each group. metacells would save calculation resources for correlation
wgcna_obj <- MetacellsByGroups(
  seurat_obj = wgcna_obj,
  group.by = c("seurat_clusters", "condition"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'seurat_clusters' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
wgcna_obj <- NormalizeMetacells(wgcna_obj)
metacell_obj <- GetMetacellObject(wgcna_obj)
wgcna_obj <- SetDatExpr(
  wgcna_obj,
  group_name = "2", # the name of the group of interest in the group.by column
  group.by='seurat_clusters', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)

# Test different soft powers:
wgcna_obj <- TestSoftPowers(
  wgcna_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(wgcna_obj)
gc()
# assemble with patchwork
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(wgcna_obj)
head(power_table)

# construct co-expression network:
wgcna_obj <- ConstructNetwork(
  wgcna_obj,
  tom_name = '2' # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(wgcna_obj, main='cluster 2 hdWGCNA Dendrogram')

# need to run ScaleData first or else harmony throws an error:
#seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

rm(plot_list,power_table,wtfemale.data)
# compute all MEs in the full single-cell dataset
wgcna_obj <- ModuleEigengenes(
  wgcna_obj,
  group.by.vars="condition"
)

# harmonized module eigengenes:
hMEs <- GetMEs(wgcna_obj)

# module eigengenes:
MEs <- GetMEs(wgcna_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
wgcna_obj <- ModuleConnectivity(
  wgcna_obj,
  group.by = 'seurat_clusters', group_name = '2'
)

# plot genes ranked by kME for each module
p <- PlotKMEs(wgcna_obj, ncol=3)
p

# get the module assignment table:
# the module table list the most coexpressed genes calcualted from cluster 2 in each module sequenced by kme
modules <- GetModules(wgcna_obj) %>% subset(module != 'grey')

#obtain the most 100 coexpressed genes in each module analyzed from cluster 2
most_coexpr100 <- GetHubGenes(wgcna_obj, n_hubs = 50)
most_coexpr100$genes <- substr(most_coexpr100$gene_name, 8, nchar(most_coexpr100$gene_name))

#build the connect network

data_expr <- GetAssayData(metacell_obj,slot = "data")
subset_expr <- data_expr[most_coexpr100$gene_name,]
subset_expr <-as.data.frame(subset_expr)

write.csv(subset_expr,"mostcoexpressgenefrom3modules.csv")

subset_expr <- read.csv("mostcoexpressgenefrom3modules.csv", row.names = 1)

subset_expr_t <- t(subset_expr)

class(subset_expr_t)

# Pick soft-thresholding power
powers <- c(1:20)
sft <- pickSoftThreshold(subset_expr, powerVector = powers, verbose = 5)

# Construct adjacency matrix
adj_matrix_a <- adjacency(subset_expr_t, power = 8)  # Use the optimal power from pickSoftThreshold

# Calculate topological overlap matrix (TOM)
TOM_a <- TOMsimilarity(adj_matrix_a)

# Convert TOM into a dissimilarity measure for hierarchical clustering
dissTOM_a <- 1 - TOM_a
# Cluster genes into modules
gene_tree <- hclust(as.dist(dissTOM_a), method = "average")

# Dynamic tree cut to detect modules
dynamic_modules <- cutreeDynamic(dendro = gene_tree, method = "tree", minClusterSize = 30)

library(pheatmap)
# Plot heatmap of the topological overlap matrix (TOM)
pheatmap(TOM_a, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, show_colnames = FALSE,
         color = colorRampPalette(c("white", "blue"))(50))

# Create the gene-module association
moduleColors <- labels2colors(dynamic_modules)

# Create a module network graph
module_network <- moduleEigengenes(subset_expr_t, moduleColors)

# Plot the module network
plotDendroAndColors(gene_tree, moduleColors, "Module colors", main = "Gene Dendrogram and Module Colors")

# If you prefer to visualize the network in R, you can use igraph
library(igraph)
gene_network <- graph_from_adjacency_matrix(TOM_a, mode = "undirected", weighted = TRUE, diag = FALSE)
plot(gene_network)

# get hub genes for the following analysis
hub_df <- GetHubGenes(wgcna_obj, n_hubs = 10)

head(hub_df)

saveRDS(wgcna_obj, file='hdWGCNA_object.rds')

#extract gene expression matrix of hub_df genes for the thrown coexpression modules
hub_exr <- GetAssayData(metacell_obj,slot = "data")[hub_df$gene_name,]
hub_exr <-as.data.frame(hub_exr)
turquoise <- hub_exr[1:10,]
blue <- hub_exr[11:20,]
brown <- hub_exr[21:30,]

rownames(turquoise)
turquoise$genes <- rownames(turquoise)
# Reshape to long format
turquoise_long <- turquoise %>%
  pivot_longer(cols = -genes, names_to = "Sample", values_to = "Expression")

#color the plot by genes
ggplot(turquoise_long, aes(x = Sample, y = Expression, group = genes, color = genes)) +
  geom_line(size = 0.2) +
  geom_point(size = 0.1) +
  theme_minimal() +
  labs(title = "Gene Expression Across Samples", x = "Sample", y = "Expression Level")


blue$genes <- rownames(blue)
# Reshape to long format
blue_long <- blue %>%
  pivot_longer(cols = -genes, names_to = "Sample", values_to = "Expression")
blue_long <- blue_long %>%
  mutate(Group = case_when(
    grepl("5x", Sample) ~ "5x",
    grepl("wt", Sample) ~ "wt"
  ))
#color the plot by Groups to indicate which group the sample come from
ggplot(blue_long, aes(x = Sample, y = Expression, group = genes, color = Group)) +
  geom_line(size = 0.2) +
  geom_point(size = 0.1) +
  theme_minimal() +
  labs(title = "Gene Expression Across Samples", x = "Sample", y = "Expression Level")


brown$genes <- rownames(brown)
# Reshape to long format
brown_long <- brown %>%
  pivot_longer(cols = -genes, names_to = "Sample", values_to = "Expression")
brown_long <- brown_long %>%
  mutate(Group = case_when(
    grepl("5x", Sample) ~ "5x",
    grepl("wt", Sample) ~ "wt"
  ))
ggplot(brown_long, aes(x = Sample, y = Expression, group = genes, color = Group)) +
  geom_line(size = 0.2) +
  geom_point(size = 0.1) +
  theme_minimal() +
  labs(title = "Gene Expression Across Samples", x = "Sample", y = "Expression Level")

#The previous are from the metacell dataset 3205 metacell, now let's see what it looks like in the raw dataset 11861 cells
#subset the matrix include the hubgene
hub_exr_raw <- GetAssayData(wgcna_obj,slot = "data")[hub_df$gene_name,]
hub_exr_raw <-as.data.frame(hub_exr_raw)
brown_raw <- hub_exr[21:25,]
brown_raw$genes <- rownames(brown_raw)
# Reshape to long format
brown_long_raw <- brown_raw %>%
  pivot_longer(cols = -genes, names_to = "Sample", values_to = "Expression")
brown_long_raw <- brown_long_raw %>%
  mutate(Group = case_when(
    grepl("5x", Sample) ~ "5x",
    grepl("wt", Sample) ~ "wt"
  ))
#add a column including cluster and group infor
brown_long_raw <- brown_long_raw %>%
  mutate(clusterandsample = substr(Sample, 1, 4))
ggplot(brown_long_raw, aes(x = Sample, y = Expression, group = genes, color = clusterandsample)) +
  geom_line(size = 0.2) +
  geom_point(size = 0.1) +
  theme_minimal() +
  labs(title = "Gene Expression Across Samples", x = "Sample", y = "Expression Level")+
  guides(color = guide_legend(override.aes = list(size = 3)))
#the results show that those coexpressed genes change according to clusters.
rownames(brown_raw)
#featrueplot show these genes are widely expressed and all higher in cluster 2 wt
plot <- FeaturePlot(Integrated_seurat,
            reduction = "harmonyUMAP",
            features = "GRCh38-TMSB4X",
            cols = c("gray","red"),
            label = TRUE,
            split.by = "condition",
            combine = FALSE)
plot[[1]] <- plot[[1]] + theme(legend.position = "right")
patchwork::wrap_plots(plot)

