setwd("/Volumes/files/EvoCell_lab/Xue/scRNA-seq_/reference")
alldata <- readRDS("seurat_with_celltype_annotation_update.rds")

library(monocle3)
library(SingleCellExperiment)
library(Seurat)
# I will subset stele cells for monole3 analysis 
# subset data 
k <- c("Initials","Phloem","Procambium","Xylem")
sub <- subset(alldata,idents=k)
data <- GetAssayData(sub,assay="RNA",slot="counts")
cell_metadata <- sub@meta.data
gene_metadata<- sub@assays$RNA@meta.features
gene_metadata$gene_short_name <- rownames(gene_metadata)
# generate monole3 object
cds <- new_cell_data_set(data,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
)
 library(tidyverse)
hvg <- sub@assays$RNA@var.features
cds <-  preprocess_cds(cds,num_dim = 30,method="PCA",norm_method="log",use_genes=hvg)
## Reduce the dimensions using UMAP
cds <- reduce_dimension(cds,umap.fast_sgd = TRUE, reduction_method = "UMAP", preprocess_method = "PCA")
# integrate seurat umap with cds
sobj_embed<- Seurat::Embeddings(sub, reduction = "umap")
 cds_embed <- cds@int_colData$reducedDims$UMAP
  cds@int_colData$reducedDims$UMAP <- sobj_embed[rownames(cds_embed),]
  plot_cells(cds, reduction_method = "UMAP", color_cells_by = "celltype",group_label_size=6)
  #cluster the cells 
cds <-cluster_cells(cds,reduction_method = "UMAP", cluster_method = "leiden",random_seed=1314)
plot_cells(cds, color_cells_by = "partition",group_label_size=6)
## order cells in pseudotime along a trajectory
cds <- learn_graph(cds)
plot_cells(cds,trajectory_graph_segment_size = 1.5,
           group_label_size = 6,group_cells_by="partition",
           color_cells_by = "celltype", label_groups_by_cluster = T,
           label_leaves = T, label_branch_points = T
)

library(ggplot2)
cds <- order_cells(cds)
h <- plot_cells(cds,
           group_label_size = 6,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 2
)
ggsave("Trajectory_pseudotime.pdf",plot=h,width=8,height=6)
############################################################################
# DEG analysis 
dea_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 7) %>%
    dplyr::filter(q_value < 0.05)
    write.csv(dea_res,"Trajectory_genes_stele.csv")
# select gene of interest
# customize the figure 
# change the gene name
which(row.names(cds)=="gene")
row.names(cds[1610]) <- "geneA"

h <- plot_genes_in_pseudotime(cds['geneA', ],label_by_short_name=FALSE,
                             
                              min_expr = 0.5, ncol = 2
)

df <- data.frame(tra1 = h$data$pseudotime,geneA_Expression = h$data$adjusted_expression,celltype=h$data$celltype)
df$celltype <- as.character(df$celltype)
df$celltype <- "vasculature"
library(reshape2)
df_melted <- melt(df, id.vars = c("geneA_Expression","celltype"), variable.name = "Trajectory", value.name = "Pseudotime")

c <- ggplot(df_melted, aes(x = Pseudotime, y = log(geneA_Expression + 1), group = Trajectory)) +
    geom_smooth(aes(linetype = Trajectory), method = "loess", se = FALSE, size = 1.5, show.legend = FALSE) +
    geom_point(aes(color = celltype), size = 2) +
    xlab("Pseudotime") +
    ylab("log counts") +
    theme_bw()

    c+theme(
    axis.text = element_text(colour = "black",size = 20, face = "bold"),
    axis.title = element_text(size =20, face = "bold"),
    legend.title = element_blank(),
    axis.text.x = element_text(size = 20, angle = 35, hjust = 1, vjust = 1),  # Set the size and angle of x-axis text
    plot.title = element_blank(),
    legend.position = "none")
ggsave("Trajectory_pseudotime_stele_geneA.pdf",plot=h,width=8,height=6)


# you can also look at the tendency for each cell type

c <- ggplot(df_melted, aes(x = Pseudotime, y = log(geneA_Expression + 1), color=celltype,group = celltype)) +
    geom_smooth(aes(linetype = celltype), method = "loess", se = FALSE, size = 1.5) +
    xlab("Pseudotime") +
    ylab("log counts") +
    theme_bw()
c+theme(
    axis.text = element_text(colour = "black",size = 20, face = "bold"),
    axis.title = element_text(size =20, face = "bold"),
    legend.title = element_blank(), legend.text = element_text(size =10, face = "bold"),
    axis.text.x = element_text(size = 20, angle = 35, hjust = 1, vjust = 1),  # Set the size and angle of x-axis text
    plot.title = element_blank(),
    legend.position = "right")

    # see if we have cell ourliers and remove them 
    m <- df_melted[df_melted$celltype=="Initials",]

    hist(m$Pseudotime)
m[which(m$Pseudotime >25),]
 v <- c(525,1074,1077)
 df_melted <- df_melted[!row.names(df_melted)%in%v,]
 c <- ggplot(df_melted, aes(x = Pseudotime, y = log(geneA_Expression + 1), color=celltype,group = celltype)) +
    geom_smooth(aes(linetype = celltype), method = "loess", se = FALSE, size = 1.5) +
    xlab("Pseudotime") +
    ylab("log counts") +
    theme_bw()
c+theme(
    axis.text = element_text(colour = "black",size = 20, face = "bold"),
    axis.title = element_text(size =20, face = "bold"),
    legend.title = element_blank(), legend.text = element_text(size =10, face = "bold"),
    axis.text.x = element_text(size = 20, angle = 35, hjust = 1, vjust = 1),  # Set the size and angle of x-axis text
    plot.title = element_blank(),
    legend.position = "right")

