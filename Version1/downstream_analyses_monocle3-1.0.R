#!/usr/bin/env Rscript

#loadlibrary
library(tidyverse)
library(Seurat)
library(monocle)
library(monocle3)
library(ggplot2)
library(gridExtra)
#set parameters
task_name="revision3"
directory="/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/cod/revised_analyses//"
#read seurat object
integrated_seurat<-readRDS(file="/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/hemoglobin_containing_dataset.rds")
#save current cluster into matadata
integrated_seurat$active_cluster<-Idents(integrated_seurat)
#Extract Seurat's log-transformed values
expression_matrix <- Seurat::GetAssayData(integrated_seurat, assay = "RNA", slot = "data")
#Extract Seurat meta_data
meta_data <- integrated_seurat@meta.data
#Extract gene names from Seurat object SCT slot to make CDS
seurat_genes <- data.frame(gene_short_name = rownames(integrated_seurat[["RNA"]]), row.names = rownames(integrated_seurat[["RNA"]]))
#make CDS object
new_cds <- monocle3::new_cell_data_set(expression_data = expression_matrix, cell_metadata = meta_data, gene_metadata = seurat_genes)
#make SCE object for data transfer
sce <- as.SingleCellExperiment(integrated_seurat, assay = "RNA")
#transfer reduction
SingleCellExperiment::reducedDims(new_cds) <- SingleCellExperiment::reducedDims(sce)
#transfer assay
SummarizedExperiment::assays(new_cds) <- SummarizedExperiment::assays(sce)
#transfer geneaname
SummarizedExperiment::rowData(new_cds) <- SummarizedExperiment::rowData(sce)
SummarizedExperiment::rowData(new_cds)$gene_short_name <-  row.names(new_cds)
#trasfer geneloading
new_cds@preprocess_aux$gene_loadings <- integrated_seurat@reductions[["pca"]]@feature.loadings.projected
#transfer cluster
new_cds@clusters[["UMAP"]]$clusters <- integrated_seurat@meta.data$integrated_snn_res.0.8
#The next two commands are run in order to allow "order_cells" to be run in monocle3
rownames(new_cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
colnames(SingleCellExperiment::reducedDims(new_cds)[["UMAP"]]) <- NULL
#transfer umap partition
new_cds@clusters$UMAP$partitions <- integrated_seurat@meta.data$active_cluster
#recluster the dataset
new_cds_leiden_01 = cluster_cells(new_cds,reduction_method = "UMAP",cluster_method = "leiden", resolution = 0.005)
new_cds_leiden_01 <- learn_graph(new_cds_leiden_01)
#select root nodes manually
new_cds_leiden_01 <- order_cells(new_cds_leiden_01)

new_cds_leiden_01<-readRDS("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/cod/revised_analyses_March/new_cds_leiden_01.rds")
#print scatterplot-pseudotime

png(filename = paste0("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/cod/revised_analyses//pseudotime_leiden_01.2.png"), width=600, height=500)
plot_cells(new_cds_leiden_01,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           show_trajectory_graph = FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           #trajectory_graph_segment_size = 3,
           graph_label_size=3)
dev.off()

#figure 4A
pseudo_plot<-plot_cells(new_cds_leiden_01,
                        color_cells_by = "pseudotime",
                        label_cell_groups=FALSE,
                        show_trajectory_graph = FALSE,
                        label_leaves=FALSE,
                        label_branch_points=FALSE,
                        label_roots = FALSE,
                        #trajectory_graph_segment_size = 3,
                        graph_label_size=3)
cluster_plot <- DimPlot(object = integrated_seurat, reduction = "umap", label = FALSE, cols = c("cyan","blue","hotpink", "red"))

png(filename = paste0("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/cod/revised_analyses//pseudotime_leiden_01.4withTree.png"), width=1100, height=500)
print(CombinePlots(plots = list(cluster_plot,pseudo_plot), ncol = 2))
dev.off()

#extract expression table
expression_table<-GetAssayData(object = integrated_seurat, assay = "RNA")%>%data.frame%>%rownames_to_column()
#select genes of interest
selected_gene_table<-expression_table%>%filter(rowname %in% c("HBA1","HBB","HBG2","SLC4A1"))%>%column_to_rownames()%>%t%>%data.frame()%>%rownames_to_column
#save cluster 
integrated_seurat$active_cluster<-Idents(integrated_seurat)
#create table of cluster for later annotation
info_table<-integrated_seurat$active_cluster%>%data.frame(active_cluster=.)%>%rownames_to_column()
#exteact pseudotime parameter from cds
hemoglobin_genes <- c("HBA1","HBB","HBG2","SLC4A1")
hemoglobin_cds <- new_cds_leiden_01[rowData(new_cds_leiden_01)$gene_short_name %in% hemoglobin_genes,]
pT_leiden<-data.frame(pseudotime_coordinate=pseudotime(hemoglobin_cds))
#combine tables
expression_table_with_pseudotime_leiden<-pT_leiden%>%rownames_to_column%>%left_join(selected_gene_table,by="rowname")%>%left_join(info_table,by="rowname")
expression_table_with_pseudotime_leiden_mod1<-expression_table_with_pseudotime_leiden%>%separate("active_cluster", into = c("time","marker"),remove = FALSE)
#figure 4B
plot_pseudo_1<-ggplot(expression_table_with_pseudotime_leiden_mod1, aes(x=pseudotime_coordinate,y=HBA1,color=active_cluster))+
  geom_point(size=0.25)+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(values=setNames(cols_pll,cluster_labels))+
  theme_classic()+
  geom_smooth(data = expression_table_with_pseudotime_leiden_mod1,
              aes(x=pseudotime_coordinate,y=HBA1),
              method = "loess",
              inherit.aes = FALSE)
plot_pseudo_facet<-ggplot(expression_table_with_pseudotime_leiden_mod1, aes(x=pseudotime_coordinate,y=HBA1,color=active_cluster))+
  geom_point(size=0.25)+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(values=setNames(cols_pll,cluster_labels))+
  theme_classic()+
  geom_smooth(data = expression_table_with_pseudotime_leiden_mod1,
              aes(x=pseudotime_coordinate,y=HBA1),
              method = "loess",
              inherit.aes = FALSE)+
  facet_grid(cols = vars(marker))

png(filename = paste0("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/cod/revised_analyses//pseudoleiden_HBA1_version3.png"), width=1100, height=250)
grid.arrange(
  grobs = list(plot_pseudo_1, plot_pseudo_facet),
  widths = c(1, 1.6)
)
dev.off()

#figure 4C
plot_pseudo_1<-ggplot(expression_table_with_pseudotime_leiden_mod1, aes(x=pseudotime_coordinate,y=SLC4A1,color=active_cluster))+
  scale_color_manual(values=setNames(cols_pll,cluster_labels))+
  geom_point(size=0.25)+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_classic()+
  geom_smooth(data = expression_table_with_pseudotime_leiden_mod1,
              aes(x=pseudotime_coordinate,y=SLC4A1),
              method = "loess",
              inherit.aes = FALSE)
plot_pseudo_facet<-ggplot(expression_table_with_pseudotime_leiden_mod1, aes(x=pseudotime_coordinate,y=SLC4A1,color=active_cluster))+
  scale_color_manual(values=setNames(cols_pll,cluster_labels))+
  geom_point(size=0.25)+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_classic()+
  geom_smooth(data = expression_table_with_pseudotime_leiden_mod1,
              aes(x=pseudotime_coordinate,y=SLC4A1),
              method = "loess",
              inherit.aes = FALSE)+
  facet_grid(cols = vars(marker))

png(filename = paste0("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/cod/revised_analyses//pseudoleiden_SLC4A1_version3.png"), width=1100, height=250)
grid.arrange(
  grobs = list(plot_pseudo_1, plot_pseudo_facet),
  widths = c(1, 1.6)
)
dev.off()
