#!/usr/bin/env Rscript
#loadlibrary
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(velocyto.R)
library(biomaRt)
library(ggplot2)
library(gridExtra)
library(grid)
library(sctransform)
library(DoubletFinder)
library(fastcluster)
library(future)
library(ComplexHeatmap)
library(circlize)
#set parameters
task_name="revision1c"
directory="/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/cod/revised_analyses/"
#load dataset
integrated_seurat<-readRDS(file="/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/hemoglobin_containing_dataset.rds")
#figure 1C
#clsuter overview
name = paste0(directory,task_name)
png(filename = paste0(name,"_clustermap_umap.png"), width=600, height=500)
p <- DimPlot(object = integrated_seurat,cols = c("cyan","blue","hotpink", "red"), reduction = "umap")
print(CombinePlots(plots = list(p), ncol = 2))
dev.off()
#figure 1D
#cell counting
#extract dataset metadata
integrated_seurat$active_cluster<-Idents(integrated_seurat)
dataset_metadata<-integrated_seurat@meta.data
#cell counting by sample (time)
dataset_metadata%>%group_by(active_cluster, time)%>%tally()%>%mutate(percentage=n/sum(n)*100)->stat_metadata
#create figure
p<-ggplot(stat_metadata%>%mutate_at(vars(percentage), round, digits =2),aes(x=time,y=percentage,fill=time))+
  geom_bar(stat="identity",position=position_dodge())+
  #scale_fill_manual(values=ucols)+
  ylim(0,60)+
  geom_text(aes(label=n), vjust=1.6, color="white", size=3.5,position = position_dodge(0.9))+
  geom_text(aes(label=percentage), vjust=-0.3, color="black", size=3.5,position = position_dodge(0.9))+
  labs(title= "grouped by time",
       subtitle = "n=11200")+
  theme(axis.title.x = element_blank())+
  facet_grid(cols=(vars(active_cluster)))
#print figure
png(paste0(directory,"/figure1c.png"), width=900, height=500,res = 150)
grid.arrange(
  grobs = list(p)
)
dev.off()

#gene expression overview
#select data slot
DefaultAssay(integrated_seurat)<-"RNA"
#define filling pattern
cols_pll = c("cyan","blue","hotpink", "red")
#figure 1E
#scatter plot
name = paste0(directory,task_name)
png(filename = paste0(name,"_fidure1d_7markers.png"), width=1200, height=1200)
print(FeaturePlot(object = integrated_seurat, features = c("CD34", "HBB","HBG1","HBG2","HBA1","SLC4A1","GYPA"), min.cutoff = "q9", cols = c("lightgrey", "blue"),reduction = "umap", pt.size = 0.5))
dev.off()
#figure 1F
#violin plot
name = paste0(directory,task_name)
png(filename = paste0(name,"_fidure1e_7markers.png"), width=1200, height=1200)
print(VlnPlot(object = integrated_seurat, features = c("CD34", "HBB","HBG1","HBG2","HBA1","SLC4A1","GYPA"),cols = cols_pll, pt.size = 0.5))
dev.off()

