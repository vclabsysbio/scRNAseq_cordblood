#loadlibrary
library(tidyverse) #manipulate table
library(Seurat) #process single cell data**main package** 
library(SeuratWrappers) #import loom file
library(velocyto.R) #RNA trajectory
library(biomaRt) #fetch gene info
library(ggplot2) #visualisation
library(gridExtra) #combine plot
library(grid) #combine plot2
library(sctransform) #normalise single cell data 
library(DoubletFinder) #eliminate doublet
library(fastcluster) #efficient clustering
library(future) #paralell processing
library(ComplexHeatmap) #heatmap creation
library(circlize) #color pallete
library(magrittr) #extract data from table
library(uwot) #need version 0.1.5
#set parameter
task_name="revision1"
directory="/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/cod/revised_analyses/"
#import 
integrated_seurat<-readRDS("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/TRUEn5000_cordblood_INTEGRATED_RNASCTnormalisedforanalysis.rds")
#A-define main clusters according to hemoglobin expression
integrated_seurat<-integrated_seurat%>%RenameIdents("0"="early_HBG","1"="late_HBB","2"="late_HBG","3"="early_HBB")
#modify cluster order
levels(integrated_seurat)<-c("early_HBG","late_HBG","early_HBB","late_HBB","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
#B-define HSC clusters according to CD34 expression
integrated_seurat<-integrated_seurat%>%RenameIdents("6"="HSC", "9"="HSC")
#modify cluster order
levels(integrated_seurat)<-c("early_HBG","late_HBG","early_HBB","late_HBB","HSC","4","5","7","8","10","11","12","13","14","15","16","17")
#C-define Megakaryocyte cluster according to PF4 expression
integrated_seurat<-integrated_seurat%>%RenameIdents("11"="Megakaryocyte")
#modify cluster order
levels(integrated_seurat)<-c("early_HBG","late_HBG","early_HBB","late_HBB","HSC","Megakaryocyte","4","5","7","8","10","12","13","14","15","16","17")
#D- define Myeloid Progenitor according to MPO, CD63, FCER1A expression
integrated_seurat<-integrated_seurat%>%RenameIdents("17"="Myeloid Progenitor","5"="Myeloid Progenitor","13"="Myeloid Progenitor","14"="Myeloid Progenitor")
#modify cluster order
levels(integrated_seurat)<-c("early_HBG","late_HBG","early_HBB","late_HBB","HSC","Megakaryocyte","Myeloid Progenitor","4","7","8","10","12","15","16")
#E-filter out non-RBC clusters
heamoglobin_expressed_dataset<-subset(integrated_seurat,idents=c("early_HBG","late_HBG","early_HBB","late_HBB","4","7","8","10","12","15","16"))
#recluster the dataset
DefaultAssay(heamoglobin_expressed_dataset)<-"integrated"
heamoglobin_expressed_dataset <- FindVariableFeatures(object = heamoglobin_expressed_dataset, dims = 1:30, assay = "integrated")
heamoglobin_expressed_dataset <- RunPCA(heamoglobin_expressed_dataset, dims = 1:30, assay = "integrated")
heamoglobin_expressed_dataset_before <- FindNeighbors(object = heamoglobin_expressed_dataset, dims = 1:30, assay = "integrated")
heamoglobin_expressed_dataset <- FindClusters(object = heamoglobin_expressed_dataset_before, resolution = 0.6, algorithm =1)
heamoglobin_expressed_dataset<- RunUMAP(heamoglobin_expressed_dataset, dims = 1:30, assay = "integrated")
heamoglobin_expressed_dataset<-NormalizeData(heamoglobin_expressed_dataset, assay = "RNA")

#define clusters
heamoglobin_expressed_dataset<-heamoglobin_expressed_dataset%>%RenameIdents("0"="late_HBB","6"="early_HBB","7"="early_HBB","2"="late_HBG", "5"="late_HBG", "8"="early_HBG", "3"="early_HBG", "1"="early_HBG","4"="early_HBG")
#modify cluster order
levels(heamoglobin_expressed_dataset)<-c("early_HBG","late_HBG","early_HBB","late_HBB","9","10","11","12")
#filter outgroup clusters
hemoglobin_containing_dataset<-subset(heamoglobin_expressed_dataset,idents=c("early_HBG","late_HBG","early_HBB","late_HBB"))
DimPlot(hemoglobin_containing_dataset)

saveRDS(hemoglobin_containing_dataset,"/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/cod/revised_analyses/hemoglobin_containing_dataset.rds")