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
library(future)
future::plan(strategy = 'multicore', workers = 12)
options(future.globals.maxSize = 10 * 1024 ^ 3)
options(warn = -1)
#define function for repeated activity
#read .loom file
import_loom<-function(){
  print("read .loom file...")
  CB_D4.ldat <- ReadVelocity(file = "/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/D4.loom")
  CB_D6.ldat <- ReadVelocity(file = "/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/D6.loom")
  CB_D8.ldat <- ReadVelocity(file = "/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/D8.loom")
  loom_list<-list(CB_D4.ldat,CB_D6.ldat,CB_D8.ldat)
  return(loom_list)
}
#load cellrangerv3 files
import_cellranger_v3<-function(){
  print("load cellrangerv3 files...")
  D4.data<-Read10X(data.dir="/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/run_count_CB_D4/outs/filtered_feature_bc_matrix/")
  D6.data<-Read10X(data.dir="/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/run_count_CB_D6/outs/filtered_feature_bc_matrix/")
  D8.data<-Read10X(data.dir="/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/run_count_CB1/outs/filtered_feature_bc_matrix/")
  cellranger_list<-list(D4.data,D6.data,D8.data)
  return(cellranger_list)
}
#create seurat object
create_seurat_object<-function(cellranger_list){
  print("create seurat object from cellranger...")
  seurat_list<-list()
  for(i in 1:length(cellranger_list)){
    seurat_list[[i]]<-CreateSeuratObject(counts = cellranger_list[[i]], project = "cordblood_v3", min.cells = 3, min.features = 200)
  }
  return(seurat_list)
}
#convesion of .loom to assay object
loom_to_seurat<-function(loomobject,seuratobject,id){
  print("begin .loom conversion process...")
  for (i in names(x = loomobject)) {
    print("Store assay in a new variable...")
    assay <- loomobject[[i]]
    
    print("Rename cell names in loom file to match cell names in Seurat object...")
    colnames(assay) <- gsub(paste0('possorted_genome_bam_',id,':'), '', colnames(assay))
    colnames(assay) <- gsub('x', '', colnames(assay))
    print("Subset to filtered cells in Seurat object...")
    assay <- assay[,colnames(seuratobject)]
    
    ### Add assay to Seurat object
    seuratobject[[i]] <- CreateAssayObject(counts = assay)
    
  }
  return(seuratobject)
}
#conversion process
loom_conversion<-function(loom_list,seurat_list){
  for(i in 1:length(seurat_list)){
    loom_id<-c("5MTBG","MZJBB","JKMX4")
    seurat_list[[i]]<-loom_to_seurat(loomobject = loom_list[[i]],seuratobject = seurat_list[[i]],id=loom_id[i])
  }
  return(seurat_list)
}
#count and subset percent mt
count_and_subset_mt<-function(seurat_list){
  print("count mitogene percentage...")
  for(i in 1:length(seurat_list)){
    seurat_list[[i]][["percent.mt"]]<-PercentageFeatureSet(seurat_list[[i]], pattern = "^MT-")
    print("subsetting...")
    seurat_list[[i]]<-subset(seurat_list[[i]], percent.mt<=10)
  }
  return(seurat_list)
}
#time assignment
time_assignment<-function(seurat_list){
  print("assign time...")
  time<-c("D4","D6","D8")
  for(i in 1:length(seurat_list)){
    seurat_list[[i]]$time<-time[i]
  }
  return(seurat_list)
}
#apply SCtransform to the list
sctransform_processing<-function(seurat_list,return.only.var.genes = FALSE){
  print("run SCtransform")
  for (i in 1:length(seurat_list)){
    seurat_list[[i]] <- SCTransform(seurat_list[[i]], verbose = FALSE, return.only.var.genes = return.only.var.genes)
  }
  return(seurat_list)
}
#dimensional redution
dimensionalreduction<-function(seurat_list){
  if(seurat_list%>%is.list){
    for (i in 1:length(seurat_list)) {
      seurat_list[[i]] <- RunPCA(seurat_list[[i]])
      seurat_list[[i]] <- FindNeighbors(object = seurat_list[[i]], dims = 1:30)
      seurat_list[[i]] <- FindClusters(object = seurat_list[[i]], resolution = 0.8)
      seurat_list[[i]] <- RunTSNE(seurat_list[[i]], dims = 1:30)
      seurat_list[[i]] <- RunUMAP(seurat_list[[i]], dims = 1:30)
    }
  }
  else{
    seurat_list <- RunPCA(seurat_list)
    seurat_list <- FindNeighbors(object = seurat_list, dims = 1:30)
    seurat_list <- FindClusters(object = seurat_list, resolution = 0.8)
    seurat_list <- RunTSNE(seurat_list, dims = 1:30)
    seurat_list <- RunUMAP(seurat_list, dims = 1:30)
  }
  return(seurat_list)
}
#doublet finding protocal
find_doublet<-function(seurat_list){
  for (i in 1:length(seurat_list)) {
    print(paste0("start ",i))
    print("pK Identification (no ground-truth)")
    sweep.res.list_temp <- paramSweep_v3(seurat_list[[i]], PCs = 1:30, sct = TRUE,num.cores = 12)
    sweep.stats_temp <- summarizeSweep(sweep.res.list_temp, GT = FALSE)
    pK_object <- find.pK(sweep.stats_temp)
    pK_value<-pK_object%>%filter(MeanBC==max(MeanBC))%>%dplyr::select(pK)%>%unlist%>%as.character()%>%as.numeric()
    print("Homotypic Doublet Proportion Estimate")
    annotations <- seurat_list[[i]]@meta.data$seurat_clusters
    annotations <- as.character(annotations)
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(0.039*length(seurat_list[[i]]%>%colnames))  ## Assuming 7.5% doublet formation rate - tailor for your dataset!!##!!
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    ##pK<0.1
    print("Run DoubletFinder with varying classification stringencies")
    seurat_list[[i]]  <- doubletFinder_v3(seurat_list[[i]] , PCs = 1:30, pN = 0.25, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    seurat_list[[i]]  <- doubletFinder_v3(seurat_list[[i]] , PCs = 1:30, pN = 0.25, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_",pK_value,"_",nExp_poi), sct = TRUE)
    print(paste0("finish ",i))
  }
  return(seurat_list)
}
#add metadata to clean seurat
add_metadata_to_clean<-function(serat_doublet,seurat_clean){
  for (i in 1:length(serat_doublet)) {
    seurat_clean[[i]]$DF.classifications<-serat_doublet[[i]]@meta.data[,ncol(serat_doublet[[i]]@meta.data)]
  }
  return(seurat_clean)
}
#subset doublet out from the data
eliminate_doublet<-function(seurat_list){
  for (i in 1:length(seurat_list)) {
    seurat_list[[1]]<-subset(seurat_list[[1]], subset = DF.classifications == 'Singlet')
  }
  return(seurat_list)
}
#integration preparation
integration_preparation<-function(seurat_preintegration,nfeatures=3000){
  seurat_features <- SelectIntegrationFeatures(object.list = seurat_preintegration, nfeatures = nfeatures)
  seurat_prepSCT <- PrepSCTIntegration(object.list = seurat_preintegration, anchor.features = seurat_features, verbose = FALSE)
  integration_factor<-list(seurat_features,seurat_prepSCT)
  return(integration_factor)
}
#integrate
integrate_seurat<-function(integration_factor){
  seurat_anchors <- FindIntegrationAnchors(object.list = integration_factor[[2]], normalization.method = "SCT", anchor.features = integration_factor[[1]], verbose = FALSE)
  seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT", verbose = FALSE)
  return(seurat_integrated)
}
#run pipeline
run_title="TRUEn5000_"
set.seed(365)
loom_list<-import_loom()
cellranger_list<-import_cellranger_v3()
seurat_list<-create_seurat_object(cellranger_list)
seurat_list<-loom_conversion(loom_list,seurat_list)
seurat_list<-count_and_subset_mt(seurat_list)
seurat_list<-time_assignment(seurat_list)
print("save the clean seurat set")
saveRDS(seurat_list, file=paste0("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/",run_title,"cordblood_preprocess_CLEAN.rds"))
print("clear workspace")
rm(cellranger_list,loom_list)
#SCtrasform
cordblood_preprocess_CLEAN=seurat_list
rm(seurat_list)
cordblood_preprocess_SCFALSE<-sctransform_processing(seurat_list=cordblood_preprocess_CLEAN,return.only.var.genes = TRUE)
rm(cordblood_preprocess_CLEAN)
cordblood_preprocess_SCFALSE_DIMREDUCT<-dimensionalreduction(cordblood_preprocess_SCFALSE)
rm(cordblood_preprocess_SCFALSE)
print("save transformed seurat set with dimensional reduction")
saveRDS(cordblood_preprocess_SCFALSE_DIMREDUCT, file=paste0("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/",run_title,"cordblood_preprocess_dimreduct.rds"))
cordblood_SCFALSE_doubletmarked<-find_doublet(cordblood_preprocess_SCFALSE_DIMREDUCT)
rm(cordblood_preprocess_SCFALSE_DIMREDUCT)
saveRDS(cordblood_SCFALSE_doubletmarked, file=paste0("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/",run_title,"cordblood_doubletmarked.rds"))
cordblood_preprocess_CLEAN<-readRDS(file=paste0("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/",run_title,"cordblood_preprocess_CLEAN.rds"))
seurat_SCFALSE_CLEAN_DOUBLET<-add_metadata_to_clean(cordblood_SCFALSE_doubletmarked,cordblood_preprocess_CLEAN)
seurat_SCFALSE_CLEAN_SINGLET<-eliminate_doublet(seurat_SCFALSE_CLEAN_DOUBLET)
rm(seurat_SCFALSE_CLEAN_DOUBLET, cordblood_SCFALSE_doubletmarked, cordblood_preprocess_CLEAN)
cordblood_prepintegrate_SCFALSE<-sctransform_processing(seurat_list=seurat_SCFALSE_CLEAN_SINGLET,return.only.var.genes = TRUE)
rm(seurat_SCFALSE_CLEAN_SINGLET)
cordblood_prepintegrate_SCFALSE_DIMREDUCT<-dimensionalreduction(cordblood_prepintegrate_SCFALSE)
print("save transformed seurat set with dimensional reduction")
saveRDS(cordblood_prepintegrate_SCFALSE_DIMREDUCT, file=paste0("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/",run_title,"cordblood_prepintegrate_DIMREDUCT.rds"))
print("clear workspace")
rm(cordblood_prepintegrate_SCFALSE)
cordblood_prepintegrate_SCFALSE_DIMREDUCT_PREPARED<-integration_preparation(cordblood_prepintegrate_SCFALSE_DIMREDUCT,nfeatures=5000)
cordblood_INTEGRATED_SCFALSE<-integrate_seurat(cordblood_prepintegrate_SCFALSE_DIMREDUCT_PREPARED)
cordblood_INTEGRATED_SCFALSE_DIMREDUCT<-dimensionalreduction(cordblood_INTEGRATED_SCFALSE)
saveRDS(cordblood_INTEGRATED_SCFALSE_DIMREDUCT, file=paste0("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/",run_title,"cordblood_INTEGRATED.rds"))
rm(cordblood_INTEGRATED_SCFALSE,cordblood_prepintegrate_SCFALSE_DIMREDUCT_PREPARED,cordblood_prepintegrate_SCFALSE_DIMREDUCT)
#print preview images
name = paste0("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/",run_title)
seurat_integrated=cordblood_INTEGRATED_SCFALSE_DIMREDUCT
png(filename = paste0(name,"_clustermap_umap.png"), width=500, height=500)
p <- DimPlot(object = seurat_integrated, reduction = "umap", label = TRUE)
CombinePlots(plots = list(p), ncol = 2)
dev.off()
png(filename = paste0(name,"_overvieweapression_umap.png"), width=1200, height=900)
FeaturePlot(object = seurat_integrated, features = c("CD34", "HBB","HBG1","HBG2","HBA1","SLC4A1","GYPA","KLF1", "BCL11A", "GATA1", "ZFPM1","MYB", "NFE4","ZBTB7A", "NR2C1", "NR2C2", "IKZF1", "MPO"), min.cutoff = "q9", cols = c("lightgrey", "blue"),reduction = "umap", pt.size = 0.5)
dev.off()
png(filename = paste0(name,"_clustermap_tsne.png"), width=500, height=500)
p <- DimPlot(object = seurat_integrated, reduction = "tsne", label = TRUE)
CombinePlots(plots = list(p), ncol = 2)
dev.off()
png(filename = paste0(name,"_overvieweapression_tsne.png"), width=1200, height=900)
FeaturePlot(object = seurat_integrated,  features = c("CD34", "HBB","HBG1","HBG2","HBA1","SLC4A1","GYPA","KLF1", "BCL11A", "GATA1", "ZFPM1","MYB", "NFE4","ZBTB7A", "NR2C1", "NR2C2", "IKZF1", "MPO"), min.cutoff = "q9", cols = c("lightgrey", "blue"),reduction = "tsne", pt.size = 0.5)
dev.off()
png(filename = paste0(name,"_clustermap_pca.png"), width=500, height=500)
p <- DimPlot(object = seurat_integrated, reduction = "pca", label = TRUE)
CombinePlots(plots = list(p), ncol = 2)
dev.off()
png(filename = paste0(name,"_overvieweapression_pca.png"), width=1200, height=900)
FeaturePlot(object = seurat_integrated,  features = c("CD34", "HBB","HBG1","HBG2","HBA1","SLC4A1","GYPA","KLF1", "BCL11A", "GATA1", "ZFPM1","MYB", "NFE4","ZBTB7A", "NR2C1", "NR2C2", "IKZF1", "MPO"), min.cutoff = "q9", cols = c("lightgrey", "blue"),reduction = "pca", pt.size = 0.5)
dev.off()


