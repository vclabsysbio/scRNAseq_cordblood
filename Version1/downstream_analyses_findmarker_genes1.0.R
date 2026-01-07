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
library(magrittr)

#define funtion to convert all markers table to tf only table
select_tf<-function(marker_table){
  marker_table%>%
    rownames_to_column()%>%
    left_join(gene_symbol_table, by=c("rowname"="external_gene_name"))%>%
    left_join(human_tf_db,by=("ensembl_gene_id"))%>%
    filter(Is_TF=="Yes")%>%
    dplyr::select(rowname,p_val,avg_logFC,pct.1,pct.2,p_val_adj)
}

select_tf_onegroup<-function(marker_table){
  marker_table%>%
    rownames_to_column()%>%
    left_join(gene_symbol_table, by=c("gene"="external_gene_name"))%>%
    left_join(human_tf_db,by=("ensembl_gene_id"))%>%
    filter(Is_TF=="Yes")%>%
    dplyr::select(gene,p.value,average_log_foldchange,percentage_of_cells_in_interested_group,percentage_of_other_cells,adjusted_p.value)
}
library(future)
plan("multiprocess", workers = 30)
options(future.globals.maxSize = 10 * 1024 ^ 3)
#set parameter
task_name="hemoglobin_containing_cells"
directory="/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/cod/analyses_hemoglobin_containing_dataset/"
print(directory)

#import file
integrated_seurat<-readRDS(file="/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/hemoglobin_containing_dataset.rds")

#Create TF list
human_tf_db<-read_tsv("/home/pakkanan/DatabaseExtract_v_1.01.txt", col_types = list(col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character(),col_character()))
colnames(human_tf_db)<-c("ensembl_gene_id","HGNC_symbol","DBD","Is_TF","TF_assessment","Binding_mode","Motif_status","Final_Notes","Final_Comments","Interpro_IDs","EntrezGene_ID","EntrezGene_Description","PDB_ID","TF_tested_by_HT-SELEX","TF _tested_by_PBM","Conditional_Binding_Requirements","Original_Comments","Vaquerizas_2009_classification","CisBP considers_it_a_TF","TFCat_classification","Is_a_GO_TF","Initial_assessment","Curator_1","Curator_2","TFclass_considers it_a_TF","Go_Evidence","Pfam_Domains_By_ENSP_ID")
#biomart_setup<-useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
#gene_symbol_table<-getBM(c("ensembl_gene_id","external_gene_name", "description"), filters="biotype",values="protein_coding", mart=biomart_setup, checkFilters = FALSE)
gene_symbol_table<-read_tsv("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/gene_symbol_table.tsv")


#select dataslot
plan("multiprocess", workers = 30)
DefaultAssay(integrated_seurat)<-"RNA"
#find markers
all_markers<-FindAllMarkers(integrated_seurat, min.pct = 0.1,assay = "RNA",slot = "data",logfc.threshold = 0.05)
#filter tf
tf_markers<-all_markers%>%rownames_to_column()%>%left_join(gene_symbol_table, by=c("rowname"="external_gene_name"))%>%left_join(human_tf_db,by=("ensembl_gene_id"))%>%filter(Is_TF=="Yes")

#write table#write table
write_tsv(all_markers,path = paste0(directory,task_name,"_all_markers_RNA.tsv"))
write_tsv(tf_markers,path = paste0(directory,task_name,"_tf_markers_RNA.tsv"))

#filter by cluster
late_HBG_DEG<-all_markers%>%filter(cluster=="late_HBG")
early_HBG_DEG<-all_markers%>%filter(cluster=="early_HBG")
late_HBB_DEG<-all_markers%>%filter(cluster=="late_HBB")
early_HBB_DEG<-all_markers%>%filter(cluster=="early_HBB")

#find DE genes by comparing clusters of interest 
marker_late_vs_early_HBG<-FindMarkers(integrated_seurat,ident.1 = "late_HBG" ,ident.2 = "early_HBG", min.pct = 0.1,assay = "RNA",slot = "data",logfc.threshold = 0.05)
marker_late_vs_early_HBB<-FindMarkers(integrated_seurat,ident.1 = "late_HBB" ,ident.2 = "early_HBB", min.pct = 0.1,assay = "RNA",slot = "data",logfc.threshold = 0.05)
marker_late_HBB_vs_HBG<-FindMarkers(integrated_seurat,ident.1 = "late_HBB" ,ident.2 = "late_HBG", min.pct = 0.1,assay = "RNA",slot = "data",logfc.threshold = 0.05)
marker_early_HBB_vs_HBG<-FindMarkers(integrated_seurat,ident.1 = "early_HBB" ,ident.2 = "early_HBG", min.pct = 0.1,assay = "RNA",slot = "data",logfc.threshold = 0.05)

#filter for TF
marker_late_vs_early_HBGtf<-marker_late_vs_early_HBG%>%select_tf()
marker_late_vs_early_HBBtf<-marker_late_vs_early_HBB%>%select_tf()
marker_late_HBB_vs_HBGtf<-marker_late_HBB_vs_HBG%>%select_tf()
marker_early_HBB_vs_HBGtf<-marker_early_HBB_vs_HBG%>%select_tf()
marker_late_HBG_DEGtf<-late_HBG_DEG%>%select_tf_onegroup()
marker_late_HBB_DEGtf<-late_HBB_DEG%>%select_tf_onegroup()
marker_early_HBG_DEGtf<-early_HBG_DEG%>%select_tf_onegroup()
marker_early_HBB_DEGtf<-early_HBB_DEG%>%select_tf_onegroup()

#write tables
write_tsv(marker_late_vs_early_HBG,path = paste0(directory,task_name,"_marker_late_vs_early_HBG.tsv"))
write_tsv(marker_late_vs_early_HBB,path = paste0(directory,task_name,"_marker_late_vs_early_HBB.tsv"))
write_tsv(marker_late_HBB_vs_HBG,path = paste0(directory,task_name,"_marker_late_HBB_vs_HBG.tsv"))
write_tsv(marker_early_HBB_vs_HBG,path = paste0(directory,task_name,"_marker_early_HBB_HBG.tsv"))
write_tsv(late_HBG_DEG,path = paste0(directory,task_name,"_marker_late_HBG_DEG.tsv"))
write_tsv(late_HBB_DEG,path = paste0(directory,task_name,"_marker_late_HBB_DEG.tsv"))
write_tsv(early_HBG_DEG,path = paste0(directory,task_name,"_marker_early_HBG_DEG.tsv"))
write_tsv(early_HBB_DEG,path = paste0(directory,task_name,"_marker_early_HBB_DEG.tsv"))
write_tsv(marker_late_vs_early_HBGtf,path = paste0(directory,task_name,"_marker_late_vs_early_HBGtf.tsv"))
write_tsv(marker_late_vs_early_HBBtf,path = paste0(directory,task_name,"_marker_late_vs_early_HBBtf.tsv"))
write_tsv(marker_late_HBB_vs_HBGtf,path = paste0(directory,task_name,"_marker_late_HBB_vs_HBGtf.tsv"))
write_tsv(marker_early_HBB_vs_HBGtf,path = paste0(directory,task_name,"_marker_early_HBB_HBGtf.tsv"))
write_tsv(marker_late_HBG_DEGtf,path = paste0(directory,task_name,"_marker_late_HBG_DEGtf.tsv"))
write_tsv(marker_late_HBB_DEGtf,path = paste0(directory,task_name,"_marker_late_HBB_DEGtf.tsv"))
write_tsv(marker_early_HBG_DEGtf,path = paste0(directory,task_name,"_marker_early_HBG_DEGtf.tsv"))
write_tsv(marker_early_HBB_DEGtf,path = paste0(directory,task_name,"_marker_early_HBB_DEGtf.tsv"))

#figure2C
#create average expression table
integrated_average_expression_table<-AverageExpression(integrated_seurat, assays = "RNA")
integrated_average_expression_table_RNA<-integrated_average_expression_table[["RNA"]]%>%rownames_to_column()

#load marker analyses result table
early_HBB_vs_HBG_table<-read_tsv(file = "./analyses_hemoglobin_containing_dataset/hemoglobin_containing_cells_marker_early_HBB_vs_HBG.tsv")
early_HBB_vs_HBGtf_table<-read_tsv(file = "./analyses_hemoglobin_containing_dataset/hemoglobin_containing_cells_marker_early_HBB_vs_HBGtf.tsv")
late_HBB_vs_HBG_table<-read_tsv(file = "./analyses_hemoglobin_containing_dataset/hemoglobin_containing_cells_marker_late_HBB_vs_HBG.tsv")
late_HBB_vs_HBGtf_table<-read_tsv(file = "./analyses_hemoglobin_containing_dataset/hemoglobin_containing_cells_marker_late_HBB_vs_HBGtf.tsv")

#load reference gene table from 10x data
gene_ref1<-read_tsv("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/run_count_CB_D4/outs/filtered_feature_bc_matrix/features.tsv.gz",col_names = F)
gene_ref2<-read_tsv("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/run_count_CB_D6//outs/filtered_feature_bc_matrix/features.tsv.gz",col_names = F)
gene_ref3<-read_tsv("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/run_count_CB1///outs/filtered_feature_bc_matrix/features.tsv.gz",col_names = F)
full_gene_ref<-gene_ref2%>%full_join(gene_ref1)%>%full_join(gene_ref3)%>%arrange(X1)
full_gene_ref->gene_ref

#extract top gene names
#DE genes
early_HBB_vs_HBG_top10_query<-early_HBB_vs_HBG_table%>%left_join(gene_ref,by=c("gene_id"="X2"))%>%filter(p_val_adj<0.05)%>%arrange(avg_logFC)%>%top_n(10, wt=avg_logFC)%>%extract2("gene_id")
early_HBB_vs_HBG_bottom10_query<-early_HBB_vs_HBG_table%>%left_join(gene_ref,by=c("gene_id"="X2"))%>%filter(p_val_adj<0.05)%>%arrange(avg_logFC)%>%top_n(-10, wt=avg_logFC)%>%extract2("gene_id")
late_HBB_vs_HBG_top13_query<-late_HBB_vs_HBG_table%>%left_join(gene_ref,by=c("gene_id"="X2"))%>%filter(p_val_adj<0.05)%>%arrange(avg_logFC)%>%top_n(13, wt=avg_logFC)%>%extract2("gene_id")
late_HBB_vs_HBG_lower_than_HBA1_query<-late_HBB_vs_HBG_table%>%left_join(gene_ref,by=c("gene_id"="X2"))%>%filter(p_val_adj<0.05)%>%arrange(avg_logFC)%>%filter(avg_logFC<late_HBB_vs_HBG_table%>%filter(gene_id=="HBA1")%>%extract2(3))%>%extract2("gene_id")
#TFs
early_HBB_vs_HBGtf_top10_query<-early_HBB_vs_HBGtf_table%>%left_join(gene_ref,by=c("rowname"="X2"))%>%filter(p_val_adj<0.05)%>%arrange(avg_logFC)%>%top_n(10, wt=avg_logFC)%>%extract2("rowname")
early_HBB_vs_HBGtf_bottom10_query<-early_HBB_vs_HBGtf_table%>%left_join(gene_ref,by=c("rowname"="X2"))%>%filter(p_val_adj<0.05)%>%arrange(avg_logFC)%>%top_n(-10, wt=avg_logFC)%>%extract2("rowname")
late_HBB_vs_HBGtf_top10_query<-late_HBB_vs_HBGtf_table%>%left_join(gene_ref,by=c("rowname"="X2"))%>%filter(p_val_adj<0.05)%>%arrange(avg_logFC)%>%top_n(10, wt=avg_logFC)%>%extract2("rowname")
late_HBB_vs_HBGtf_bottom10_query<-late_HBB_vs_HBGtf_table%>%left_join(gene_ref,by=c("rowname"="X2"))%>%filter(p_val_adj<0.05)%>%arrange(avg_logFC)%>%top_n(-10, wt=avg_logFC)%>%extract2("rowname")

#extract expression table according to gene names
#DE genes
early_HBBhigherthanHBG_table<-integrated_average_expression_table_RNA%>%filter(rowname %in% early_HBB_vs_HBG_top10_query)%>%cbind.data.frame(group="early_HBBhigherthanHBG_allgene")
early_HBGhigherthanHBB_table<-integrated_average_expression_table_RNA%>%filter(rowname %in% early_HBB_vs_HBG_bottom10_query)%>%cbind.data.frame(group="early_HBGhigherthanHBB_allgene")
late_HBBhigherthanHBG_table<-integrated_average_expression_table_RNA%>%filter(rowname %in% late_HBB_vs_HBG_top13_query)%>%cbind.data.frame(group="late_HBBhigherthanHBG_allgene")
late_HBGhigherthanHBB_table<-integrated_average_expression_table_RNA%>%filter(rowname %in% late_HBB_vs_HBG_lower_than_HBA1_query)%>%cbind.data.frame(group="late_HBGhigherthanHBB_allgene")
#TFs
early_HBBhigherthanHBG_table_tf<-integrated_average_expression_table_RNA%>%filter(rowname %in% early_HBB_vs_HBGtf_top10_query)%>%cbind.data.frame(group="early_HBBhigherthanHBG_TF")
early_HBGhigherthanHBB_table_tf<-integrated_average_expression_table_RNA%>%filter(rowname %in% early_HBB_vs_HBGtf_bottom10_query)%>%cbind.data.frame(group="early_HBGhigherthanHBB_TF")
late_HBBhigherthanHBG_table_tf<-integrated_average_expression_table_RNA%>%filter(rowname %in% late_HBB_vs_HBGtf_top10_query)%>%cbind.data.frame(group="late_HBBhigherthanHBG_TF")
late_HBGhigherthanHBB_table_tf<-integrated_average_expression_table_RNA%>%filter(rowname %in% late_HBB_vs_HBGtf_bottom10_query)%>%cbind.data.frame(group="late_HBGhigherthanHBB_TF")
#create matrix for heatmap
#DE genes
early_HBBhigherthanHBG_matrix<-early_HBBhigherthanHBG_table%>%dplyr::select(-group)%>%column_to_rownames()%>%as.matrix()%>%t%>%scale()%>%t%>%as.matrix()
early_HBBhigherthanHBG_hclust<-early_HBBhigherthanHBG_matrix%>%hclust.vector(method="ward",metric="euclidean")
early_HBBhigherthanHBG_table_sorted<-early_HBBhigherthanHBG_table[early_HBBhigherthanHBG_hclust$order,]%>%cbind.data.frame(idx=seq(1:nrow(.)))%>%remove_rownames()
late_HBBhigherthanHBG_matrix<-late_HBBhigherthanHBG_table%>%dplyr::select(-group)%>%column_to_rownames()%>%as.matrix()%>%t%>%scale()%>%t%>%as.matrix()

late_HBBhigherthanHBG_hclust<-late_HBBhigherthanHBG_matrix%>%hclust.vector(method="ward",metric="euclidean")
late_HBBhigherthanHBG_table_sorted<-late_HBBhigherthanHBG_table[late_HBBhigherthanHBG_hclust$order,]%>%cbind.data.frame(idx=seq(1:nrow(.)))%>%remove_rownames()
#---------------#
early_HBGhigherthanHBB_matrix<-early_HBGhigherthanHBB_table%>%dplyr::select(-group)%>%column_to_rownames()%>%as.matrix()%>%t%>%scale()%>%t%>%as.matrix()
early_HBGhigherthanHBB_hclust<-early_HBGhigherthanHBB_matrix%>%hclust.vector(method="ward",metric="euclidean")
early_HBGhigherthanHBB_table_sorted<-early_HBGhigherthanHBB_table[early_HBGhigherthanHBB_hclust$order,]%>%cbind.data.frame(idx=seq(1:nrow(.)))%>%remove_rownames()

late_HBGhigherthanHBB_matrix<-late_HBGhigherthanHBB_table%>%dplyr::select(-group)%>%column_to_rownames()%>%as.matrix()%>%t%>%scale()%>%t%>%as.matrix()
late_HBGhigherthanHBB_hclust<-late_HBGhigherthanHBB_matrix%>%hclust.vector(method="ward",metric="euclidean")
late_HBGhigherthanHBB_table_sorted<-late_HBGhigherthanHBB_table[late_HBGhigherthanHBB_hclust$order,]%>%cbind.data.frame(idx=seq(1:nrow(.)))%>%remove_rownames()
#TFs
early_HBBhigherthanHBG_matrix_tf<-early_HBBhigherthanHBG_table_tf%>%dplyr::select(-group)%>%column_to_rownames()%>%as.matrix()%>%t%>%scale()%>%t%>%as.matrix()
early_HBBhigherthanHBG_hclust_tf<-early_HBBhigherthanHBG_matrix_tf%>%hclust.vector(method="ward",metric="euclidean")
early_HBBhigherthanHBG_table_sorted_tf<-early_HBBhigherthanHBG_table_tf[early_HBBhigherthanHBG_hclust_tf$order,]%>%cbind.data.frame(idx=seq(1:nrow(.)))%>%remove_rownames()

late_HBBhigherthanHBG_matrix_tf<-late_HBBhigherthanHBG_table_tf%>%dplyr::select(-group)%>%column_to_rownames()%>%as.matrix()%>%t%>%scale()%>%t%>%as.matrix()
late_HBBhigherthanHBG_hclust_tf<-late_HBBhigherthanHBG_matrix_tf%>%hclust.vector(method="ward",metric="euclidean")
late_HBBhigherthanHBG_table_sorted_tf<-late_HBBhigherthanHBG_table_tf[late_HBBhigherthanHBG_hclust_tf$order,]%>%cbind.data.frame(idx=seq(1:nrow(.)))%>%remove_rownames()
#--------------#
early_HBGhigherthanHBB_matrix_tf<-early_HBGhigherthanHBB_table_tf%>%dplyr::select(-group)%>%column_to_rownames()%>%as.matrix()%>%t%>%scale()%>%t%>%as.matrix()
early_HBGhigherthanHBB_hclust_tf<-early_HBGhigherthanHBB_matrix_tf%>%hclust.vector(method="ward",metric="euclidean")
early_HBGhigherthanHBB_table_sorted_tf<-early_HBGhigherthanHBB_table_tf[early_HBGhigherthanHBB_hclust_tf$order,]%>%cbind.data.frame(idx=seq(1:nrow(.)))%>%remove_rownames()

late_HBGhigherthanHBB_matrix_tf<-late_HBGhigherthanHBB_table_tf%>%dplyr::select(-group)%>%column_to_rownames()%>%as.matrix()%>%t%>%scale()%>%t%>%as.matrix()
late_HBGhigherthanHBB_hclust_tf<-late_HBGhigherthanHBB_matrix_tf%>%hclust.vector(method="ward",metric="euclidean")
late_HBGhigherthanHBB_table_sorted_tf<-late_HBGhigherthanHBB_table_tf[late_HBGhigherthanHBB_hclust_tf$order,]%>%cbind.data.frame(idx=seq(1:nrow(.)))%>%remove_rownames()

#define function for heatmap ##lazy heatmap
make_lazyheatmap2<-function(matrix_or_interest,heatmap_name,name){Heatmap(matrix_or_interest, 
                                                                          #heatmap internal name
                                                                          name=name,
                                                                          #column ordering and dendogram config
                                                                          show_column_dend = FALSE,
                                                                          cluster_columns = FALSE,
                                                                          #row ordering and dendrogram config
                                                                          show_row_dend = FALSE, 
                                                                          cluster_rows = FALSE,
                                                                          #title config
                                                                          column_title = heatmap_name,
                                                                          #row and column name config
                                                                          show_row_names = TRUE, 
                                                                          show_column_names = TRUE,
                                                                          column_names_side = "top",
                                                                          row_names_side = "left", 
                                                                          row_names_gp = gpar(fontsize = 12),
                                                                          #color scale config
                                                                          col = colorRamp2(c(-1.5,0,1.5), c("blue","black","yellow")),
                                                                          #borderconfig
                                                                          #rect_gp = gpar(col = "white", lwd = 2),
                                                                          #split=5
                                                                          #set legend
                                                                          
                                                                          #split=splitter,
                                                                          use_raster = FALSE,
                                                                          #annotation
                                                                          #top_annotation = column_annotation_object
                                                                          column_names_rot = 45,
                                                                          width = unit(40, "mm")
)
  
}

#join all tables into one and sort the data
#DEG
all_gene_table_sorted<-early_HBBhigherthanHBG_table_sorted%>%rbind.data.frame(early_HBGhigherthanHBB_table_sorted)%>%rbind.data.frame(late_HBBhigherthanHBG_table_sorted)%>%rbind.data.frame(late_HBGhigherthanHBB_table_sorted)%>%separate("group",c("time","direction","type"))
#TF
all_gene_table_sorted_tf<-early_HBBhigherthanHBG_table_sorted_tf%>%rbind.data.frame(early_HBGhigherthanHBB_table_sorted_tf)%>%rbind.data.frame(late_HBBhigherthanHBG_table_sorted_tf)%>%rbind.data.frame(late_HBGhigherthanHBB_table_sorted_tf)%>%separate("group",c("time","direction","type"))

#print heatmap
#figure 2B
temp_heatmap<-make_lazyheatmap2(matrix_or_interest = all_gene_table_sorted%>%filter(time=="early")%>%dplyr::select(-direction,-time,-type,-idx)%>%column_to_rownames()%>%as.matrix()%>%t%>%scale()%>%t%>%as.matrix(), "Relative average transcription levels of\n differentially expressed genes (DEGs)\n between the 'early' HBB- or HBG-high groups",name="z-score")
name = paste0(directory,task_name)
png(filename = paste0(name,"_heatmap_DEG_RNA.png"), width=360, height=420)
draw(temp_heatmap+Heatmap(all_gene_table_sorted%>%filter(time=="early")%>%dplyr::select(direction)%>%mutate_if(is.character,str_replace_all, pattern = "HBBhigherthanHBG", replacement = "HBB > HBG")%>%mutate_if(is.character,pattern = "HBGhigherthanHBB", replacement = "HBG > HBB"), name="direction",col=c("black","grey"), width = unit(5, "mm")),
     heatmap_legend_side="right",
     padding = unit(c(5, 3, 3, 3), "mm")
)
dev.off()
#figure 2C
temp_heatmap<-make_lazyheatmap2(matrix_or_interest = all_gene_table_sorted%>%filter(time=="late")%>%dplyr::select(-direction,-time,-type,-idx)%>%column_to_rownames()%>%as.matrix()%>%t%>%scale()%>%t%>%as.matrix(), "Relative average transcription levels of\n differentially expressed genes (DEGs)\n between the 'late' HBB- or HBG-high groups",name="z-score")
name = paste0(directory,task_name)
png(filename = paste0(name,"_heatmap_DEG_RNA_late.png"), width=360, height=490)
draw(temp_heatmap+Heatmap(all_gene_table_sorted%>%filter(time=="late")%>%dplyr::select(direction)%>%mutate_if(is.character,str_replace_all, pattern = "HBBhigherthanHBG", replacement = "HBB > HBG")%>%mutate_if(is.character,   str_replace_all, pattern = "HBGhigherthanHBB", replacement = "HBG > HBB"), name="direction",col=c("black","grey"), width = unit(5, "mm")),
     heatmap_legend_side="right",
     padding = unit(c(5, 3, 3, 3), "mm")
)
dev.off()
#figure 2D
temp_heatmap<-make_lazyheatmap2(matrix_or_interest = all_gene_table_sorted_tf%>%filter(time=="early")%>%dplyr::select(-direction,-time,-type,-idx)%>%column_to_rownames()%>%as.matrix()%>%t%>%scale()%>%t%>%as.matrix(), "Relative average transcription levels of\n transciption factors (TFs)\n between the 'early' HBB- or HBG-high groups\n (hard cutoff)",name="z-score")
name = paste0(directory,task_name)
png(filename = paste0(name,"_heatmap_TF_RNA.png"), width=360, height=430)
draw(temp_heatmap+Heatmap(all_gene_table_sorted_tf%>%filter(time=="early")%>%dplyr::select(direction)%>%mutate_if(is.character,str_replace_all, pattern = "HBBhigherthanHBG", replacement = "HBB > HBG")%>%mutate_if(is.character,pattern = "HBGhigherthanHBB", replacement = "HBG > HBB"), name="direction",col=c("black","grey"), width = unit(5, "mm")),
     heatmap_legend_side="right",
     padding = unit(c(5, 3, 3, 3), "mm")
)
dev.off()
#figure 2E
temp_heatmap<-make_lazyheatmap2(matrix_or_interest = all_gene_table_sorted_tf%>%filter(time=="late")%>%dplyr::select(-direction,-time,-type,-idx)%>%column_to_rownames()%>%as.matrix()%>%t%>%scale()%>%t%>%as.matrix(), "Relative average transcription levels of\n transciption factors (TFs)\n between the 'late' HBB- or HBG-high groups\n (hard cutoff)",name="z-score")
name = paste0(directory,task_name)
png(filename = paste0(name,"_heatmap_TF_RNA_late.png"), width=360, height=440)
draw(temp_heatmap+Heatmap(all_gene_table_sorted_tf%>%filter(time=="late")%>%dplyr::select(direction)%>%mutate_if(is.character,str_replace_all, pattern = "HBBhigherthanHBG", replacement = "HBB > HBG")%>%mutate_if(is.character,str_replace_all, pattern = "HBGhigherthanHBB", replacement = "HBG > HBB"), name="direction",col=c("black","grey"), width = unit(5, "mm")),
     heatmap_legend_side="right",
     padding = unit(c(5, 3, 3, 3), "mm")
)
dev.off()

#mark selected cluster
#figure 2A
#extract_color from seurat
dataset_umap<-UMAPPlot(integrated_seurat)
dataset_umap_build<-ggplot_build(dataset_umap)
Idents(integrated_seurat)<-integrated_seurat$time
dataset_umap_time<-UMAPPlot(integrated_seurat)
dataset_umap_data <- dataset_umap_build$data[[1]]
dataset_umap_data <-  dataset_umap_data[order(dataset_umap_data$group), ]
ucols <- unique(dataset_umap_data$colour)
names(ucols) <- unique(dataset_umap_data$group)
dataset_umap_time_build<-ggplot_build(dataset_umap_time)
dataset_umap_time_data<-dataset_umap_time_build$data[[1]]
dataset_umap_time_data<-dataset_umap_time_data%>%dplyr::select(x,y,group,colour)
dataset_umap_time_data<-dataset_umap_time_data%>%mutate(time=cut(group, breaks = c(0,1,2,3),labels=c("D4","D6","D8")))
dataset_umap_data<-dataset_umap_data%>%mutate(active_cluster=cut(group, breaks = c(0,1,2,3,4),labels=c("early_HBG","late_HBG","early_HBB", "late_HBB")))
dataset_umap_table<-dataset_umap_data%>%left_join(dataset_umap_time_data,by=c("x","y"))
#define labels
cols_pll = c("cyan","blue","hotpink", "red")
cluster_labels=c("early_HBG","late_HBG","early_HBB", "late_HBB")
#print figure
umap_overall<-ggplot(dataset_umap_table, aes(x=x,y=y, color=active_cluster))+
  geom_point(size=0.001)+
  labs(title= "Scatterplot of dataset",
       subtitle = "n=11200",
       x="UMAP_1",
       y="UMAP_2")+
  #facet_grid(cols=(vars(time)),rows=vars(active_cluster))+
  guides(colour = guide_legend(override.aes =list(size=3)))+
  scale_color_manual(values=setNames(cols_pll,cluster_labels))+
  theme_classic()

png("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/cod/revised_analyses//revision3_cluster_umap.png", width=700, height=500,res = 150)
grid.arrange(
  grobs = list(umap_overall)
)
dev.off()
#define labels
cols_pll2.1 = c("cyan","grey90","hotpink", "grey90")
cluster_labels2.1=c("early_HBG","late_HBG","early_HBB", "late_HBB")

umap_overall<-ggplot(dataset_umap_table, aes(x=x,y=y, color=active_cluster))+
  geom_point(size=0.001)+
  labs(title= "Scatterplot of Dataset",
       subtitle = "Early",
       x="UMAP_1",
       y="UMAP_2")+
  #facet_grid(cols=(vars(time)),rows=vars(active_cluster))+
  guides(colour = guide_legend(override.aes =list(size=3)))+
  scale_color_manual(values=setNames(cols_pll2.1,cluster_labels2.1))+
  theme_classic()

png("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/cod/revised_analyses//revision3_figure2.1a.png", width=700, height=500,res = 150)
grid.arrange(
  grobs = list(umap_overall)
)
dev.off()
#define labels
cols_pll2.2 = c("grey90","blue","grey90", "red")
cluster_labels2.2=c("early_HBG","late_HBG","early_HBB", "late_HBB")

umap_overall<-ggplot(dataset_umap_table, aes(x=x,y=y, color=active_cluster))+
  geom_point(size=0.001)+
  labs(title= "Scatterplot of Dataset",
       subtitle = "Late",
       x="UMAP_1",
       y="UMAP_2")+
  #facet_grid(cols=(vars(time)),rows=vars(active_cluster))+
  guides(colour = guide_legend(override.aes =list(size=3)))+
  scale_color_manual(values=setNames(cols_pll2.2,cluster_labels2.2))+
  theme_classic()

png("/home/shared_data_vclab_2/SCRNAseq_cordblood/SCRNAseq_cordblood_pakkanan/cod/revised_analyses//revision3_figure2.2a.png", width=700, height=500,res = 150)
grid.arrange(
  grobs = list(umap_overall)
)
dev.off()


