rm(list = ls())
source("0_ScFunctions.R")
SI_merged_prefix <- "SI_merged_CD4w_removed/"
SI_Enterocyte_monocle_prefix <- "SI_Enterocyte_monocle/"
SI_CD_Enterocyte_monocle_prefix <- "SI_Enterocyte_monocle/CD_group/"
SI_HfiD_Enterocyte_monocle_prefix <- "SI_Enterocyte_monocle/HfiD_group/"

SI_merged <- readRDS(file = paste0(SI_merged_prefix,"Harmony_SI_merged_final_final_annotated.rds"))
SI_merged_Enterocyte <- subset(SI_merged, My_annotation %in% c("Enterocyte Progenitor",
                                                               "Villus Tip Enterocyte",
                                                               "Immune Related Enterocyte",
                                                               "Goblet Cell"))
SI_merged_Enterocyte_CD <- subset(SI_merged_Enterocyte, orig.ident %in% c("CD_SI_HFD0w", "CD_SI_HFD4w",  "CD_SI_HFD8w"))
SI_merged_Enterocyte_HfiD <- subset(SI_merged_Enterocyte, orig.ident %in% c("HfiD_SI_HFD0w", "HfiD_SI_HFD4w", "HfiD_SI_HFD8w"))


###########################CD group Enterocyte##################################

SI_merged_Enterocyte_CD$My_annotation <- Idents(SI_merged_Enterocyte_CD)
Idents(SI_merged_Enterocyte_CD) <- SI_merged_Enterocyte_CD$orig.ident
SI_merged_Enterocyte_CD <- subset(x=SI_merged_Enterocyte_CD, downsample=1000)
Idents(SI_merged_Enterocyte_CD) <- SI_merged_Enterocyte_CD$My_annotation

data_CD <- as(as.matrix(SI_merged_Enterocyte_CD@assays$SCT@counts), 'sparseMatrix')
pd_CD <- new('AnnotatedDataFrame', data = SI_merged_Enterocyte_CD@meta.data)
fData_CD <- data.frame(gene_short_name = row.names(data_CD), row.names = row.names(data_CD))
fd_CD <- new('AnnotatedDataFrame', data = fData_CD)

myCDs_CD <- newCellDataSet(data_CD,
                        phenoData = pd_CD,
                        featureData = fd_CD,
                        expressionFamily = negbinomial.size())

myCDs_CD <- estimateSizeFactors(myCDs_CD)
myCDs_CD <- estimateDispersions(myCDs_CD)                        
myCDs_CD <- detectGenes(myCDs_CD, min_expr = 0.1)
myCDs_expressed_genes_CD <- row.names(myCDs_CD)
diff_test_res_CD <- differentialGeneTest(myCDs_CD[myCDs_expressed_genes_CD,],fullModelFormulaStr = "~My_annotation",cores = 4,verbose = T)

diff_test_res_CD[1:4,1:4]
write.csv(diff_test_res_CD, file=paste0(SI_CD_Enterocyte_monocle_prefix,"CD_group_diff_test_result_", Sys.Date(),".csv"),sep="")
#diff_test_res_CD <- read.csv(file=paste0(AdMonocle_prefix,"CD_group_diff_test_result.csv",sep=""))

ordering_genes_CD <- row.names(subset(diff_test_res_CD, qval < 0.01))
ordering_genes_CD[1:5]
myCDs_CD <- setOrderingFilter(myCDs_CD, ordering_genes_CD)
plot_ordering_genes(myCDs_CD)
ggsave(paste0(SI_CD_Enterocyte_monocle_prefix,"CD_group_Ordering.genes_", Sys.Date(), ".pdf"), height = 7, width = 8)

myCDs_CD <- reduceDimension(myCDs_CD, max_components = 2, method = 'DDRTree',verbose = T)
myCDs_CD <- orderCells(myCDs_CD)

saveRDS(myCDs_CD, file = paste0(SI_CD_Enterocyte_monocle_prefix, "CD_group_myCDs_",Sys.Date(),".rds"))
#myCDs_CD <- readRDS(paste0(SI_CD_Enterocyte_monocle_prefix, "CD_group_myCDs.rds"))

#plot heatmap of CD group
diff_test_res_CD <- diff_test_res_CD[order(diff_test_res_CD$qval, decreasing = F),]
sig_gene_names_CD <- row.names(diff_test_res_CD[1:700,])
#nrow(subset(diff_test_res_CD, qval < 0.0001))

bins <- monocle3_scale_to_100(myCDs_CD[sig_gene_names_CD,])

p1 <- plot_pseudotime_heatmap(myCDs_CD[sig_gene_names_CD,],
                              num_clusters = 3,
                              cores = 6,
                              show_rownames = T, return_heatmap = T,
                              add_annotation_col = bins,use_gene_short_name = T)
pdf(paste0(SI_CD_Enterocyte_monocle_prefix,"CD_group_monocel_heatmap_",Sys.Date(),".pdf"), width = 5, height = 15)
p1
dev.off()

df <- data.frame((myCDs_CD[sig_gene_names_CD,]@phenoData@data))
df <- df[,c("Pseudotime", "My_annotation","orig.ident")]

#"#A6CEE3", "#A2D48E",
#"#66A5CC","#47A93A",
#"#267CB6","#20854EFF"
color.celtype <- data.frame(celltype = c("Enterocyte Progenitor",
                                         "Villus Tip Enterocyte",
                                         "Immune Related Enterocyte",
                                         "Goblet Cell"), 
                            color = c(hughie_color[1],hughie_color[4],hughie_color[5],hughie_color[7]))  
mWAT_SAMcolors <- data.frame(orig.ident = c("CD_SI_HFD0w","CD_SI_HFD4w","CD_SI_HFD8w"),
                             color = c("#A6CEE3","#66A5CC","#267CB6"))

df <- merge(df,color.celtype,by.x = "My_annotation", by.y = "celltype",drop = F)
df <- merge(df,mWAT_SAMcolors,by.x = "orig.ident", by.y = "orig.ident",drop = F)
df <- df[order(df$Pseudotime, decreasing = F),]

write.csv(df, file=paste0(SI_CD_Enterocyte_monocle_prefix,"CD_group_colnames_",Sys.Date(),".csv",sep=""))

matrix <- as.matrix(df[,3])
pdf(paste0(SI_CD_Enterocyte_monocle_prefix,"CD_group_monocel_heatmap_celltype_",Sys.Date(),".pdf"))
image(matrix, col = df$color.x)
dev.off()

pdf(paste0(SI_CD_Enterocyte_monocle_prefix,"CD_group_monocel_heatmap_timepoint_",Sys.Date(),".pdf"))

image(matrix, col = df$color.y)
dev.off()


#######monocle trajectory
data <- GetAssayData(SI_merged_Enterocyte_CD, assay = 'SCT', slot = 'counts')
cell_metadata <- SI_merged_Enterocyte_CD@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
monocle3_CDs <- new_cell_data_set(data,
                                  cell_metadata = cell_metadata,
                                  gene_metadata = gene_annotation)
#preprocess_CDs函数相当于seurat中NormalizeData ScaleData RunPCA
monocle3_CDs <- preprocess_cds(monocle3_CDs, num_dim = 50)
#umap降维
monocle3_CDs <- reduce_dimension(monocle3_CDs, preprocess_method = "PCA")
p1 <- plot_cells(monocle3_CDs, reduction_method="UMAP", color_cells_by="My_annotation"
                 #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
)+ggtitle('CDs.umap')
##从seurat导入整合过的umap坐标
monocle3_CDs.embed <- monocle3_CDs@int_colData$reducedDims$UMAP
int.embed <- Embeddings(SI_merged_Enterocyte_CD, reduction = "umap")
int.embed <- int.embed[rownames(monocle3_CDs.embed),]
monocle3_CDs@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(monocle3_CDs, reduction_method="UMAP", color_cells_by="My_annotation"
                 #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
)+ggtitle('int.umap')
p1 + p2
ggsave(paste0(SI_CD_Enterocyte_monocle_prefix,"Monocle3_umap_embed_", Sys.Date(), ".pdf"), height = 7, width = 16)

## Monocle3聚类分区
monocle3_CDs <- cluster_cells(monocle3_CDs, resolution = 0.01)
p3 <- plot_cells(monocle3_CDs, show_trajectory_graph = FALSE)+ggtitle("Label by clusterID")
p4 <- plot_cells(monocle3_CDs, color_cells_by = "partition", show_trajectory_graph = FALSE)+ggtitle("Label by partitionID")
p3 + p4
ggsave(paste0(SI_CD_Enterocyte_monocle_prefix,"Monocle3_umap_cluster_", Sys.Date(), ".pdf"), height = 7, width = 16)


## 识别轨迹
monocle3_CDs <- learn_graph(monocle3_CDs)
p5 <- plot_cells(monocle3_CDs, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                 label_branch_points = FALSE)
p5
ggsave(paste0(SI_CD_Enterocyte_monocle_prefix,"Monocle3_umap_trajectory_", Sys.Date(), ".pdf"), width = 5, height = 5)


##细胞按拟时排序
# CDs <- order_cells(CDs) 存在bug，使用辅助线选择root细胞
#p   geom_vline(xintercept = seq(-7,-6,0.25))   geom_hline(yintercept = seq(0,1,0.25))
#embed <- data.frame(Embeddings(all, reduction = "umap"))
#embed <- subset(embed, umap_1 > -6.75 & umap_1 < -6.5 & umap_2 > 0.24 & umap_2 < 0.25)
#root.cell <- rownames(embed)
#monocle3_CDs <- order_cells(monocle3_CDs, root_cells = root.cell)
monocle3_CDs <- order_cells(monocle3_CDs)
p6 <- plot_cells(monocle3_CDs, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                 label_leaves = F,  label_branch_points = FALSE,graph_label_size = 3,
                 cell_size = 1, alpha = 0.5, trajectory_graph_color = "black",trajectory_graph_segment_size = 1)
p6
ggsave(paste0(SI_CD_Enterocyte_monocle_prefix,"Monocle3_umap_trajectory_line_", Sys.Date(), ".pdf"), width = 6.5, height = 5)








###########################HfiD group Enterocyte##################################

SI_merged_Enterocyte_HfiD$My_annotation <- Idents(SI_merged_Enterocyte_HfiD)
Idents(SI_merged_Enterocyte_HfiD) <- SI_merged_Enterocyte_HfiD$orig.ident
SI_merged_Enterocyte_HfiD <- subset(x=SI_merged_Enterocyte_HfiD, downsample=1000)
Idents(SI_merged_Enterocyte_HfiD) <- SI_merged_Enterocyte_HfiD$My_annotation

data_HfiD <- as(as.matrix(SI_merged_Enterocyte_HfiD@assays$SCT@counts), 'sparseMatrix')
pd_HfiD <- new('AnnotatedDataFrame', data = SI_merged_Enterocyte_HfiD@meta.data)
fData_HfiD <- data.frame(gene_short_name = row.names(data_HfiD), row.names = row.names(data_HfiD))
fd_HfiD <- new('AnnotatedDataFrame', data = fData_HfiD)

myCDs_HfiD <- newCellDataSet(data_HfiD,
                           phenoData = pd_HfiD,
                           featureData = fd_HfiD,
                           expressionFamily = negbinomial.size())

myCDs_HfiD <- estimateSizeFactors(myCDs_HfiD)
myCDs_HfiD <- estimateDispersions(myCDs_HfiD)                        
myCDs_HfiD <- detectGenes(myCDs_HfiD, min_expr = 0.1)
myCDs_expressed_genes_HfiD <- row.names(myCDs_HfiD)
diff_test_res_HfiD <- differentialGeneTest(myCDs_HfiD[myCDs_expressed_genes_HfiD,],fullModelFormulaStr = "~My_annotation",cores = 4,verbose = T)

diff_test_res_HfiD[1:4,1:4]
write.csv(diff_test_res_HfiD, file=paste0(SI_HfiD_Enterocyte_monocle_prefix,"HfiD_group_diff_test_result_", Sys.Date(),".csv"),sep="")
#diff_test_res_HfiD <- read.csv(file=paste0(AdMonocle_prefix,"HfiD_group_diff_test_result.csv",sep=""))

ordering_genes_HfiD <- row.names(subset(diff_test_res_HfiD, qval < 0.01))
ordering_genes_HfiD[1:5]
myCDs_HfiD <- setOrderingFilter(myCDs_HfiD, ordering_genes_HfiD)
plot_ordering_genes(myCDs_HfiD)
ggsave(paste0(SI_HfiD_Enterocyte_monocle_prefix,"HfiD_group_Ordering.genes_", Sys.Date(), ".pdf"), height = 7, width = 8)

myCDs_HfiD <- reduceDimension(myCDs_HfiD, max_components = 2, method = 'DDRTree',verbose = T)
myCDs_HfiD <- orderCells(myCDs_HfiD)

saveRDS(myCDs_HfiD, file = paste0(SI_HfiD_Enterocyte_monocle_prefix, "HfiD_group_myHfiDs_",Sys.Date(),".rds"))
#myCDs_HfiD <- readRDS(paste0(SI_HfiD_Enterocyte_monocle_prefix, "HfiD_group_myHfiDs.rds"))

#plot heatmap of HfiD group
diff_test_res_HfiD <- diff_test_res_HfiD[order(diff_test_res_HfiD$qval, decreasing = F),]
sig_gene_names_HfiD <- row.names(diff_test_res_HfiD[1:300,])
#nrow(subset(diff_test_res_HfiD, qval < 0.0001))

bins <- monocle3_scale_to_100(myCDs_HfiD[sig_gene_names_HfiD,])

p1 <- plot_pseudotime_heatmap(myCDs_HfiD[sig_gene_names_HfiD,],
                              num_clusters = 3,
                              cores = 6,
                              show_rownames = T, return_heatmap = T,
                              add_annotation_col = bins,use_gene_short_name = T)
pdf(paste0(SI_HfiD_Enterocyte_monocle_prefix,"HfiD_group_monocel_heatmap_",Sys.Date(),".pdf"), width = 5, height = 15)
p1
dev.off()

df <- data.frame((myCDs_HfiD[sig_gene_names_HfiD,]@phenoData@data))
df <- df[,c("Pseudotime", "My_annotation","orig.ident")]

#"#A6CEE3", "#A2D48E",
#"#66A5CC","#47A93A",
#"#267CB6","#20854EFF"
color.celtype <- data.frame(celltype = c("Enterocyte Progenitor",
                                         "Villus Tip Enterocyte",
                                         "Immune Related Enterocyte",
                                         "Goblet Cell"), 
                            color = c(hughie_color[1],hughie_color[4],hughie_color[5],hughie_color[7]))  
mWAT_SAMcolors <- data.frame(orig.ident = c("HfiD_SI_HFD0w","HfiD_SI_HFD4w","HfiD_SI_HFD8w"),
                             color = c("#A2D48E","#47A93A","#20854EFF"))

df <- merge(df,color.celtype,by.x = "My_annotation", by.y = "celltype",drop = F)
df <- merge(df,mWAT_SAMcolors,by.x = "orig.ident", by.y = "orig.ident",drop = F)
df <- df[order(df$Pseudotime, decreasing = F),]

write.csv(df, file=paste0(SI_HfiD_Enterocyte_monocle_prefix,"HfiD_group_colnames_",Sys.Date(),".csv",sep=""))

matrix <- as.matrix(df[,3])
pdf(paste0(SI_HfiD_Enterocyte_monocle_prefix,"HfiD_group_monocel_heatmap_celltype_",Sys.Date(),".pdf"))
image(matrix, col = df$color.x)
dev.off()

pdf(paste0(SI_HfiD_Enterocyte_monocle_prefix,"HfiD_group_monocel_heatmap_timepoint_",Sys.Date(),".pdf"))

image(matrix, col = df$color.y)
dev.off()


#######monocle trajectory
data <- GetAssayData(SI_merged_Enterocyte_HfiD, assay = 'SCT', slot = 'counts')
cell_metadata <- SI_merged_Enterocyte_HfiD@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
monocle3_HfiDs <- new_cell_data_set(data,
                                  cell_metadata = cell_metadata,
                                  gene_metadata = gene_annotation)
#preprocess_CDs函数相当于seurat中NormalizeData ScaleData RunPCA
monocle3_HfiDs <- preprocess_cds(monocle3_HfiDs, num_dim = 50)
#umap降维
monocle3_HfiDs <- reduce_dimension(monocle3_HfiDs, preprocess_method = "PCA")
p1 <- plot_cells(monocle3_HfiDs, reduction_method="UMAP", color_cells_by="My_annotation"
                 #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
)+ggtitle('HfiDs.umap')
##从seurat导入整合过的umap坐标
monocle3_HfiDs.embed <- monocle3_HfiDs@int_colData$reducedDims$UMAP
int.embed <- Embeddings(SI_merged_Enterocyte_HfiD, reduction = "umap")
int.embed <- int.embed[rownames(monocle3_HfiDs.embed),]
monocle3_HfiDs@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(monocle3_HfiDs, reduction_method="UMAP", color_cells_by="My_annotation"
                 #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
)+ggtitle('int.umap')
p1 + p2
ggsave(paste0(SI_HfiD_Enterocyte_monocle_prefix,"Monocle3_umap_embed_", Sys.Date(), ".pdf"), height = 7, width = 16)

## Monocle3聚类分区
monocle3_HfiDs <- cluster_cells(monocle3_HfiDs, resolution = 0.01)
p3 <- plot_cells(monocle3_HfiDs, show_trajectory_graph = FALSE)+ggtitle("Label by clusterID")
p4 <- plot_cells(monocle3_HfiDs, color_cells_by = "partition", show_trajectory_graph = FALSE)+ggtitle("Label by partitionID")
p3 + p4
ggsave(paste0(SI_HfiD_Enterocyte_monocle_prefix,"Monocle3_umap_cluster_", Sys.Date(), ".pdf"), height = 7, width = 16)


## 识别轨迹
monocle3_HfiDs <- learn_graph(monocle3_HfiDs)
p5 <- plot_cells(monocle3_HfiDs, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                 label_branch_points = FALSE)
p5
ggsave(paste0(SI_HfiD_Enterocyte_monocle_prefix,"Monocle3_umap_trajectory_", Sys.Date(), ".pdf"), width = 5, height = 5)


##细胞按拟时排序
# HfiDs <- order_cells(HfiDs) 存在bug，使用辅助线选择root细胞
#p   geom_vline(xintercept = seq(-7,-6,0.25))   geom_hline(yintercept = seq(0,1,0.25))
#embed <- data.frame(Embeddings(all, reduction = "umap"))
#embed <- subset(embed, umap_1 > -6.75 & umap_1 < -6.5 & umap_2 > 0.24 & umap_2 < 0.25)
#root.cell <- rownames(embed)
#monocle3_HfiDs <- order_cells(monocle3_HfiDs, root_cells = root.cell)
monocle3_HfiDs <- order_cells(monocle3_HfiDs)
p6 <- plot_cells(monocle3_HfiDs, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                 label_leaves = F,  label_branch_points = FALSE,graph_label_size = 3,
                 cell_size = 1, alpha = 0.5, trajectory_graph_color = "black",trajectory_graph_segment_size = 1)
p6
ggsave(paste0(SI_HfiD_Enterocyte_monocle_prefix,"Monocle3_umap_trajectory_line_", Sys.Date(), ".pdf"), width = 6.5, height = 5)