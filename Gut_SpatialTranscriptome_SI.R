rm(list = ls())
source("0_ScFunctions.R")
SI_merged_prefix <- "SI_merged_CD4w_removed/"

####################### 0. Pre-processing dataset###############################
#Load Spatial datasets
CD_SI_HFD0w <- readRDS("../raw/0w/sample_CD_SI_bin50_seurat.rds")
CD_SI_HFD4w <- readRDS("../raw/4w/sample_SI_CD4W_bin50_seurat.rds")
CD_SI_HFD8w <- readRDS("../raw/8w/sample_SI_CD8W_bin50_seurat.rds")
HfiD_SI_HFD0w <- readRDS("../raw/0w/sample_HfiD_SI_bin50_seurat.rds")
HfiD_SI_HFD4w <- readRDS("../raw/4w/sample_SI_HfiD4W_bin50_seurat.rds")
HfiD_SI_HFD8w <- readRDS("../raw/8w/sample_SI_HfiD8W_bin50_seurat.rds")

CD_SI_HFD0w$orig.ident <- "CD_SI_HFD0w"
CD_SI_HFD4w$orig.ident <- "CD_SI_HFD4w"
CD_SI_HFD8w$orig.ident <- "CD_SI_HFD8w"
HfiD_SI_HFD0w$orig.ident <- "HfiD_SI_HFD0w"
HfiD_SI_HFD4w$orig.ident <- "HfiD_SI_HFD4w"
HfiD_SI_HFD8w$orig.ident <- "HfiD_SI_HFD8w"

CD_SI_HFD4w <- subset(CD_SI_HFD4w, subset = x < 8750 | y < 20300)
CD_SI_HFD4w <- subset(CD_SI_HFD4w, subset = x < 9050 | y < 20000)
CD_SI_HFD4w <- subset(CD_SI_HFD4w, subset = x < 9200 | y < 19850)
CD_SI_HFD4w <- subset(CD_SI_HFD4w, subset = x < 9800 | y < 19550)
CD_SI_HFD4w <- subset(CD_SI_HFD4w, subset = x < 10300 | y < 19300)
CD_SI_HFD4w <- subset(CD_SI_HFD4w, subset = x < 10550 | y < 18950)
CD_SI_HFD4w <- subset(CD_SI_HFD4w, subset = x < 10900 | y < 18800)
CD_SI_HFD4w <- subset(CD_SI_HFD4w, subset = x < 11100 | y < 18750)

CD_SI_HFD0w <- PercentageFeatureSet(CD_SI_HFD0w, "^mt-", col.name = "percent_mito")
CD_SI_HFD4w <- PercentageFeatureSet(CD_SI_HFD4w, "^mt-", col.name = "percent_mito")
CD_SI_HFD8w <- PercentageFeatureSet(CD_SI_HFD8w, "^mt-", col.name = "percent_mito")
HfiD_SI_HFD0w <- PercentageFeatureSet(HfiD_SI_HFD0w, "^mt-", col.name = "percent_mito")
HfiD_SI_HFD4w <- PercentageFeatureSet(HfiD_SI_HFD4w, "^mt-", col.name = "percent_mito")
HfiD_SI_HFD8w <- PercentageFeatureSet(HfiD_SI_HFD8w, "^mt-", col.name = "percent_mito")

if(F){
  SI_merged <- merge(CD_SI_HFD0w, c(CD_SI_HFD4w, CD_SI_HFD8w, HfiD_SI_HFD0w,HfiD_SI_HFD4w,HfiD_SI_HFD8w))
  VlnPlot(SI_merged, features = c("nFeature_Spatial","nCount_Spatial","percent_mito"), pt.size = 0, group.by = "orig.ident", cols = hughie_color)
  ggsave(paste0(SI_merged_prefix,"All_merged_QC_before_",Sys.Date(),".pdf"), width = 12, height = 8)
  table(SI_merged$orig.ident)
  SI_merged <- subset(SI_merged, subset = nFeature_Spatial > 100 & nCount_Spatial > 100 &percent_mito < 10)
  VlnPlot(SI_merged, features = c("nFeature_Spatial","nCount_Spatial","percent_mito"), pt.size = 0, group.by = "orig.ident", cols = hughie_color)
  ggsave(paste0(SI_merged_prefix,"All_merged_QC_after_",Sys.Date(),".pdf"), width = 12, height = 8)
}


########################## 1. Integrate SI CD & HFiD############################
library(harmony)

mWAT_SAMcolors <- c("#A6CEE3", "#66A5CC","#267CB6",
                    "#A2D48E","#47A93A","#20854EFF")

if(file.exists(paste0(SI_merged_prefix, "Harmony_SI_merged_final.rds"))){
  SI_merged <- readRDS(paste0(SI_merged_prefix, "Harmony_SI_merged_final.rds"))
} else {
  SI_merged <- merge(CD_SI_HFD0w, c(CD_SI_HFD4w, CD_SI_HFD8w, HfiD_SI_HFD0w,HfiD_SI_HFD4w,HfiD_SI_HFD8w))
  SI_merged <- PercentageFeatureSet(SI_merged, "^mt-", col.name = "percent_mito")
  SI_merged <- subset(SI_merged, subset = nFeature_Spatial > 100 & nCount_Spatial > 100 &percent_mito < 10)
  DefaultAssay(SI_merged) <- "SCT"
  names(SI_merged@images) <- c("CD_SI_HFD0w","CD_SI_HFD4w","CD_SI_HFD8w",
                               "HfiD_SI_HFD0w","HfiD_SI_HFD4w","HfiD_SI_HFD8w")
  VariableFeatures(SI_merged) <- c(VariableFeatures(CD_SI_HFD0w), 
                                   VariableFeatures(CD_SI_HFD4w), 
                                   VariableFeatures(CD_SI_HFD8w),
                                   VariableFeatures(HfiD_SI_HFD0w),
                                   VariableFeatures(HfiD_SI_HFD4w),
                                   VariableFeatures(HfiD_SI_HFD8w))
  rm(CD_SI_HFD0w,CD_SI_HFD4w,CD_SI_HFD8w,HfiD_SI_HFD0w,HfiD_SI_HFD4w,HfiD_SI_HFD8w)
  gc()
  
  SI_merged <- ScaleData(SI_merged)
  SI_merged <- RunPCA(SI_merged, verbose = T)
  SI_merged <- RunHarmony(SI_merged, group.by.vars = "orig.ident")
  SI_merged <- RunUMAP(SI_merged, reduction = "harmony", dims = 1:30)
  SI_merged <- FindNeighbors(SI_merged, reduction = "harmony", dims = 1:30)
  SI_merged <- FindClusters(SI_merged, resolution = seq(from=0, by=0.1, length=10))
  
  
  clustree(SI_merged)
  ggsave(paste0(SI_merged_prefix,"Harmony_clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
  
  Idents(SI_merged) <- SI_merged$SCT_snn_res.0.4

  saveRDS(SI_merged, file = paste0(SI_merged_prefix,"Harmony_SI_merged_final.rds"))
  
  my_cols <- c('0' = hughie_color[1], '1' = hughie_color[2], '2' = hughie_color[3],
               '3' = hughie_color[4], '4' = hughie_color[5], '5' = hughie_color[6],
               '6' = hughie_color[7], '7' = hughie_color[8], '8' = hughie_color[9],
               '9' = hughie_color[10], '10' = hughie_color[11], '11' = hughie_color[12],
               '12' = hughie_color[13], '13' = hughie_color[14], 
               '14' = hughie_color[15], '15' = hughie_color[16],'16' = hughie_color[17])
  
  #detemine of res, select res
  SI_merged$SCT_snn_res.0.4 <- factor(SI_merged$SCT_snn_res.0.4, levels = 0:9)
  SI_merged$SCT_snn_res.0.7 <- factor(SI_merged$SCT_snn_res.0.7, levels = 0:14)
  Idents(SI_merged) <- SI_merged$SCT_snn_res.0.4
  Idents(SI_merged) <- SI_merged$SCT_snn_res.0.7
  
  pdf(paste0(SI_merged_prefix,"Harmony_UMAP_res0.7_SpatialMaps_",Sys.Date(),".pdf"), width = 12, height = 8)
  SpatialDimPlot(SI_merged, label = TRUE, label.size = 3, repel = T,ncol = 3,pt.size.factor = 75, cols = my_cols)
  dev.off()
  

  # Visualization of the clusters landscape
  pdf(paste0(SI_merged_prefix,"Harmony_UMAP-res0.7-ByClusters_",Sys.Date(), ".pdf"), width = 16, height = 8)
  DimPlot(SI_merged, reduction = "umap", group.by = c("ident", "orig.ident"), label = T, repel = T, cols = hughie_color)
  dev.off()
  #ggsave(paste0(SI_merged_prefix,"UMAP-res0.2-ByClusters_",Sys.Date(), ".pdf"), width = 16, height = 8)
  
  DimPlot(SI_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = T, repel = T, cols = hughie_color)
  ggsave(paste0(SI_merged_prefix,"Harmony_UMAP-res0.7-Byorig.ident_",Sys.Date(), ".pdf"), width = 14, height = 8)
  
  p1 <- DimPlot(SI_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = TRUE, repel = T, cols = hughie_color)
  p2 <- SpatialDimPlot(SI_merged, label = TRUE, label.size = 3, repel = T,ncol = 3,pt.size.factor = 75, cols = my_cols)
  pdf(paste0(SI_merged_prefix,"Harmony_UMAP-res0.7-Dim+SpatialMaps_",Sys.Date(),".pdf"), width = 24, height = 8)
  p1 + p2
  dev.off()
  
  #Find markers
  SI_merged <- PrepSCTFindMarkers(SI_merged, assay = "SCT", verbose = TRUE)
  SI_merged_markers <- FindAllMarkers(SI_merged, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2, verbose = T)
  write.csv(SI_merged_markers, file = paste0(SI_merged_prefix,"Harmony_Res0.7_markers_",Sys.Date(), ".csv"), quote = F)
  SI_merged_TopMarkers <- SI_merged_markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC) %>% as.data.frame()
  
  #Plot
  VlnPlot(SI_merged, features = SI_merged_TopMarkers$gene, stack=T, flip=T) + NoLegend()
  ggsave(paste0(SI_merged_prefix, "Harmony_SI_All_markers_Vlnplot_",Sys.Date(), ".pdf"), width = 12, height = 30)
  
  DotPlot(SI_merged, features = unique(SI_merged_TopMarkers$gene)) + 
    theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1), axis.text.y = element_text(size=16)) 
  ggsave(paste0(SI_merged_prefix, "Harmony_SI_All_markers_dotplot_",Sys.Date(), ".pdf"), width = 20, height = 12)
  
  DoHeatmap(SI_merged, features = SI_merged_markers$gene)
  ggsave(paste0(SI_merged_prefix, "Harmony_SI_All_markers_heatmap_",Sys.Date(), ".pdf"), width = 15, height = 15)
  
  Plot_Cell_compoistion(Seurat_Obj = SI_merged, OutPrefix = SI_merged_prefix, ColorUse1 = hughie_color, ColorUse2 = mWAT_SAMcolors)
  DEG_enrichment(Seruat_DEG_file = paste0(SI_merged_prefix,"Harmony_Res0.7_markers_",Sys.Date(), ".csv"), showCategoryNum = 15, filterLevel = 4)
  #PlotObjMetrices(Seurat_Obj = SI_merged, OutPrefix= paste0(SI_merged_prefix,"Plot.obj_"), ColorUse = mWAT_SAMcolors)
  
  saveRDS(SI_merged, file = paste0(SI_merged_prefix,"Harmony_SI_merged_final.rds"))
}


################################ 2. Annotation##################################
# library("SCINA")
# exp_SI_merged <- as.matrix(SI_merged@assays$SCT@counts)
# sig_raw <- read.csv("CellMarker_from_ref.csv",header = T)
# sig_raw$Tissue <- gsub(" ", "_", sig_raw$Tissue)
# sig_raw$Cell.Type <- gsub(" ", "_", sig_raw$Cell.Type)
# sig_raw$Tissue.Cell.type <- paste0(sig_raw$Tissue,"_", sig_raw$Cell.Type)
# 
# # for (i in 1:nrow(sig_raw)) {
# #   Cell.Marker <- strsplit(sig_raw$Cell.Marker[i],", ")
# #   sig_SI_merged <- Cell.Marker
# #   sig_SI_merged <- list(sig_raw$Cell.Marker[1], sig_raw$Cell.Marker[2],sig_raw$Cell.Marker[3])
# # }
# 
# 
# sig_SI_merged <- strsplit(sig_raw$Cell.Marker,", ")
# names(sig_SI_merged) <- sig_raw$Tissue.Cell.type
# for (i in 1:nrow(sig_raw)) {
#   if(length(sig_SI_merged[[i]]) >49){
#     sig_SI_merged[[i]] <- sample(sig_SI_merged[[i]], 50, replace = F)
#   }
# }
# results_SI_merged = SCINA(exp_SI_merged, sig_SI_merged, max_iter = 120, convergence_n = 12, 
#                           convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap = FALSE)
# SI_merged$Cell_annotation_ref <- results_SI_merged$cell_labels
# auto_anno <- table(SI_merged$SCT_snn_res.0.7, SI_merged$Cell_annotation_ref)
# write.csv(auto_anno, paste0(SI_merged_prefix, "Auto_annotation_res0.7_ref_", Sys.Date(), ".csv"))

############Manual annotation############

SI_markers <- c("Alpi","Sis","Apoa1","Apoa4",                     #Enterocyte               0,2,3,4,6,7
                "Cps1","Rpsa","Dmbt1",                            #Enterocyte Progenitor    5,8
                "Muc2","Tff3","Pigr",                             #Goblet Cell              1
                "Guca2a","Cpe",                                   #Enteroendocrine Cell     9
                "Lyz1","Defa17","Defa24", "Defa30",               #Paneth Cell              13
                "Acta2","Des","Mylk","Myl9",                      #Myofibroblasts           10,11,14
                "Igkc","Jchain","Iglc1"                           #Plasma B cell            12
)

p1 <- VlnPlot(SI_merged, features = SI_markers, stack=T, flip=T)
#p2 <- VlnPlot(SI_merged, features = SI_markers, stack=T, flip=T, split.by = "orig.ident", cols = hughie_color)
p1
ggsave(paste0(SI_merged_prefix, "SI_markers_VlnPlot_",Sys.Date(), ".pdf"), width = 20, height = 12)

SI_merged <- RenameIdents(SI_merged, 
                          '0' = "Enterocyte-3",
                          '1' = "Goblet Cell",
                          '2' = "Enterocyte-2",
                          '3' = "Enterocyte-4",
                          '4' = "Enterocyte-2",
                          '5' = "Enterocyte Progenitor",
                          '6' = "Enterocyte-1",
                          '7' = "Enterocyte-1",
                          '8' = "Enterocyte Progenitor",
                          '9' = "Enteroendocrine Cell",
                          '10' = "Myofibroblast",
                          '11' = "Myofibroblast",
                          '12' = "Plasma B Cell",
                          '13' = "Paneth Cell",
                          '14' = "Myofibroblast")



SI_merged$My_annotation <- Idents(SI_merged)
SI_merged$My_annotation <- factor(SI_merged$My_annotation, levels = c("Enterocyte-1",
                                                                      "Enterocyte-2",
                                                                      "Enterocyte-3", 
                                                                      "Enterocyte-4",
                                                                      "Enterocyte Progenitor",
                                                                      "Enteroendocrine Cell",
                                                                      "Goblet Cell",
                                                                      "Paneth Cell",
                                                                      "Myofibroblast",
                                                                      "Plasma B Cell"))
Idents(SI_merged) <- SI_merged$My_annotation
saveRDS(SI_merged, file = paste0(SI_merged_prefix,"Harmony_SI_merged_final_annotated.rds"))


pdf(paste0(SI_merged_prefix,"Harmony_UMAP-res0.7-ByClusters_annotated_",Sys.Date(),".pdf"), width = 18, height = 8)
p1 <- DimPlot(SI_merged, reduction = "umap", group.by = c("ident"), label = T, repel = T, cols = hughie_color)
p2 <- DimPlot(SI_merged, reduction = "umap", group.by = c("orig.ident"), label = T, repel = T, cols = mWAT_SAMcolors)
p1 + p2
dev.off()


DimPlot(SI_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = T, repel = T, cols = hughie_color)
ggsave(paste0(SI_merged_prefix,"Harmony_UMAP-res0.7-Byorig.ident_annotated_",Sys.Date(),".pdf"), width = 14, height = 8)


my_cols_annotated <- c("Enterocyte-1" = hughie_color[1],
                       "Enterocyte-2" = hughie_color[2],
                       "Enterocyte-3" = hughie_color[3],
                       "Enterocyte-4" = hughie_color[4],
                       "Enterocyte Progenitor" = hughie_color[5],
                       "Goblet Cell" = hughie_color[6],
                       "Enteroendocrine Cell" = hughie_color[7],
                       "Paneth Cell" = hughie_color[8],
                       "Myofibroblast" = hughie_color[9],
                       "Plasma B Cell" = hughie_color[10])

p1 <- DimPlot(SI_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = TRUE, repel = T, cols = hughie_color)
p2 <- SpatialDimPlot(SI_merged, label = TRUE, label.size = 3, repel = T, ncol = 3, pt.size.factor = 75,cols = my_cols_annotated)
pdf(paste0(SI_merged_prefix,"Harmony_UMAP-res0.7-Dim+SpatialMaps_annotated_",Sys.Date(),".pdf"), width = 24, height = 8)
p1 + p2
dev.off()

pdf(paste0(SI_merged_prefix,"Harmony_UMAP_res0.7_SpatialMaps_annotated_",Sys.Date(),".pdf"), width = 20, height = 8)
p2
dev.off()

##################Find markers#####################
SI_merged <- PrepSCTFindMarkers(SI_merged, assay = "SCT", verbose = TRUE)
SI_merged_markers <- FindAllMarkers(SI_merged, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2, verbose = T)
write.csv(SI_merged_markers, file = paste0(SI_merged_prefix,"Harmony_Res0.7_markers_annotated_",Sys.Date(), ".csv"), quote = F)
SI_merged_TopMarkers <- SI_merged_markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC) %>% as.data.frame()

####################Plot###################
VlnPlot(SI_merged, features = SI_merged_TopMarkers$gene, stack=T, flip=T) + NoLegend()
ggsave(paste0(SI_merged_prefix, "Harmony_SI_All_markers_Vlnplot_annotated_",Sys.Date(), ".pdf"), width = 12, height = 30)

DotPlot(SI_merged, features = unique(SI_merged_TopMarkers$gene)) + 
  theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1), axis.text.y = element_text(size=16)) 
ggsave(paste0(SI_merged_prefix, "Harmony_SI_All_markers_dotplot_annotated_",Sys.Date(), ".pdf"), width = 20, height = 12)

DoHeatmap(SI_merged, features = SI_merged_markers$gene)
ggsave(paste0(SI_merged_prefix, "Harmony_SI_All_markers_heatmap_annotated_",Sys.Date(), ".pdf"), width = 15, height = 15)

Plot_Cell_compoistion(Seurat_Obj = SI_merged, OutPrefix = SI_merged_prefix, ColorUse1 = hughie_color, ColorUse2 = mWAT_SAMcolors)
DEG_enrichment(Seruat_DEG_file = paste0(SI_merged_prefix,"Harmony_Res0.7_markers_annotated_",Sys.Date(), ".csv"), showCategoryNum = 15, filterLevel = 4)

saveRDS(SI_merged, file = paste0(SI_merged_prefix,"Harmony_SI_merged_final_annotated.rds"))
#SI_merged <- readRDS(file = paste0(SI_merged_prefix,"Harmony_SI_merged_final_annotated.rds"))

#######################Final Markers###################

SI_markers <- c("Alpi","Sis","Apoa1","Apoa4",                     #Enterocyte               0,2,3,4,6,7
                "Rpsa","Dmbt1",                                   #Enterocyte Progenitor    5,8
                "Guca2a",                                         #Enteroendocrine Cell     9
                "Muc2","Tff3","Pigr",                             #Goblet Cell              1
                "Lyz1","Defa17","Defa24", "Defa30",               #Paneth Cell              13
                "Acta2","Des","Mylk","Myl9",                      #Myofibroblasts           10,11,14
                "Igkc","Jchain","Iglc1"                           #Plasma B cell            12
)

p1 <- VlnPlot(SI_merged, features = SI_markers, stack=T, flip=T)
p1
ggsave(paste0(SI_merged_prefix, "Harmony_SI_markers_VlnPlot_annotated_",Sys.Date(), ".pdf"), width = 12, height = 12)

SI_markers <- c("Alpi","Apoa1","Apoa4",                     #Enterocyte               0,2,3,4,6,7
                "Rpsa","Dmbt1",                                   #Enterocyte Progenitor    5,8
                "Guca2a",                                         #Enteroendocrine Cell     9
                "Muc2","Pigr",                             #Goblet Cell              1
                "Lyz1","Defa17","Defa24", "Defa30",               #Paneth Cell              13
                "Acta2","Des","Mylk","Myl9",                      #Myofibroblasts           10,11,14
                "Igkc","Jchain","Iglc1"                           #Plasma B cell            12
)


FeatureHeatmapPlot(Seurat_Obj = SI_merged,OutPrefix = SI_merged_prefix,
                   features = SI_markers,cols = my_cols_annotated, group.by = "My_annotation")


DotPlot(SI_merged, features = unique(SI_markers)) + 
  theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1), axis.text.y = element_text(size=16)) 
ggsave(paste0(SI_merged_prefix, "Haromony_SI_markers_dotplot_",Sys.Date(), ".pdf"), width = 12, height = 8)


SpatialFeaturePlot(subset(SI_merged, orig.ident %in% "HfiD_SI_HFD8w"), 
                   features = c("Alpi","Apoa1","Apoa4"), 
                   ncol = 3, 
                   alpha = c(0.1, 1), 
                   crop = T,
                   pt.size.factor = 75)
ggsave(paste0(SI_merged_prefix,"Harmony_SI_Enterocyte_markers_SpatialFeatures_",Sys.Date(), ".pdf"), width = 10, height = 5)


saveRDS(SI_merged, file = paste0(SI_merged_prefix,"Harmony_SI_merged_final_annotated.rds"))



#######################  3. Enterocytes re-annatation ##########################
SI_merged_Enterocyte_prefix <- "SI_Enterocyte_CD4w_removed/"

SI_merged_Enterocyte <- subset(SI_merged, idents = c("Enterocyte-1","Enterocyte-2","Enterocyte-3", "Enterocyte-4"))
Idents(SI_merged_Enterocyte) <- factor(Idents(SI_merged_Enterocyte), levels = c("Enterocyte-1", 
                                                                                "Enterocyte-2",
                                                                                "Enterocyte-3", 
                                                                                "Enterocyte-4"))
SI_merged_Enterocyte$My_annotation <- Idents(SI_merged_Enterocyte)

Enterocyte_markers <- c("Reg1", "Reg3b", "Reg3g","Nlrp6","Il18","Ccl25","Lypd8",#Antimicrobial Program, bottom
                                                                                #absorb middle
                        "Apobec1","Apob","Apoa4","Apoa1","Npc1l1",              #Apolipoproteins Cholesterol
                        "Slc15a1",                                              #Peptides
                        "Slc5a1","Slc2a5","Slc2a2",                             #Carbohydrates
                        "Slc7a9","Slc7a8","Slc7a7",                             #Amino acids
                        "Ada","Egfr","Klf4","Enpp3","Nt5e","Slc28a2"            #shedding top Purine metabolism
)

Enterocyte_markers <- c("Reg1", "Reg3g","Nlrp6","Il18","Ccl25",                 #Antimicrobial Program, bottom
                        "Apobec1","Npc1l1",                                     #Apolipoproteins Cholesterol
                        "Slc15a1",                                              #Peptides
                        "Slc7a9","Slc7a8","Slc7a7",                             #Amino acids
                        "Klf4","Enpp3","Ada","Nt5e","Slc28a2",                   #shedding top Purine metabolism
                        "Cd44","Cd74"
)


p1 <- VlnPlot(SI_merged_Enterocyte, features = Enterocyte_markers, stack=T, flip=T)
p1
ggsave(paste0(SI_merged_Enterocyte_prefix, "Harmony_SI_markers_VlnPlot_annotated_",Sys.Date(), ".pdf"), width = 12, height = 12)


DotPlot(SI_merged_Enterocyte, features = Enterocyte_markers) + 
  theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1), axis.text.y = element_text(size=16)) 
ggsave(paste0(SI_merged_Enterocyte_prefix, "Harmony_SI_All_markers_dotplot_annotated_",Sys.Date(), ".pdf"), width = 20, height = 12)

DoHeatmap(SI_merged_Enterocyte, features = Enterocyte_markers)
ggsave(paste0(SI_merged_Enterocyte_prefix, "Harmony_SI_All_markers_heatmap_annotated_",Sys.Date(), ".pdf"), width = 15, height = 15)

FeatureHeatmapPlot(Seurat_Obj = SI_merged_Enterocyte,OutPrefix = SI_merged_Enterocyte_prefix,
                   features = Enterocyte_markers,cols = hughie_color[1:4], group.by = "My_annotation")


Selected_enrichment_Plot(GO_enrichment_file = paste0(SI_merged_Enterocyte_prefix,"Harmony_Res0.7_markers_annotated_2024-10-01_Selected_GO-BP_list_2024-10-01.csv"),
                         OutPrefix = SI_merged_Enterocyte_prefix, 
                         cluster_order = c("Enterocyte Progenitor","Enterocyte-1","Enterocyte-2","Enterocyte-3","Enterocyte-4"),
                         levels = c("RNA splicing",
                                    "regulation of RNA splicing",
                                    "cytoplasmic translation",
                                    "positive regulation of translation",
                                    "glutathione metabolic process",
                                    "glutathione biosynthetic process",
                                    "mitochondrial ATP synthesis coupled electron transport",
                                    "mitochondrial transport",
                                    "intestinal absorption",
                                    "organic anion transport",
                                    "negative regulation of cell-cell adhesion",
                                    "negative regulation of cell adhesion"))

###re-name enterocyte
SI_merged <- readRDS(file = paste0(SI_merged_prefix,"Harmony_SI_merged_final_annotated.rds"))
SI_merged <- RenameIdents(SI_merged, 
                          'Enterocyte-1' = "Crypt Mitochondrial Enterocyte",
                          'Enterocyte-2' = "Immune Related Enterocyte",
                          'Enterocyte-3' = "Villus Absorption Enterocyte",
                          'Enterocyte-4' = "Villus Tip Enterocyte")
SI_merged$My_annotation <- Idents(SI_merged)
SI_merged$My_annotation <- factor(SI_merged$My_annotation, levels = c("Enterocyte Progenitor",
                                                                      "Crypt Mitochondrial Enterocyte",
                                                                      "Villus Absorption Enterocyte",
                                                                      "Villus Tip Enterocyte", 
                                                                      "Immune Related Enterocyte",
                                                                      "Enteroendocrine Cell",
                                                                      "Goblet Cell",
                                                                      "Paneth Cell",
                                                                      "Myofibroblast",
                                                                      "Plasma B Cell"))
Idents(SI_merged) <- SI_merged$My_annotation
saveRDS(SI_merged, file = paste0(SI_merged_prefix,"Harmony_SI_merged_final_final_annotated.rds"))

########re-plot after re-name
pdf(paste0(SI_merged_prefix,"Harmony_UMAP-res0.7-ByClusters_annotated_final_",Sys.Date(),".pdf"), width = 18, height = 8)
p1 <- DimPlot(SI_merged, reduction = "umap", group.by = c("ident"), label = T, repel = T, cols = hughie_color)
p2 <- DimPlot(SI_merged, reduction = "umap", group.by = c("orig.ident"), label = T, repel = T, cols = mWAT_SAMcolors)
p1 + p2
dev.off()

DimPlot(SI_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = T, repel = T, cols = hughie_color)
ggsave(paste0(SI_merged_prefix,"Harmony_UMAP-res0.7-Byorig.ident_annotated_final_",Sys.Date(),".pdf"), width = 14, height = 8)

my_cols_annotated <- c("Enterocyte Progenitor" = hughie_color[1],
                       "Crypt Mitochondrial Enterocyte" = hughie_color[2],
                       "Villus Absorption Enterocyte" = hughie_color[3],
                       "Villus Tip Enterocyte" = hughie_color[4],
                       "Immune Related Enterocyte" = hughie_color[5],
                       "Enteroendocrine Cell" = hughie_color[6],
                       "Goblet Cell" = hughie_color[7],
                       "Paneth Cell" = hughie_color[8],
                       "Myofibroblast" = hughie_color[9],
                       "Plasma B Cell" = hughie_color[10])

p1 <- DimPlot(SI_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = TRUE, repel = T, cols = hughie_color)
p2 <- SpatialDimPlot(SI_merged, label = TRUE, label.size = 3, repel = T, ncol = 3, pt.size.factor = 75,cols = my_cols_annotated)
pdf(paste0(SI_merged_prefix,"Harmony_UMAP-res0.7-Dim+SpatialMaps_annotated_final_",Sys.Date(),".pdf"), width = 24, height = 8)
p1 + p2
dev.off()

pdf(paste0(SI_merged_prefix,"Harmony_UMAP_res0.7_SpatialMaps_annotated_final_",Sys.Date(),".pdf"), width = 20, height = 8)
p2
dev.off()

#######Find markers
SI_merged <- PrepSCTFindMarkers(SI_merged, assay = "SCT", verbose = TRUE)
SI_merged_markers <- FindAllMarkers(SI_merged, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2, verbose = T)
write.csv(SI_merged_markers, file = paste0(SI_merged_prefix,"Harmony_Res0.7_markers_annotated_final_",Sys.Date(), ".csv"), quote = F)
SI_merged_TopMarkers <- SI_merged_markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC) %>% as.data.frame()

#######Plot
VlnPlot(SI_merged, features = SI_merged_TopMarkers$gene, stack=T, flip=T) + NoLegend()
ggsave(paste0(SI_merged_prefix, "Harmony_SI_All_markers_Vlnplot_annotated_final_",Sys.Date(), ".pdf"), width = 12, height = 30)

DotPlot(SI_merged, features = unique(SI_merged_TopMarkers$gene)) + 
  theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1), axis.text.y = element_text(size=16)) 
ggsave(paste0(SI_merged_prefix, "Harmony_SI_All_markers_dotplot_annotated_final_",Sys.Date(), ".pdf"), width = 20, height = 12)

DoHeatmap(SI_merged, features = SI_merged_markers$gene)
ggsave(paste0(SI_merged_prefix, "Harmony_SI_All_markers_heatmap_annotated_final_",Sys.Date(), ".pdf"), width = 15, height = 15)

Plot_Cell_compoistion(Seurat_Obj = SI_merged, OutPrefix = SI_merged_prefix, ColorUse1 = hughie_color, ColorUse2 = mWAT_SAMcolors)
Plot_Cell_compoistion(Seurat_Obj = subset(SI_merged, idents = c("Villus Tip Enterocyte",
                                                                "Immune Related Enterocyte",
                                                                "Goblet Cell")), 
                      OutPrefix = paste0(SI_merged_prefix,"Diff_"), 
                      ColorUse1 = c(hughie_color[4],hughie_color[5],hughie_color[7]), 
                      ColorUse2 = mWAT_SAMcolors)

DEG_enrichment(Seruat_DEG_file = paste0(SI_merged_prefix,"Harmony_Res0.7_markers_annotated_final_",Sys.Date(), ".csv"), showCategoryNum = 15, filterLevel = 4)

SI_markers_final <- c(
  "Rpsa","Dmbt1",                                   #Enterocyte Progenitor    5,8
  "Alpi","Apoa1","Apoa4",                     #Enterocyte               0,2,3,4,6,7
  "Neurod1","Neurog3",                                         #Enteroendocrine Cell     9
  "Muc2","Pigr",                             #Goblet Cell              1
  "Lyz1","Defa17","Defa24", "Defa30",               #Paneth Cell              13
  "Acta2","Des","Mylk","Myl9",                      #Myofibroblasts           10,11,14
  "Igkc","Jchain","Iglc1"                           #Plasma B cell            12
)
                
                
FeatureHeatmapPlot(Seurat_Obj = SI_merged,OutPrefix = SI_merged_prefix,
                   features = SI_markers_final,cols = hughie_color[1:10], group.by = "My_annotation")


saveRDS(SI_merged, file = paste0(SI_merged_prefix,"Harmony_SI_merged_final_final_annotated.rds"))

#re plot Enterocytes marker after re-name
SI_merged_Enterocyte <- subset(SI_merged, idents = c("Crypt Mitochondrial Enterocyte",
                                                     "Villus Absorption Enterocyte",
                                                     "Villus Tip Enterocyte", 
                                                     "Immune Related Enterocyte"))

FeatureHeatmapPlot(Seurat_Obj = SI_merged_Enterocyte,OutPrefix = SI_merged_Enterocyte_prefix,
                   features = Enterocyte_markers,cols = hughie_color[2:5], group.by = "My_annotation")


Selected_enrichment_Plot(GO_enrichment_file = paste0(SI_merged_Enterocyte_prefix,"Harmony_Res0.7_markers_annotated_2024-10-01_Selected_GO-BP_list_2024-10-01.csv"),
                         OutPrefix = SI_merged_Enterocyte_prefix, 
                         cluster_order = c("Enterocyte Progenitor",
                                           "Crypt Mitochondrial Enterocyte",
                                           "Villus Absorption Enterocyte",
                                           "Villus Tip Enterocyte", 
                                           "Immune Related Enterocyte"),
                         levels = c("RNA splicing",
                                    "regulation of RNA splicing",
                                    "cytoplasmic translation",
                                    "positive regulation of translation",
                                    "glutathione metabolic process",
                                    "glutathione biosynthetic process",
                                    "mitochondrial ATP synthesis coupled electron transport",
                                    "mitochondrial transport",
                                    "intestinal absorption",
                                    "organic anion transport",
                                    "negative regulation of cell-cell adhesion",
                                    "negative regulation of cell adhesion",
                                    "B cell homeostasis",
                                    "lymphocyte homeostasis"))



####################### 4. Enterocytes Merge with Cell 2018 ####################
SI_merged_Enterocyte_overlap_prefix <- "SI_Enterocyte_CD4w_removed/Andreas_Cell_2018_"

SI_merged_Enterocyte <- subset(SI_merged, idents = c("Enterocyte Progenitor",
                                                     "Crypt Mitochondrial Enterocyte",
                                                     "Villus Absorption Enterocyte",
                                                     "Villus Tip Enterocyte", 
                                                     "Immune Related Enterocyte"))
SI_merged_Enterocyte <- subset(SI_merged_Enterocyte, downsample = 400)
Enterocyte_cell_2018 <- read.table(file = paste0(SI_merged_Enterocyte_prefix, "Cell_2018/table_B_scRNAseq_UMI_counts.tsv"), row.names = 1,header = T)
Enterocyte_cell_2018_meta <- read.table(file = paste0(SI_merged_Enterocyte_prefix, "Cell_2018/table_C_scRNAseq_tsne_coordinates_zones.tsv"), row.names = 1,header = T)
Enterocyte_cell_2018 <- CreateSeuratObject(counts = Enterocyte_cell_2018, project = "Cell_2018", meta.data = Enterocyte_cell_2018_meta)
Enterocyte_cell_2018 <- NormalizeData(Enterocyte_cell_2018)

ObjList = c(SI_merged_Enterocyte, Enterocyte_cell_2018)
features <- SelectIntegrationFeatures(object.list = ObjList, nfeatures = 2000)
anchors <-  FindIntegrationAnchors(object.list = ObjList, anchor.features = features)
  
combined <- IntegrateData(anchorset = anchors)
#houyu pepline
Enterocyte_overlap_merged <- ScaleData(combined) %>% RunPCA(npcs = 20) %>% RunUMAP(reduction = "pca", dims = 1:20) %>% 
    FindNeighbors(reduction = "pca", dims = 1:20) %>% FindClusters(resolution = seq(from=0, by=0.05, length=9))

Enterocyte_overlap_merged$orig.ident <- gsub("^.{11}$", "Spatial_Enterocytes", Enterocyte_overlap_merged$orig.ident)
Enterocyte_overlap_merged$orig.ident <- gsub("^.{13}$", "Spatial_Enterocytes", Enterocyte_overlap_merged$orig.ident)
Enterocyte_overlap_merged$orig.ident <- gsub("^.{9}$", "Andreas_Cell_2018", Enterocyte_overlap_merged$orig.ident)
Enterocyte_overlap_merged$orig.ident <- as.factor(Enterocyte_overlap_merged$orig.ident)

  
Idents(Enterocyte_overlap_merged) <- Enterocyte_overlap_merged$integrated_snn_res.0.35
(p00 <- DimPlot(Enterocyte_overlap_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = hughie_color, label.color = "white"))
(p01 <- DimPlot(Enterocyte_overlap_merged, reduction = "umap", group.by = "orig.ident", label = F, repel = T, label.box = T, order = T, cols = c("#DC0000FF","#1F78B4")))
Plot_Cell_compoistion(Seurat_Obj = Enterocyte_overlap_merged, 
                      OutPrefix = SI_merged_Enterocyte_overlap_prefix, 
                      ColorUse1 = rev(hughie_color[1:7]), 
                      ColorUse2 = c("#DC0000FF","#1F78B4"))

Enterocyte_overlap_merged$zone <- factor(Enterocyte_overlap_merged$zone, levels = c("Crypt","V1","V2","V3","V4","V5","V6"))
Idents(Enterocyte_overlap_merged) <- Enterocyte_overlap_merged$zone
(p10 <- DimPlot(Enterocyte_overlap_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = c("#3E499A","#5C4A9B","#8F61A7","#C888AE","#DBA990","#EDD368","#F1E60D")))

Idents(Enterocyte_overlap_merged) <- Enterocyte_overlap_merged$My_annotation
(p20 <- DimPlot(Enterocyte_overlap_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = c("#0072B5","#BC3C29","#E18727","#20854E","#6F99AD")))

(p00+p01)/(p10)/(p20)
ggsave(paste0(SI_merged_Enterocyte_overlap_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 10, height = 15)

saveRDS(Enterocyte_overlap_merged, file = paste0(SI_merged_Enterocyte_overlap_prefix,"Seurat.rds"))



############################### 5. Scvelo#######################################
SI_merged <- readRDS(file = paste0(SI_merged_prefix,"Harmony_SI_merged_final_final_annotated.rds"))
library(reticulate)

conda_list()
use_condaenv("r-velo", required = TRUE)
scv <- import("scvelo")
SI_CD_Enterocyte_velocyto_prefix <- "SI_Enterocyte_velocity/CD_group/"
SI_HfiD_Enterocyte_velocyto_prefix <- "SI_Enterocyte_velocity/HfiD_group/"


########################Prepare emat and nmat############
ldat.CD.SI.0w <- read.loom.matrices.spatial(file = "../PZ_analysis/raw/0w/sample_CD_SI_bin50_seurat.loom")
emat.CD.SI.0w <- ldat.CD.SI.0w$spliced
#colnames(emat.CD.SI.0w) <- gsub("mWAT_CD2:", "CD_HFD0w_", colnames(emat.CD.0w))
#colnames(emat.CD.SI.0w) <- gsub("x", "-1", colnames(emat.CD.0w))

ldat.CD.SI.4w <- read.loom.matrices.spatial(file = "../PZ_analysis/raw/4w/sample_SI_CD4W_bin50_seurat.loom")
emat.CD.SI.4w <- ldat.CD.SI.4w$spliced

ldat.CD.SI.8w <- read.loom.matrices.spatial(file = "../PZ_analysis/raw/8w/sample_SI_CD8W_bin50_seurat.loom")
emat.CD.SI.8w <- ldat.CD.SI.8w$spliced

ldat.HfiD.SI.0w <- read.loom.matrices.spatial(file = "../PZ_analysis/raw/0w/sample_HfiD_SI_bin50_seurat.loom")
emat.HfiD.SI.0w <- ldat.HfiD.SI.0w$spliced

ldat.HfiD.SI.4w <- read.loom.matrices.spatial(file = "../PZ_analysis/raw/4w/sample_SI_HfiD4W_bin50_seurat.loom")
emat.HfiD.SI.4w <- ldat.HfiD.SI.4w$spliced

ldat.HfiD.SI.8w <- read.loom.matrices.spatial(file = "../PZ_analysis/raw/8w/sample_SI_HfiD8W_bin50_seurat.loom")
emat.HfiD.SI.8w <- ldat.HfiD.SI.8w$spliced

nmat.CD.SI.0w <- ldat.CD.SI.0w$unspliced
nmat.CD.SI.4w <- ldat.CD.SI.4w$unspliced
nmat.CD.SI.8w <- ldat.CD.SI.8w$unspliced
nmat.HfiD.SI.0w <- ldat.HfiD.SI.0w$unspliced
nmat.HfiD.SI.4w <- ldat.HfiD.SI.4w$unspliced
nmat.HfiD.SI.8w <- ldat.HfiD.SI.8w$unspliced

SI_merged_Enterocyte <- subset(SI_merged, My_annotation %in% c("Enterocyte Progenitor",
                                                               "Crypt Mitochondrial Enterocyte",
                                                               "Villus Absorption Enterocyte",
                                                               "Villus Tip Enterocyte", 
                                                               "Immune Related Enterocyte",
                                                               "Goblet Cell"))
Idents(SI_merged_Enterocyte)

cell.chose <- as.data.frame(SI_merged_Enterocyte$My_annotation)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(SI_merged_Enterocyte@assays$SCT@meta.features)
feature.chose$feature <- row.names(feature.chose)

mat.feature.chose <- intersect(row.names(emat.CD.SI.0w),row.names(emat.CD.SI.4w))
mat.feature.chose <- intersect(mat.feature.chose, row.names(emat.CD.SI.8w))
mat.feature.chose <- intersect(mat.feature.chose, row.names(feature.chose))

emat.CD.SI.0w <- emat.CD.SI.0w[rownames(emat.CD.SI.0w) %in% mat.feature.chose, ]
emat.CD.SI.4w <- emat.CD.SI.4w[rownames(emat.CD.SI.4w) %in% mat.feature.chose, ]
emat.CD.SI.8w <- emat.CD.SI.8w[rownames(emat.CD.SI.8w) %in% mat.feature.chose, ]

emat <- cbind(emat.CD.SI.0w,emat.CD.SI.4w,emat.CD.SI.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]

nmat.CD.SI.0w <- nmat.CD.SI.0w[rownames(nmat.CD.SI.0w) %in% mat.feature.chose, ]
nmat.CD.SI.4w <- nmat.CD.SI.4w[rownames(nmat.CD.SI.4w) %in% mat.feature.chose, ]
nmat.CD.SI.8w <- nmat.CD.SI.8w[rownames(nmat.CD.SI.8w) %in% mat.feature.chose, ]

nmat <- cbind(nmat.CD.SI.0w,nmat.CD.SI.4w,nmat.CD.SI.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% row.names(cell.chose)]

cell.type <- cell.chose[colnames(emat),1]

saveRDS(emat, paste0("SI_Enterocyte_velocity/raw/emat.CD.rds"))
saveRDS(nmat, paste0("SI_Enterocyte_velocity/raw/nmat.CD.rds"))

#确保细胞数和emat，nmat一致
SI_merged_Enterocyte <- subset(SI_merged_Enterocyte, cells = colnames(emat))
SI_merged_Enterocyte <- subset(SI_merged_Enterocyte, features = rownames(emat))
saveRDS(SI_merged_Enterocyte, paste0(SI_CD_Enterocyte_velocyto_prefix,"SI_merged_Enterocyte.rds"))





mat.feature.chose <- intersect(row.names(emat.HfiD.SI.0w),row.names(emat.HfiD.SI.4w))
mat.feature.chose <- intersect(mat.feature.chose, row.names(emat.HfiD.SI.8w))
mat.feature.chose <- intersect(mat.feature.chose, row.names(feature.chose))

emat.HfiD.SI.0w <- emat.HfiD.SI.0w[rownames(emat.HfiD.SI.0w) %in% mat.feature.chose, ]
emat.HfiD.SI.4w <- emat.HfiD.SI.4w[rownames(emat.HfiD.SI.4w) %in% mat.feature.chose, ]
emat.HfiD.SI.8w <- emat.HfiD.SI.8w[rownames(emat.HfiD.SI.8w) %in% mat.feature.chose, ]

emat <- cbind(emat.HfiD.SI.0w,emat.HfiD.SI.4w,emat.HfiD.SI.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]

nmat.HfiD.SI.0w <- nmat.HfiD.SI.0w[rownames(nmat.HfiD.SI.0w) %in% mat.feature.chose, ]
nmat.HfiD.SI.4w <- nmat.HfiD.SI.4w[rownames(nmat.HfiD.SI.4w) %in% mat.feature.chose, ]
nmat.HfiD.SI.8w <- nmat.HfiD.SI.8w[rownames(nmat.HfiD.SI.8w) %in% mat.feature.chose, ]

nmat <- cbind(nmat.HfiD.SI.0w,nmat.HfiD.SI.4w,nmat.HfiD.SI.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% row.names(cell.chose)]

cell.type <- cell.chose[colnames(emat),1]

saveRDS(emat, paste0("SI_Enterocyte_velocity/raw/emat.HfiD.rds"))
saveRDS(nmat, paste0("SI_Enterocyte_velocity/raw/nmat.HfiD.rds"))





SI_merged_Enterocyte_CD <- subset(SI_merged_Enterocyte, orig.ident %in% c("CD_SI_HFD0w", "CD_SI_HFD4w",  "CD_SI_HFD8w"))
SI_merged_Enterocyte_HfiD <- subset(SI_merged_Enterocyte, orig.ident %in% c("HfiD_SI_HFD0w", "HfiD_SI_HFD4w", "HfiD_SI_HFD8w"))


####################CD group scvelo##################################
Idents(SI_merged_Enterocyte_CD)
DimPlot(SI_merged_Enterocyte_CD, reduction = "umap", label = F, repel = T, label.box = T, order = T, cols = hughie_color)
ggsave(paste0(SI_CD_Enterocyte_velocyto_prefix,"UMAP_Cluster_",Sys.Date(),".pdf"), width = 8, height = 5)

cell.chose <- as.data.frame(SI_merged_Enterocyte_CD$My_annotation)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(SI_merged_Enterocyte_CD@assays$SCT@meta.features)
feature.chose$feature <- row.names(feature.chose)

mat.feature.chose <- intersect(row.names(emat.CD.SI.0w),row.names(emat.CD.SI.4w))
mat.feature.chose <- intersect(mat.feature.chose, row.names(emat.CD.SI.8w))
mat.feature.chose <- intersect(mat.feature.chose, row.names(feature.chose))

emat.CD.SI.0w <- emat.CD.SI.0w[rownames(emat.CD.SI.0w) %in% mat.feature.chose, ]
emat.CD.SI.4w <- emat.CD.SI.4w[rownames(emat.CD.SI.4w) %in% mat.feature.chose, ]
emat.CD.SI.8w <- emat.CD.SI.8w[rownames(emat.CD.SI.8w) %in% mat.feature.chose, ]

emat <- cbind(emat.CD.SI.0w,emat.CD.SI.4w,emat.CD.SI.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]

nmat.CD.SI.0w <- nmat.CD.SI.0w[rownames(nmat.CD.SI.0w) %in% mat.feature.chose, ]
nmat.CD.SI.4w <- nmat.CD.SI.4w[rownames(nmat.CD.SI.4w) %in% mat.feature.chose, ]
nmat.CD.SI.8w <- nmat.CD.SI.8w[rownames(nmat.CD.SI.8w) %in% mat.feature.chose, ]

nmat <- cbind(nmat.CD.SI.0w,nmat.CD.SI.4w,nmat.CD.SI.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% row.names(cell.chose)]

cell.type <- cell.chose[colnames(emat),1]

saveRDS(emat, paste0(SI_CD_Enterocyte_velocyto_prefix,"emat.CD.rds"))
saveRDS(nmat, paste0(SI_CD_Enterocyte_velocyto_prefix,"nmat.CD.rds"))

#确保细胞数和emat，nmat一致
SI_merged_Enterocyte_CD <- subset(SI_merged_Enterocyte_CD, cells = colnames(emat))
SI_merged_Enterocyte_CD <- subset(SI_merged_Enterocyte_CD, features = rownames(emat))
saveRDS(SI_merged_Enterocyte_CD, paste0(SI_CD_Enterocyte_velocyto_prefix,"SI_merged_Enterocyte_CD.rds"))


#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(SI_merged_Enterocyte_CD@reductions$umap@cell.embeddings)))

#获取每个细胞在图上的坐标位置
emb <- SI_merged_Enterocyte_CD@reductions$umap@cell.embeddings
pca <- SI_merged_Enterocyte_CD@reductions$pca@cell.embeddings

#Generate_scvelo_InputFiles(SI_merged_Enterocyte_CD,OutPrefix = paste0(SI_CD_Enterocyte_velocyto_prefix,"Scvelo_files_"))

counts_matrix <- GetAssayData(SI_merged_Enterocyte_CD, assay='SCT', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(SI_merged_Enterocyte_CD@meta.data)
#rownames(dfobs) <- cells
dfvar <- data.frame(genes_attributes)
rownames(dfvar) <- dfvar$gene
adata <- ad$AnnData(
  X=t(counts_matrix),
  obs=dfobs,
  var=dfvar,
  layers=list('spliced'=t(emat), 'unspliced'=t(nmat)),
  obsm=list('X_umap'=emb, 'X_pca'=pca[,1:2]) 
)


## get embedding
emb <- as.matrix(adata$obsm['X_umap'])
#clusters <- as.factor(adata$obs$Myannotation_Ad_ASPC)
#rownames(emb) <- names(clusters) <- as.matrix(adata$obs_names$values)

rm(emat.CD.SI.0w,emat.CD.SI.4w,emat.CD.SI.8w,nmat.CD.SI.0w,nmat.CD.SI.4w,nmat.CD.SI.8w,SI_merged,SI_merged_Enterocyte)
gc()

## run scvelo dynamic model
scv$pp$filter_genes(adata) ## filter
scv$pp$moments(adata) ## normalize and compute moments
scv$tl$recover_dynamics(adata) ## model

## plot (creates pop up window)
scv$tl$velocity(adata, mode='dynamical')
scv$tl$velocity_graph(adata)
pdf(paste0(SI_CD_Enterocyte_velocyto_prefix,"Scvelo_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='My_annotation',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300)
dev.off()

pdf(paste0(SI_CD_Enterocyte_velocyto_prefix,"Scvelo_PCA_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='My_annotation',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300) ## other embedding
dev.off()

##save
adata$write(paste0(SI_CD_Enterocyte_velocyto_prefix,'Scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(mWAT_Ad_ASPC_CD_velocyto_prefix,'scvelo_anndata.h5ad'))

####################HfiD group scvelo##################################
Idents(SI_merged_Enterocyte_HfiD)
DimPlot(SI_merged_Enterocyte_HfiD, reduction = "umap", label = F, repel = T, label.box = T, order = T, cols = hughie_color)
ggsave(paste0(SI_HfiD_Enterocyte_velocyto_prefix,"UMAP_Cluster_",Sys.Date(),".pdf"), width = 8, height = 5)

cell.chose <- as.data.frame(SI_merged_Enterocyte_HfiD$My_annotation)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(SI_merged_Enterocyte_HfiD@assays$SCT@meta.features)
feature.chose$feature <- row.names(feature.chose)

mat.feature.chose <- intersect(row.names(emat.HfiD.SI.0w),row.names(emat.HfiD.SI.4w))
mat.feature.chose <- intersect(mat.feature.chose, row.names(emat.HfiD.SI.8w))
mat.feature.chose <- intersect(mat.feature.chose, row.names(feature.chose))

emat.HfiD.SI.0w <- emat.HfiD.SI.0w[rownames(emat.HfiD.SI.0w) %in% mat.feature.chose, ]
emat.HfiD.SI.4w <- emat.HfiD.SI.4w[rownames(emat.HfiD.SI.4w) %in% mat.feature.chose, ]
emat.HfiD.SI.8w <- emat.HfiD.SI.8w[rownames(emat.HfiD.SI.8w) %in% mat.feature.chose, ]

emat <- cbind(emat.HfiD.SI.0w,emat.HfiD.SI.4w,emat.HfiD.SI.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]

nmat.HfiD.SI.0w <- nmat.HfiD.SI.0w[rownames(nmat.HfiD.SI.0w) %in% mat.feature.chose, ]
nmat.HfiD.SI.4w <- nmat.HfiD.SI.4w[rownames(nmat.HfiD.SI.4w) %in% mat.feature.chose, ]
nmat.HfiD.SI.8w <- nmat.HfiD.SI.8w[rownames(nmat.HfiD.SI.8w) %in% mat.feature.chose, ]

nmat <- cbind(nmat.HfiD.SI.0w,nmat.HfiD.SI.4w,nmat.HfiD.SI.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% row.names(cell.chose)]

cell.type <- cell.chose[colnames(emat),1]

saveRDS(emat, paste0(SI_HfiD_Enterocyte_velocyto_prefix,"emat.HfiD.rds"))
saveRDS(nmat, paste0(SI_HfiD_Enterocyte_velocyto_prefix,"nmat.HfiD.rds"))

#确保细胞数和emat，nmat一致
SI_merged_Enterocyte_HfiD <- subset(SI_merged_Enterocyte_HfiD, cells = colnames(emat))
SI_merged_Enterocyte_HfiD <- subset(SI_merged_Enterocyte_HfiD, features = rownames(emat))
saveRDS(SI_merged_Enterocyte_HfiD, paste0(SI_HfiD_Enterocyte_velocyto_prefix,"SI_merged_Enterocyte_HfiD.rds"))

#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(SI_merged_Enterocyte_HfiD@reductions$umap@cell.embeddings)))

#获取每个细胞在图上的坐标位置
emb <- SI_merged_Enterocyte_HfiD@reductions$umap@cell.embeddings
pca <- SI_merged_Enterocyte_HfiD@reductions$pca@cell.embeddings

#Generate_scvelo_InputFiles(SI_merged_Enterocyte_HfiD,OutPrefix = paste0(SI_HfiD_Enterocyte_velocyto_prefix,"Scvelo_files_"))

counts_matrix <- GetAssayData(SI_merged_Enterocyte_HfiD, assay='SCT', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(SI_merged_Enterocyte_HfiD@meta.data)
#rownames(dfobs) <- cells
dfvar <- data.frame(genes_attributes)
rownames(dfvar) <- dfvar$gene
adata <- ad$AnnData(
  X=t(counts_matrix),
  obs=dfobs,
  var=dfvar,
  layers=list('spliced'=t(emat), 'unspliced'=t(nmat)),
  obsm=list('X_umap'=emb, 'X_pca'=pca[,1:2]) 
)


## get embedding
emb <- as.matrix(adata$obsm['X_umap'])
#clusters <- as.factor(adata$obs$Myannotation_Ad_ASPC)
#rownames(emb) <- names(clusters) <- as.matrix(adata$obs_names$values)


## run scvelo dynamic model
scv$pp$filter_genes(adata) ## filter
scv$pp$moments(adata) ## normalize and compute moments
scv$tl$recover_dynamics(adata) ## model

## plot (creates pop up window)
scv$tl$velocity(adata, mode='dynamical')
scv$tl$velocity_graph(adata)
pdf(paste0(SI_HfiD_Enterocyte_velocyto_prefix,"Scvelo_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='My_annotation',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300)
dev.off()

pdf(paste0(SI_HfiD_Enterocyte_velocyto_prefix,"Scvelo_PCA_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='My_annotation',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300) ## other embedding
dev.off()

##save
adata$write(paste0(SI_HfiD_Enterocyte_velocyto_prefix,'Scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(mWAT_Ad_ASPC_HfiD_velocyto_prefix,'scvelo_anndata.h5ad'))