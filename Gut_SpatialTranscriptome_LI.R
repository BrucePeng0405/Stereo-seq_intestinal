rm(list = ls())
source("0_ScFunctions.R")
LI_merged_prefix <- "1.LI_merged/"


####################### 0. Pre-processing dataset###############################
#Load Spatial datasets
CD_LI_HFD0w <- readRDS("../raw/0w/sample_CD_LI_bin50_seurat.rds")
CD_LI_HFD4w <- readRDS("../raw/4w/sample_LI_CD4W_bin50_seurat.rds")
CD_LI_HFD8w <- readRDS("../raw/8w/sample_LI_CD8W_bin50_seurat.rds")
HfiD_LI_HFD0w <- readRDS("../raw/0w/sample_HfiD_LI_bin50_seurat.rds")
HfiD_LI_HFD4w <- readRDS("../raw/4w/sample_LI_HfiD4W_bin50_seurat.rds")
HfiD_LI_HFD8w <- readRDS("../raw/8w/sample_LI_HfiD8W_bin50_seurat.rds")

CD_LI_HFD0w$orig.ident <- "CD_LI_HFD0w"
CD_LI_HFD4w$orig.ident <- "CD_LI_HFD4w"
CD_LI_HFD8w$orig.ident <- "CD_LI_HFD8w"
HfiD_LI_HFD0w$orig.ident <- "HfiD_LI_HFD0w"
HfiD_LI_HFD4w$orig.ident <- "HfiD_LI_HFD4w"
HfiD_LI_HFD8w$orig.ident <- "HfiD_LI_HFD8w"

CD_LI_HFD0w <- subset(CD_LI_HFD0w, subset = !((x > 21450 & x < 21650) & y == 11450))
CD_LI_HFD0w <- subset(CD_LI_HFD0w, subset = !((x > 21400 & x < 21750) & y == 11500))
CD_LI_HFD0w <- subset(CD_LI_HFD0w, subset = !((x > 21350 & x < 21800) & y == 11550))
CD_LI_HFD0w <- subset(CD_LI_HFD0w, subset = !((x > 21300 & x < 21850) & (y > 11550 & y < 11800)))
CD_LI_HFD0w <- subset(CD_LI_HFD0w, subset = !((x > 21350 & x < 21850) & (y > 11750 & y < 11900)))
CD_LI_HFD0w <- subset(CD_LI_HFD0w, subset = !((x > 21400 & x < 21850) & y == 11900))
CD_LI_HFD0w <- subset(CD_LI_HFD0w, subset = !((x > 21550 & x < 21800) & y == 11950))

CD_LI_HFD0w <- PercentageFeatureSet(CD_LI_HFD0w, "^mt-", col.name = "percent_mito")
CD_LI_HFD4w <- PercentageFeatureSet(CD_LI_HFD4w, "^mt-", col.name = "percent_mito")
CD_LI_HFD8w <- PercentageFeatureSet(CD_LI_HFD8w, "^mt-", col.name = "percent_mito")
HfiD_LI_HFD0w <- PercentageFeatureSet(HfiD_LI_HFD0w, "^mt-", col.name = "percent_mito")
HfiD_LI_HFD4w <- PercentageFeatureSet(HfiD_LI_HFD4w, "^mt-", col.name = "percent_mito")
HfiD_LI_HFD8w <- PercentageFeatureSet(HfiD_LI_HFD8w, "^mt-", col.name = "percent_mito")

if(F){
  LI_merged <- merge(CD_LI_HFD0w, c(CD_LI_HFD4w, CD_LI_HFD8w, HfiD_LI_HFD0w,HfiD_LI_HFD4w,HfiD_LI_HFD8w))
  VlnPlot(LI_merged, features = c("nFeature_Spatial","nCount_Spatial","percent_mito"), pt.size = 0, group.by = "orig.ident", cols = hughie_color)
  ggsave(paste0(LI_merged_prefix,"All_merged_QC_before_",Sys.Date(),".pdf"), width = 12, height = 8)
  table(LI_merged$orig.ident)
  LI_merged <- subset(LI_merged, subset = nFeature_Spatial > 100 & nCount_Spatial > 100 &percent_mito < 10)
  VlnPlot(LI_merged, features = c("nFeature_Spatial","nCount_Spatial","percent_mito"), pt.size = 0, group.by = "orig.ident", cols = hughie_color)
  ggsave(paste0(LI_merged_prefix,"All_merged_QC_after_",Sys.Date(),".pdf"), width = 12, height = 8)
}


########################## 1. Integrate LI CD & HFiD############################
library(harmony)

mWAT_SAMcolors <- c("#A6CEE3", "#66A5CC","#267CB6",
                    "#A2D48E","#47A93A","#20854EFF")

if(file.exists(paste0(LI_merged_prefix, "Harmony_LI_merged.rds"))){
  LI_merged <- readRDS(paste0(LI_merged_prefix, "Harmony_LI_merged_final.rds"))
} else {
  LI_merged <- merge(CD_LI_HFD0w, c(CD_LI_HFD4w, CD_LI_HFD8w, HfiD_LI_HFD0w,HfiD_LI_HFD4w,HfiD_LI_HFD8w))
  LI_merged <- PercentageFeatureSet(LI_merged, "^mt-", col.name = "percent_mito")
  LI_merged <- subset(LI_merged, subset = nFeature_Spatial > 100 & nCount_Spatial > 100 &percent_mito < 10)
  DefaultAssay(LI_merged) <- "SCT"
  names(LI_merged@images) <- c("CD_LI_HFD0w","CD_LI_HFD4w","CD_LI_HFD8w",
                               "HfiD_LI_HFD0w","HfiD_LI_HFD4w","HfiD_LI_HFD8w")
  VariableFeatures(LI_merged) <- c(VariableFeatures(CD_LI_HFD0w), 
                                   VariableFeatures(CD_LI_HFD4w), 
                                   VariableFeatures(CD_LI_HFD8w),
                                   VariableFeatures(HfiD_LI_HFD0w),
                                   VariableFeatures(HfiD_LI_HFD4w),
                                   VariableFeatures(HfiD_LI_HFD8w))
  rm(CD_LI_HFD0w,CD_LI_HFD4w,CD_LI_HFD8w,HfiD_LI_HFD0w,HfiD_LI_HFD4w,HfiD_LI_HFD8w)
  gc()

  my_cols <- c('0' = hughie_color[1], '1' = hughie_color[2], '2' = hughie_color[3],
               '3' = hughie_color[4], '4' = hughie_color[5], '5' = hughie_color[6],
               '6' = hughie_color[7], '7' = hughie_color[8], '8' = hughie_color[9],
               '9' = hughie_color[10], '10' = hughie_color[11], '11' = hughie_color[12],
               '12' = hughie_color[13], '13' = hughie_color[14], 
               '14' = hughie_color[15], '15' = hughie_color[16],'16' = hughie_color[17],
               '17' = hughie_color[18],'18' = hughie_color[19],'19' = hughie_color[20],
               '20' = hughie_color[21],'21' = hughie_color[22])
  
  LI_merged <- ScaleData(LI_merged)
  LI_merged <- RunPCA(LI_merged, verbose = T)
  LI_merged <- RunUMAP(LI_merged, dims = 1:30)
  p1 <- DimPlot(LI_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = TRUE, repel = T, cols = hughie_color)
  p2 <- SpatialDimPlot(LI_merged, label = TRUE, label.size = 3, repel = T,ncol = 3,pt.size.factor = 75, cols = my_cols)
  pdf(paste0(LI_merged_prefix,"Before_harmony_UMAP-Dim+SpatialMaps_",Sys.Date(),".pdf"), width = 24, height = 8)
  p1 + p2
  dev.off()
  
  LI_merged <- RunHarmony(LI_merged, group.by.vars = "orig.ident")
  LI_merged <- RunUMAP(LI_merged, reduction = "harmony", dims = 1:30)
  LI_merged <- FindNeighbors(LI_merged, reduction = "harmony", dims = 1:30)
  LI_merged <- FindClusters(LI_merged, resolution = seq(from=0, by=0.1, length=10))
  
  
  clustree(LI_merged)
  ggsave(paste0(LI_merged_prefix,"Harmony_clustree_",Sys.Date(),".pdf"), width = 12, height = 10)



  
  #detemine of res, select res  
  Idents(LI_merged) <- LI_merged$SCT_snn_res.0.9
  pdf(paste0(LI_merged_prefix,"Harmony_UMAP_res0.9_SpatialMaps_",Sys.Date(),".pdf"), width = 12, height = 8)
  SpatialDimPlot(LI_merged, label = TRUE, label.size = 3, repel = T,ncol = 3,pt.size.factor = 75, cols = my_cols)
  dev.off()


  Idents(LI_merged) <- LI_merged$SCT_snn_res.0.6

  saveRDS(LI_merged, file = paste0(LI_merged_prefix,"Harmony_LI_merged.rds"))
  LI_merged <- readRDS(file = paste0(LI_merged_prefix,"Harmony_LI_merged.rds"))

  #Visualization of the clusters landscape
  pdf(paste0(LI_merged_prefix,"Harmony_UMAP-res0.6-ByClusters_",Sys.Date(), ".pdf"), width = 16, height = 8)
  DimPlot(LI_merged, reduction = "umap", group.by = c("ident", "orig.ident"), label = T, repel = T, cols = hughie_color)
  dev.off()
  #ggsave(paste0(LI_merged_prefix,"UMAP-res0.2-ByClusters_",Sys.Date(), ".pdf"), width = 16, height = 8)
  
  DimPlot(LI_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = T, repel = T, cols = hughie_color)
  ggsave(paste0(LI_merged_prefix,"Harmony_UMAP-res0.6-Byorig.ident_",Sys.Date(), ".pdf"), width = 14, height = 8)
  
  p1 <- DimPlot(LI_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = TRUE, repel = T, cols = hughie_color)
  p2 <- SpatialDimPlot(LI_merged, label = TRUE, label.size = 3, repel = T,ncol = 3,pt.size.factor = 75, cols = my_cols)
  pdf(paste0(LI_merged_prefix,"Harmony_UMAP-res0.6-Dim+SpatialMaps_",Sys.Date(),".pdf"), width = 24, height = 8)
  p1 + p2
  dev.off()
  
  #Find markers
  LI_merged <- PrepSCTFindMarkers(LI_merged, assay = "SCT", verbose = TRUE)
  LI_merged_markers <- FindAllMarkers(LI_merged, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2, verbose = T)
  write.csv(LI_merged_markers, file = paste0(LI_merged_prefix,"Harmony_Res0.6_markers_",Sys.Date(), ".csv"), quote = F)
  #LI_merged_TopMarkers <- LI_merged_markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC) %>% as.data.frame()
  
  ####################Plot###################
  # VlnPlot(LI_merged, features = LI_merged_TopMarkers$gene, stack=T, flip=T) + NoLegend()
  # ggsave(paste0(LI_merged_prefix, "Harmony_LI_All_markers_Vlnplot_",Sys.Date(), ".pdf"), width = 12, height = 30)
  # 
  # DotPlot(LI_merged, features = unique(LI_merged_TopMarkers$gene)) + 
  #   theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1), axis.text.y = element_text(size=16)) 
  # ggsave(paste0(LI_merged_prefix, "Harmony_LI_All_markers_dotplot_",Sys.Date(), ".pdf"), width = 20, height = 12)
  
  # DoHeatmap(LI_merged, features = LI_merged_markers$gene)
  # ggsave(paste0(LI_merged_prefix, "Harmony_LI_All_markers_heatmap_",Sys.Date(), ".pdf"), width = 15, height = 15)
  
  Plot_Cell_compoistion(Seurat_Obj = LI_merged, OutPrefix = LI_merged_prefix, ColorUse1 = hughie_color, ColorUse2 = mWAT_SAMcolors)
  DEG_enrichment(Seruat_DEG_file = paste0(LI_merged_prefix,"Harmony_Res0.6_markers_",Sys.Date(), ".csv"), showCategoryNum = 15, filterLevel = 4)
  #PlotObjMetrices(Seurat_Obj = LI_merged, OutPrefix= paste0(LI_merged_prefix,"Plot.obj_"), ColorUse = mWAT_SAMcolors)
  
  saveRDS(LI_merged, file = paste0(LI_merged_prefix,"Harmony_LI_merged_final.rds"))
  #LI_merged <- readRDS(paste0(LI_merged_prefix,"Harmony_LI_merged_final.rds"))
}


################################ 2. Annotation##################################
LI_merged_annotation_prefix <- "2.LI_annotation/"
# library("SCINA")
# exp_LI_merged <- as.matrix(LI_merged@assays$SCT@counts)
# sig_raw <- read.csv("CellMarker_from_ref.csv",header = T)
# sig_raw$Tissue <- gsub(" ", "_", sig_raw$Tissue)
# sig_raw$Cell.Type <- gsub(" ", "_", sig_raw$Cell.Type)
# sig_raw$Tissue.Cell.type <- paste0(sig_raw$Tissue,"_", sig_raw$Cell.Type)
# 
# # for (i in 1:nrow(sig_raw)) {
# #   Cell.Marker <- strsplit(sig_raw$Cell.Marker[i],", ")
# #   sig_LI_merged <- Cell.Marker
# #   sig_LI_merged <- list(sig_raw$Cell.Marker[1], sig_raw$Cell.Marker[2],sig_raw$Cell.Marker[3])
# # }
# 
# 
# sig_LI_merged <- strsplit(sig_raw$Cell.Marker,", ")
# names(sig_LI_merged) <- sig_raw$Tissue.Cell.type
# for (i in 1:nrow(sig_raw)) {
#   if(length(sig_LI_merged[[i]]) >49){
#     sig_LI_merged[[i]] <- sample(sig_LI_merged[[i]], 50, replace = F)
#   }
# }
# results_LI_merged = SCINA(exp_LI_merged, sig_LI_merged, max_iter = 120, convergence_n = 12, 
#                           convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap = FALSE)
# LI_merged$Cell_annotation_ref <- results_LI_merged$cell_labels
# auto_anno <- table(LI_merged$SCT_snn_res.0.7, LI_merged$Cell_annotation_ref)
# write.csv(auto_anno, paste0(LI_merged_annotation_prefix, "Auto_annotation_res0.6_ref_", Sys.Date(), ".csv"))

############Manual annotation############

JEC_markers <- c("Fabp2", #Enterocyte 2
                 "Zg16", "Clca1", "Tff3", "Spink4", #Goblet cells 0,1,5,8,9,11,14
                 "Lyz2", "Fkbp1a", #PCs & PLCs  7,10,12,14,15
                 "Cpe", #Enteroendocrine cells 4,7,10,12,14,15
                 #Progenitor cells
                 "Ybx1", #Transient-amplifying cells
                 "Rgmb" #Stem cells 10
                 #Tuft cells
                 #secretory progenitors
                 #enterocyte progenitors
                 )

NM_markers <- c("Slc12a2", #ISC 3
                 #Enterocyte progenitors
                 "Fabp1", "Aldob","Apoa4", #Enterocyte  
                 "Tff3","Muc2", #Goblet Progenitor
                 "Agr2","Klf4"#Goblet cell
                 )

JEM_markers <- c("Atp5a1","Selenop", #Enterocyte  2,6,    4,15 Progenitor?
                "Slc12a2",#Stem cell
                "Agr2","Dmbt1", #Goblet cell 1,3,6,11,12,14,15
                 #Paneth cell
                "Tac1","Cpe"#EEC
                # Tuft cell
)

MB_markers <- c("Fabp2","Fabp1","S100g", # Enterocyte  2
                "Guca2a",#BEST4_epithelial
                "Sst",#Enteroendocrine cell 10
                "Tff3","Clca1"#Goblet cell 0,1,5,8,11,14
)


LI_markers <- c("Atp5a1","Fabp1","Pigr",            #Enterocyte 1,2,6   Progenitor 15 
                "Ybx1",                             #Transient-amplifying cells 3 (Stem cell) 
                "Tff3","Clca1",                     #Goblet cell 0,5,8,9,11
                "Lyz2",                             #PCs  12
                "Cpe",                              #Enteroendocrine cells 14
                "Acta2","Des","Mylk","Myl9"         #Myofibroblast 4,7,10,13
)

LI_markers_final <- c(
  "Mapk3","Itgb1",                    #Progenitor 15 
  "Atp5a1","Fabp1","Pigr",            #Enterocyte 1,2,6   
  #"Reg3b","Lypd8",                    #Antimicrobial Program, bottom
  #"Apoa4",                            #Apolipoproteins Cholesterol
  #"Atp2a3",                          #Metal ion
  #"Slc5a1",                           #Carbohydrates
  "Ybx1",                             #Transient-amplifying cells 3 (Stem cell) 
  "Tff3","Clca1",                     #Goblet cell 0,5,8,9,11
  "Lyz2",                             #PCs  12
  "Cpe",                              #Enteroendocrine cells 14
  "Acta2","Des","Mylk","Myl9"         #Myofibroblast 4,7,10,13
)

LI_markers_final <- c(
  "Dmbt1","Rpsa",                          #Enterocyte Progenitor    1
  "Apoa1","Apoa4",                #Enterocyte               0,2&9,6
  "Clca1","Klf4",                          #Goblet Cell              5,8
  "Muc2","Spink4",           #Plasma Goblet Cell       3,11  
  "Lyz1",                                  #Paneth Cell              12  
  "Neurod1",                               #Enteroendocrine Cell     15
  "Des"              #Myofibroblasts           4,7,10,13,14
)

p1 <- VlnPlot(LI_merged, features = LI_markers_final, stack=T, flip=T)
#p2 <- VlnPlot(LI_merged, features = LI_markers, stack=T, flip=T, split.by = "orig.ident", cols = hughie_color)
p1
ggsave(paste0(LI_merged_annotation_prefix, "SI_markers_VlnPlot_",Sys.Date(), ".pdf"), width = 20, height = 12)

FeatureHeatmapPlot(Seurat_Obj = LI_merged,OutPrefix = LI_merged_annotation_prefix,
                   features = LI_markers_final,cols = hughie_color[1:16], group.by = "SCT_snn_res.0.6")


LI_merged <- RenameIdents(LI_merged, 
                          '0' = "Goblet Cell-3",
                          '1' = "Enterocyte Progenitor",
                          '2' = "Enterocyte-2",
                          '3' = "Plasma Goblet Cell",
                          '4' = "Myofibroblast",
                          '5' = "Goblet Cell-2",
                          '6' = "Enterocyte-1",
                          '7' = "Myofibroblast",
                          '8' = "Goblet Cell-1",
                          '9' = "Enterocyte-2",
                          '10' = "Myofibroblast",
                          '11' = "Goblet Cell-3",
                          '12' = "Paneth Cell",
                          '13' = "Myofibroblast",
                          '14' = "Myofibroblast",
                          '15' = "Enteroendocrine Cell")

# "Enterocyte-2" = "Ribosome-Enriched Enterocyte",
# "Enterocyte-1" = "Absorptive Enterocyte",
# "Goblet Cell-3" = "Mitochondrial Goblet Cell",
# "Goblet Cell-1" = "Lipid Metabolic Goblet Cell",
# "Goblet Cell-2" = "Barrier Protective Goblet Cell"

LI_merged$My_annotation <- Idents(LI_merged)
LI_merged$My_annotation <- factor(LI_merged$My_annotation, levels = c("Enterocyte Progenitor",
                                                                      "Enterocyte-1",
                                                                      "Enterocyte-2",
                                                                      "Goblet Cell-3",
                                                                      "Goblet Cell-1",
                                                                      "Goblet Cell-2",
                                                                      "Plasma Goblet Cell",
                                                                      "Paneth Cell",
                                                                      "Enteroendocrine Cell",
                                                                      "Myofibroblast"))
Idents(LI_merged) <- LI_merged$My_annotation
#saveRDS(LI_merged, file = paste0(LI_merged_annotation_prefix,"Harmony_LI_merged_final_annotated.rds"))


pdf(paste0(LI_merged_annotation_prefix,"Harmony_UMAP-res0.6-ByClusters_annotated_",Sys.Date(),".pdf"), width = 18, height = 8)
p1 <- DimPlot(LI_merged, reduction = "umap", group.by = c("ident"), label = T, repel = T, cols = hughie_color)
p2 <- DimPlot(LI_merged, reduction = "umap", group.by = c("orig.ident"), label = F, repel = T, cols = mWAT_SAMcolors)
p1 + p2
dev.off()


p3 <- DimPlot(LI_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = T, repel = T, cols = hughie_color)
pdf(paste0(LI_merged_annotation_prefix,"Harmony_UMAP-res0.6-Byorig.ident_annotated_",Sys.Date(),".pdf"), width = 14, height = 8)
p3
dev.off()

my_cols_annotated <- c(
  "Enterocyte Progenitor"= hughie_color[1],
  "Enterocyte-1"= hughie_color[2],
  "Enterocyte-2"= hughie_color[3],
  "Goblet Cell-3"= hughie_color[4],
  "Goblet Cell-1"= hughie_color[5],
  "Goblet Cell-2"= hughie_color[6],
  "Plasma Goblet Cell"= hughie_color[7],
  "Paneth Cell"= hughie_color[8],
  "Enteroendocrine Cell"= hughie_color[9],
  "Myofibroblast"= hughie_color[10]
  )
  

p1 <- DimPlot(LI_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = TRUE, repel = T, cols = hughie_color)
p2 <- SpatialDimPlot(LI_merged, label = F, label.size = 3, repel = T, ncol = 3, pt.size.factor = 75,cols = my_cols_annotated)
pdf(paste0(LI_merged_annotation_prefix,"Harmony_UMAP-res0.6-Dim+SpatialMaps_annotated_",Sys.Date(),".pdf"), width = 24, height = 8)
p1 + p2
dev.off()

pdf(paste0(LI_merged_annotation_prefix,"Harmony_UMAP_res0.6_SpatialMaps_annotated_",Sys.Date(),".pdf"), width = 20, height = 8)
p2
dev.off()

##################Find markers#####################
LI_merged <- PrepSCTFindMarkers(LI_merged, assay = "SCT", verbose = TRUE)
LI_merged_markers <- FindAllMarkers(LI_merged, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2, verbose = T)
write.csv(LI_merged_markers, file = paste0(LI_merged_annotation_prefix,"Harmony_Res0.6_markers_annotated_",Sys.Date(), ".csv"), quote = F)
LI_merged_TopMarkers <- LI_merged_markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC) %>% as.data.frame()

####################Plot###################
VlnPlot(LI_merged, features = LI_merged_TopMarkers$gene, stack=T, flip=T) + NoLegend()
ggsave(paste0(LI_merged_annotation_prefix, "Harmony_LI_All_markers_Vlnplot_annotated_",Sys.Date(), ".pdf"), width = 12, height = 30)

DotPlot(LI_merged, features = unique(LI_merged_TopMarkers$gene)) + 
  theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1), axis.text.y = element_text(size=16)) 
ggsave(paste0(LI_merged_annotation_prefix, "Harmony_LI_All_markers_dotplot_annotated_",Sys.Date(), ".pdf"), width = 20, height = 12)

DoHeatmap(LI_merged, features = LI_merged_markers$gene)
ggsave(paste0(LI_merged_annotation_prefix, "Harmony_LI_All_markers_heatmap_annotated_",Sys.Date(), ".pdf"), width = 15, height = 15)

Plot_Cell_compoistion(Seurat_Obj = LI_merged, OutPrefix = LI_merged_annotation_prefix, ColorUse1 = hughie_color, ColorUse2 = mWAT_SAMcolors)
DEG_enrichment(Seruat_DEG_file = paste0(LI_merged_annotation_prefix,"Harmony_Res0.6_markers_annotated_",Sys.Date(), ".csv"), showCategoryNum = 15, filterLevel = 4)

saveRDS(LI_merged, file = paste0(LI_merged_annotation_prefix,"Harmony_LI_merged_final_annotated.rds"))
#LI_merged <- readRDS(file = paste0(LI_merged_annotation_prefix,"Harmony_LI_merged_final_annotated.rds"))

#######################Final Markers###################

LI_markers_final <- c(
  "Dmbt1","Rpsa",                          #Enterocyte Progenitor    1
  "Apoa1","Apoa4",                 #Enterocyte               0&11,2&9,6
  "Clca1",                          #Goblet Cell              5,8
  "Muc2",           #Plasma Goblet Cell       3  
  "Lyz1",                                  #Paneth Cell              12  
  "Neurod1",                               #Enteroendocrine Cell     15
  "Acta2","Des","Mylk","Myl9"              #Myofibroblasts           4,7,10,13,14
)

p1 <- VlnPlot(LI_merged, features = LI_markers_final, stack=T, flip=T)
p1
ggsave(paste0(LI_merged_annotation_prefix, "Harmony_LI_markers_VlnPlot_annotated_",Sys.Date(), ".pdf"), width = 12, height = 12)

FeatureHeatmapPlot(Seurat_Obj = LI_merged,OutPrefix = LI_merged_annotation_prefix,
                   features = LI_markers_final,cols = my_cols_annotated, group.by = "My_annotation")


DotPlot(LI_merged, features = unique(LI_markers_final)) + 
  theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1), axis.text.y = element_text(size=16)) 
ggsave(paste0(LI_merged_annotation_prefix, "Haromony_LI_markers_dotplot_",Sys.Date(), ".pdf"), width = 12, height = 8)


SpatialFeaturePlot(subset(LI_merged, orig.ident %in% "HfiD_LI_HFD8w"), 
                   features = c("Alpi","Apoa1","Apoa4"), 
                   ncol = 3, 
                   alpha = c(0.1, 1), 
                   crop = T,
                   pt.size.factor = 75)
ggsave(paste0(LI_merged_annotation_prefix,"Harmony_LI_Enterocyte_markers_SpatialFeatures_",Sys.Date(), ".pdf"), width = 10, height = 5)


#saveRDS(LI_merged, file = paste0(LI_merged_annotation_prefix,"Harmony_LI_merged_final_annotated.rds"))
LI_merged <- readRDS(file = paste0(LI_merged_annotation_prefix,"Harmony_LI_merged_final_annotated.rds"))


#######################  3. Enterocytes re-annotation ##########################
LI_merged_Enterocyte_prefix <- "3.Enterocyte_annotation/"

LI_merged_Enterocyte <- subset(LI_merged, idents = c("Enterocyte-1","Enterocyte-2","Goblet Cell-3",
                                                     "Goblet Cell-1",
                                                     "Goblet Cell-2"))
Idents(LI_merged_Enterocyte) <- factor(Idents(LI_merged_Enterocyte), levels = c("Enterocyte-2",
                                                                                "Enterocyte-1", 
                                                                                "Goblet Cell-3",
                                                                                "Goblet Cell-1",
                                                                                "Goblet Cell-2"))
LI_merged_Enterocyte$My_annotation <- Idents(LI_merged_Enterocyte)


Enterocyte_markers <- c("Reg3b", "Reg3g","Ccl25",#Antimicrobial Program, bottom
                        #absorb middle
                        "Apoa4","Apoa1",              #Apolipoproteins Cholesterol
                        "Atp2a3",                                          #Metal ion
                        "Slc5a1",                             #Carbohydrates
                        "Slc7a8",                             #Amino acids
                        "Slc28a2",            #shedding top Purine metabolism
                        "Cldn7","Cldn3",       #Tight junction
                        "Actn4","Tmsb4x",              #Actin filament 
                        "Ehhadh","Lpin2"               #Fat metabolic process
                        )
Enterocyte_markers <- c("Rpl41","Rps29","Rps20",      #Translation
                        "Apoa4","Apoa1",              #Apolipoproteins Cholesterol
                        "Atp2a3",                     #Metal ion
                        "Slc7a8",                     #Amino acids
                        "Slc28a2",                    #shedding top Purine metabolism
                        "Atp1a1",            #ATP metabolism
                        "Reg3b",               #Antimicro
                        #"Actn4","Tmsb4x",              #Actin filament 
                        "Ehhadh","Lpin2",               #Fat metabolic process
                        "Pof1b","Nedd4l"                #Cell junction 
)

p1 <- VlnPlot(LI_merged_Enterocyte, features = Enterocyte_markers, stack=T, flip=T)
p1
ggsave(paste0(LI_merged_Enterocyte_prefix, "Harmony_LI_markers_VlnPlot_annotated_",Sys.Date(), ".pdf"), width = 12, height = 12)

FeatureHeatmapPlot(Seurat_Obj = LI_merged_Enterocyte,OutPrefix = LI_merged_Enterocyte_prefix,
                   features = Enterocyte_markers,cols = hughie_color[1:5], group.by = "My_annotation")

DotPlot(LI_merged_Enterocyte, features = Enterocyte_markers) + 
  theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1), axis.text.y = element_text(size=16)) 
ggsave(paste0(LI_merged_Enterocyte_prefix, "Harmony_LI_All_markers_dotplot_annotated_",Sys.Date(), ".pdf"), width = 20, height = 12)

DoHeatmap(LI_merged_Enterocyte, features = Enterocyte_markers)
ggsave(paste0(LI_merged_Enterocyte_prefix, "Harmony_LI_All_markers_heatmap_annotated_",Sys.Date(), ".pdf"), width = 15, height = 15)


###re-name enterocyte
#LI_merged <- readRDS(file = paste0(LI_merged_prefix,"Harmony_LI_merged_final_annotated.rds"))
LI_merged <- RenameIdents(LI_merged, 
                          "Enterocyte-2" = "Ribosome-Enriched Enterocyte",
                          "Enterocyte-1" = "Absorptive Enterocyte",
                          "Goblet Cell-3" = "Mitochondrial Goblet Cell",
                          "Goblet Cell-1" = "Lipid Metabolic Goblet Cell",
                          "Goblet Cell-2" = "Barrier Protective Goblet Cell"
                          )

LI_merged$My_annotation <- Idents(LI_merged)
LI_merged$My_annotation <- factor(LI_merged$My_annotation, levels = c(
  "Enterocyte Progenitor",
  "Absorptive Enterocyte",  
  "Ribosome-Enriched Enterocyte",
  "Mitochondrial Goblet Cell",
  "Lipid Metabolic Goblet Cell",
  "Barrier Protective Goblet Cell",
  "Plasma Goblet Cell",
  "Paneth Cell",
  "Enteroendocrine Cell",
  "Myofibroblast"))

Idents(LI_merged) <- LI_merged$My_annotation
# saveRDS(LI_merged, file = paste0(LI_merged_prefix,"Harmony_LI_merged_final_final_annotated.rds"))

########re-plot after re-name
pdf(paste0(LI_merged_prefix,"Harmony_UMAP-res0.6-ByClusters_annotated_final_",Sys.Date(),".pdf"), width = 18, height = 8)
p1 <- DimPlot(LI_merged, reduction = "umap", group.by = c("ident"), label = T, repel = T, cols = hughie_color)
p2 <- DimPlot(LI_merged, reduction = "umap", group.by = c("orig.ident"), label = T, repel = T, cols = mWAT_SAMcolors)
p1 + p2
dev.off()

pdf(paste0(LI_merged_prefix,"Harmony_UMAP-res0.6-Byorig.ident_annotated_final_",Sys.Date(),".pdf"), width = 14, height = 8)
DimPlot(LI_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = T, repel = T, cols = hughie_color)
dev.off()

my_cols_annotated <- c(
  "Enterocyte Progenitor"= hughie_color[1],
  "Absorptive Enterocyte"= hughie_color[2],
  "Ribosome-Enriched Enterocyte"= hughie_color[3],
  "Mitochondrial Goblet Cell"=hughie_color[4],
  "Lipid Metabolic Goblet Cell"= hughie_color[5],
  "Barrier Protective Goblet Cell"= hughie_color[6],
  "Plasma Goblet Cell"= hughie_color[7],
  "Paneth Cell"= hughie_color[8],
  "Enteroendocrine Cell"= hughie_color[9],
  "Myofibroblast"= hughie_color[10]
)

p1 <- DimPlot(LI_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = TRUE, repel = T, cols = hughie_color)
p2 <- SpatialDimPlot(LI_merged, label = F, label.size = 3, repel = T, ncol = 3, pt.size.factor = 75,cols = my_cols_annotated)
pdf(paste0(LI_merged_prefix,"Harmony_UMAP-res0.6-Dim+SpatialMaps_annotated_final_",Sys.Date(),".pdf"), width = 24, height = 8)
p1 + p2
dev.off()

pdf(paste0(LI_merged_prefix,"Harmony_UMAP_res0.6_SpatialMaps_annotated_final_",Sys.Date(),".pdf"), width = 20, height = 8)
p2
dev.off()

#######Find markers
LI_merged <- PrepSCTFindMarkers(LI_merged, assay = "SCT", verbose = TRUE)
LI_merged_markers <- FindAllMarkers(LI_merged, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2, verbose = T)
write.csv(LI_merged_markers, file = paste0(LI_merged_prefix,"Harmony_Res0.6_markers_annotated_final_",Sys.Date(), ".csv"), quote = F)
LI_merged_TopMarkers <- LI_merged_markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC) %>% as.data.frame()

#######Plot
VlnPlot(LI_merged, features = LI_merged_TopMarkers$gene, stack=T, flip=T) + NoLegend()
ggsave(paste0(LI_merged_prefix, "Harmony_LI_All_markers_Vlnplot_annotated_final_",Sys.Date(), ".pdf"), width = 12, height = 30)

DotPlot(LI_merged, features = unique(LI_merged_TopMarkers$gene)) + 
  theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1), axis.text.y = element_text(size=16)) 
ggsave(paste0(LI_merged_prefix, "Harmony_LI_All_markers_dotplot_annotated_final_",Sys.Date(), ".pdf"), width = 20, height = 12)

DoHeatmap(LI_merged, features = LI_merged_markers$gene)
ggsave(paste0(LI_merged_prefix, "Harmony_LI_All_markers_heatmap_annotated_final_",Sys.Date(), ".pdf"), width = 15, height = 15)

Plot_Cell_compoistion(Seurat_Obj = LI_merged, OutPrefix = LI_merged_prefix, ColorUse1 = hughie_color, ColorUse2 = mWAT_SAMcolors)
Plot_Cell_compoistion(Seurat_Obj = subset(LI_merged, idents = c("Ribosome-Enriched Enterocyte",
                                                                "Absorptive Enterocyte",
                                                                "Mitochondrial Goblet Cell",
                                                                "Lipid Metabolic Goblet Cell",
                                                                "Barrier Protective Goblet Cell",
                                                                "Plasma Goblet Cell")), 
                      OutPrefix = paste0(LI_merged_prefix,"Diff_"), 
                      ColorUse1 = c(hughie_color[2],hughie_color[3],hughie_color[4],hughie_color[5],hughie_color[6],hughie_color[7]), 
                      ColorUse2 = mWAT_SAMcolors)

DEG_enrichment(Seruat_DEG_file = paste0(LI_merged_prefix,"Harmony_Res0.6_markers_annotated_final_",Sys.Date(), ".csv"), showCategoryNum = 15, filterLevel = 4)

# First Batch
# LI_markers_final <- c(
#   "Mapk3","Itgb1",                    #Progenitor 15 
#   "Atp5a1","Fabp1","Pigr",            #Enterocyte 1,2,6   
#   #"Reg3b","Lypd8",                    #Antimicrobial Program, bottom
#   #"Apoa4",                            #Apolipoproteins Cholesterol
#   #"Atp2a3",                          #Metal ion
#   #"Slc5a1",                           #Carbohydrates
#   "Ybx1",                             #Transient-amplifying cells 3 (Stem cell) 
#   "Tff3","Clca1",                     #Goblet cell 0,5,8,9,11
#   "Lyz2",                             #PCs  12
#   "Cpe",                              #Enteroendocrine cells 14
#   "Acta2","Des","Mylk","Myl9"         #Myofibroblast 4,7,10,13
# )

LI_markers_final <- c(
  "Dll1",                         #Enterocyte Progenitor    1
  "Apoa1","Apoa4",                #Enterocyte               0,2&9,6
  "Clca1","Klf4",                          #Goblet Cell              5,8
  "Muc2","Spink4",           #Plasma Goblet Cell       3,11  
  "Lyz1",                                  #Paneth Cell              12  
  "Neurod1",                               #Enteroendocrine Cell     15
  "Des"              #Myofibroblasts           4,7,10,13,14
)

FeatureHeatmapPlot(Seurat_Obj = LI_merged,OutPrefix = LI_merged_prefix,
                   features = LI_markers_final,cols = hughie_color[1:10], group.by = "My_annotation")


#re plot Enterocytes marker after re-name
LI_merged_Enterocyte <- subset(LI_merged, idents = c("Absorptive Enterocyte",
                                                     "Ribosome-Enriched Enterocyte",
                                                     
                                                     "Mitochondrial Goblet Cell",
                                                     "Lipid Metabolic Goblet Cell",
                                                     "Barrier Protective Goblet Cell",
                                                     "Plasma Goblet Cell"))

Enterocyte_markers <- c("Apoa4","Apoa1",              #Apolipoproteins Cholesterol
                        "Atp2a3",                     #Metal ion
                        "Slc7a8",                     #Amino acids
                        "Slc28a2",                    #shedding top Purine metabolism
                        "Rps20","Rpl41","Rps29",      #Translation
                        
                        "Atp1a1",            #ATP metabolism
                        #"Reg3b",               #Antimicro
                        #"Actn4","Tmsb4x",              #Actin filament 
                        "Acsl5","Acss1",               #Fat metabolic process
                        "Pof1b","Nedd4l",                #Cell junction
                        "Igkc","Jchain"
)

FeatureHeatmapPlot(Seurat_Obj = LI_merged_Enterocyte,OutPrefix = LI_merged_Enterocyte_prefix,
                   features = Enterocyte_markers,cols = hughie_color[1:6], group.by = "My_annotation")


 Selected_enrichment_Plot(GO_enrichment_file = paste0(LI_merged_Enterocyte_prefix,"Harmony_Res0.6_markers_annotated_final_2024-11-13_Selected_GO-BP_list_2024-11-13.csv"),
                          OutPrefix = LI_merged_Enterocyte_prefix,
                          cluster_order = c("Absorptive Enterocyte",
                                            "Ribosome-Enriched Enterocyte",
                                            
                                            "Mitochondrial Goblet Cell",
                                            "Lipid Metabolic Goblet Cell",
                                            "Barrier Protective Goblet Cell",
                                            "Plasma Goblet Cell"),
                          levels = c("intestinal cholesterol absorption",
                                     "intestinal lipid absorption",
                                     "intestinal absorption",
                                     "cytoplasmic translation",
                                     
                                     "cellular respiration",
                                     "ATP biosynthetic process",
                                     "mitochondrial respiratory chain complex assembly",
                                     "lipid catabolic process",
                                     "fatty acid metabolic process",
                                     "fatty acid catabolic process",
                                     "lipid oxidation",
                                     "cell-cell junction assembly",
                                     "bicellular tight junction assembly",
                                     "tight junction assembly",
                                     "protein glycosylation",
                                     "protein folding"))

 saveRDS(LI_merged, file = paste0(LI_merged_prefix,"Harmony_LI_merged_final_final_annotated.rds"))
 #LI_merged <- readRDS(file = paste0(LI_merged_prefix,"Harmony_LI_merged_final_final_annotated.rds"))
 

############################### 4. Scvelo#######################################
LI_merged <- readRDS(file = paste0(LI_merged_prefix,"Harmony_LI_merged_final_final_annotated.rds"))

LI_velocyto_prefix <- "4.Scvelo_velocity/"
########################Prepare emat and nmat############
ldat.CD.LI.0w <- read.loom.matrices.spatial(file = "../PZ_analysis/raw/0w/sample_CD_LI_bin50_seurat.loom")
emat.CD.LI.0w <- ldat.CD.LI.0w$spliced
#colnames(emat.CD.LI.0w) <- gsub("mWAT_CD2:", "CD_HFD0w_", colnames(emat.CD.0w))
#colnames(emat.CD.LI.0w) <- gsub("x", "-1", colnames(emat.CD.0w))

ldat.CD.LI.4w <- read.loom.matrices.spatial(file = "../PZ_analysis/raw/4w/sample_LI_CD4W_bin50_seurat.loom")
emat.CD.LI.4w <- ldat.CD.LI.4w$spliced

ldat.CD.LI.8w <- read.loom.matrices.spatial(file = "../PZ_analysis/raw/8w/sample_LI_CD8W_bin50_seurat.loom")
emat.CD.LI.8w <- ldat.CD.LI.8w$spliced

ldat.HfiD.LI.0w <- read.loom.matrices.spatial(file = "../PZ_analysis/raw/0w/sample_HfiD_LI_bin50_seurat.loom")
emat.HfiD.LI.0w <- ldat.HfiD.LI.0w$spliced

ldat.HfiD.LI.4w <- read.loom.matrices.spatial(file = "../PZ_analysis/raw/4w/sample_LI_HfiD4W_bin50_seurat.loom")
emat.HfiD.LI.4w <- ldat.HfiD.LI.4w$spliced

ldat.HfiD.LI.8w <- read.loom.matrices.spatial(file = "../PZ_analysis/raw/8w/sample_LI_HfiD8W_bin50_seurat.loom")
emat.HfiD.LI.8w <- ldat.HfiD.LI.8w$spliced

nmat.CD.LI.0w <- ldat.CD.LI.0w$unspliced
nmat.CD.LI.4w <- ldat.CD.LI.4w$unspliced
nmat.CD.LI.8w <- ldat.CD.LI.8w$unspliced
nmat.HfiD.LI.0w <- ldat.HfiD.LI.0w$unspliced
nmat.HfiD.LI.4w <- ldat.HfiD.LI.4w$unspliced
nmat.HfiD.LI.8w <- ldat.HfiD.LI.8w$unspliced

cell.chose <- as.data.frame(LI_merged$My_annotation)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- row.names(LI_merged@assays$SCT@meta.features)
feature.chose$feature <- row.names(feature.chose)

mat.feature.chose <- intersect(row.names(emat.CD.LI.0w),row.names(emat.CD.LI.4w))
mat.feature.chose <- intersect(mat.feature.chose, row.names(emat.CD.LI.8w))
mat.feature.chose <- intersect(mat.feature.chose, feature.chose)

emat.CD.LI.0w <- emat.CD.LI.0w[rownames(emat.CD.LI.0w) %in% mat.feature.chose, ]
emat.CD.LI.4w <- emat.CD.LI.4w[rownames(emat.CD.LI.4w) %in% mat.feature.chose, ]
emat.CD.LI.8w <- emat.CD.LI.8w[rownames(emat.CD.LI.8w) %in% mat.feature.chose, ]

emat <- cbind(emat.CD.LI.0w,emat.CD.LI.4w,emat.CD.LI.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]

nmat.CD.LI.0w <- nmat.CD.LI.0w[rownames(nmat.CD.LI.0w) %in% mat.feature.chose, ]
nmat.CD.LI.4w <- nmat.CD.LI.4w[rownames(nmat.CD.LI.4w) %in% mat.feature.chose, ]
nmat.CD.LI.8w <- nmat.CD.LI.8w[rownames(nmat.CD.LI.8w) %in% mat.feature.chose, ]

nmat <- cbind(nmat.CD.LI.0w,nmat.CD.LI.4w,nmat.CD.LI.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% row.names(cell.chose)]

cell.type <- cell.chose[colnames(emat),1]

saveRDS(emat, paste0(LI_velocyto_prefix,"emat.CD.rds"))
saveRDS(nmat, paste0(LI_velocyto_prefix,"nmat.CD.rds"))

#确保细胞数和emat，nmat一致
LI_merged <- subset(LI_merged, cells = colnames(emat))
LI_merged <- subset(LI_merged, features = rownames(emat))
saveRDS(LI_merged, paste0(LI_velocyto_prefix,"LI_merged_CD_scvelo.rds"))




LI_merged <- readRDS(file = paste0(LI_merged_prefix,"Harmony_LI_merged_final_final_annotated.rds"))

cell.chose <- as.data.frame(LI_merged$My_annotation)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- row.names(LI_merged@assays$SCT@meta.features)
feature.chose$feature <- row.names(feature.chose)

mat.feature.chose <- intersect(row.names(emat.HfiD.LI.0w),row.names(emat.HfiD.LI.4w))
mat.feature.chose <- intersect(mat.feature.chose, row.names(emat.HfiD.LI.8w))
mat.feature.chose <- intersect(mat.feature.chose, feature.chose)

emat.HfiD.LI.0w <- emat.HfiD.LI.0w[rownames(emat.HfiD.LI.0w) %in% mat.feature.chose, ]
emat.HfiD.LI.4w <- emat.HfiD.LI.4w[rownames(emat.HfiD.LI.4w) %in% mat.feature.chose, ]
emat.HfiD.LI.8w <- emat.HfiD.LI.8w[rownames(emat.HfiD.LI.8w) %in% mat.feature.chose, ]

emat <- cbind(emat.HfiD.LI.0w,emat.HfiD.LI.4w,emat.HfiD.LI.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]

nmat.HfiD.LI.0w <- nmat.HfiD.LI.0w[rownames(nmat.HfiD.LI.0w) %in% mat.feature.chose, ]
nmat.HfiD.LI.4w <- nmat.HfiD.LI.4w[rownames(nmat.HfiD.LI.4w) %in% mat.feature.chose, ]
nmat.HfiD.LI.8w <- nmat.HfiD.LI.8w[rownames(nmat.HfiD.LI.8w) %in% mat.feature.chose, ]

nmat <- cbind(nmat.HfiD.LI.0w,nmat.HfiD.LI.4w,nmat.HfiD.LI.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% row.names(cell.chose)]

cell.type <- cell.chose[colnames(emat),1]

saveRDS(emat, paste0(LI_velocyto_prefix,"emat.HfiD.rds"))
saveRDS(nmat, paste0(LI_velocyto_prefix,"nmat.HfiD.rds"))

#确保细胞数和emat，nmat一致
LI_merged <- subset(LI_merged, cells = colnames(emat))
LI_merged <- subset(LI_merged, features = rownames(emat))
saveRDS(LI_merged, paste0(LI_velocyto_prefix,"LI_merged_HfiD_scvelo.rds"))


# LI_merged_Enterocyte_CD <- subset(LI_merged_Enterocyte, orig.ident %in% c("CD_LI_HFD0w", "CD_LI_HFD4w",  "CD_LI_HFD8w"))
# LI_merged_Enterocyte_HfiD <- subset(LI_merged_Enterocyte, orig.ident %in% c("HfiD_LI_HFD0w", "HfiD_LI_HFD4w", "HfiD_LI_HFD8w"))




##################### 5. Monocle################
LI_merged_prefix <- "1.LI_merged/"
LI_monocle_prefix <- "5.Monocle/New/"

LI_merged <- readRDS(file = paste0(LI_merged_prefix,"Harmony_LI_merged_final_final_annotated.rds"))
    LI_merged_monocle <- subset(LI_merged, My_annotation %in% c(
                                                                "Enterocyte Progenitor",
                                                                "Ribosome-Enriched Enterocyte",
                                                                "Absorptive Enterocyte",
                                                                "Mitochondrial Goblet Cell",
                                                                "Lipid Metabolic Goblet Cell",
                                                                "Barrier Protective Goblet Cell",
                                                                "Plasma Goblet Cell"
                                                                ))
    # LI_merged_monocle <- subset(LI_merged_monocle, orig.ident %in% c("CD_LI_HFD4w","CD_LI_HFD8w",
    #                                                                  "HfiD_LI_HFD4w","HfiD_LI_HFD8w"))

separate_timepoint(Timepoint = "HFD0w")
separate_timepoint(Timepoint = "HFD4w")
separate_timepoint(Timepoint = "HFD8w")

separate_timepoint(Timepoint = "HFD4w8w")


separate_timepoint <- function(Timepoint = Timepoint) {
  
  #CD group Enterocyte
  LI_CD_monocle_prefix <- paste0(LI_monocle_prefix,"CD_group/",Timepoint,"/EP_six_diff/")
  LI_merged_monocle_CD <- subset(LI_merged_monocle, orig.ident %in% paste0("CD_LI_",Timepoint))
  
  LI_merged_monocle_CD$My_annotation <- Idents(LI_merged_monocle_CD)
  Idents(LI_merged_monocle_CD) <- LI_merged_monocle_CD$orig.ident
  #LI_merged_monocle_CD <- subset(x=LI_merged_monocle_CD, downsample=1000)
  Idents(LI_merged_monocle_CD) <- LI_merged_monocle_CD$My_annotation
  
  LI_merged_monocle_CD <- RunPCA(LI_merged_monocle_CD, verbose = T)
  LI_merged_monocle_CD <- RunUMAP(LI_merged_monocle_CD, dims = 1:3)
  
  p1 <- DimPlot(LI_merged_monocle_CD, reduction = "umap", cols = c(
    hughie_color[1],
    hughie_color[2],
    hughie_color[3],
    hughie_color[4],
    hughie_color[5],
    hughie_color[6],
    hughie_color[7]
  ))
  p1
  ggsave(paste0(LI_CD_monocle_prefix,"Monocle_UMAP-Bycluster_",Sys.Date(), ".pdf"), width = 14, height = 8)
  
  data_CD <- as(as.matrix(LI_merged_monocle_CD@assays$SCT@counts), 'sparseMatrix')
  pd_CD <- new('AnnotatedDataFrame', data = LI_merged_monocle_CD@meta.data)
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
  write.csv(diff_test_res_CD, file=paste0(LI_CD_monocle_prefix,"CD_group_diff_test_result_", Sys.Date(),".csv"),sep="")
  #diff_test_res_CD <- read.csv(file=paste0(AdMonocle_prefix,"CD_group_diff_test_result.csv",sep=""))
  
  ordering_genes_CD <- row.names(subset(diff_test_res_CD, qval < 0.01))
  ordering_genes_CD[1:5]
  myCDs_CD <- setOrderingFilter(myCDs_CD, ordering_genes_CD)
  plot_ordering_genes(myCDs_CD)
  ggsave(paste0(LI_CD_monocle_prefix,"CD_group_Ordering.genes_", Sys.Date(), ".pdf"), height = 7, width = 8)
  
  myCDs_CD <- reduceDimension(myCDs_CD, max_components = 2, method = 'DDRTree',verbose = T)
  myCDs_CD <- orderCells(myCDs_CD)
  
  saveRDS(myCDs_CD, file = paste0(LI_CD_monocle_prefix, "CD_group_myCDs_",Sys.Date(),".rds"))
  #myCDs_CD <- readRDS(paste0(LI_CD_monocle_prefix, "CD_group_myCDs.rds"))
  
  #######monocle trajectory
  data <- GetAssayData(LI_merged_monocle_CD, assay = 'SCT', slot = 'counts')
  cell_metadata <- LI_merged_monocle_CD@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  monocle3_CDs <- new_cell_data_set(data,
                                    cell_metadata = cell_metadata,
                                    gene_metadata = gene_annotation)
  #preprocess_CDs函数相当于seurat中NormalizeData ScaleData RunPCA
  monocle3_CDs <- preprocess_cds(monocle3_CDs, num_dim = 100)
  #umap降维
  monocle3_CDs <- reduce_dimension(monocle3_CDs, preprocess_method = "PCA")
  p1 <- plot_cells(monocle3_CDs, reduction_method="UMAP", color_cells_by="My_annotation"
                   #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
  )+ggtitle('CDs.umap')
  ##从seurat导入整合过的umap坐标
  monocle3_CDs.embed <- monocle3_CDs@int_colData$reducedDims$UMAP
  int.embed <- Embeddings(LI_merged_monocle_CD, reduction = "umap")
  int.embed <- int.embed[rownames(monocle3_CDs.embed),]
  monocle3_CDs@int_colData$reducedDims$UMAP <- int.embed
  p2 <- plot_cells(monocle3_CDs, reduction_method="UMAP", color_cells_by="My_annotation"
                   #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
  )+ggtitle('int.umap')
  p1 + p2
  ggsave(paste0(LI_CD_monocle_prefix,"Monocle3_umap_embed_", Sys.Date(), ".pdf"), height = 7, width = 16)
  
  ## Monocle3聚类分区
  monocle3_CDs <- cluster_cells(monocle3_CDs, resolution = 1e-3)
  p3 <- plot_cells(monocle3_CDs, show_trajectory_graph = FALSE)+ggtitle("Label by clusterID")
  p4 <- plot_cells(monocle3_CDs, color_cells_by = "partition", show_trajectory_graph = FALSE)+ggtitle("Label by partitionID")
  p3 + p4
  ggsave(paste0(LI_CD_monocle_prefix,"Monocle3_umap_cluster_", Sys.Date(), ".pdf"), height = 7, width = 16)
  
  
  ## 识别轨迹
  monocle3_CDs <- learn_graph(monocle3_CDs)
  p5 <- plot_cells(monocle3_CDs, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                   label_branch_points = FALSE)
  p5
  ggsave(paste0(LI_CD_monocle_prefix,"Monocle3_umap_trajectory_", Sys.Date(), ".pdf"), width = 5, height = 5)
  
  
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
  ggsave(paste0(LI_CD_monocle_prefix,"Monocle3_umap_trajectory_line_", Sys.Date(), ".pdf"), width = 6.5, height = 5)

  
  #HfiD group Enterocyte
  LI_HfiD_monocle_prefix <- paste0(LI_monocle_prefix,"HfiD_group/",Timepoint,"/EP_six_diff/")
  LI_merged_monocle_HfiD <- subset(LI_merged_monocle, orig.ident %in% paste0("HfiD_LI_",Timepoint))
  
  LI_merged_monocle_HfiD$My_annotation <- Idents(LI_merged_monocle_HfiD)
  Idents(LI_merged_monocle_HfiD) <- LI_merged_monocle_HfiD$orig.ident
  #LI_merged_monocle_HfiD <- subset(x=LI_merged_monocle_HfiD, downsample=1000)
  Idents(LI_merged_monocle_HfiD) <- LI_merged_monocle_HfiD$My_annotation
  
  LI_merged_monocle_HfiD <- RunPCA(LI_merged_monocle_HfiD, verbose = T)
  LI_merged_monocle_HfiD <- RunUMAP(LI_merged_monocle_HfiD, dims = 1:3)

  p1 <- DimPlot(LI_merged_monocle_HfiD, reduction = "umap", cols = c(
    hughie_color[1],
    hughie_color[2],
    hughie_color[3],
    hughie_color[4],
    hughie_color[5],
    hughie_color[6],
    hughie_color[7]
  ))
  p1 
  ggsave(paste0(LI_HfiD_monocle_prefix,"Monocle_UMAP-Bycluster_",Sys.Date(), ".pdf"), width = 14, height = 8)
  
  data_HfiD <- as(as.matrix(LI_merged_monocle_HfiD@assays$SCT@counts), 'sparseMatrix')
  pd_HfiD <- new('AnnotatedDataFrame', data = LI_merged_monocle_HfiD@meta.data)
  fData_HfiD <- data.frame(gene_short_name = row.names(data_HfiD), row.names = row.names(data_HfiD))
  fd_HfiD <- new('AnnotatedDataFrame', data = fData_HfiD)
  
  myHfiDs_HfiD <- newCellDataSet(data_HfiD,
                             phenoData = pd_HfiD,
                             featureData = fd_HfiD,
                             expressionFamily = negbinomial.size())
  
  myHfiDs_HfiD <- estimateSizeFactors(myHfiDs_HfiD)
  myHfiDs_HfiD <- estimateDispersions(myHfiDs_HfiD)                        
  myHfiDs_HfiD <- detectGenes(myHfiDs_HfiD, min_expr = 0.1)
  myHfiDs_expressed_genes_HfiD <- row.names(myHfiDs_HfiD)
  diff_test_res_HfiD <- differentialGeneTest(myHfiDs_HfiD[myHfiDs_expressed_genes_HfiD,],fullModelFormulaStr = "~My_annotation",cores = 4,verbose = T)
  
  diff_test_res_HfiD[1:4,1:4]
  write.csv(diff_test_res_HfiD, file=paste0(LI_HfiD_monocle_prefix,"HfiD_group_diff_test_result_", Sys.Date(),".csv"),sep="")
  #diff_test_res_HfiD <- read.csv(file=paste0(LI_HfiD_monocle_prefix,"HfiD_group_diff_test_result_", Sys.Date(),".csv"),sep="")
  
  ordering_genes_HfiD <- row.names(subset(diff_test_res_HfiD, qval < 0.01))
  ordering_genes_HfiD[1:5]
  myHfiDs_HfiD <- setOrderingFilter(myHfiDs_HfiD, ordering_genes_HfiD)
  plot_ordering_genes(myHfiDs_HfiD)
  ggsave(paste0(LI_HfiD_monocle_prefix,"HfiD_group_Ordering.genes_", Sys.Date(), ".pdf"), height = 7, width = 8)
  
  myHfiDs_HfiD <- reduceDimension(myHfiDs_HfiD, max_components = 2, method = 'DDRTree',verbose = T)
  myHfiDs_HfiD <- orderCells(myHfiDs_HfiD)
  
  saveRDS(myHfiDs_HfiD, file = paste0(LI_HfiD_monocle_prefix, "HfiD_group_myHfiDs_",Sys.Date(),".rds"))
  #myHfiDs_HfiD <- readRDS(paste0(LI_HfiD_monocle_prefix, "HfiD_group_myHfiDs_2024-11-13.rds"))
  
  #######monocle trajectory
  data <- GetAssayData(LI_merged_monocle_HfiD, assay = 'SCT', slot = 'counts')
  cell_metadata <- LI_merged_monocle_HfiD@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  monocle3_HfiDs <- new_cell_data_set(data,
                                    cell_metadata = cell_metadata,
                                    gene_metadata = gene_annotation)
  #preprocess_HfiDs函数相当于seurat中NormalizeData ScaleData RunPCA
  monocle3_HfiDs <- preprocess_cds(monocle3_HfiDs, num_dim = 100)
  #umap降维
  monocle3_HfiDs <- reduce_dimension(monocle3_HfiDs, preprocess_method = "PCA")
  p1 <- plot_cells(monocle3_HfiDs, reduction_method="UMAP", color_cells_by="My_annotation"
                   
  )+ggtitle('HfiDs.umap')
  p1
  
  ##从seurat导入整合过的umap坐标
  monocle3_HfiDs.embed <- monocle3_HfiDs@int_colData$reducedDims$UMAP
  int.embed <- Embeddings(LI_merged_monocle_HfiD, reduction = "umap")
  int.embed <- int.embed[rownames(monocle3_HfiDs.embed),]
  monocle3_HfiDs@int_colData$reducedDims$UMAP <- int.embed
  p2 <- plot_cells(monocle3_HfiDs, reduction_method="UMAP", color_cells_by="My_annotation"
                   #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
  )+ggtitle('int.umap')
  p1 + p2
  ggsave(paste0(LI_HfiD_monocle_prefix,"Monocle3_umap_embed_", Sys.Date(), ".pdf"), height = 7, width = 16)
  
  ## Monocle3聚类分区
  monocle3_HfiDs <- cluster_cells(monocle3_HfiDs,resolution = 1e-3)
  p3 <- plot_cells(monocle3_HfiDs, show_trajectory_graph = FALSE)+ggtitle("Label by clusterID")
  p4 <- plot_cells(monocle3_HfiDs, color_cells_by = "partition", show_trajectory_graph = FALSE)+ggtitle("Label by partitionID")
  p3 + p4
  ggsave(paste0(LI_HfiD_monocle_prefix,"Monocle3_umap_cluster_", Sys.Date(), ".pdf"), height = 7, width = 16)
  
  
  ## 识别轨迹
  monocle3_HfiDs <- learn_graph(monocle3_HfiDs)
  p5 <- plot_cells(monocle3_HfiDs, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                   label_branch_points = FALSE)
  p5
  ggsave(paste0(LI_HfiD_monocle_prefix,"Monocle3_umap_trajectory_", Sys.Date(), ".pdf"), width = 5, height = 5)
  
  
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
  ggsave(paste0(LI_HfiD_monocle_prefix,"Monocle3_umap_trajectory_line_", Sys.Date(), ".pdf"), width = 6.5, height = 5)
  
}
