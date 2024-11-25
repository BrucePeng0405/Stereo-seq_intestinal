#===============================================================================
# This script stores functions to analyze sc/snRNA-seq dataset
# Version 1.7 at 2024-04-22
# Created by Houyu Zhang
# Issue report on houyuzhang@stu.pku.edu.cn
# Copyright (c) 2023 __CarlosLab@PKU__. All rights reserved.
#===============================================================================
#Global setup
R.utils::setOption("clusterProfiler.download.method",'auto')
#General packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BiocManager))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(venn))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsankey))
suppressPackageStartupMessages(library(ggalluvial))
suppressPackageStartupMessages(library(ggh4x))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(EnhancedVolcano))
#Single cell packages
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratWrappers))
suppressPackageStartupMessages(library(SeuratDisk)) #remotes::install_github("mojaveazure/seurat-disk")
suppressPackageStartupMessages(library(sctransform))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(SingleR)) #BiocManager::install("SingleR")
suppressPackageStartupMessages(library(velocyto.R)) #devtools::install_github("velocyto-team/velocyto.R"), install boost first!
suppressPackageStartupMessages(library(monocle3)) #devtools::install_github('cole-trapnell-lab/monocle3')
suppressPackageStartupMessages(library(DoubletFinder)) #remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
suppressPackageStartupMessages(library(clustree))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(CellChat)) #devtools::install_github("sqjin/CellChat") also, circlize
suppressPackageStartupMessages(library(NMF))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(clusterProfiler)) #BiocManager::install("clusterProfiler")
suppressPackageStartupMessages(library(ReactomePA)) #BiocManager::install("ReactomePA")
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(org.Mm.eg.db)) #BiocManager::install("org.Mm.eg.db")
suppressPackageStartupMessages(library(Orthology.eg.db))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(org.Rn.eg.db))
suppressPackageStartupMessages(library(SCENIC)) #devtools::install_github("aertslab/SCENIC")
suppressPackageStartupMessages(library(cerebroApp))
suppressPackageStartupMessages(library(UCell))
suppressPackageStartupMessages(library(wpa))
suppressPackageStartupMessages(library(Scillus))
suppressPackageStartupMessages(library(scCustomize))
suppressPackageStartupMessages(library(scRNAtoolVis)) # devtools::install_github('junjunlab/scRNAtoolVis')
suppressPackageStartupMessages(library(monocle))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(loomR)) # Install devtools from CRAN install.packages("devtools") # Use devtools to install hdf5r and loomR from GitHub devtools::install_github(repo = "hhoeflin/hdf5r") devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

#===============================================================================
# Download some utilities 
#===============================================================================
if(F){#SCENIC mouse motif database (# mc9nr: Motif collection version 9: 24k motifs)
  dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
               "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
  dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
  for(featherURL in dbFiles){download.file(featherURL, destfile=basename(featherURL))}
}
#===============================================================================
# Pre-set color schemes for snRNA-seq datasets
#===============================================================================
# display.brewer.all()
# brewer.pal(n = 12, name = 'Paired')
hughie_color <- c("#0072B5FF","#BC3C29FF", "#E18727FF", "#20854EFF", "#6F99ADFF", "#7876B1FF","#EE4C97FF", #nemj scheme
                  "#E64B35FF", "#4DBBD5FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", #npg scheme
                  "#3B4992FF", "#631879FF", "#5F559BFF", "#A20056FF", "#808180FF", "#1B1919FF") #AAAS scheme

mWAT_colors <- c( "#AAD0AC","#EED0E0","#EFE2AA","#83B4EF", "#DBC9B3","#8ECFF8","#7EC0C3","#EBAEA9","#BFA6C9",
                    "#F5E0BA","#AED0DF","#89B780","#F5D8D0","#CB95BB","#808180FF")

mWAT_SAMcolors <- c("#A6CEE3", "#66A5CC","#267CB6",
                    "#A2D48E","#47A93A","#20854EFF")

PairedColor <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")

c("#A6CEE3", "#66A5CC", "#267CB6","#5CA3A2", "#A2D48E", "#83C767", "#47A93A", "#729E4E","#D09B82", "#F47878",
  "#E93B3C", "#E83C2D", "#F48A54", "#FDB45D","#FE9628", "#F98416", "#E09C7B", "#C4ABD2","#9774B6", "#6A3D9A")

#===============================================================================
# List for some generally used markers
#===============================================================================
PlotGenes <- c("Pdgfra","Nova1","Dcn", # ASPC
               "Plin1","Acvr1c","Ucp1","Cidea","Pank1",# Adipocyte-1
               "Cyp2e1","Adrb3","Tshr", "Lep","Adipoq",# Adipocyte-2
               "Cyyr1","Flt1","Jam2","Esam","Vwf","Cdh5","Kdr","Icam2","Pecam1", # Endothelial Cell
               "Pkhd1","Dcdc2a","Abcc9","Sgip1",# Epithelial Cell
               "Steap4","Rgs5","Abcc9", # Pericyte
               "Unc5c","Myh11","Acta2","Tbx18", # Smooth muscle Cell
               "Muc16", "Pkhd1l1","Msln","Upk3b", # Mesothelial Cell
               "Dock2","Ptprc","Mctp1","Mrc1","Adgre1",# Macrophage
               "Prkcq","Themis",# T Cell
               "Pax5","Bcl11a", # B Cell
               "Fn1","Retnla", "Plac8","Clec4e","Ccr2","Ly6c2","Ear2","Cx3xr1","Csf3r","Fcar","Sell","Fcer1a","Fcgr3a","Hes4",#Monocytes
               "Clec9a","Xcr1","Ccr7","Pclaf","Stmn1","Cd209","Sirpa", "Ear2","Irf8","Dpp4","Cadm1","Cd1c","Il7r","Lamp3","Clec10a","Cd1c", #DC
               "Btk","Tec","Mafb","Cybb","Flt3","Cpa3","Csf3r","Ucp1","Il7r",
               "Sox10","Ncam1","Egr2", #Schwann
               "Slc34a1","Acsm2",# Kidney PT Cell
               "Slc12a3","Slc14a2",# Kidney DT Cell
               "Ctnnd2","Tmprss2","Upk3a" # Urothelial Cell
)

ASPC_markes <- list()
ASPC_markes[["WAT_FAPs"]] <- c("Dpp4","Anxa3","Pi16","Cd55","Wnt2","Ly6c1","Fn1","Il33","Cd34")
ASPC_markes[["BAT_FAPs"]] <- c("Col5a3","Cxcl14","Bmper","Dpp4","Pi16","Clec11a","Gdf10","Trpv1","Vipr2","Gli1","Myh11","")
ASPC_markes[["WAT_Preadipocytes"]] <- c("Pparg","Col4a2","Icam1","Dlk1","Cd36","Fabp4")
ASPC_markes[["BAT_Preadipocytes"]] <- c("","","","","","","","","","","","") #Cd142 is F3
ASPC_markes[["WAT_Aregs"]] <- c("Clec11a","Fmo2","Cd36","Lgals3","Ccl6")

# unname(unlist(ASPC_markes))

#These markers are integrated from the cell metabolism review
White_adipocytes_markes <- list()
White_adipocytes_markes[["LGA/Adipo-PLIN/mAd1/hAd4"]] <- c("Acaca","Acly","Fasn","Dgat1","Dgat2",#Lipogenic adipocyte 
                                                           "Plin1","Plin4","Lipe", "Ebf2","Pck1", "Gria4") #hAd4 
White_adipocytes_markes[["LSA/Adipo-SAA/mAd3,4/hAd3"]] <- c("Car3","Apoe","Cd36",#Lipid-scavenging adipocyte
                                                            "Mgll","Agpat2", "Apoe","Cacna1a","Fgf1","Lep", "Cav2") #hAd3
White_adipocytes_markes[["SLSA/mAd5,6/hAd2,5"]] <- c("Apoe","Cd36","Hif1a","Ddr2", # Stressed #Lipid-scavenging adipocyte
                                                     "Acss2","Prune2","Mt2","Casp4","Bcl3", "Nav2","Fabp4","Atp1b3","Pgap1") #hAd2,5
White_adipocytes_markes[["mAd2/hAd1,7"]] <- c("Sorbs2","Crim1", "Aff3","Wdpcp","Thsd7b","Agmo")#hAd1,7

Patrickl_Seale_review <- c("Adipoq","Fabp4","Pparg","Plin1","Acvr1c", #pan markers
                           "Cebpa","Cebpb",
                           "Ucp1","Dio2","Cidea","Ppargc1a","Ppara","Cox7a1","Cox8b","Mecom","Prdm16","Ebf2", #brown and beige 
                           "Lhx8","Zic1","Mpzl2","Pdk4","Epsti1",# brown-only
                           "Tbx1","Cited1","Shox2","Tnfrsf9","Tmem26","Slc36a2","P2rx5",#beige-only, Cd137 (Tnfrsf9, Pat2 (Slc36a2
                           "Lep","Retn","Agt" #White-only Pref1 (Dlk1), Cd29 (Itgb1)
                           )
Prdm_Genes <- c("Prdm1","Prdm2","Mecom","Prdm4","Prdm5","Prdm6","Prdm7","Prdm8","Prdm9","Prdm10",
                "Prdm11","Prdm12","Prdm13","Prdm14","Prdm15","Prdm16","Prdm17")
#===============================================================================
# Pathway genelist
#===============================================================================
#The overall process of triglyceride biosynthesis consists of four biochemical pathways: fatty acyl-CoA biosynthesis, \
#conversion of fatty acyl-CoA to phosphatidic acid, conversion of phosphatidic acid to diacylglycerol, and conversion of diacylglycerol to triacylglycerol.
triglyceride_biosynthetic <- c("Mogat1","Lpin1","Dgat1","Lpin2","Gpam","Gk","Lpin3","Dgat2","Gpat2","Agmo","Gk2","Mogat2","Gykl1") #R-MMU-75109.1
triglyceride_catabolic <- c("Pnpla5","Fabp7","Fabp2","Gpd2","Fabp9","Fabp12","Fabp5","Fabp3","Fabp1","Fabp4") #R-MMU-163560.1

Lipogeneic_genes <- c("Fasn","Acly","Acaca","Elovl6","Scd","Gpat","Agpat2","Lipin",#Enzymes
                      "Srebp1","Srebp1c","Srebp1a","Chrebp","Lxra","Usf1", "Igf2r","Pparg" ) 
Lipolytic_genes <-c("Abcg1","Apoe","Cd36")

Nuclear_receptors <- c("Esr1","Esr2","Esrra","Esrrb","Esrrg","Ar","Nr3c1","Nr3c2",#Estrogen Receptor-like
                       "Thra","Thrb","Ppara","Ppard","Pparg","Nr1d1","Nr1d2","Vdr","Nr1i2","Nr1i3", #Thyroid Hormone Receptor-like
                       "Rxra","Rxrb","Rxrg","Rora","Rorb","Rorc",#Retinoid X Receptor-like
                       "Nr1h4","Nr1h5","Nr1h3","Nr1h2",#Liver X receptor-like receptors
                       "Hnf4a","Nr2c1","Nr2c2", #Hepatocyte nuclear factor-4 receptors/Testicular receptors
                       "Nr2f1","Nr2f2","Nr2f6", #COUP-TF-like receptors
                       "Nr4a1","Nr4a2","Nr4a3", #Nerve growth factor IB-like receptors
                       "Nr5a1","Nr5a2","Nr6a1", #Fushi tarazu F1-like receptors/Germ cell nuclear factor receptors 
                       "Erbb2","Tshr","Tgfb1","Lepr"
)

metabolic_pathways <- c("adaptive thermogenesis","antibiotic metabolic process","biosynthetic process","catabolic process","catabolite repression",
                        "cellular metabolic process","collagen metabolic process","demethylation","epinephrine metabolic process","futile creatine cycle",
                        "glycosylation","hormone metabolic process","hydrogen metabolic process","iron-sulfur-molybdenum cofactor metabolic process",
                        "methylation","negative regulation of metabolic process","nicotine metabolic process","nitrogen compound metabolic process",
                        "organic substance metabolic process","phosphinothricin metabolic process","pigment metabolic process",
                        "positive regulation of metabolic process","primary metabolic process","prosthetic group metabolic process",
                        "regulation of metabolic process","respiratory burst","secondary metabolic process","small molecule metabolic process",
                        "tartrate metabolic process","tyrocidine metabolic process")

#===============================================================================
# Seurat command cheetsheet
#===============================================================================
if(F){ 
  # Get cell and feature names, and total numbers
  rownames(All_combined);colnames(All_combined)
  Cells(All_combined)
  ncol(All_combined);nrow(All_combined)
  # Get cell identity classes
  Idents(All_combined);levels(All_combined)
  # Subset Seurat object based on identity class, also see ?SubsetData
  subset(All_combined, idents = "B cells")
  subset(All_combined, idents = c("CD4 T cells", "CD8 T cells"), invert = TRUE)
  # Subset on the expression level of a gene/feature
  subset(All_combined, subset = Ucp1 > 3)
  # Subset on a combination of criteria
  subset(All_combined, subset = Ucp1 > 3 & PC1 > 5)
  subset(All_combined, subset = Ucp1 > 3, idents = "Adipocytes")
  # Downsample the number of cells per identity class
  subset(All_combined, downsample = 100)
  
  # Retrieve or set data in an expression matrix ('counts', 'data', and 'scale.data')
  GetAssayData(All_combined, slot = "scale.data")
  All_combined <- SetAssayData(All_combined, slot = "scale.data", new.data = new.data)
  # Get cell embeddings and feature loadings
  Embeddings(All_combined, reduction = "pca")
  Loadings(All_combined, reduction = "pca")
  # FetchData can pull anything from expression matrices, cell embeddings, or metadata
  FetchData(All_combined, vars = c("PC_1", "Myannotation", "Ucp1"))
  FetchData(All_combined, vars = c("rna_CD3E", "adt_CD3"))  # Pull feature expression from both assays by using keys
  FeatureScatter(All_combined, feature1 = "Ucp1", feature2 = "Cidea")
  FeatureScatter(All_combined, feature1 = "Ucp1", feature2 = "CD3D")
  CellScatter(All_combined, cell1 = "GCCAGCATCCAAGAGG-1_1", cell2 = "ATCATTCGTAGTCTTG-1_1")
  VariableFeaturePlot(All_combined)
  RidgePlot(All_combined, feature = c("Ucp1", "Cidea", "Aldh1a1"))

  plot <- DimPlot(All_combined, reduction = "umap",label = T, repel = T, label.box =T,  cols = rAT_Ad_colorMaps, order = T)
  HoverLocator(plot, information = FetchData(All_combined, vars = c("Myannotation", "Ucp1")))
  select.cells <- CellSelector(plot = plot)
  Idents(All_combined, cells = select.cells) <- "NewCells"
  LabelPoints(plot, points = TopCells(All_combined[["pca"]]), repel = TRUE)  # Label points on a ggplot object
  # Plot data from multiple assays using keys
  FeatureScatter(All_combined, feature1 = "rna_CD3E", feature2 = "adt_CD3")
  
  #remove cells by identity
  Obj <- Obj[,!colnames(Obj) %in% WhichCells(Obj, idents = c("3"))]
  #remove cells by selected ID
  Idents(Obj) <- Obj$integrated_snn_res.0.2
  select.cells <- CellSelector(plot = DimPlot(rATM_CE_comb, reduction = "umap"))
  select.cells <- c(CellSelector(plot = DimPlot(rATM_CE_comb, reduction = "umap")), select.cells)
  Obj <- Obj[,!colnames(Obj) %in% select.cells]
  Obj <- ScaleData(Obj) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = seq(from=0, by=0.05, length=10))
  #Remove genes
  mito.genes <- grep(pattern = "^MT-", rownames(Obj@assays$RNA), value = TRUE)
  Obj <- subset(Obj, features = setdiff(rownames(Obj@assays$RNA), mito.genes))
}

#===============================================================================
# Custome functions for sc/snRNA-seq analysis 
#===============================================================================
if(F){ #------------------- 1 Trajectory analysis ---------------------------
  cds <- as.cell_data_set(mPRAT_Ad_merged,assay ="RNA")
  rowData(cds)$gene_short_name <- rownames(rowData(cds))
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- order_cells(cds)
  
  pdf(paste0(mPRAT_Ad_prefix,"Sudotime_",Sys.Date(),".pdf"), width = 5, height = 4)
  plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=F, cell_size = 0.35, label_roots = F, label_leaves=F, 
             label_branch_points=F, graph_label_size=3, trajectory_graph_segment_size = 1, trajectory_graph_color = "black")
  dev.off()
  
  mPRAT_Ad_merged <- AddMetaData(object = mPRAT_Ad_merged, metadata = cds@principal_graph_aux@listData$UMAP$pseudotime, col.name = "pseudotime")
  modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 20)
  write.csv(modulated_genes, file = paste0(mPRAT_Ad_prefix, "modulated_genes_",Sys.Date(),".csv"), quote = F)
  
  # modulated_genes <- read.csv("mPRAT_Adipocytes_modulated_genes_2023-06-08.csv", row.names = 1)
  genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.1))
  genes <- setdiff(genes, grep("^mt-",genes,value = T))
  pt.matrix <- normalized_counts(cds, norm_method = "log")[match(genes, rownames(rowData(cds))), order(pseudotime(cds))]
  metainfor <- colData(cds)[order(pseudotime(cds)),]
  metainfor$pseudotime <- pseudotime(cds)[order(pseudotime(cds))]
  pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
  pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
  rownames(pt.matrix) <- genes
  
  pdf(paste0(mPRAT_Ad_prefix,"Sudotime_heatmap3_",Sys.Date(),".pdf"), width = 8, height = 8)
  ha = columnAnnotation(CellType = metainfor$Myannotation_Ad,col = list(CellType = mPRAT_Ad_colors))
  p <- Heatmap(pt.matrix, name  = "z-score", col = colorRamp2(seq(from=-2,to=2, length=11), rev(brewer.pal(11, "Spectral"))),
               show_row_names = TRUE, show_column_names = FALSE, row_names_gp = gpar(fontsize = 6), raster_device = "png",
               row_title_rot = 0, km = 3, cluster_rows = T, cluster_row_slices = FALSE, cluster_columns = FALSE, bottom_annotation = ha)
  draw(p)
  dev.off()
  saveRDS(p, paste0(mPRAT_Ad_prefix,"Sudotime_heatmap_res_",Sys.Date(),".rds"))
}
if(F){ #------------------- 2 label transfer from Time-series mPRAT dataset to rat ---------------------------
  RatTransfer_prefix <- "RmPRAT_Adipocyte_Transfer_"
  mPRAT_Ad_merged <- subset(mPRAT_Ad_merged, subset = orig.ident == "mPRAT-2mo")
  mRPATanchors <- FindTransferAnchors(reference = mPRAT_Ad_merged, query = RmPRAT_Ad_merged, dims = 1:15, reference.reduction = "pca")
  predictions <- TransferData(anchorset = mRPATanchors, refdata = mPRAT_Ad_merged$Myannotation_Ad, dims = 1:15)
  RmPRAT_Ad_merged <- AddMetaData(RmPRAT_Ad_merged, metadata = predictions)
  
  mRPAT_Ad_merged <- RunUMAP(mPRAT_Ad_merged, dims = 1:15, reduction = "pca", return.model = TRUE)
  RmPRAT_Ad_merged <- MapQuery(anchorset = mRPATanchors, reference = mPRAT_Ad_merged, refdata = list(celltype = "Myannotation_Ad"),  
                               reference.reduction = "pca", reduction.model = "umap", query = RmPRAT_Ad_merged)
  RmPRAT_Ad_merged$predicted.id <- factor(RmPRAT_Ad_merged$predicted.id, levels = c("mPRAT-ad1","mPRAT-ad2","mPRAT-ad3","mPRAT-ad4"))
  
  celltypes <- c("mPRAT-ad1","mPRAT-ad2","mPRAT-ad3","mPRAT-ad4")
  colors <- c("#BC3C29FF","#FB9A99","#B2DF8A","#20854EFF")
  names(colors) <- celltypes
  
  Plot_Cell_compoistion(Seurat_Obj = RmPRAT_Ad_merged, OutPrefix = paste0(RatTransfer_prefix,"NewLabel_"), ColorUse1 = colors, ColorUse2 = RmPRAT_SAMcolors)
  
  (p0 <- DimPlot(mRPAT_Ad_merged, reduction = "umap", group.by = "Myannotation_Ad", split.by = "orig.ident", label = TRUE, label.size = 3, repel = TRUE, cols = colors) + NoLegend())
  (p1 <- DimPlot(RmPRAT_Ad_merged, reduction = "umap", group.by = "Myannotation_Ad", split.by = "orig.ident", label = TRUE, label.size = 3, repel = TRUE, cols = RmPRAT_Ad_colors) + NoLegend())
  (p2 <- DimPlot(RmPRAT_Ad_merged, reduction = "umap", group.by = "predicted.id", split.by = "orig.ident", label = TRUE, label.size = 3, repel = TRUE, cols = colors) + NoLegend())
  # (p3 <- DimPlot(RmPRAT_Ad_merged, reduction = "ref.umap", split.by = "orig.ident", label = TRUE, label.size = 3, repel = TRUE, cols = colors) + NoLegend())
  p1/p2/p0
  ggsave(paste0(RatTransfer_prefix, "UMAP_", Sys.Date(), ".pdf"), width = 7, height = 16)
  
  subset(RmPRAT_Ad_merged, subset = orig.ident == "RmPRAT-PV")
  
  subset(RmPRAT_Ad_merged, subset = orig.ident == "RmPRAT-PU")@meta.data[,c("Myannotation_Ad","predicted.celltype")] %>% 
    dplyr::count(Myannotation_Ad, predicted.celltype) %>% filter(Myannotation_Ad != "NA") %>% 
    ggplot(aes(axis1 = Myannotation_Ad, axis2 = predicted.celltype, y = n)) +
    geom_alluvium(aes(fill = Myannotation_Ad), alpha=0.5) +
    geom_stratum(aes(fill = Myannotation_Ad)) + 
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_manual(values = c("#BC3C29FF","#EE4C97FF","#47A93A","#20854EFF")) + 
    theme_void() + theme(legend.position = "none")
  ggsave(paste0(mRPAT_Ad_prefix,"SankeyPlot_RmPRAT-PU_",Sys.Date(),".pdf"), width = 4, height = 8)
  
  CompositionFreq <- table(RmPRAT_Ad_merged$predicted.id, RmPRAT_Ad_merged$orig.ident)%>% as.data.frame() 
  mat <- CompositionFreq %>% group_by(Var2) %>% mutate(Freq = Freq/sum(Freq) * 100)
  mat$Var2 <- factor(mat$Var2, levels = levels(mat$Var2))
  
  ggplot(mat, aes(x="", y=Freq, fill = Var1, label = paste0(round(Freq, 1), "%"))) + 
    geom_bar(stat="identity", width = 1, color = "black") +
    coord_polar("y", start = 0) +
    geom_text(position=position_stack(vjust=0.5), color="white",size=8) + 
    scale_fill_manual(values = c("#BC3C29FF","#FB9A99","#B2DF8A","#20854EFF")) + 
    facet_wrap(~Var2) +
    theme_void() + 
    labs(x="", y="Percentage of each cell type") +
    theme(
      axis.title.x = element_text(color="black", size=16, face="bold"),
      axis.text.x = element_blank(),
      axis.title.y = element_text(color="black", size=16, face="bold"),
      axis.text.y = element_text(color="black", size=16, face="bold")
    )
  
  ggsave(paste0(RatTransfer_prefix, "cellNumbersPerCluster_Pie_pvpu_",Sys.Date(),".pdf"), width = 8, height = 4)
}

if(F){ #------------------- 3 Use CellCycleScoring function to perform cell cycle scoring ------------------- 
  DefaultAssay(iWAT_Adipocyte) <- "RNA"
  iWAT_Adipocyte <- CellCycleScoring(iWAT_Adipocyte, s.features = read.csv("../../Mouse_S-phase_genes.csv")$x, g2m.features = read.csv("../../Mouse_G2M-phase_genes.csv")$x, set.ident = TRUE)
  table(iWAT_Adipocyte$Phase)
  
  Idents(iWAT_Adipocyte) <- iWAT_Adipocyte$Myannotation_Ad
  DimPlot(iWAT_Adipocyte, reduction = "umap", group.by = c("Phase"), cols = hughie_color, split.by = "orig.ident", ncol = 3)
  ggsave(paste0(iWAT_merged_prefix, "UMAP-ByCellPhase_",Sys.Date(),".pdf"), width = 12, height = 10)
  Plot_CellCycle(iWAT_Adipocyte, OutPrefix = "iWAT_15CDvsHFD_CellCycle_", cellRanks = iWAT_AdCellRanks, ColorUse = hughie_color)
  Idents(iWAT_Adipocyte) <- iWAT_Adipocyte$Myannotation
}

if(F){ #------------------ Compare DEGs among groups ------------------
  Idents(iWAT_merged) <- iWAT_merged$Myannotation
  Sample_names <- unique(iWAT_merged$orig.ident)
  cell_names <- as.vector(unique(Idents(iWAT_merged)))
  iWAT_merged$CT_SAMPLE <- paste(Idents(iWAT_merged), iWAT_merged$orig.ident, sep = "_")
  iWAT_merged$celltype <- Idents(iWAT_merged)
  Idents(iWAT_merged) <- "CT_SAMPLE"
  
  iWAT_merged_ad <- subset(iWAT_merged, idents = c("Adipocytes_iWAT_KO_15CD","Adipocytes_iWAT_WT_15CD","Adipocytes_iWAT_WT_HFD","Adipocytes_iWAT_KO_HFD"))
  compare_comb <- combn(unique(iWAT_merged_ad$CT_SAMPLE),2)
  compare_comb[1,2] <- "Adipocytes_iWAT_KO_HFD"; compare_comb[2,2] <- "Adipocytes_iWAT_KO_15CD"
  compare_comb[1,5] <- "Adipocytes_iWAT_WT_HFD"; compare_comb[2,5] <- "Adipocytes_iWAT_WT_15CD"
  
  DEG_genes <- list()
  for (i in c(1,2,5,6)){
    Seurat_obj_DEG <- FindMarkers(iWAT_merged_ad, ident.1 = compare_comb[1,i], ident.2 = compare_comb[2,i], min.pct = 0.1, logfc.threshold = 0.25, assay = "RNA")
    Seurat_obj_DEG$change <- ifelse(Seurat_obj_DEG$p_val_adj < 0.05 & abs(Seurat_obj_DEG$avg_log2FC) >= 0.20, 
                                    ifelse(Seurat_obj_DEG$avg_log2FC > 0.20 ,'Up','Down'),'Stable')
    write.csv(Seurat_obj_DEG, paste0(iWAT_merged_prefix,compare_comb[1,i],"vs",compare_comb[2,i][1],"_DEG_List_",Sys.Date(),".csv"))
    
    Seurat_obj_DEG_up <- row.names(Seurat_obj_DEG)[Seurat_obj_DEG$change == "Up"] 
    Seurat_obj_DEG_down <- row.names(Seurat_obj_DEG)[Seurat_obj_DEG$change == "Down"] 
    DEG_genes[[paste0(compare_comb[1,i],"vs",compare_comb[2,i],"-Up")]] <- Seurat_obj_DEG_up
    DEG_genes[[paste0(compare_comb[1,i],"vs",compare_comb[2,i],"-Down")]] <-Seurat_obj_DEG_down
  }
  
  names(DEG_genes) <- gsub("Adipocytes_iWAT_","",names(DEG_genes))
  Result_suffix <- paste0(iWAT_merged_prefix,"Adipocytes_comparison_",Sys.Date())
  compGO <- compareCluster(geneCluster = DEG_genes, fun= "enrichGO", pvalueCutoff = 0.05, keyType = "SYMBOL",
                           pAdjustMethod = "BH", OrgDb = org.Mm.eg.db, readable = T, ont = "BP")
  write.csv(compGO@compareClusterResult, file = paste0(Result_suffix,"_GO-BP_list.csv"))
  compGO_filtered <- gofilter(compGO, level = 4)
  dotplot(compGO_filtered, showCategory = 20, title = paste0(Result_suffix,"\n - GO-BP Enrichment Analysis")) +  
    scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(Result_suffix,"_GO-BP_dotplot_Filter4","_20-Items.pdf"), width = 12, height = 18)
  
  #==== Pair-wise enrichment analysis ====
  if(F){
    DefaultAssay(iWAT_Adipocyte) <- "RNA"
    PairWise_DEG_enrichment(Seurat_obj = iWAT_Adipocyte, OutPrefix = "_iWAT_Adipocyte_", PlotDEGHeatmap = F,
                            CellTypeInclude = as.vector(unique(iWAT_Adipocyte@active.ident)))
    DEGs <- read.csv("2022-10-23_iWAT_Adipocyte_All DEG_List.csv")
    DEGs <- strsplit(DEGs$x, split = " ")
    GeneID_list <- list()
    GeneID_list[[DEGs[[1]][1]]] <- strsplit(DEGs[[1]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[3]][1]]] <- strsplit(DEGs[[3]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[9]][1]]] <- strsplit(DEGs[[9]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[7]][1]]] <- strsplit(DEGs[[7]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[5]][1]]] <- strsplit(DEGs[[5]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[2]][1]]] <- strsplit(DEGs[[2]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[4]][1]]] <- strsplit(DEGs[[4]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[10]][1]]] <- strsplit(DEGs[[10]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[8]][1]]] <- strsplit(DEGs[[8]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[6]][1]]] <- strsplit(DEGs[[6]][2], split = ",")[[1]]
    
    compGO <- compareCluster(geneCluster = GeneID_list, fun= "enrichGO", pvalueCutoff = 0.05, keyType = "SYMBOL",
                             pAdjustMethod = "BH", OrgDb = org.Mm.eg.db, readable = T, ont = "BP")
    
    write.csv(compGO@compareClusterResult, file = paste0(Sys.Date(), iWAT_Adipocyte_prefix,"_Adipocytes_1-5_GO-BP_list.csv"))
    compGO_filtered <- gofilter(compGO, level = 4)
    dotplot(compGO_filtered, showCategory = 10, title = paste0("GO Pathways up/down regulated of Adipocytes_1-5 after ALK7 KO")) + 
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0(Sys.Date(), iWAT_Adipocyte_prefix,"Adipocytes_1-5_GO-BP_dotplot.pdf"), width = 14, height = 14)
    
    #Heatmap
    Sample_names <- unique(iWAT_Adipocyte$orig.ident)
    cell_names <- as.vector(unique(Idents(iWAT_Adipocyte)))
    iWAT_Adipocyte$CT_SAMPLE <- paste(Idents(iWAT_Adipocyte), iWAT_Adipocyte$orig.ident, sep = "_")
    iWAT_Adipocyte$celltype <- Idents(iWAT_Adipocyte)
    Idents(iWAT_Adipocyte) <- "CT_SAMPLE"
    iWAT_Adipocyte_averge <- AverageExpression(iWAT_Adipocyte)
    PlotMatrix <- iWAT_Adipocyte_averge$RNA %>% as.data.frame()
    colnames(PlotMatrix) <- gsub("_iWAT|_HFD","",colnames(PlotMatrix))
    PlotMatrix <- PlotMatrix[,c(6,8:9,7,10,1,2,5,4,3)]
    
    PathwaysTOPlot <- c("lipid modification","lipid oxidation","lipid catabolic process", "lipid storage","regulation of lipid localization",
                        "response to transforming growth factor beta","cold-induced thermogenesis",
                        "temperature homeostasis","regulation of membrane potential")
    
    for (Pathway in PathwaysTOPlot){
      GenesToPlot <- unique(strsplit(paste0(compGO@compareClusterResult[compGO@compareClusterResult$Description %in% Pathway,]$geneID, 
                                            collapse = "/"), split = "/")[[1]])
      pdf(paste0(Sys.Date(), iWAT_Adipocyte_prefix, Pathway, "_DEG_heatmap.pdf"), width = 12, height = 12)
      p <- pheatmap(PlotMatrix[GenesToPlot,], 
                    show_colnames = T, show_rownames = T, scale = "row",border_color = "black",
                    cluster_rows = T, cluster_cols = F, main=paste0('Gene expression in (', Pathway, ") Pathway"))
      plot(p)
      dev.off()
    }
  }
  #Plot Alk7 pathways
  if(F){
    GenesToPlot <- c("Inhbb","Gdf3","Smad2","Smad3","Adrb1","Adrb2","Adrb3","Rgs2",
                     "Ppargc1a","Ucp1","Ucp2","Ucp3","Pparg","Cebpa","Lipe","Pnpla2",
                     "Klf15", "Prodh","Gpt","Bcat2", #aa degrading enzymes
                     "Lep","Fabp4","Zfp423","Cd24a","Acvr1b","Adipoq", "Cycs")
    
    pdf(paste0(iWAT_Adipocyte_prefix, "SelectedPaperGenes_heatmap",Sys.Date(),".pdf"), width = 10, height = 10)
    pheatmap(PlotMatrix[GenesToPlot,],show_colnames = T, show_rownames = T, scale = "row", border_color = "black",cluster_rows = F, cluster_cols = F)
    dev.off()
  }
  
  #==== Step 4. Pair-wise enrichment analysis ====
  if(F){
    DefaultAssay(iWAT_ASPC) <- "RNA"
    PairWise_DEG_enrichment(Seurat_obj = iWAT_ASPC, OutPrefix = "_iWAT_ASPC_", PlotDEGHeatmap = F, RunGO = F,
                            CellTypeInclude = as.vector(unique(iWAT_ASPC@active.ident)))
    DEGs <- read.csv("2022-10-24_iWAT_ASPC_All DEG_List.csv")
    DEGs <- strsplit(DEGs$x, split = " ")
    GeneID_list <- list()
    GeneID_list[[DEGs[[1]][1]]] <- strsplit(DEGs[[1]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[3]][1]]] <- strsplit(DEGs[[3]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[5]][1]]] <- strsplit(DEGs[[5]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[7]][1]]] <- strsplit(DEGs[[7]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[2]][1]]] <- strsplit(DEGs[[2]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[4]][1]]] <- strsplit(DEGs[[4]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[6]][1]]] <- strsplit(DEGs[[6]][2], split = ",")[[1]]
    GeneID_list[[DEGs[[8]][1]]] <- strsplit(DEGs[[8]][2], split = ",")[[1]]
    
    compGO <- compareCluster(geneCluster = GeneID_list, fun= "enrichGO", pvalueCutoff = 0.05, keyType = "SYMBOL",
                             pAdjustMethod = "BH", OrgDb = org.Mm.eg.db, readable = T, ont = "BP")
    
    write.csv(compGO@compareClusterResult, file = paste0(Sys.Date(), iWAT_ASPC_prefix,"1-4_GO-BP_list.csv"))
    compGO_filtered <- gofilter(compGO, level = 5)
    dotplot(compGO_filtered, showCategory = 20, title = paste0("GO Pathways up/down regulated of ASPC1-4 after ALK7 KO")) + 
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0(Sys.Date(), iWAT_ASPC_prefix,"1-4_GO-BP_dotplot.pdf"), width = 14, height = 14)
    
    #Heatmap
    Sample_names <- unique(iWAT_ASPC$orig.ident)
    cell_names <- as.vector(unique(Idents(iWAT_ASPC)))
    iWAT_ASPC$CT_SAMPLE <- paste(Idents(iWAT_ASPC), iWAT_ASPC$orig.ident, sep = "_")
    iWAT_ASPC$celltype <- Idents(iWAT_ASPC)
    Idents(iWAT_ASPC) <- "CT_SAMPLE"
    iWAT_Adipocyte_averge <- AverageExpression(iWAT_ASPC)
    PlotMatrix <- iWAT_Adipocyte_averge$RNA %>% as.data.frame()
    colnames(PlotMatrix) <- gsub("_iWAT|_HFD","",colnames(PlotMatrix))
    PlotMatrix <- PlotMatrix[,c(6,5,8,7,1:4)]
    
    PathwaysTOPlot <- c("extracellular matrix organization","epithelial cell proliferation","fibroblast proliferation",
                        "fat cell differentiation","negative regulation of lipid metabolic process",
                        "negative regulation of locomotion","cell-cell signaling by wnt","response to transforming growth factor beta")
    
    for (Pathway in PathwaysTOPlot){
      GenesToPlot <- unique(strsplit(paste0(compGO@compareClusterResult[compGO@compareClusterResult$Description %in% Pathway,]$geneID, 
                                            collapse = "/"), split = "/")[[1]])
      pdf(paste0(Sys.Date(), iWAT_ASPC_prefix,Pathway, "_DEG_heatmap.pdf"), width = 12, height = 12)
      p <- pheatmap(PlotMatrix[GenesToPlot,], 
                    show_colnames = T, show_rownames = T, scale = "row",border_color = "black",
                    cluster_rows = T, cluster_cols = F, main=paste0('Gene expression in (', Pathway, ") Pathway"))
      plot(p)
      dev.off()
    }
  }
}

##' A warpper to run Seurat preprocessing
##' @param Seurat_Obj A Seurat object
##' @param NormalizationMethod Choose either VST or SCT
##' @param removeDoublet Logic value to act on doublet remove
##' @param DoubletPercent Percentage of doublet in your dataset (8% for 20K cells loaded in 10X)
##' @param species Specify species(mouse, rat, human) you researched (default: mouse)
##' @param Filter_nFeature_RNA threshold of minimal and maximal genes detected in each cell
##' @param Filter_nCount_RNA threshold of minimal and maximal UMI counts detected in each cell
##' @param MitoPercent threshold to remove cells based on mitochondria reads (default: 15%)
##' @param OutPrefix Prefix for output file
##' @examples
##'
##' RunSeurat(ObjList = c(rATM1mo_Obj, rATM2mo_Obj, rATM6mo_Obj), OutPrefix = rATM_prefix, NormalizationMethod = "VST",
##'           Filter_nFeature_RNA = c(500, 6000), Filter_nCount_RNA = c(1000, 50000), MitoPercent = 15)

PlotFeature <- function(ObjList, Prefix){
  pdf(paste0(Prefix,".pdf"), width = 9, height = 6)
  for(i in 1:length(ObjList)){p <- VlnPlot(ObjList[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3); plot(p)}
  dev.off()
}

RunSeurat_single <- function(Obj, OutPrefix, NormalizationMethod = "VST", removeDoublet = T, species = "mouse",
                             Filter_nFeature_RNA = c(500, 6000), Filter_nCount_RNA = c(1000, 50000), MitoPercent = 15, DoubletPercent = 0.08){
  if(species == "mouse"){mtPatterm <- "^mt-"
  } else if(species == "rat"){mtPatterm <- "^Mt-"
  } else if(species == "human"){mtPatterm <- "^MT-"
  }
  
  Obj[["percent.mt"]] <- PercentageFeatureSet(Obj, pattern = mtPatterm)
  Obj <- subset(Obj, subset = nFeature_RNA > Filter_nFeature_RNA[1] & nFeature_RNA < Filter_nFeature_RNA[2] 
                & nCount_RNA > Filter_nCount_RNA[1] & nCount_RNA < Filter_nCount_RNA[2] & percent.mt < MitoPercent)
  Obj <- NormalizeData(Obj, normalization.method = "LogNormalize", scale.factor = 10000)
  Obj <- FindVariableFeatures(Obj, selection.method = "vst", nfeatures = 2000)
  Obj <- ScaleData(Obj) %>% RunPCA(dims = 1:15) %>% RunUMAP(dims = 1:15) %>% 
    FindNeighbors(reduction = "pca", dims = 1:15) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  
  pdf(paste0(OutPrefix,"UMAP-ByDoublets_",Sys.Date(),".pdf"), width = 8, height = 8)
  ## pK Identification
  sweep.res.list <- paramSweep(Obj, PCs = 1:15, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  optimal.pk <- as.numeric(as.vector(bcmvn[which(bcmvn$BCmetric == max(bcmvn$BCmetric)),]$pK))
  
  ## Homotypic Doublet Proportion Estimate
  annotations <- Obj@meta.data$RNA_snn_res.0.1
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(DoubletPercent*nrow(Obj@meta.data)) ## Assuming 8% doublet formation rate based 10X guideline when load 20,000 cells
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies
  Obj <- doubletFinder(Obj, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  # Obj <- doubletFinder(Obj, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_",optimal.pk,"_",nExp_poi), sct = FALSE)
  Idents(Obj) <- Obj@meta.data[,paste0("DF.classifications_0.25_",optimal.pk,"_",nExp_poi)]
  p <- DimPlot(Obj, reduction = "umap", label = T, repel = T, label.box =T, cols = hughie_color, label.color = "white");plot(p)
  Idents(Obj) <- Obj$seurat_clusters
  toRemove <- rownames(Obj@meta.data[Obj@meta.data[,ncol(Obj@meta.data)] == "Doublet",])
  Obj <- Obj[,!colnames(Obj) %in% toRemove]
  cat(paste0(length(toRemove)," doublets were removed !!!\n"))
  dev.off()
  
  Obj <- FindVariableFeatures(Obj, selection.method = "vst", nfeatures = 2000)
  Obj <- ScaleData(Obj) %>% RunPCA(dims = 1:15) %>% RunUMAP(dims = 1:15) %>% 
    FindNeighbors(reduction = "pca", dims = 1:15) %>% FindClusters(resolution = seq(from=0, by=0.05, length=3))  

  saveRDS(Obj, file = paste0(OutPrefix,"Seurat.rds"))
  
  # Visualization of the clusters landscape
  DimPlot(Obj, reduction = "umap", group.by = c("ident","orig.ident"), label = T, repel = T, label.box =T, cols = hughie_color, label.color = "white")
  ggsave(paste0(OutPrefix,"UMAP-ByClusters_",Sys.Date(),".pdf"), width = 16, height = 8)
  
  DimPlot(Obj, reduction = "umap", split.by = "orig.ident", ncol = 2, label = T, repel = T, cols = hughie_color)
  ggsave(paste0(OutPrefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 12, height = 12)
  
  Plot_Cell_compoistion(Seurat_Obj = Obj, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = OutPrefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = Obj, ColorUse = hughie_color, OutPrefix = OutPrefix)
  return(Obj)
}

RunSeurat <- function(ObjList, OutPrefix, NormalizationMethod = "VST", removeDoublet = T, species = "mouse",k_weight = 100,
                      Filter_nFeature_RNA = c(500, 6000), Filter_nCount_RNA = c(1000, 50000), MitoPercent = 15, DoubletPercent = 0.08){

  if(species == "mouse"){mtPatterm <- "^mt-"
  } else if(species == "rat"){mtPatterm <- "^Mt-"
  } else if(species == "human"){mtPatterm <- "^MT-"
  }
  
  # Plot data metrics before filtering
  PlotFeature(ObjList, Prefix = paste0(OutPrefix,"BeforeFilter_",Sys.Date()))
  # future::plan("multisession", workers = 20) # do parallel

  if(NormalizationMethod == "VST"){
    ObjList <- lapply(X = ObjList, FUN = function(Obj) {
      #1. Low-quality cells or empty droplets will often have very few genes. 
      #2. Cell doublets or multiplets may exhibit an aberrantly high gene count. 
      #3. Low-quality / dying cells often exhibit extensive mitochondrial contamination
      Obj[["percent.mt"]] <- PercentageFeatureSet(Obj, pattern = mtPatterm)
      Obj <- subset(Obj, subset = nFeature_RNA > Filter_nFeature_RNA[1] & nFeature_RNA < Filter_nFeature_RNA[2] 
                    & nCount_RNA > Filter_nCount_RNA[1] & nCount_RNA < Filter_nCount_RNA[2] & percent.mt < MitoPercent)
      Obj <- NormalizeData(Obj, normalization.method = "LogNormalize", scale.factor = 10000)
      Obj <- FindVariableFeatures(Obj, selection.method = "vst", nfeatures = 2000)
      Obj <- ScaleData(Obj) %>% RunPCA(dims = 1:15) %>% RunUMAP(dims = 1:15) %>% 
        FindNeighbors(reduction = "pca", dims = 1:15) %>% FindClusters(resolution = seq(from=0, by=0.05, length=3))
    })
    
    pdf(paste0(OutPrefix,"UMAP-ByDoublets_",Sys.Date(),".pdf"), width = 8, height = 8)
    if(removeDoublet){
      for (i in 1:length(ObjList)){
        ## pK Identification
        sweep.res.list <- paramSweep(ObjList[[i]], PCs = 1:15, sct = F)
        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        optimal.pk <- as.numeric(as.vector(bcmvn[which(bcmvn$BCmetric == max(bcmvn$BCmetric)),]$pK))
        
        ## Homotypic Doublet Proportion Estimate
        annotations <- ObjList[[i]]@meta.data$RNA_snn_res.0.1
        homotypic.prop <- modelHomotypic(annotations)
        nExp_poi <- round(DoubletPercent*nrow(ObjList[[i]]@meta.data)) ## Assuming 8% doublet formation rate based 10X guideline when load 20,000 cells
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        
        ## Run DoubletFinder with varying classification stringencies
        ObjList[[i]] <- doubletFinder(ObjList[[i]], PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
        # ObjList[[i]] <- doubletFinder(ObjList[[i]], PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_",optimal.pk,"_",nExp_poi), sct = FALSE)
        Idents(ObjList[[i]]) <- ObjList[[i]]@meta.data[,paste0("DF.classifications_0.25_",optimal.pk,"_",nExp_poi)]
        p <- DimPlot(ObjList[[i]], reduction = "umap", label = T, repel = T, label.box =T, cols = hughie_color, label.color = "white");plot(p)
        Idents(ObjList[[i]]) <- ObjList[[i]]$seurat_clusters
        toRemove <- rownames(ObjList[[i]]@meta.data[ObjList[[i]]@meta.data[,ncol(ObjList[[i]]@meta.data)] == "Doublet",])
        ObjList[[i]] <- ObjList[[i]][,!colnames(ObjList[[i]]) %in% toRemove]
        cat(paste0(length(toRemove)," doublets were removed !!!\n"))
      }
    }
    dev.off()

    ObjList <- lapply(X = ObjList, FUN = function(Obj) {
      Obj <- FindVariableFeatures(Obj, selection.method = "vst", nfeatures = 2000)
      Obj <- ScaleData(Obj)
    })
    
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = ObjList, nfeatures = 2000)
    anchors <-  FindIntegrationAnchors(object.list = ObjList, anchor.features = features)
    combined <- IntegrateData(anchorset = anchors, k.weight = k_weight) 
    DefaultAssay(combined) <- "integrated"
  }
  
  #Using SCT for normalization
  if(NormalizationMethod == "SCT"){
    ObjList <- lapply(X = ObjList, FUN = function(Obj) {
      Obj[["percent.mt"]] <- PercentageFeatureSet(Obj, pattern = "^mt-")
      Obj <- subset(Obj, subset = nFeature_RNA > Filter_nFeature_RNA[1] & nFeature_RNA < Filter_nFeature_RNA[2] & 
                      nCount_RNA > Filter_nCount_RNA[1] & nCount_RNA < Filter_nCount_RNA[2] & percent.mt < MitoPercent)      
      Obj <- NormalizeData(Obj, normalization.method = "LogNormalize", scale.factor = 10000)
      Obj <- SCTransform(Obj, method = "glmGamPoi", vst.flavor="v2", vars.to.regress = "percent.mt", verbose = T)
      Obj <- RunPCA(Obj, dims = 1:15) %>% RunUMAP(dims = 1:15)
    })
    
    if(removeDoublet){
      for (i in 1:length(ObjList)){
        ## pK Identification
        sweep.res.list <- paramSweep(ObjList[[i]], PCs = 1:15, sct = T)
        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        optimal.pk <- as.numeric(as.vector(bcmvn[which(bcmvn$BCmetric == max(bcmvn$BCmetric)),]$pK))
        
        ## Homotypic Doublet Proportion Estimate
        annotations <- ObjList[[i]]@meta.data$SCT_snn_res.0.1
        homotypic.prop <- modelHomotypic(annotations)
        nExp_poi <- round(DoubletPercent*nrow(ObjList[[i]]@meta.data)) ## Assuming 8% doublet formation rate based 10X guideline when load 20,000 cells
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        
        ## Run DoubletFinder with varying classification stringencies
        ObjList[[i]] <- doubletFinder(ObjList[[i]], PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
        # ObjList[[i]] <- doubletFinder(ObjList[[i]], PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_",optimal.pk,"_",nExp_poi), sct = T)
        toRemove <- rownames(ObjList[[i]]@meta.data[ObjList[[i]]@meta.data[,ncol(ObjList[[i]]@meta.data)] == "Doublet",])
        ObjList[[i]] <- ObjList[[i]][,!colnames(ObjList[[i]]) %in% toRemove]
        cat(paste0(length(toRemove)," doublets were removed !!!\n"))
      }
    }
    
    ObjList <- lapply(X = ObjList, FUN = function(Obj) {
      Obj <- FindVariableFeatures(Obj, selection.method = "vst", nfeatures = 2000)
      Obj <- ScaleData(Obj)
    })
    features <- SelectIntegrationFeatures(object.list = ObjList, nfeatures = 2000)
    ObjList <-  PrepSCTIntegration(object.list = ObjList, anchor.features = features)
    anchors <-  FindIntegrationAnchors(object.list = ObjList, anchor.features = features, normalization.method = "SCT")
    combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
}
  
  # Plot data metrics after filtering
  PlotFeature(ObjList, Prefix = paste0(OutPrefix,"AfterFilter_",Sys.Date()))

  if(NormalizationMethod == "SCT"){pc = 15}else{pc = 30}
  combined <- ScaleData(combined) %>% RunPCA(npcs = pc) %>% RunUMAP(reduction = "pca", dims = 1:pc) %>% 
    FindNeighbors(reduction = "pca", dims = 1:pc) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  clustree(combined)
  ggsave(paste0(OutPrefix,"clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
  
  Idents(combined) <- combined$integrated_snn_res.0.1
  pdf(paste0(OutPrefix,"AllFeatureQC_",Sys.Date(),".pdf"), width = 12, height = 12)
  p <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1); plot(p); dev.off()
  saveRDS(combined, file = paste0(OutPrefix,"Seurat.rds"))
  
  # Visualization of the clusters landscape
  DimPlot(combined, reduction = "umap", group.by = c("ident","orig.ident"), label = T, repel = T, label.box =T, cols = hughie_color, label.color = "white")
  ggsave(paste0(OutPrefix,"UMAP-ByClusters_",Sys.Date(),".pdf"), width = 16, height = 8)
  
  #DimPlot(combined, reduction = "umap", split.by = "orig.ident", ncol = 2, label = T, repel = T, cols = hughie_color)
  #ggsave(paste0(OutPrefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 12, height = 12)
  
  #Plot_Cell_compoistion(Seurat_Obj = combined, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = OutPrefix)
  #MarkersPlot(Seurat_Obj = combined, ColorUse = hughie_color, OutPrefix = OutPrefix, NMarkers = 10)
  
  # future::plan("multisession", workers = 4)
  return(combined)
}

##' Function to caculate cell type module score based on marker genes
##' @param Seurat_Obj A Seurat object
##' @param Seruat_DEG_file DEG file for reference cell types, returned by the FindAllMarkers function
##' @param NtopGene Number of top genes (ranked by log2 FC) used as a module (Default: 50)
##' @param Ylimits Range for Y-axis on plotting
##' @param convertGeneName convert gene names when comparing different species, H2M indicate human to mouse, etc
##' @param ColorUse Color vector 
##' @param OutPrefix Prefix for output file
##' @examples
##'
##' Calculate_moduleScore(Seurat_Obj = subset(iBAT_Ad_combined, subset = orig.ident == "iBAT-2mo"), ColorUse = iBAT_Ad_colors, NtopGene = 50, Ylimits = c(0, 1),
##'                       Seruat_DEG_file = "../1. rATM_126mo/rATM_Adipocytes_2mo_Allmarkers_2023-07-05.csv", OutPrefix = paste0(iBAT_Ad_prefix, "2mo_mPRAT_"))
Calculate_moduleScore <- function(Seurat_Obj, Seruat_DEG_file, NtopGene = 50, OutPrefix, Ylimits = c(-2,2), convertGeneName = FALSE, ColorUse = hughie_color){
  DEG_matrix <- read.table(Seruat_DEG_file, sep = ",", header = T, row.names = 1)
  DEG_matrix1 <- DEG_matrix %>% group_by(cluster) %>% slice_max(n = NtopGene, order_by = avg_log2FC)
  DEG_list <- split(DEG_matrix1$gene, f = DEG_matrix1$cluster)
  
  if(convertGeneName == "H2M"){DEG_list <- lapply(DEG_list, HumanGeneToMouseGene)}
  if(convertGeneName == "M2H"){DEG_list <- lapply(DEG_list, MouseGeneToHumanGene)}
  
  pdf(paste0(OutPrefix,"UMAP-ModuleScore_",NtopGene,"TopMarkers_",Sys.Date(),".pdf"), width = 8, height = 10)
  for(CellTypes in names(DEG_list)){
    Seurat_Obj <- AddModuleScore_UCell(Seurat_Obj, features = list(DEG_list[[CellTypes]]), name = paste0(make.names(CellTypes),"_score"), assay = "RNA")
    p <- FeaturePlot(Seurat_Obj, features = paste0("signature_1",make.names(CellTypes),"_score"), min.cutoff = "q10", order = T, pt.size=0.3) + 
      scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
    plot(p)
  }
  p1 <- VlnPlot(Seurat_Obj, features = paste0("signature_1",make.names(names(DEG_list)),"_score"), split.by = "orig.ident", stack = T, flip = T,
                fill.by = "ident", cols = hughie_color) + FontSize(x.text = 16, y.text = 16) + scale_y_continuous(limits = Ylimits)
  p2 <- VlnPlot(Seurat_Obj, features = paste0("signature_1",make.names(names(DEG_list)),"_score"), stack = T, flip = T,
                fill.by = "ident", cols = hughie_color) + FontSize(x.text = 16, y.text = 16) + scale_y_continuous(limits = Ylimits)
  plot(p1); plot(p2)
  dev.off()
  
  Seurat_Obj$activeId <- Idents(Seurat_Obj)
  sigMat <- Seurat_Obj@meta.data[,c("activeId",paste0("signature_1",make.names(names(DEG_list)),"_score"))]
  sigMat <- reshape2::melt(sigMat,c("activeId"))
  sigMat$variable <- gsub("signature_1","",sigMat$variable)
  
  p3 <- ggplot(sigMat, aes(x = activeId, y = value)) + 
    geom_violin(aes(fill = activeId),width=1, linewidth =0.1) +
    geom_boxplot(color="black", alpha=0.2, width=0.1, linewidth =0.1, outlier.shape = NA) +
    stat_compare_means(method = "t.test", comparisons = combn(as.vector(unique(sigMat$activeId)), 2, simplify = FALSE)) +
    # stat_compare_means(label.y = 0.6) +
    facet_wrap(~variable, ncol=6, scales = "free_x") +
    scale_fill_manual(values = ColorUse) +
    expand_limits(y=0) +
    labs(x = "", y = "Module score") +
    theme_bw() +
    theme(
      plot.title = element_text(color="black", size=10, face="bold"),
      axis.title.x = element_text(color="black", size=10, face="bold"),
      axis.text.x = element_text(color="black", size=10, face="bold", angle = 45, hjust = 1),
      axis.title.y = element_text(color="black", size=10, face="bold"),
      axis.text.y = element_text(color="black", size=10, face="bold"),
      legend.title = element_text(color="black", size=10, face="bold"), 
      legend.text = element_text(color="black", size=10, face="bold"),
      strip.text.x = element_text(size = 8, face="bold")
    )
  plot(p3)
  ggsave(paste0(OutPrefix,"UMAP-ModuleScore_",NtopGene,"TopMarkers_Integrated",Sys.Date(),".pdf"), width = 14, height = 10)

  return(Seurat_Obj)
}

##' A wrapper to do DEG and plotting
##' @param Seurat_Obj A Seurat object
##' @param NMarkers Number of markers for violin and heatmap
##' @param ColorUse List of colors for plot
##' @param OutPrefix Prefix for output file
##' @examples
##'
##' MarkersPlot(Seurat_Obj = iWAT_combined, ColorUse = iWAT_HFD_color_maps, OutPrefix = "_iWAT_combined_",NMarkers = 20)
MarkersPlot <- function(Seurat_Obj, OutPrefix, NMarkers = 20, ColorUse = hughie_color, logfc.threshold = 0.25, min.pct = 0.25, only.pos = TRUE){
  
  markers <- FindAllMarkers(Seurat_Obj, only.pos = only.pos, logfc.threshold = logfc.threshold, min.pct = min.pct) %>% 
    filter(p_val < 0.05) %>% filter(!str_detect(gene,"^mt-"))
  write.csv(markers, file = paste0(OutPrefix, "Allmarkers_",Sys.Date(),".csv"), quote = F)
  
  TopMarkersList <- markers %>% group_by(cluster) %>% slice_max(n = NMarkers, order_by = avg_log2FC) %>% 
    arrange(match(cluster, levels(Idents(Seurat_Obj)))) %>% as.data.frame()
  
  # Plot all markers
  VlnPlot(Seurat_Obj, features = unique(TopMarkersList$gene), stack=T, flip=T, cols = ColorUse, fill.by = "ident") + NoLegend()
  ggsave(paste0(OutPrefix,"Allmarkers_Vlnplot_",Sys.Date(),".pdf"), width = 12, height = 45)
  
  # DotPlot(Seurat_Obj, features = unique(TopMarkersList$gene), split.by = "orig.ident", cols = hughie_color, scale = T) + RotatedAxis()
  # ggsave(paste0(OutPrefix,  "Allmarkers_Dotplot_",Sys.Date(),".pdf"), width = 16, height = 10)
  
  DoHeatmap(Seurat_Obj, features = unique(TopMarkersList$gene), group.colors = ColorUse)
  ggsave(paste0(OutPrefix, "Allmarkers_Heatmap_",Sys.Date(),".pdf"), width = 16, height = 20)
  return(TopMarkersList)
}

##' A wrapper to plot selected features heatmap
##' @param Seurat_Obj A Seurat object
##' @param features Selected genes
##' @param ColorUse List of colors for plot
##' @param OutPrefix Prefix for output file
FeatureHeatmapPlot <- function(Seurat_Obj, OutPrefix, features,cols, group.by){
  #browser()

  mycol2<- colorRamp2(c(-2, 0, 2), c("#0da9ce", "white", "#e74a32"))
  
  #各亚群平均表达量提取
  aver_dt<- AverageExpression(Seurat_Obj,
                              features= features,
                              group.by = group.by,
                              slot= 'data')
  aver_dt<- as.data.frame(aver_dt$SCT)
  
  #归一化：
  aver_dtt<- t(scale(t(aver_dt)))
  
  #添加行列注释：
  cell_anno <- data.frame(cell_anno = colnames(aver_dt))
  
  cols<- cols
  #cols<- c("#782AB6", "#3283FE", "#16FF32", "#85660D", "#565656", "#F7E1A0", 
  #         "#1C8356", "#AA0DFE")
  
  names(cols) <- cell_anno$cell_anno

  #列注释：
  cell<- data.frame(colnames(aver_dtt))
  colnames(cell) <- 'cell'
  col_anno<- HeatmapAnnotation(df = cell,
                               show_annotation_name= F,
                               gp= gpar(col = 'white', lwd = 2),
                               col= list(cell = cols))
  
  #热图绘制：
  pdf(paste0(OutPrefix, "SelectedMarkers_Heatmap_",Sys.Date(),".pdf"), width = 8, height = 12)
  print(Heatmap(aver_dtt,
          name= 'Expression',
          col= mycol2,
          cluster_columns= F,
          cluster_rows= F,
          show_column_names = T,
          #        column_names_side= c('top'),
          #        column_names_rot= 60,
          row_names_gp= gpar(fontsize = 12, fontface = 'italic'),
          rect_gp= gpar(col = "white", lwd = 1.5),
          #top_annotation= col_anno
  )) #+ row_anno
  dev.off()
 
}


##' Function to plot cell quality metrice from a Seurat Obj
##' @param Seurat_Obj A Seurat object
##' @param ColorUse Color vector 
##' @param OutPrefix Prefix for output file
PlotObjMetrices <- function(Seurat_Obj, OutPrefix, ColorUse = hughie_color){
  met <- Seurat_Obj@meta.data[c("orig.ident","nCount_RNA","percent.mt","nFeature_RNA")] %>% reshape2::melt(c("orig.ident"))
  met$variable <- factor(met$variable, levels = c("nCount_RNA","nFeature_RNA","percent.mt" ))
  
  p <- ggplot(met, aes(x = orig.ident, y = value)) + 
    geom_violin(aes(fill = orig.ident), width=1, linewidth=0.1) +
    geom_boxplot(color="black", alpha=0.1, width=0.1, outlier.colour = NA, linewidth=0.1) +
    facet_wrap(~variable, ncol=6, scales = "free") +
    # facetted_pos_scales(y = list(variable == "nCount_RNA" ~ scale_y_continuous(limits = c(0, 25000),  breaks = seq(0, 25000, 5000)),
    #                              variable == "nFeature_RNA" ~ scale_y_continuous(limits = c(0, 5000),  breaks = seq(0, 5000, 1000)),
    #                              variable == "percent.mt" ~ scale_y_continuous(limits = c(0, 10),  breaks = seq(0, 10, 2)))) +
    scale_fill_manual(values = ColorUse) +
    labs(x = "", y = "") + theme_bw() +
    theme(
      axis.title.x = element_text(color="black", size=8, face="bold"),
      axis.text.x = element_text(color="black", size=8, face="bold", angle = 45, hjust = 1),
      axis.title.y = element_text(color="black", size=8, face="bold"),
      axis.text.y = element_text(color="black", size=8, face="bold"),
      strip.text.x = element_text(size = 8, face="bold")
    )
  plot(p)
  ggsave(paste0(OutPrefix,"snRNA-seq_metrics.pdf"), width = 8, height = 4)
}

##' Function to plot cell quality metrice from raw Cellranger outputs
##' @param dir directory contain all cellranger outputs
##' @examples
##'
##' PlotAllMetrics(dir = "../_rawdata/")
##' PlotAllMetrics(dir = "../../3_ALK7_project/_rawdata/")
##' PlotAllMetrics(dir = "../../4_Gut_mWAT snRNA-seq/_rawdata/")
##' PlotAllMetrics(dir = "../_rawdata/_Arch/")
PlotAllMetrics <- function(dir){
  Res <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(Res) <- c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt" )
  files <- list.files(pattern = "CellBender_feature_bc_matrix_filtered.h5", path = dir, recursive = T)
  for (file in files){
    data <- ReadCB_h5(paste0(dir,"/",file), use.names = TRUE)
    data_Obj <- CreateSeuratObject(counts = data, project = gsub("CellBender_feature_bc_matrix_filtered.h5|/","",file), min.cells = 10)
    data_Obj[["percent.mt"]] <- PercentageFeatureSet(data_Obj, pattern = "^Mt-|^mt-")
    Res <- rbind(Res,data_Obj@meta.data)
  }
  
  Res2 <- reshape2::melt(Res, "orig.ident")
  
  pdf(paste0(dir,"_snRNA-seq_metrics.pdf"), width = 13, height = 13)
  p <- ggviolin(Res2, x = "orig.ident", y = "value", fill = "orig.ident",
                palette = hughie_color, size=0.1, draw_quantiles = T, 
                add = c("boxplot"), add.params = list(fill = "white")) +
    labs(x = "", y = "") + 
    facet_wrap(~variable, ncol=1, scales = "free_y") +
    theme_bw() +
    theme(
      axis.text.x = element_text(color="black", size=11, face="bold", angle = 60,hjust = 1),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      axis.text.y = element_text(color="black", size=11, face="bold"),
      legend.title = element_text(color="black", size=14, face="bold"), 
      legend.text = element_text(color="black", size=12, face="bold"),
      strip.text.x = element_text(size = 12, face="bold")
    )
  plot(p)
  dev.off()
}

##' Read h5 files
##' @param filename h5 filename
##' @examples
##'
ReadCB_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
  genomes <- names(x = infile)
  output <- list()
  if (hdf5r::existsGroup(infile, 'matrix')) {
    # cellranger version 3
    message('CellRanger version 3+ format H5')
    if (use.names) {
      feature_slot <- 'features/name'
    } else {
      feature_slot <- 'features/id'
    }
  } else {
    message('CellRanger version 2 format H5')
    if (use.names) {
      feature_slot <- 'gene_names'
    } else {
      feature_slot <- 'genes'
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    features <- infile[[paste0(genome, '/', feature_slot)]][]
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    sparse.mat <- sparseMatrix(
      i = indices[] + 1,
      p = indptr[],
      x = as.numeric(x = counts[]),
      dims = shp[],
      repr="T"
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = 'CsparseMatrix')
    # Split v3 multimodal
    if (infile$exists(name = paste0(genome, '/features'))) {
      types <- infile[[paste0(genome, '/features/feature_type')]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(
          X = types.unique,
          FUN = function(x) {
            return(sparse.mat[which(x = types == x), ])
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(x = output) == 1) {
    return(output[[genome]])
  } else{
    return(output)
  }
}

##' Print cell composition
GetCellComposition <- function(Seurat_Obj){
  CompositionFreq <- table(Seurat_Obj$orig.ident, Seurat_Obj@active.ident) %>% as.data.frame() 
  mat <- CompositionFreq %>% group_by(Var1) %>% mutate(Freq = Freq/sum(Freq) * 100) 
  mat$Var2 <- factor(mat$Var2, levels = rev(levels(mat$Var2)))
  colnames(mat) <- c("Sample","CellType","Perc(%)")
  print(as.data.frame(mat),digits=2)
  # gridExtra::grid.table(mat)
}

##' Plot multiple marker genes expression
##' @param Seurat_Obj A Seurat object
##' @param ColorUse1 List of colors for cell types
##' @param ColorUse2 List of colors for sample types
##' @param OutPrefix Prefix for output file
##' @examples
##'
##' Plot_Cell_compoistion(Seurat_Obj = iWAT_combined, OutPrefix = "_iWAT_combined_", ColorUse = iWAT_HFD_color_maps)
Plot_Cell_compoistion <- function(Seurat_Obj, ColorUse1, ColorUse2, OutPrefix, PlotBySample = F){
  
  CompositionFreq <- table(Seurat_Obj$orig.ident, Seurat_Obj@active.ident) %>% as.data.frame() 
  
  mat <- CompositionFreq %>% group_by(Var1) %>% mutate(Freq = Freq/sum(Freq) * 100) 
  mat$Var2 <- factor(mat$Var2, levels = levels(mat$Var2))
  write.csv(mat,paste0(OutPrefix, "cellNumbersPercentage_",Sys.Date(),".csv"))
  p2 <- ggplot(mat,aes(x=Var1, y=Freq, fill = Var2, label = paste0(round(Freq, 1), "%"))) + 
    geom_bar(stat="identity", width = 0.7, color = "black",
             position = position_stack(reverse = T)) +
    geom_text(position=position_stack(vjust=0.5, reverse = T), color="white") + 
    # scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = ColorUse1) + 
    coord_flip() +  scale_x_discrete(limits = rev(levels(mat$Var1))) + theme_bw() + 
    labs(x="", y="Percentage of cell types") +
    theme(
      axis.title.x = element_text(color="black", size=16, face="bold"),
      axis.text.x = element_text(color="black", size=16, face="bold"),
      axis.title.y = element_text(color="black", size=16, face="bold"),
      axis.text.y = element_text(color="black", size=16, face="bold")
    )
  plot(p2)
  ggsave(paste0(OutPrefix, "cellNumbersPerCluster_filled_",Sys.Date(),".pdf"), width = 10, height = 5)
  #browser()
  p4<-ggplot(mat,aes(x="", y=Freq, fill = Var2, label = paste0(round(Freq, 1), "%"))) + 
    geom_bar(stat="identity", width = 1, color = "black") +
    coord_polar("y", start = 0) +
    geom_text(position=position_stack(vjust=0.5), color="white",size=5) + 
    scale_fill_manual(values = ColorUse1) + 
    facet_wrap(~Var1) +
    theme_void() + 
    labs(x="", y="Percentage of cell types") +
    theme(
      axis.title.x = element_text(color="black", size=16, face="bold"),
      axis.text.x=element_blank(),
      axis.title.y = element_text(color="black", size=16, face="bold"),
      axis.text.y = element_text(color="black", size=16, face="bold")
    )
  plot(p4)
  ggsave(paste0(OutPrefix, "cellNumbersPerCluster_Pie_",Sys.Date(),".pdf"), width = 8, height = 4)
  
  NCellType <- length(unique(Seurat_Obj@active.ident))
  Quan <- table(Seurat_Obj$orig.ident,Seurat_Obj@active.ident) %>% as.data.frame.matrix()
  Quan$Sum <- rowSums(Quan)
  Quan <- Quan*100/Quan$Sum
  Quan <- Quan[,-ncol(Quan)]
  Quan$Samples <- rownames(Quan)
  Quan <- reshape2::melt(Quan,c("Samples"))
  Quan$Samples <- factor( Quan$Samples, levels = unique(Quan$Samples))
  p3 <- ggplot(Quan,aes(x=variable, y = value, width=0.7)) +
    geom_bar(aes(fill=Samples),color="black",stat="identity", position=position_dodge()) +
    scale_fill_manual(values = ColorUse2) +
    # geom_text(aes(label = round(value,2)), size = 5) +
    theme_bw() +
    labs(x="", y="Percentage of cell types") +
    theme(
      axis.title.x = element_text(color="black", size=16, face="bold"),
      axis.text.x = element_text(color="black", size=14, face="bold", angle = 45, hjust = 1),
      axis.title.y = element_text(color="black", size=16, face="bold"),
      axis.text.y = element_text(color="black", size=14, face="bold")
    )
  plot(p3)
  ggsave(paste0(OutPrefix,"cellNumbersPerCluster_dodged_",Sys.Date(),".pdf"), width = 8, height = 6)
  
  if(PlotBySample == T){
    Quan$Samples <- factor(Quan$Samples, levels = unique(Quan$Samples))
    p3 <- ggplot(Quan,aes(x=Samples, y = value, width=0.7)) +
      geom_bar(aes(fill=variable),color="black",stat="identity", position=position_dodge()) +
      scale_fill_manual(values = ColorUse1) +
      theme_bw() +
      labs(x="", y="Percentage of cell types") +
      theme(
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.text.x = element_text(color="black", size=14, face="bold", angle = 45, hjust = 1),
        axis.title.y = element_text(color="black", size=16, face="bold"),
        axis.text.y = element_text(color="black", size=14, face="bold")
      )
    plot(p3)
    ggsave(paste0(OutPrefix,"cellNumbersPerCluster_dodged_BySample_",Sys.Date(),".pdf"), width = 8, height = 6)
  }
}

convert_human_seu_to_mouse <- function(seu, ...) {
  new_rownames <- convert_symbols_by_species(src_genes = rownames(seu), src_species = "human")
  # seu_slots <- c("counts", "data", "scale.data", "meta.features")
  if(seu@active.assay == "RNA"){
    rownames(seu@assays$RNA@meta.features) <- new_rownames
    rownames(seu@assays$RNA@counts) <- new_rownames
    rownames(seu@assays$RNA@data) <- new_rownames
  }
  if(seu@active.assay == "integrated"){
    rownames(seu@assays$integrated@meta.features) <- new_rownames
    rownames(seu@assays$integrated@data) <- new_rownames
    rownames(seu@assays$integrated@scale.data) <- new_rownames
  }
  return(seu)
}

convert_symbols_by_species <- function(src_genes, src_species) {
  if(file.exists("../human_to_mouse_homologs.rds")){
    human_to_mouse_homologs <- readRDS("../human_to_mouse_homologs.rds")
  } else {
    #Generate human_mouse_homologs
    human_mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
    mouse_mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
    human_to_mouse_homologs = biomaRt::getLDS(attributes = c("hgnc_symbol","entrezgene_id","ensembl_gene_id"),
                                              # filters = "hgnc_symbol", values = rownames(human_gene_transcript_seu), filters = "entrezgene_id", values = hs_entrez,
                                              mart = human_mart, attributesL = c("mgi_symbol","entrezgene_id","ensembl_gene_id"), martL = mouse_mart)
    human_to_mouse_homologs <- human_to_mouse_homologs[c("HGNC.symbol", "MGI.symbol")]
    # write_rds(human_to_mouse_homologs,"human_to_mouse_homologs.rds")
  }
  
  if (src_species == "human") {
    dest_species <- "mouse"
    dest_symbols <- src_genes %>% tibble::enframe("gene_index", "HGNC.symbol") %>%
      dplyr::left_join(human_to_mouse_homologs, by = "HGNC.symbol") %>% dplyr::distinct(HGNC.symbol, .keep_all = TRUE) %>%
      dplyr::mutate(MGI.symbol = dplyr::case_when(is.na(MGI.symbol) ~ stringr::str_to_sentence(HGNC.symbol), TRUE ~ MGI.symbol)) %>% 
      dplyr::select(-gene_index) %>% identity()
  } else if (src_species == "mouse") {
    dest_species <- "human"
    dest_symbols <- src_genes %>% tibble::enframe("gene_index", "MGI.symbol") %>%
      dplyr::left_join(human_to_mouse_homologs, by = "MGI.symbol") %>% dplyr::distinct(MGI.symbol, .keep_all = TRUE) %>%
      dplyr::mutate(HGNC.symbol = dplyr::case_when( is.na(HGNC.symbol) ~ stringr::str_to_upper(MGI.symbol),TRUE ~ HGNC.symbol)) %>%
      dplyr::select(-gene_index) %>% identity()
    # dplyr::mutate(HGNC.symbol = make.unique(HGNC.symbol)) %>%
  }
  return(make.unique(dest_symbols[[2]]))
}


##' Transform the Human and Mouse gene names
##' @param GenesNames A vector of Human gene names
##' @examples
##'
##' write.csv(convertHumanGeneList(cc.genes$g2m.genes),"Mouse_G2M-phase_genes.csv")
##' write.csv(convertHumanGeneList(cc.genes$s.genes),"Mouse_S-phase_genes.csv")
convertHumanGeneList <- function(GenesToBeConverted){
  require("biomaRt")
  human = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  mouse = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
  # human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  # mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = GenesToBeConverted , mart = human, 
                   attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  mousex <- unique(genesV2[, 2]) # unique the 2nd column values
  if(length(GenesToBeConverted) != length(mousex)){
    genes_not_trans <- setdiff(GenesToBeConverted, genesV2$HGNC.symbol)
    print("These genes could not be translated:")
    print(genes_not_trans)
    print(paste("A total number of ",length(genes_not_trans),"genes could not be translated!"), sep=" ")
  }else{
    print("All genes were translated successfully!")
  }
  return(mousex)
}
minMax <- function(x) {(x - min(x)) / (max(x) - min(x))}
##' @examples
##' MouseGeneToHumanGene(head(keys(org.Mm.eg.db, "SYMBOL"), 40))
##' HumanGeneToMouseGene(head(keys(org.Hs.eg.db, "SYMBOL"), 40))
MouseGeneToHumanGene <- function(mousegenes){
  gns <- mapIds(org.Mm.eg.db, mousegenes, "ENTREZID", "SYMBOL")
  mapped <- select(Orthology.eg.db, gns, "Homo_sapiens","Mus_musculus")
  naind <- is.na(mapped$Homo_sapiens)
  hsymb <- mapIds(org.Hs.eg.db, as.character(mapped$Homo_sapiens[!naind]), "SYMBOL", "ENTREZID")
  out <- data.frame(Mouse_symbol = mousegenes, mapped, Human_symbol = NA)
  out$Human_symbol[!naind] <- hsymb
  na.omit(out$Human_symbol) %>% as.vector()
}

MouseGeneToHumanGene2 <- function(mousegenes){
  gns <- mapIds(org.Mm.eg.db, mousegenes, "ENTREZID", "SYMBOL")
  mapped <- select(Orthology.eg.db, gns, "Homo_sapiens","Mus_musculus")
  naind <- is.na(mapped$Homo_sapiens)
  hsymb <- mapIds(org.Hs.eg.db, as.character(mapped$Homo_sapiens[!naind]), "SYMBOL", "ENTREZID")
  out <- data.frame(Mouse_symbol = mousegenes, mapped, Human_symbol = NA)
  out$Human_symbol[!naind] <- hsymb
  na.omit(out)
}

HumanGeneToMouseGene <- function(Humangenes){
  gns <- mapIds(org.Hs.eg.db, Humangenes, "ENTREZID", "SYMBOL")
  mapped <- select(Orthology.eg.db, gns, "Mus_musculus","Homo_sapiens")
  naind <- is.na(mapped$Mus_musculus)
  mmymb <- mapIds(org.Mm.eg.db, as.character(mapped$Mus_musculus[!naind]), "SYMBOL", "ENTREZID")
  out <- data.frame(Human_symbol = Humangenes, mapped, Mouse_symbol = NA)
  out$Mouse_symbol[!naind] <- mmymb
  na.omit(out$Mouse_symbol) %>% as.vector()
}

RatGeneToMouseGene <- function(Ratgenes){
  gns <- mapIds(org.Rn.eg.db, Ratgenes, "ENTREZID", "SYMBOL")
  mapped <- select(Orthology.eg.db, gns, "Mus_musculus","Rattus_norvegicus")
  naind <- is.na(mapped$Mus_musculus)
  mmymb <- mapIds(org.Mm.eg.db, as.character(mapped$Mus_musculus[!naind]), "SYMBOL", "ENTREZID")
  out <- data.frame(Rat_symbol = Ratgenes, mapped, Mouse_symbol = NA)
  out$Mouse_symbol[!naind] <- mmymb
  na.omit(out$Mouse_symbol) %>% as.vector()
}

##' Enrichment analysis for markers and plot
##' @param Seruat_DEG_file A Seurat DEG file from FindAllMarkers()
##' @param showCategoryNum How many GO/KEGG items to show?
##' @param cellRanks Vector ranking cell types to plot 
##' @param filterLevel GO levels for plot
##' @param species Specify species for analysis (mosue, rat, human)
##' @examples
##'
##' DEG_enrichment(Seruat_DEG_file = "2022-10-05_LI_merged_res0.3_markers.csv", showCategoryNum = 10, filterLevel = 4, cellRanks = LI_ranks)
DEG_enrichment <- function(Seruat_DEG_file, showCategoryNum = 10, cellRanks = NULL, filterLevel = 3, species = "mouse"){
  #browser()
  DEG_matrix <- read.table(Seruat_DEG_file, sep = ",", header = T, row.names = 1)
  DEG_list <- split(rownames(DEG_matrix), f = DEG_matrix$cluster)
  Result_suffix <- tools::file_path_sans_ext(Seruat_DEG_file)
  
  if(species == "mouse"){
    annoDB <- org.Mm.eg.db
    organism <- "mmu"
  } else if(species == "rat"){
    annoDB <- org.Rn.eg.db
    organism <- "rno"
  } else if(species == "human"){
    annoDB <- org.Hs.eg.db
    organism <- "hsa"
  }
  
  #Transfer gene symbol into entrez id
  GeneID_list_ID <- DEG_list %>% map(~{
    gene.df <- AnnotationDbi::select(annoDB, keys = .x, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
    gene <- gene.df$ENTREZID
    gene <- gene[which(!is.na(gene))]
    gene <- unique(gene)
    return(gene)
  })
  
  #GO Analysis
  compGO <- compareCluster(geneCluster = DEG_list, fun= "enrichGO", pvalueCutoff = 0.05, keyType = "SYMBOL",
                           pAdjustMethod = "BH", OrgDb = annoDB, readable = T, ont = "BP")
  saveRDS(compGO, paste0(Result_suffix,"_comp_",Sys.Date(),".rds") )
  write.csv(compGO@compareClusterResult, file = paste0(Result_suffix,"_GO-BP_list_",Sys.Date(),".csv"))
  
  if(!is.null(cellRanks)){compGO@compareClusterResult$Cluster <- factor(compGO@compareClusterResult$Cluster, levels = cellRanks)}
  compGO_filtered <- gofilter(compGO, level = filterLevel)
  dotplot(compGO_filtered, showCategory = showCategoryNum,includeAll=T,title = paste0(Result_suffix,"\n - GO-BP Enrichment Analysis")) +  
    scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(Result_suffix,"_GO-BP_dotplot_Filter",filterLevel,"_",showCategoryNum,"-Items_",Sys.Date(),".pdf"), width = 14, height = 30)
  
  # KEGG analysis
  compKEGG <- compareCluster(geneCluster = GeneID_list_ID, fun = "enrichKEGG", pvalueCutoff = 0.05, pAdjustMethod = "BH", organism = organism)
  write.csv(compKEGG@compareClusterResult, file = paste0(Result_suffix,"_KEGG_list_",Sys.Date(),".csv"))
  if(!is.null(cellRanks)){compKEGG@compareClusterResult$Cluster <- factor(compKEGG@compareClusterResult$Cluster, levels = cellRanks)}
  
  dotplot(compKEGG, showCategory = showCategoryNum, includeAll=T, title = paste0(Result_suffix,"\n - KEGG Pathway Enrichment Analysis")) +  
    scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(Result_suffix,"_KEGG_dotplot_",showCategoryNum,"-Items_",Sys.Date(),".pdf"), width = 14, height = 30)
}

##' Enrichment analysis for markers in different conditions and plot
##' @param Seurat_Obj A Seurat object
##' @param CellTypeInclude Specify cell type that included for DEG analysis
##' @param PlotDEGHeatmap Whether to plot DEG for each cell type using heatmap?
##' @param OutPrefix Prefix for output file
##' @examples
##'
##' PairWise_DEG_enrichment(Seurat_obj = iWAT_combined, OutPrefix = "_iWAT_combined_", PlotDEGHeatmap = F,RunGO = F,
##'              CellTypeInclude = as.vector(unique(iWAT_combined$Myannotation)[c(5,6)]))
PairWise_DEG_enrichment <- function(Seurat_obj, CellTypeInclude, PlotDEGHeatmap = F, RunGO = F, OutPrefix){
  
  Sample_names <- unique(Seurat_obj$orig.ident)
  cell_names <- as.vector(unique(Idents(Seurat_obj)))
  Seurat_obj$CT_SAMPLE <- paste(Idents(Seurat_obj), Seurat_obj$orig.ident, sep = "_")
  Seurat_obj$celltype <- Idents(Seurat_obj)
  Idents(Seurat_obj) <- "CT_SAMPLE"
  
  #List to store all pair-wise DEGs
  DEG_genes <- list()
  
  for (cell_name in intersect(cell_names, CellTypeInclude)){
    cat("Processing",cell_name,"...\n")
    A_comp <- paste0(cell_name,"_", Sample_names[1])
    B_comp <- paste0(cell_name,"_", Sample_names[2])
    
    Seurat_obj_DEG <- FindMarkers(Seurat_obj, ident.1 = A_comp, ident.2 = B_comp, min.pct = 0.1, logfc.threshold = 0.25) 
    Seurat_obj_DEG$change <- ifelse(Seurat_obj_DEG$p_val_adj < 0.05 & abs(Seurat_obj_DEG$avg_log2FC) >= 0.20, ifelse(Seurat_obj_DEG$avg_log2FC > 0.20 ,'Up','Down'),'Stable')
    
    write.csv(Seurat_obj_DEG, paste0(OutPrefix,cell_name,"_",Sample_names[1],"-vs-",Sample_names[2],"_DEG_List_",Sys.Date(),".csv"))
    
    Seurat_obj_DEG_up <- row.names(Seurat_obj_DEG[Seurat_obj_DEG$avg_log2FC > 0 & Seurat_obj_DEG$p_val_adj < 0.05,]) 
    Seurat_obj_DEG_down <- row.names(Seurat_obj_DEG[Seurat_obj_DEG$avg_log2FC < 0 & Seurat_obj_DEG$p_val_adj < 0.05,]) 
    
    DEG_genes[[paste0(A_comp,"-Up")]] <- Seurat_obj_DEG_up
    DEG_genes[[paste0(A_comp,"-Down")]] <-Seurat_obj_DEG_down
    
    Seurat_obj@active.ident <- as.factor(Seurat_obj$CT_SAMPLE)
    Seurat_obj_averge <- AverageExpression(Seurat_obj)
    
    ## Extract the average value of the expression of each sample and draw heatmap
    if(PlotDEGHeatmap){
      Plot_expr <- Seurat_obj_averge$RNA[c(Seurat_obj_DEG_up,Seurat_obj_DEG_down),c(A_comp,B_comp)]
      pdf(paste0(OutPrefix,cell_name,"_",Sample_names[1],"-vs-",Sample_names[2],"_DEG_heatmap_",Sys.Date(),".pdf"), width = 5, height = 10)
      p <- pheatmap(Plot_expr, show_colnames = T, show_rownames = T, scale = "row")
      plot(p)
      dev.off()
    }
    if(RunGO){
      ## GO enrichment analysis of up- and down-regulated genes
      Go_BP_up <- enrichGO(gene = Seurat_obj_DEG_up, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = 'BH', qvalueCutoff = 0.05)
      Go_BP_down <- enrichGO(gene = Seurat_obj_DEG_down, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = 'BH',  qvalueCutoff = 0.05)
      GO_BP_All <- enrichGO(gene = c(Seurat_obj_DEG_up,Seurat_obj_DEG_down), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = 'BH',  qvalueCutoff = 0.05)
      
      write.csv(Go_BP_up@result, paste0(OutPrefix,"_GO-BP analysis of ",cell_name," up DEG_GO_List_",Sys.Date(),".csv"))
      write.csv(Go_BP_down@result, paste0(OutPrefix,"_GO-BP analysis of ",cell_name," down DEG_GO_List_",Sys.Date(),".csv"))
      write.csv(GO_BP_All@result, paste0(OutPrefix,"_GO-BP analysis of ",cell_name," All DEG_GO_List_",Sys.Date(),".csv"))
      
      pdf(paste0(OutPrefix,"_GO-BP analysis of ", "_",cell_name, Sample_names[1]," and ",Sample_names[2], " DEG_barplot_",Sys.Date(),".pdf"),width = 18, height = 12)
      p1 <- barplot(Go_BP_up, showCategory = 10, title = paste0("GO-BP analysis of ",cell_name, "_",Sample_names[1]," and",Sample_names[2], " up DEG")) +  
        scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      p2 <- barplot(Go_BP_down, showCategory = 10, title  = paste0("GO-BP analysis of ",cell_name, "_", Sample_names[1]," and",Sample_names[2], " down DEG"))+  
        scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      p3 <- barplot(GO_BP_All, showCategory = 10, title  = paste0("GO-BP analysis of ",cell_name, "_", Sample_names[1]," and",Sample_names[2], " All DEG"))+  
        scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      plot(p1/p2/p3)
      dev.off()
    }
  }
  
  write.csv(sapply(names(DEG_genes), function(x) paste(x,paste(DEG_genes[[x]],collapse=","))), paste0(OutPrefix, "_AllDEG_List_",Sys.Date(),".csv"))

  AllSamples <- unique(Seurat_obj$CT_SAMPLE)[grep(paste0(CellTypeInclude,collapse = "|"),unique(Seurat_obj$CT_SAMPLE))]
  AllSamples <- AllSamples[order(AllSamples)]
  
  #Plot heatmap for all DEGs on all cell types
  Allgenes <- unique(unlist(DEG_genes[names(DEG_genes)], use.names = F))
  if("RNA" %in% names(Seurat_obj@assays)){Plot_All_expr <- Seurat_obj_averge$RNA[Allgenes, AllSamples]}
  if("Spatial" %in% names(Seurat_obj@assays)){Plot_All_expr <- Seurat_obj_averge$Spatial[Allgenes, AllSamples]}
  if(!all(is.finite(rowSums(Plot_All_expr)))){ Plot_All_expr <- Plot_All_expr[is.finite(rowSums(Plot_All_expr)),]}
  pdf(paste0(OutPrefix, "_DEG_heatmap_",Sys.Date(),".pdf"), width = 15, height = 15)
  p <- pheatmap(Plot_All_expr, show_colnames = T, show_rownames = T, scale = "row",cluster_rows = T, cluster_cols = F)
  plot(p)
  dev.off()
}

##' Enrichment plot for selected terms
##' @param GO_enrichment_file A csv file containing selected GO terms
##' @param OutPrefix OutPrefix
Selected_enrichment_Plot <- function(GO_enrichment_file, OutPrefix, levels = levels, cluster_order =cluster_order){
  #browser()
  df <- read.table(GO_enrichment_file,sep = ",", header = T, row.names = 1)
  df$GeneRatio <- parse_ratio(df$GeneRatio)
  df$Description <- factor(df$Description,levels = rev(levels))
  df$Cluster <- factor(df$Cluster, levels = cluster_order)
  
  ggplot(df, aes(x=Cluster, y = Description)) +
    geom_point(aes(colour = p.adjust, size = GeneRatio))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
    scale_colour_gradient2(low = "#e74a32",mid = "orange",high =  "yellow", midpoint = 0.025)+
    labs(x = NULL, u = NULL) + 
    guides(size = guide_legend(order = 1)) +
    theme(legend.direction = "horizontal", legend.position = "top") +
    scale_y_discrete(position = "right")
  ggsave(paste0(OutPrefix,"Selectedmarkers__GO_dotplot_",Sys.Date(),".pdf"), width = 8, height = 12)
  
}
  


##' Inference cell-cell communication on each sample using CellChat
##' @param Seurat_Obj A Seurat object
##' @param Seurat_MetaInfo Meta information of your Seurat object
##' @param ColorUse List of colors for plot
##' @param OutPrefix Prefix for output file
##' @examples
##' Run_CellChat(SeuratObj = iWAT_KO_HFD_data, Seurat_MetaInfo = iWAT_KO_HFD_meta, OutPrefix = "_iWAT_KO_HFD_CellChat_", ColorUse = iWAT_HFD_maps)
Run_CellChat <- function(SeuratObj, Seurat_MetaInfo, OutPrefix, ColorUse, species = "mouse"){
  
  if(file.exists(paste0(OutPrefix,".rds"))){
    cellchat <- readRDS(paste0(OutPrefix,".rds"))
  } else {
    cellchat <- createCellChat(object = SeuratObj, meta = Seurat_MetaInfo, group.by = "Anno")
    if(species == "mouse"){CellChatDB <- CellChatDB.mouse} else if(species == "human"){CellChatDB <- CellChatDB.human}
    
    # suppressPackageStartupMessages(library(openxlsx))
    # wb <- createWorkbook(creator = "AUTO")
    # addWorksheet(wb, sheetName = "Mouse_CellChat_interactions"); writeData(wb, sheet = "Mouse_CellChat_interactions", x = CellChatDB$interaction)
    # addWorksheet(wb, sheetName = "Mouse_CellChat_complex"); writeData(wb, sheet = "Mouse_CellChat_complex", x = CellChatDB$complex)
    # addWorksheet(wb, sheetName = "Mouse_CellChat_cofactor"); writeData(wb, sheet = "Mouse_CellChat_cofactor", x = CellChatDB$cofactor)
    # addWorksheet(wb, sheetName = "Mouse_CellChat_geneInfo"); writeData(wb, sheet = "Mouse_CellChat_geneInfo", x = CellChatDB$geneInfo)
    # saveWorkbook(wb, file = paste0("Mouse_CellChat.xlsx"), overwrite = TRUE)
    
    CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    # showDatabaseCategory(CellChatDB)
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan("multisession", workers = 20) # do parallel
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat) # project gene expression data onto PPI network (optional)
    cellchat <- projectData(cellchat, PPI.mouse)
    cellchat <- computeCommunProb(cellchat, type="triMean")  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10) #Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat, thresh = 0.05) #Calculate the aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    saveRDS(cellchat, paste0(OutPrefix,".rds"))
    future::plan("multisession", workers = 2) #Finish parallel
  }
  # cellchat <- updateCellChat(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  
  #Plot interactions and Weights
  pdf(paste0(OutPrefix, "_Interactions.pdf"), width = 16, height = 8)
  par(mfrow = c(1,3), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions", vertex.label.cex = 1.5,color.use = ColorUse)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 1.5,color.use = ColorUse)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F,  title.name = "Interaction strength", vertex.label.cex = 1.5,color.use = ColorUse)
  dev.off()
  
  pdf(paste0(OutPrefix, "_Interactions_Seperate.pdf"), width = 12, height = 12)
  mat <- cellchat@net$weight
  par(mfrow = c(3,4), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, label.edge = F, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = ColorUse)
  }
  dev.off()
  
  SigPathways <- cellchat@netP$pathways

  pdf(paste0(OutPrefix, "_PathwayInteractions_circle.pdf"), width = 12, height = 12)
  par(mfrow = c(4,4), xpd=TRUE)
  for (pathway in SigPathways){netVisual_aggregate(cellchat, signaling = pathway, layout = "circle", color.use = ColorUse)}
  dev.off()

  # pdf(paste0(OutPrefix, "_PathwayInteractions_chord_",Sys.Date(),".pdf"), width = 16, height = 16)
  # par(mfrow = c(4,4), xpd=TRUE)
  # for (pathway in SigPathways){netVisual_aggregate(cellchat, signaling = pathway, layout = "chord", color.use = ColorUse)}
  # dev.off()
  
  pdf(paste0(OutPrefix, "_PathwayInteractions_netAnalysis_contribution.pdf"), width = 10, height = 10)
  for (pathway in SigPathways){p <- netAnalysis_contribution(cellchat, signaling = pathway);plot(p)}
  dev.off()
  
  # show all the significant interactions (L-R pairs)
  pdf(paste0(OutPrefix,"_PathwayInteractions_netVisual_bubble.pdf"), width = 12, height = 10)
  p <- netVisual_bubble(cellchat, remove.isolate = FALSE, color.text.use = ColorUse);plot(p)
  dev.off()
  
  pdf(paste0(OutPrefix,"_PathwayInteractions_PathwayGeneExpression.pdf"), width = 12, height = 12)
  for (pathway in SigPathways){p <- plotGeneExpression(cellchat, signaling = pathway, color.use = ColorUse, enriched.only = F);plot(p)}
  dev.off()

  # Compute the network centrality scores
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
  pdf(paste0(OutPrefix,"_PathwayInteractions_netAnalysis_signalingRole_network.pdf"), width = 8, height = 6)
  for (pathway in SigPathways){netAnalysis_signalingRole_network(cellchat, signaling = pathway, width = 8, height = 2.5, font.size = 10, color.use = ColorUse)}
  dev.off()

  pdf(paste0(OutPrefix,"_PathwayInteractions_netAnalysis_signalingRole_scatter.pdf"), width = 6, height = 6)
  p <- netAnalysis_signalingRole_scatter(cellchat,color.use = ColorUse);plot(p)
  dev.off()
  
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  pdf(paste0(OutPrefix,"_PathwayInteractions_netAnalysis_signalingRole_heatmap.pdf"), width = 14, height = 10)
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.use = ColorUse, width = 10, height = 18, color.heatmap = "OrRd")
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.use = ColorUse, width = 10, height = 18, color.heatmap = "OrRd")
  draw(ht1 + ht2)
  dev.off()

  pdf(paste0(OutPrefix, "_PathwayInteractions_Pattern_outgoing.pdf"), width = 8, height = 6)
  # selectK(cellchat, pattern = "outgoing")
  nPatterns = 3
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  p1 <- netAnalysis_river(cellchat, pattern = "outgoing",color.use = ColorUse)
  p2 <- netAnalysis_dot(cellchat, pattern = "outgoing",color.use = ColorUse)
  plot(p1+p2)
  dev.off()
  
  pdf(paste0(OutPrefix, "_PathwayInteractions_Pattern_incoming.pdf"), width = 8, height = 6)
  # selectK(cellchat, pattern = "incoming")
  nPatterns = 3
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
  p1 <- netAnalysis_river(cellchat, pattern = "incoming",color.use = ColorUse)
  p2 <- netAnalysis_dot(cellchat, pattern = "incoming",color.use = ColorUse)
  plot(p1+p2)
  dev.off()
}

##' Comparison analysis of multiple datasets
##' @param CellChat_Obj_list A list of Seurat object
##' @param Seurat_MetaInfo Meta information of your Seurat object
##' @param ColorUse List of colors for plot
##' @param OutPrefix Prefix for output file
##' @examples
##' 
##' Run_CellChat_Compare(CellChat_Obj_list = iWAT_object.list, OutPrefix = "_iWAT_CellChat_Merged_", ColorUse = iWAT_HFD_maps)
Run_CellChat_Compare <- function(CellChat_Obj_list, OutPrefix, ColorUse, ColorUse_SAM, runShort = F){
  
  if(runShort){ #Run individual function
    SigPathways <- c()
    for (i in 1:length(CellChat_Obj_list)){SigPathways <- unique(c(SigPathways,CellChat_Obj_list[[i]]@netP$pathways))}
    cellchat <- mergeCellChat(CellChat_Obj_list, add.names = names(CellChat_Obj_list),cell.prefix = TRUE)
    
    pdf(paste0(OutPrefix, "_netVisual_circle_inte.pdf"), width = 10, height = 8)
    weight.max <- getMaxWeight(CellChat_Obj_list, attribute = c("idents","count"))
    par(mfrow = c(2,2), xpd=TRUE)
    for (i in 1:length(CellChat_Obj_list)) {
      # if(i %in% c(1,2)){ColorUse1 = ColorUse[-5]} else{ColorUse1 = ColorUse}
      netVisual_circle(CellChat_Obj_list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], vertex.label.cex = 0, edge.label.cex = 1.2,
                       edge.width.max = 5, title.name = paste0("Number of interactions - ", names(CellChat_Obj_list)[i]), color.use = ColorUse1)
    }
    dev.off()
    
    pdf(paste0(OutPrefix,"_netAnalysis_signalingRole_heatmap_inte.pdf"), width = 24, height = 12)
    par(mfrow = c(1,6), xpd=TRUE)
    W = 6
    ht1 <- netAnalysis_signalingRole_heatmap(CellChat_Obj_list[[1]], pattern = "all", signaling = SigPathways, color.use = ColorUse,
                                             title = names(CellChat_Obj_list)[1], width = W, height = 18, color.heatmap = "OrRd")
    ht2 <- netAnalysis_signalingRole_heatmap(CellChat_Obj_list[[2]], pattern = "all", signaling = SigPathways, color.use = ColorUse,
                                             title = names(CellChat_Obj_list)[2], width = W, height = 18, color.heatmap = "OrRd")
    ht3 <- netAnalysis_signalingRole_heatmap(CellChat_Obj_list[[3]], pattern = "all", signaling = SigPathways, color.use = ColorUse,
                                             title = names(CellChat_Obj_list)[3], width = W, height = 18, color.heatmap = "OrRd")
    ht4 <- netAnalysis_signalingRole_heatmap(CellChat_Obj_list[[4]], pattern = "all", signaling = SigPathways, color.use = ColorUse,
                                             title = names(CellChat_Obj_list)[4], width = W, height = 18, color.heatmap = "OrRd")
    # ht5 <- netAnalysis_signalingRole_heatmap(CellChat_Obj_list[[5]], pattern = "all", signaling = SigPathways, color.use = ColorUse,
    #                                          title = names(CellChat_Obj_list)[5], width = W, height = 18, color.heatmap = "OrRd")
    # ht6 <- netAnalysis_signalingRole_heatmap(CellChat_Obj_list[[6]], pattern = "all", signaling = SigPathways, color.use = ColorUse,
    #                                          title = names(CellChat_Obj_list)[6], width = W, height = 18, color.heatmap = "OrRd")
    draw(ht1 + ht2 + ht3 + ht4 , ht_gap = unit(1, "cm"))
    dev.off()
    
  } else{
    ListLen <- length(CellChat_Obj_list)
    combinations <- combn(ListLen,2)
    SigPathways <- c()
    for (i in 1:ListLen){SigPathways <- unique(c(SigPathways, CellChat_Obj_list[[i]]@netP$pathways))}
    
    cellchat <- mergeCellChat(CellChat_Obj_list, add.names = names(CellChat_Obj_list),cell.prefix = TRUE)
    # Compare the total number of interactions and interaction strength
    pdf(paste0(OutPrefix, "_compareInteractions.pdf"), width = 8, height = 8)
    gg1 <- compareInteractions(cellchat, show.legend = F)
    gg2 <- compareInteractions(cellchat, show.legend = F, measure = "weight")
    plot(gg1 + gg2)
    dev.off()
    
    # Compare the number of interactions and interaction strength among different cell populations
    pdf(paste0(OutPrefix, "_netVisual_diffInteraction.pdf"), width = 16, height = 10)
    for (i in 1:ncol(combinations)){
      SampleA <- combinations[1,i]; SampleB <- combinations[2,i]
      par(mfrow = c(1,3), xpd=TRUE)
      netVisual_diffInteraction(cellchat, weight.scale = T, color.use = ColorUse, comparison = c(SampleA, SampleB), label.edge = T,
                                title.name = paste0(names(CellChat_Obj_list[SampleA]),"_vs_",names(CellChat_Obj_list[SampleB])))
      netVisual_diffInteraction(cellchat, weight.scale = T, color.use = ColorUse, comparison = c(SampleA, SampleB), label.edge = F,
                                title.name = paste0(names(CellChat_Obj_list[SampleA]),"_vs_",names(CellChat_Obj_list[SampleB])))
      netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", color.use = ColorUse,comparison =  c(SampleA, SampleB))
    }
    dev.off()
    
    pdf(paste0(OutPrefix, "_netVisual_heatmap.pdf"), width = 12, height = 6)
    for (i in 1:ncol(combinations)){
      SampleA <- combinations[1,i]; SampleB <- combinations[2,i]
      gg1 <- netVisual_heatmap(cellchat, color.use = ColorUse, comparison = c(SampleA, SampleB),
                               title.name = paste0(names(CellChat_Obj_list[SampleA]),"_vs_",names(CellChat_Obj_list[SampleB])))
      gg2 <- netVisual_heatmap(cellchat, measure = "weight", color.use = ColorUse,comparison = c(SampleA, SampleB))
      plot(gg1 + gg2)
    }
    dev.off()
    
    pdf(paste0(OutPrefix, "_netVisual_circle.pdf"), width = 10, height = 10)
    weight.max <- getMaxWeight(CellChat_Obj_list, attribute = c("idents","count"))
    par(mfrow = c(3,2), xpd=TRUE)
    for (i in 1:length(CellChat_Obj_list)) {
      netVisual_circle(CellChat_Obj_list[[i]]@net$count, weight.scale = T, label.edge= T, edge.weight.max = weight.max[2],vertex.label.cex = 1, edge.label.cex = 1.2,
                       edge.width.max = 5, title.name = paste0("Number of interactions - ", names(CellChat_Obj_list)[i]), color.use = ColorUse)
    }
    dev.off()
    
    # Identify and visualize the conserved and context-specific signaling pathways, Compare the overall information flow of each signaling pathway
    pdf(paste0(OutPrefix, "_rankNet.pdf"), width = 12, height = 6)
    for (i in 1:ncol(combinations)){
      SampleA <- combinations[1,i];SampleB <- combinations[2,i]
      gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(SampleA, SampleB))
      gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,comparison = c(SampleA, SampleB))
      plot(gg1 + gg2)
    }
    dev.off()
    
    # Compare outgoing (or incoming) signaling associated with each cell population
    pdf(paste0(OutPrefix,"_netAnalysis_signalingRole_heatmap.pdf"), width = 14, height = 12)
    for (i in 1:ListLen){
      par(mfrow = c(2,2), xpd=TRUE)
      ht <- netAnalysis_signalingRole_heatmap(CellChat_Obj_list[[i]], pattern = "all", signaling = SigPathways, color.use = ColorUse,
                                              title = names(CellChat_Obj_list)[i], width = 10, height = 18, color.heatmap = "OrRd")
      draw(ht)
    }
    dev.off()
    
    # Identify the upgulated and down-regulated signaling ligand-receptor pairs
    pdf(paste0(OutPrefix, "_netVisual_bubble.pdf"), width = 14, height = 30)
    for (i in 1:ncol(combinations)){
      SampleA <- combinations[1,i]; SampleB <- combinations[2,i]
      SigPathways <- intersect(CellChat_Obj_list[[SampleA]]@netP$pathways, CellChat_Obj_list[[SampleB]]@netP$pathways)
      p <- netVisual_bubble(cellchat, comparison = c(SampleA, SampleB), signaling = SigPathways, angle.x = 90,
                            title.name = paste0(names(CellChat_Obj_list[SampleA]),"_vs_",names(CellChat_Obj_list[SampleB])));plot(p)
    }
    dev.off()
    
    pdf(paste0(OutPrefix, "_plotGeneExpression.pdf"), width = 12, height = 12)
    cellchat@meta$orig.ident = factor(cellchat@meta$orig.ident, levels = unique(cellchat@meta$orig.ident)) # set factor level
    for (pathway in SigPathways){p <- plotGeneExpression(cellchat, signaling = pathway, split.by = "orig.ident", colors.ggplot = T, color.use = ColorUse_SAM); plot(p)}
    dev.off()
  }
}

##' Plot multiple marker genes expression
##' @param Seurat_Obj A Seurat object
##' @param annotationPath Path store markers list parsed from the CellMarker 2.0 database (http://bio-bigdata.hrbmu.edu.cn/CellMarker/)
##' @param OutPrefix Prefix for output file
##' @examples
##'
##' Traversal_Cellmarkers(Seurat_Obj = iWAT_combined, annotationPath = "../Adipocytes_Markers/", OutPrefix = "_iWAT_combined_")

Traversal_Cellmarkers <- function(Seurat_Obj, annotationPath, OutPrefix){
  
  CellMarker_annotaions <- list.files(pattern = ".txt", path = annotationPath, full.names = T)
  
  #Vlnplot
  pdf(paste0(Sys.Date(), "_",OutPrefix,"Markers_Vlnplot.pdf"), width = 60, height = 15)
  for (file in CellMarker_annotaions){
    MarkersGenes <- read.table(file, header = F, sep = "\t")
    PVlnPlot <- VlnPlot(Seurat_Obj, features = unique(MarkersGenes$V9), stack=T) +
      theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1),
            axis.text.y = element_text(size=16)) +ggtitle(file) + NoLegend()
    plot(PVlnPlot)
  }
  dev.off()
  
  #Dotplot
  pdf(paste0(Sys.Date(),"_",OutPrefix,"Markers_Dotplot.pdf"), width = 60, height = 15)
  for (file in CellMarker_annotaions){
    MarkersGenes <- read.table(file, header = F, sep = "\t")
    PDotPlot <- DotPlot(Seurat_Obj, features = unique(MarkersGenes$V9)) +
      theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1),
            axis.text.y = element_text(size=16)) + ggtitle(file) + NoLegend()
    plot(PDotPlot)
  }
  dev.off()
}

##' Plot Cell Cycle percentages across cell types and samples
##' @param Seurat_Obj A Seurat object
##' @param ColorUse List of colors for plot
##' @param cellRanks Vector ranking cell types to plot 
##' @param OutPrefix Prefix for output file
##' @examples
##'
##' Plot_CellCycle(Seurat_Obj = iWAT_combined, OutPrefix = "_iWAT_combined_", ColorUse = iWAT_HFD_color_maps, cellRanks = NULL)

Plot_CellCycle <- function(Seurat_Obj, cellRanks = NULL, ColorUse, OutPrefix){
  
  CellCycle_freq <- as.data.frame(table(Seurat_Obj$Phase, Seurat_Obj$Myannotation,Seurat_Obj$orig.ident))
  if(!is.null(cellRanks)){CellCycle_freq$Var2 <- factor(CellCycle_freq$Var2, levels = cellRanks)}
  ggplot(data = CellCycle_freq, aes(x = Var2, y=Freq, fill = Var1)) +
    geom_bar(stat="identity", position = "fill", width = 0.8, color='black') +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() + labs(y = 'Fraction of Cell Cycle in each cell type',x="") +
    scale_fill_manual(values = ColorUse) + facet_grid(Var3 ~ .) +
    theme(
      axis.text.x = element_text(color="black", size=14, face="bold", angle = 60, hjust = 1),
      axis.title.y = element_text(color="black", size=16, face="bold"),
      axis.text.y = element_text(color="black", size=14, face="bold")
    )
  ggsave(paste0(OutPrefix, "Histogram_",Sys.Date(),".pdf"), width = 8, height = 8)
}

add_cell_annotation <- function(loom, cellAnnotation){
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation))){
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation)){add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])}
  invisible(loom)
}

##' Generate scVelo input files
Generate_scvelo_InputFiles <- function(Seurat_obj, OutPrefix){
  # save metadata table:
  Seurat_obj$barcode <- colnames(Seurat_obj)
  Seurat_obj$UMAP_1 <- Seurat_obj@reductions$umap@cell.embeddings[,1]
  Seurat_obj$UMAP_2 <- Seurat_obj@reductions$umap@cell.embeddings[,2]
  write.csv(Seurat_obj@meta.data, file=paste0(OutPrefix,'scvelo_metadata.csv'), quote=F, row.names=F)
  # write expression counts matrix
  counts_matrix <- GetAssayData(Seurat_obj, assay='RNA', slot='counts')
  writeMM(counts_matrix, file=paste0(OutPrefix, 'scvelo_counts.mtx'))
  # write dimesnionality reduction matrix, in this example case pca matrix
  write.csv(Seurat_obj@reductions$pca@cell.embeddings, file=paste0(OutPrefix,'scvelo_PCA.csv'), quote=F, row.names=F)
  # write gene names
  write.table(data.frame('gene'=rownames(counts_matrix)), file=paste0(OutPrefix,'scvelo_GeneNames.csv'), quote=F, row.names=F, col.names=F)
  return()
}

Generate_cellex_InputFiles <- function(Seurat_obj, OutPrefix){
  # write expression counts matrix
  counts_matrix <- GetAssayData(Seurat_obj, assay='RNA', slot='counts') %>% as.data.frame()
  # ref <- read.table("D:/OneDrive - pku.edu.cn/2_CarlosLab/4_Omics/mm10.allgenes_54838.bed")
  # name1 <- rownames(counts_matrix) %in% ref$V7
  # counts_matrix <- counts_matrix[name1,]
  # rownames(counts_matrix) <- ref[match(rownames(counts_matrix), ref$V7),]$V4
  write.csv(counts_matrix, file=paste0(OutPrefix, 'cellex_counts.csv'), quote = F)
  
  #Write cell annotation file
  Seurat_obj$cellexanno <- Idents(Seurat_obj)
  anno <- Seurat_obj@meta.data %>% rownames_to_column()
  anno <- anno[,c("rowname","cellexanno")]
  colnames(anno) <- c("","cell_type")
  write.csv(anno, file=paste0(OutPrefix,'cellex_Metadata.csv'), quote=F, row.names=F)
}

##' Create "spliced" and "unspliced" assays in existing Seurat object from loom files
##' @param Seurat_Obj A Seurat object
##' @param Loom_Obj Merged loom object
##' @param Convert_scVelo whether convert Seurat object to h5ad file for scVelo/CellRank for velocity analysis
##' @param OutPrefix Prefix for output file
##' @examples
##'
##' bm <- Velocyto_Into_Seurat(Seurat_Obj = BA_combined, Loom_Obj = ldat, Convert_scVelo = T, OutPrefix = BA_prefix)

Velocyto_Into_Seurat <- function(Seurat_Obj, Loom_Obj, Convert_scVelo=TRUE, OutPrefix) {
  dmatrix <- function(mat,...){
    cell=10000;submar=list(); end=0; j=0
    win <- ceiling(ncol(mat)/cell)
    for(i in 1:win){
      j = j + 1; start = end + 1; end = start + cell
      if(end > ncol(mat)){end = ncol(mat)}
      matTr <- as.matrix(mat[,start:end])
      submar[[j]] = matTr
    }
    return(submar)
  }
  
  MergeVelo <- function(Seurat_Obj, ldata){
    ldat_data <- dmatrix(ldata)
    ldat_data <- do.call(cbind, ldat_data)
    colnames(ldat_data) <- gsub("x$", "",colnames(ldat_data))
    sample <- as.data.frame(limma::strsplit2(split=":",colnames(ldat_data)))
    sample$V3 <-""
    
    CellSuffix <- table(sub("[A-Z]+","",rownames(Seurat_Obj@meta.data)),Seurat_Obj@meta.data$orig.ident) %>% 
      as.data.frame() %>% filter(Freq != 0)
    for (i in 1:nrow(CellSuffix)){
      sample$V3[which(sample$V1 == CellSuffix[i,]$Var2)] <- paste0(sample$V2[which(sample$V1 == CellSuffix[i,]$Var2)],CellSuffix[i,]$Var1)
    }
    colnames(ldat_data) <- sample$V3
    index <- duplicated(rownames(ldat_data))
    ldat_data <- do.call(cbind,dmatrix(ldat_data[!index,]))
    
    ## Create a matrix of the same size
    tmp_zero <- matrix(0, nrow=dim(as.matrix(Seurat_Obj@assays$RNA@counts))[1], ncol=dim(as.matrix(Seurat_Obj@assays$RNA@counts))[2])
    colnames(tmp_zero) <- colnames(as.matrix(Seurat_Obj@assays$RNA@counts))
    rownames(tmp_zero) <- rownames(as.matrix(Seurat_Obj@assays$RNA@counts))
    
    #Match velocyto counts to Seurat RNA counts
    ldat_data <- ldat_data[,colnames(ldat_data) %in% colnames(tmp_zero)]
    ldat_data <- ldat_data[rownames(ldat_data) %in% rownames(tmp_zero),]
    
    TURE_FAlES <- !colnames(tmp_zero) %in% colnames(ldat_data)
    tmp_col <- tmp_zero[,TURE_FAlES]
    colnames(tmp_col) <- colnames(tmp_zero)[TURE_FAlES]
    if (length(colnames(tmp_zero)[TURE_FAlES])>0) {
      for (i in length(colnames(tmp_col))) {
        add_col <- matrix(0, nrow=dim(ldat_data)[1], ncol=1)
        colnames(add_col) <- colnames(tmp_col)[i]
        rownames(add_col) <- rownames(ldat_data)
        ldat_data <- cbind(ldat_data, add_col)
      }
    }
    
    TURE_FAlES <- !rownames(tmp_zero) %in% rownames(ldat_data)
    tmp_row <- tmp_zero[TURE_FAlES,]
    rownames(tmp_row) <- rownames(tmp_zero)[TURE_FAlES]
    if (length(rownames(tmp_zero)[TURE_FAlES])>0) {
      for (j in 1:length(rownames(tmp_row))) {
        add_row <- matrix(0, nrow=1, ncol=dim(ldat_data)[2])
        colnames(add_row) <- colnames(ldat_data)
        rownames(add_row) <- rownames(tmp_row)[j]
        ldat_data <- rbind(ldat_data, add_row)
      }
    }
    return(ldat_data)
  }
  
  Seurat_Obj_new <- Seurat_Obj
  Seurat_Obj_new[["spliced"]] <- CreateAssayObject(counts = MergeVelo(Seurat_Obj, ldata = Loom_Obj$spliced))
  Seurat_Obj_new[["unspliced"]] <- CreateAssayObject(counts = MergeVelo(Seurat_Obj, ldata = Loom_Obj$unspliced))
  
  if(Convert_scVelo){
    Seurat_Obj_new[["RNA"]] <- Seurat_Obj_new[["spliced"]]
    DefaultAssay(Seurat_Obj_new) <- "RNA"
    h5SeuratFile <- paste0(OutPrefix,".h5Seurat")
    SaveH5Seurat(Seurat_Obj_new, filename = h5SeuratFile)
    Convert(h5SeuratFile, dest = "h5ad")
    file.remove(h5SeuratFile)
    if(F){#Pyhton code for scVelo on h5ad file
      # import scvelo as scv
      # adata = scv.read("OutPrefix.h5ad", cache=True)
      # scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
      # scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
      # scv.tl.velocity(adata)
      # scv.tl.velocity_graph(adata)
      # scv.pl.velocity_embedding_stream(adata, basis="umap", color=["clusters"])
    }
  }
  
  # saveRDS(Seurat_Obj_new, file = paste0(OutPrefix,"WithLoom.rds"))
  return(Seurat_Obj_new)
}

get_earliest_principal_node <- function(cds, root_type="mAd-1"){
  # cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
  # cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  cell_ids <- rownames(colData(cds))
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

#Scale monocle3 elements to 100
monocle3_scale_to_100 <- function(cells_subset){
  #browser()
  df <- data.frame((cells_subset@phenoData@data))
  df <- df[,c("Pseudotime", "My_annotation","orig.ident")]
  df <- df[order(df$Pseudotime, decreasing = F),]
  df$My_annotation <- as.factor(df$My_annotation)
  df$orig.ident <- as.factor(df$orig.ident)
  #df <- df[,!(colnames(df) %in% c("Pseudotime")), drop = F]
  len <- length(df$Pseudotime)
  bin<-round(len/100)
  Celltype <- c()
  Pseudotime <- c()
  orig.ident <- c()
  value1 <- c()
  value2 <- c()
  value3 <- c()
  for(i in 0:99){
    if(i < 99){
      start <- 1+(bin*i)
      stop <- bin+(bin*i)
      freq_table1 <- table(df$My_annotation[c(start:stop)])
      value1 <- names(freq_table1)[which.max(freq_table1)]
      value2 <- median(as.numeric(as.vector(df$Pseudotime[c(start:stop)])))
      freq_table3 <- table(df$orig.ident[c(start:stop)])
      value3 <- names(freq_table3)[which.max(freq_table3)]
      #browser()
      Celltype <- c(Celltype, value1)
      Pseudotime <- c(Pseudotime, value2)
      orig.ident <- c(orig.ident,value3)
    }
    else{
      Celltype <- c(Celltype, value1)
      Pseudotime <- c(Pseudotime, value2)
      orig.ident <- c(orig.ident,value3)
    }
  }
  return(data.frame(Celltype = Celltype,Pseudotime = Pseudotime,Timepoint = orig.ident))
}

#Scale monocle3 elements to 100
monocle3_scale_to_100_ASPC <- function(cells_subset){
  #browser()
  df <- data.frame((cells_subset@phenoData@data))
  df <- df[,c("Pseudotime", "Myannotation_ASPC","timepoint")]
  df <- df[order(df$Pseudotime, decreasing = F),]
  df$Myannotation_ASPC <- as.factor(df$Myannotation_ASPC)
  df$timepoint <- as.factor(df$timepoint)
  #df <- df[,!(colnames(df) %in% c("Pseudotime")), drop = F]
  len <- length(df$Pseudotime)
  bin<-round(len/100)
  Celltype <- c()
  Pseudotime <- c()
  timepoint <- c()
  value1 <- c()
  value2 <- c()
  value3 <- c()
  for(i in 0:99){
    if(i < 99){
      start <- 1+(bin*i)
      stop <- bin+(bin*i)
      freq_table1 <- table(df$Myannotation_ASPC[c(start:stop)])
      value1 <- names(freq_table1)[which.max(freq_table1)]
      value2 <- median(as.numeric(as.vector(df$Pseudotime[c(start:stop)])))
      freq_table3 <- table(df$timepoint[c(start:stop)])
      value3 <- names(freq_table3)[which.max(freq_table3)]
      #browser()
      Celltype <- c(Celltype, value1)
      Pseudotime <- c(Pseudotime, value2)
      timepoint <- c(timepoint,value3)
    }
    else{
      Celltype <- c(Celltype, value1)
      Pseudotime <- c(Pseudotime, value2)
      timepoint <- c(timepoint,value3)
    }
  }
  return(data.frame(Celltype = Celltype,Pseudotime = Pseudotime,Timepoint = timepoint))
}

##' Calulate Velocity and plot
##' @param emat spliced
##' @param nmat unspliced
##' @param cell.dist cell distance of umap
##' @param prefix Prefix for output file
##' @param emb Umap information
##' @param cell.colors cell.colors
##' @examples
##'
##' bm <- Velocyto_Into_Seurat(Seurat_Obj = BA_combined, Loom_Obj = ldat, Convert_scVelo = T, OutPrefix = BA_prefix)

Velocity_calu_plot <- function(emat,nmat,cell.dist,prefix,emb,cell.colors){
#计算RNA速率
rvel <- gene.relative.velocity.estimates(
  emat,nmat,
  deltaT = 1, kCells = 10,
  cell.dist = cell.dist,
  fit.quantile = 0.02, n.cores = 1)
save(rvel,file = paste0(prefix, "velocity_rvel.Rdata"))

#全局速率可视化
pdf(file = paste0(prefix, "velocity.pdf"), width = 10, height = 10)

show.velocity.on.embedding.cor(
  emb,rvel,n=50,
  scale='sqrt',
  cell.colors=ac(cell.colors,alpha=1),
  cex=1.8,
  arrow.scale=2.5,
  show.grid.flow=T,
  min.grid.cell.mass=3,
  grid.n=40,
  arrow.lwd=1.1
)
dev.off()

#pca降维并绘图
pdf(file = paste0(prefix, "velocity_PCA.pdf"), width = 10, height = 10)
pca.velocity.plot(rvel,nPcs=5,
                  plot.cols=2,
                  cell.colors=ac(cell.colors,alpha=0.7),
                  cex=1,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1),
                  arrow.scale = 0.1,
                  arrow.lwd = 0.1,
                  max.arrow.size= 0.1)
dev.off()
}

# reads layers/spliced/unspliced/ambiguous from a saptial loom file
##' Read in loom matrices of spliced/unpsliced reads as
##' prepared by velocyto.py CLI
##'
##' @param file loom file name
##' @param engine Use hdf5r or h5 to import loom file
##' @return a list containing spliced, unspliced, ambiguous and spanning matrices
##' @export
read.loom.matrices.spatial <- function(file, engine='hdf5r') {
  if (engine == 'h5'){
    cat('reading loom file via h5...\n')
    f <- h5::h5file(file,mode='r');
    cells <- f["col_attrs/CellID"][];
    genes <- f["row_attrs/Gene"][];
    dl <- c(spliced="/layers/spliced",unspliced="/layers/unspliced",ambiguous="/layers/ambiguous");
    if("/layers/spanning" %in% h5::list.datasets(f)) {
      dl <- c(dl,c(spanning="/layers/spanning"))
    }
    dlist <- lapply(dl,function(path) {
      m <- as(f[path][],'dgCMatrix'); rownames(m) <- genes; colnames(m) <- cells; return(m)
    })
    h5::h5close(f)
    return(dlist)
  } else if (engine == 'hdf5r') {
    cat('reading loom file via hdf5r...\n')
    #browser()
    f <- hdf5r::H5File$new(file, mode='r')
    cells <- f[["col_attrs/CellID"]][]
    genes <- f[["row_attrs/Gene"]][]
    dl <- c(spliced="layers/spliced",
            unspliced="layers/unspliced",
            ambiguous="layers/ambiguous")
    if("layers/spanning" %in% hdf5r::list.datasets(f)) {
      dl <- c(dl, c(spanning="layers/spanning"))
    }
    dlist <- lapply(dl, function(path) {
      m <- as(t(f[["matrix"]][,]),'dgCMatrix')
      rownames(m) <- genes; colnames(m) <- cells;
      return(m)
    })
    f$close_all()
    return(dlist)
  }
  else {
    warning('Unknown engine. Use hdf5r or h5 to import loom file.')
    return(list())
  }
}


