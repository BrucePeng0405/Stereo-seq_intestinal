#===============================================================================
.libPaths()
rm(list = ls())
print("[1]")
getwd()

library(Seurat)
library(velocyto.R)
library(ggplot2)
library(reticulate)
print("Yes")
raw_data_path <- 'raw/'
conda_list()
use_condaenv("r-velo", required = TRUE)
scv <- import("scvelo")

Timepoint <- "HFD4w8w/EP_six_diff"

#hughie_color <- c("#0072B5FF","#BC3C29FF", "#E18727FF", "#20854EFF", "#6F99ADFF", "#7876B1FF","#EE4C97FF", #nemj scheme
#"#E64B35FF", "#4DBBD5FF", "#3C5488FF", "#F39B7FFF")

#separate_timepoint <- function(Timepoint = Timepoint) {
################################CD group scvelo##################################
##Output path
CD_velocity_prefix <- paste0('OUTPUT/CD_group/',Timepoint,'/')

##read raw data
LI_merged_scvelo_CD <- readRDS(paste0(raw_data_path, 'LI_merged_CD_scvelo.rds'))
emat <- readRDS(paste0(raw_data_path, 'emat.CD.rds'))
nmat <- readRDS(paste0(raw_data_path, 'nmat.CD.rds'))


#??????emat,nmat??
LI_merged_scvelo_CD <- subset(LI_merged_scvelo_CD, My_annotation %in% c(                                                                                                    "Enterocyte Progenitor",
                                                                "Ribosome-Enriched Enterocyte",
                                                                "Absorptive Enterocyte",
                                                                "Mitochondrial Enterocyte",
                                                               "Lipid Metabolic Goblet Cell",
                                                               "Barrier Protective Goblet Cell",
                                                               "Plasma Goblet Cell"
                                                                ))
LI_merged_scvelo_CD <- subset(LI_merged_scvelo_CD, orig.ident %in% c("CD_LI_HFD4w","CD_LI_HFD8w"))
#c("CD_LI_HFD0w", "CD_LI_HFD4w",  "CD_LI_HFD8w"))

LI_merged_scvelo_CD <- subset(LI_merged_scvelo_CD, cells = colnames(emat))
LI_merged_scvelo_CD <- subset(LI_merged_scvelo_CD, features = rownames(emat))

emat <- emat[, colnames(emat) %in% row.names(as.data.frame(LI_merged_scvelo_CD$My_annotation))]
nmat <- nmat[, colnames(nmat) %in% row.names(as.data.frame(LI_merged_scvelo_CD$My_annotation))]

DimPlot(LI_merged_scvelo_CD, reduction = "umap", label = F, repel = T, label.box = T, order = T, cols=c(
    "#0072B5FF",
    "#BC3C29FF", 
    "#E18727FF", 
    "#20854EFF", 
   "#6F99ADFF",
   "#7876B1FF",
   "#EE4C97FF"
))
ggsave(paste0(CD_velocity_prefix,"UMAP_Cluster_",Sys.Date(),".pdf"), width = 8, height = 5)

#????????
cell.dist <- as.dist(1 - armaCor(t(LI_merged_scvelo_CD@reductions$umap@cell.embeddings)))

#??????????????
emb <- LI_merged_scvelo_CD@reductions$umap@cell.embeddings
pca <- LI_merged_scvelo_CD@reductions$pca@cell.embeddings

#Generate_scvelo_InputFiles(LI_merged_scvelo_CD,OutPrefix = paste0(CD_velocity_prefix,"Scvelo_files_"))

counts_matrix <- GetAssayData(LI_merged_scvelo_CD, assay='SCT', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(LI_merged_scvelo_CD@meta.data)
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

## run scvelo dynamic model
scv$pp$filter_genes(adata) ## filter
scv$pp$moments(adata) ## normalize and compute moments
scv$tl$recover_dynamics(adata, n_jobs=56) ## model

## plot (creates pop up window)
scv$tl$velocity(adata, mode='dynamical')
scv$tl$velocity_graph(adata)

##save
adata$write(paste0(CD_velocity_prefix,'Scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read('Scvelo_anndata.h5ad')

svg(paste0(CD_velocity_prefix,"Scvelo_HFD4w8w_","umap_",Sys.Date(),".svg"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='My_annotation',
                                 min_mass=2.5,
                                 arrow_size=1,
                                 density=1.3,
                                 linewidth=2,
                                 size=100,
                                 palette=c(
                                     "#0072B5FF",
                                     "#BC3C29FF",
                                     "#E18727FF",
                                     "#20854EFF", 
                                     "#6F99ADFF", 
                                     "#7876B1FF",
                                     "#EE4C97FF"
                                     ),
                                 alpha=0.5,
                                 legend_loc="none",
                                 title="",
                                 add_margin=0,
                                 dpi=300)
                               
dev.off()

svg(paste0(CD_velocity_prefix,"Scvelo_HFD4w8w_","pca_",Sys.Date(),".svg"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='My_annotation',
                                 min_mass=2.5,
                                 arrow_size=1,
                                 density=1.3,
                                 linewidth=1,
                                 size=100,
                                 palette=c(
                                     "#0072B5FF",
                                     "#BC3C29FF",
                                     "#E18727FF",
                                     "#20854EFF", 
                                        "#6F99ADFF",
                                        "#7876B1FF",
                                        "#EE4C97FF"
                                     ),
                                 alpha=0.5,
                                 legend_loc="none",
                                 title="",
                                 add_margin=0,
                                 dpi=300) ## other embedding
dev.off()
                              
print(paste0("CD_group ",Timepoint," Done!"))





################################HfiD group scvelo##################################
##Output path
HfiD_velocity_prefix <- paste0('OUTPUT/HfiD_group/',Timepoint,'/')

##read raw data
LI_merged_scvelo_HfiD <- readRDS(paste0(raw_data_path, 'LI_merged_HfiD_scvelo.rds'))
emat <- readRDS(paste0(raw_data_path, 'emat.HfiD.rds'))
nmat <- readRDS(paste0(raw_data_path, 'nmat.HfiD.rds'))


#??????emat,nmat??
LI_merged_scvelo_HfiD <- subset(LI_merged_scvelo_HfiD, My_annotation %in% c(
                                                                "Enterocyte Progenitor",
                                                                "Ribosome-Enriched Enterocyte",
                                                                "Absorptive Enterocyte",
                                                                "Mitochondrial Enterocyte",
                                                               "Lipid Metabolic Goblet Cell",
                                                               "Barrier Protective Goblet Cell",
                                                               "Plasma Goblet Cell"
                                                                ))
LI_merged_scvelo_HfiD <- subset(LI_merged_scvelo_HfiD, orig.ident %in% c("HfiD_LI_HFD4w",  "HfiD_LI_HFD8w"))
#c("HfiD_LI_HFD0w", "HfiD_LI_HFD4w",  "HfiD_LI_HFD8w"))

LI_merged_scvelo_HfiD <- subset(LI_merged_scvelo_HfiD, cells = colnames(emat))
LI_merged_scvelo_HfiD <- subset(LI_merged_scvelo_HfiD, features = rownames(emat))

emat <- emat[, colnames(emat) %in% row.names(as.data.frame(LI_merged_scvelo_HfiD$My_annotation))]
nmat <- nmat[, colnames(nmat) %in% row.names(as.data.frame(LI_merged_scvelo_HfiD$My_annotation))]

DimPlot(LI_merged_scvelo_HfiD, reduction = "umap", label = F, repel = T, label.box = T, order = T, cols=c(
    "#0072B5FF",
    "#BC3C29FF", 
    "#E18727FF", 
    "#20854EFF", 
   "#6F99ADFF",
   "#7876B1FF",
   "#EE4C97FF"
))
ggsave(paste0(HfiD_velocity_prefix,"UMAP_Cluster_",Sys.Date(),".pdf"), width = 8, height = 5)

#????????
cell.dist <- as.dist(1 - armaCor(t(LI_merged_scvelo_HfiD@reductions$umap@cell.embeddings)))

#??????????????
emb <- LI_merged_scvelo_HfiD@reductions$umap@cell.embeddings
pca <- LI_merged_scvelo_HfiD@reductions$pca@cell.embeddings

#Generate_scvelo_InputFiles(LI_merged_scvelo_HfiD,OutPrefix = paste0(HfiD_velocity_prefix,"Scvelo_files_"))

counts_matrix <- GetAssayData(LI_merged_scvelo_HfiD, assay='SCT', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(LI_merged_scvelo_HfiD@meta.data)
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

## run scvelo dynamic model
scv$pp$filter_genes(adata) ## filter
scv$pp$moments(adata) ## normalize and compute moments
scv$tl$recover_dynamics(adata, n_jobs=56) ## model

## plot (creates pop up window)
scv$tl$velocity(adata, mode='dynamical')
scv$tl$velocity_graph(adata)

##save
adata$write(paste0(HfiD_velocity_prefix,'Scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(HfiD_velocity_prefix,'Scvelo_anndata.h5ad'))

svg(paste0(HfiD_velocity_prefix,"Scvelo_HFD4w8w_","umap_",Sys.Date(),".svg"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='My_annotation',
                                 min_mass=2.5,
                                 arrow_size=1,
                                 density=1.3,
                                 linewidth=2,
                                 size=100,
                                 palette=c(
                                   "#0072B5FF",
                                   "#BC3C29FF",
                                   "#E18727FF",
                                   "#20854EFF", 
                                   "#6F99ADFF", 
                                   "#7876B1FF",
                                   "#EE4C97FF"
                                 ),
                                 alpha=0.5,
                                 legend_loc="none",
                                 title="",
                                 add_margin=0,
                                 dpi=300)
#                               
dev.off()

svg(paste0(HfiD_velocity_prefix,"Scvelo_HFD4w8w_","pca_",Sys.Date(),".svg"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='My_annotation',
                                 min_mass=2.5,
                                 arrow_size=1,
                                 density=1.3,
                                 linewidth=1,
                                 size=100,
                                 palette=c(
                                   "#0072B5FF",
                                   "#BC3C29FF",
                                   "#E18727FF",
                                   "#20854EFF", 
                                   "#6F99ADFF",
                                   "#7876B1FF",
                                   "#EE4C97FF"
                                 ),
                                 alpha=0.5,
                                 legend_loc="none",
                                 title="",
                                 add_margin=0,
                                 dpi=300) ## other embedding
dev.off()

print(paste0("HfiD_group ",Timepoint," Done!"))

#}


#separate_timepoint(Timepoint = "HFD0w")
#separate_timepoint(Timepoint = "HFD4w")
#separate_timepoint(Timepoint = "HFD8w")

print("GREAT JOB")