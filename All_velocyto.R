rm(list = ls())


source("PZ_ScFunctions.R")
library(reticulate)

conda_list()
use_condaenv("r-velo", required = TRUE)
scv <- import("scvelo")
#scv$logging$print_version()
########################Prepare emat and nmat############
ldat.CD.0w <- read.loom.matrices("raw/loom/mWAT_CD.loom")
emat.CD.0w <- ldat.CD.0w$spliced
colnames(emat.CD.0w) <- gsub("mWAT_CD2:", "CD_HFD0w_", colnames(emat.CD.0w))
colnames(emat.CD.0w) <- gsub("x", "-1", colnames(emat.CD.0w))

ldat.CD.4w <- read.loom.matrices("raw/loom/mWAT_CD-4HFD.loom")
emat.CD.4w <- ldat.CD.4w$spliced
colnames(emat.CD.4w) <- gsub("mWAT_CD-4HFD:", "CD_HFD4w_", colnames(emat.CD.4w))
colnames(emat.CD.4w) <- gsub("x", "-1", colnames(emat.CD.4w))

ldat.CD.8w <- read.loom.matrices("raw/loom/mWAT_CD-8HFD.loom")
emat.CD.8w <- ldat.CD.8w$spliced
colnames(emat.CD.8w) <- gsub("mWAT_CD-8HFD:", "CD_HFD8w_", colnames(emat.CD.8w))
colnames(emat.CD.8w) <- gsub("x", "-1", colnames(emat.CD.8w))

ldat.HfiD.0w <- read.loom.matrices("raw/loom/mWAT_HFiD.loom")
emat.HfiD.0w <- ldat.HfiD.0w$spliced
colnames(emat.HfiD.0w) <- gsub("mWAT_HFiD:", "HfiD_HFD0w_", colnames(emat.HfiD.0w))
colnames(emat.HfiD.0w) <- gsub("x", "-1", colnames(emat.HfiD.0w))

ldat.HfiD.4w <- read.loom.matrices("raw/loom/mWAT_HFiD-4HFD.loom")
emat.HfiD.4w <- ldat.HfiD.4w$spliced
colnames(emat.HfiD.4w) <- gsub("mWAT_HFiD-4HFD:", "HfiD_HFD4w_", colnames(emat.HfiD.4w))
colnames(emat.HfiD.4w) <- gsub("x", "-1", colnames(emat.HfiD.4w))

ldat.HfiD.8w <- read.loom.matrices("raw/loom/mWAT_HFiD-8HFD.loom")
emat.HfiD.8w <- ldat.HfiD.8w$spliced
colnames(emat.HfiD.8w) <- gsub("mWAT_HFiD-8HFD:", "HfiD_HFD8w_", colnames(emat.HfiD.8w))
colnames(emat.HfiD.8w) <- gsub("x", "-1", colnames(emat.HfiD.8w))

nmat.CD.0w <- ldat.CD.0w$unspliced
colnames(nmat.CD.0w) <- gsub("mWAT_CD2:", "CD_HFD0w_", colnames(nmat.CD.0w))
colnames(nmat.CD.0w) <- gsub("x", "-1", colnames(nmat.CD.0w))

nmat.CD.4w <- ldat.CD.4w$unspliced
colnames(nmat.CD.4w) <- gsub("mWAT_CD-4HFD:", "CD_HFD4w_", colnames(nmat.CD.4w))
colnames(nmat.CD.4w) <- gsub("x", "-1", colnames(nmat.CD.4w))

nmat.CD.8w <- ldat.CD.8w$unspliced
colnames(nmat.CD.8w) <- gsub("mWAT_CD-8HFD:", "CD_HFD8w_", colnames(nmat.CD.8w))
colnames(nmat.CD.8w) <- gsub("x", "-1", colnames(nmat.CD.8w))

nmat.HfiD.0w <- ldat.HfiD.0w$unspliced
colnames(nmat.HfiD.0w) <- gsub("mWAT_HFiD:", "HfiD_HFD0w_", colnames(nmat.HfiD.0w))
colnames(nmat.HfiD.0w) <- gsub("x", "-1", colnames(nmat.HfiD.0w))

nmat.HfiD.4w <- ldat.HfiD.4w$unspliced
colnames(nmat.HfiD.4w) <- gsub("mWAT_HFiD-4HFD:", "HfiD_HFD4w_", colnames(nmat.HfiD.4w))
colnames(nmat.HfiD.4w) <- gsub("x", "-1", colnames(nmat.HfiD.4w))

nmat.HfiD.8w <- ldat.HfiD.8w$unspliced
colnames(nmat.HfiD.8w) <- gsub("mWAT_HFiD-8HFD:", "HfiD_HFD8w_", colnames(nmat.HfiD.8w))
colnames(nmat.HfiD.8w) <- gsub("x", "-1", colnames(nmat.HfiD.8w))

mWAT_Ad_merged <- readRDS("Adipocyte/mWAT_Adipocytes_Seurat.rds")

mWAT_ASPC_merged <- readRDS("ASPC/mWAT_ASPC_Seurat.rds")

mWAT_Ad_ASPC_merged <- merge(mWAT_Ad_merged,y=mWAT_ASPC_merged,
                             add.cell.ids = NULL,
                             project = "CD_Ad_ASPC")

mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC <- mWAT_Ad_ASPC_merged@active.ident





########################CD HfiD & Ad ASPC###########################
mWAT_Ad_ASPC_velocyto_prefix <- "Adipocyte_velocyto/mWAT_Ad_ASPC_"

mWAT_Ad_ASPC_merged <- FindVariableFeatures(mWAT_Ad_ASPC_merged,selection.method = "vst", nfeatures = 2000)

mWAT_Ad_ASPC_merged <- ScaleData(mWAT_Ad_ASPC_merged) %>% RunPCA(npcs = 20) %>% RunUMAP(dims = 1:10) %>% 
  RunTSNE(dims = 1:10) %>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))

Idents(mWAT_Ad_ASPC_merged) <- mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC
DimPlot(mWAT_Ad_ASPC_merged, reduction = "umap", label = F, repel = T, label.box = T, order = T, cols = hughie_color)
ggsave(paste0(mWAT_Ad_ASPC_velocyto_prefix,"UMAP_Cluster_",Sys.Date(),".pdf"), width = 6, height = 5)

#scvelo
cell.chose <- as.data.frame(mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(mWAT_Ad_ASPC_merged@assays$RNA@meta.features)
feature.chose$feature <- row.names(feature.chose)

emat <- cbind(emat.CD.0w,emat.CD.4w,emat.CD.8w,
              emat.HfiD.0w,emat.HfiD.4w,emat.HfiD.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
emat <- emat[rownames(emat) %in% row.names(mWAT_Ad_ASPC_merged), ]

nmat <- cbind(nmat.CD.0w,nmat.CD.4w,nmat.CD.8w,
              nmat.HfiD.0w,nmat.HfiD.4w,nmat.HfiD.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% colnames(emat)]
nmat <- nmat[rownames(nmat) %in% row.names(mWAT_Ad_ASPC_merged), ]


cell.type <- cell.chose[colnames(emat),1]


#确保细胞数和emat，nmat一致
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, cells = colnames(emat))
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, features = rownames(emat))

#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings)))

#获取每个细胞在图上的坐标位置
emb <- mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings
pca <- mWAT_Ad_ASPC_merged@reductions$pca@cell.embeddings

Generate_scvelo_InputFiles(mWAT_Ad_ASPC_merged,OutPrefix = paste0(mWAT_Ad_ASPC_velocyto_prefix,"scvelo_files_"))

counts_matrix <- GetAssayData(mWAT_Ad_ASPC_merged, assay='RNA', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(mWAT_Ad_ASPC_merged@meta.data)
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
pdf(paste0(mWAT_Ad_ASPC_velocyto_prefix,"scvelo_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300)
dev.off()

pdf(paste0(mWAT_Ad_ASPC_velocyto_prefix,"scvelo_PCA_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300) ## other embedding
dev.off()

##save
adata$write(paste0(mWAT_Ad_ASPC_velocyto_prefix,'scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(mWAT_Ad_ASPC_velocyto_prefix,'scvelo_anndata.h5ad'))





########################CD & Ad ASPC############################
mWAT_Ad_ASPC_CD_velocyto_prefix <- "Adipocyte_velocyto/mWAT_Ad_ASPC_CD_"

mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, group %in% "CD")

mWAT_Ad_ASPC_merged <- FindVariableFeatures(mWAT_Ad_ASPC_merged,selection.method = "vst", nfeatures = 2000)

mWAT_Ad_ASPC_merged <- ScaleData(mWAT_Ad_ASPC_merged) %>% RunPCA(npcs = 20) %>% RunUMAP(dims = 1:10) %>% 
  RunTSNE(dims = 1:10) %>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))

Idents(mWAT_Ad_ASPC_merged) <- mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC

DimPlot(mWAT_Ad_ASPC_merged, reduction = "umap", label = F, repel = T, label.box = T, order = T, cols = hughie_color)
ggsave(paste0(mWAT_Ad_ASPC_CD_velocyto_prefix,"UMAP_Cluster_",Sys.Date(),".pdf"), width = 6, height = 5)

#scvelo
cell.chose <- as.data.frame(mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(mWAT_Ad_ASPC_merged@assays$RNA@meta.features)
feature.chose$feature <- row.names(feature.chose)

emat <- cbind(emat.CD.0w,emat.CD.4w,emat.CD.8w,
              emat.HfiD.0w,emat.HfiD.4w,emat.HfiD.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
emat <- emat[rownames(emat) %in% row.names(mWAT_Ad_ASPC_merged), ]

nmat <- cbind(nmat.CD.0w,nmat.CD.4w,nmat.CD.8w,
              nmat.HfiD.0w,nmat.HfiD.4w,nmat.HfiD.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% colnames(emat)]
nmat <- nmat[rownames(nmat) %in% row.names(mWAT_Ad_ASPC_merged), ]


cell.type <- cell.chose[colnames(emat),1]


#确保细胞数和emat，nmat一致
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, cells = colnames(emat))
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, features = rownames(emat))

#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings)))

#获取每个细胞在图上的坐标位置
emb <- mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings
pca <- mWAT_Ad_ASPC_merged@reductions$pca@cell.embeddings

Generate_scvelo_InputFiles(mWAT_Ad_ASPC_merged,OutPrefix = paste0(mWAT_Ad_ASPC_CD_velocyto_prefix,"scvelo_files_"))

counts_matrix <- GetAssayData(mWAT_Ad_ASPC_merged, assay='RNA', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(mWAT_Ad_ASPC_merged@meta.data)
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
pdf(paste0(mWAT_Ad_ASPC_CD_velocyto_prefix,"scvelo_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300)
dev.off()

pdf(paste0(mWAT_Ad_ASPC_CD_velocyto_prefix,"scvelo_PCA_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300) ## other embedding
dev.off()

##save
adata$write(paste0(mWAT_Ad_ASPC_CD_velocyto_prefix,'scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(mWAT_Ad_ASPC_CD_velocyto_prefix,'scvelo_anndata.h5ad'))





########################HfiD & Ad ASPC############################
mWAT_Ad_ASPC_HfiD_velocyto_prefix <- "Adipocyte_velocyto/mWAT_Ad_ASPC_HfiD_"

mWAT_Ad_ASPC_merged <- merge(mWAT_Ad_merged,y=mWAT_ASPC_merged,
                             add.cell.ids = NULL,
                             project = "CD_Ad_ASPC")

mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC <- mWAT_Ad_ASPC_merged@active.ident

mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, group %in% "HfiD")

mWAT_Ad_ASPC_merged <- FindVariableFeatures(mWAT_Ad_ASPC_merged,selection.method = "vst", nfeatures = 2000)

mWAT_Ad_ASPC_merged <- ScaleData(mWAT_Ad_ASPC_merged) %>% RunPCA(npcs = 20) %>% RunUMAP(dims = 1:10) %>% 
  RunTSNE(dims = 1:10) %>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))

Idents(mWAT_Ad_ASPC_merged) <- mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC

DimPlot(mWAT_Ad_ASPC_merged, reduction = "umap", label = F, repel = T, label.box = T, order = T, cols = hughie_color)
ggsave(paste0(mWAT_Ad_ASPC_HfiD_velocyto_prefix,"UMAP_Cluster_",Sys.Date(),".pdf"), width = 6, height = 5)

#scvelo
cell.chose <- as.data.frame(mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(mWAT_Ad_ASPC_merged@assays$RNA@meta.features)
feature.chose$feature <- row.names(feature.chose)

emat <- cbind(emat.CD.0w,emat.CD.4w,emat.CD.8w,
              emat.HfiD.0w,emat.HfiD.4w,emat.HfiD.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
emat <- emat[rownames(emat) %in% row.names(mWAT_Ad_ASPC_merged), ]

nmat <- cbind(nmat.CD.0w,nmat.CD.4w,nmat.CD.8w,
              nmat.HfiD.0w,nmat.HfiD.4w,nmat.HfiD.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% colnames(emat)]
nmat <- nmat[rownames(nmat) %in% row.names(mWAT_Ad_ASPC_merged), ]


cell.type <- cell.chose[colnames(emat),1]


#确保细胞数和emat，nmat一致
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, cells = colnames(emat))
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, features = rownames(emat))

#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings)))

#获取每个细胞在图上的坐标位置
emb <- mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings
pca <- mWAT_Ad_ASPC_merged@reductions$pca@cell.embeddings

Generate_scvelo_InputFiles(mWAT_Ad_ASPC_merged,OutPrefix = paste0(mWAT_Ad_ASPC_HfiD_velocyto_prefix,"scvelo_files_"))

counts_matrix <- GetAssayData(mWAT_Ad_ASPC_merged, assay='RNA', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(mWAT_Ad_ASPC_merged@meta.data)
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
pdf(paste0(mWAT_Ad_ASPC_HfiD_velocyto_prefix,"scvelo_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300)
dev.off()

pdf(paste0(mWAT_Ad_ASPC_HfiD_velocyto_prefix,"scvelo_PCA_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300) ## other embedding
dev.off()

##save
adata$write(paste0(mWAT_Ad_ASPC_HfiD_velocyto_prefix,'scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(mWAT_Ad_ASPC_CD_velocyto_prefix,'scvelo_anndata.h5ad'))





########################CD HfiD & Ad ASPC-1###########################
mWAT_Ad_ASPC1_velocyto_prefix <- "Adipocyte_velocyto/mWAT_Ad_ASPC1_"

mWAT_Ad_ASPC_merged <- merge(mWAT_Ad_merged,y=mWAT_ASPC_merged,
                             add.cell.ids = NULL,
                             project = "CD_Ad_ASPC")
mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC <- mWAT_Ad_ASPC_merged@active.ident
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, Myannotation_Ad_ASPC %in% 
                                c("mWAT-ad1","mWAT-ad2","mWAT-ad3",
                                  "mWAT-ASPC1"))

mWAT_Ad_ASPC_merged <- FindVariableFeatures(mWAT_Ad_ASPC_merged,selection.method = "vst", nfeatures = 2000)

mWAT_Ad_ASPC_merged <- ScaleData(mWAT_Ad_ASPC_merged) %>% RunPCA(npcs = 20) %>% RunUMAP(dims = 1:10) %>% 
  RunTSNE(dims = 1:10) %>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))

Idents(mWAT_Ad_ASPC_merged) <- mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC
DimPlot(mWAT_Ad_ASPC_merged, reduction = "umap", label = F, repel = T, label.box = T, order = T, cols = hughie_color)
ggsave(paste0(mWAT_Ad_ASPC1_velocyto_prefix,"UMAP_Cluster_",Sys.Date(),".pdf"), width = 6, height = 5)

#scvelo
cell.chose <- as.data.frame(mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(mWAT_Ad_ASPC_merged@assays$RNA@meta.features)
feature.chose$feature <- row.names(feature.chose)

emat <- cbind(emat.CD.0w,emat.CD.4w,emat.CD.8w,
              emat.HfiD.0w,emat.HfiD.4w,emat.HfiD.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
emat <- emat[rownames(emat) %in% row.names(mWAT_Ad_ASPC_merged), ]

nmat <- cbind(nmat.CD.0w,nmat.CD.4w,nmat.CD.8w,
              nmat.HfiD.0w,nmat.HfiD.4w,nmat.HfiD.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% colnames(emat)]
nmat <- nmat[rownames(nmat) %in% row.names(mWAT_Ad_ASPC_merged), ]


cell.type <- cell.chose[colnames(emat),1]


#确保细胞数和emat，nmat一致
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, cells = colnames(emat))
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, features = rownames(emat))

#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings)))

#获取每个细胞在图上的坐标位置
emb <- mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings
pca <- mWAT_Ad_ASPC_merged@reductions$pca@cell.embeddings

Generate_scvelo_InputFiles(mWAT_Ad_ASPC_merged,OutPrefix = paste0(mWAT_Ad_ASPC1_velocyto_prefix,"scvelo_files_"))

counts_matrix <- GetAssayData(mWAT_Ad_ASPC_merged, assay='RNA', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(mWAT_Ad_ASPC_merged@meta.data)
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
pdf(paste0(mWAT_Ad_ASPC1_velocyto_prefix,"scvelo_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=2.5,
                                 arrow_size=1.6,
                                 density=1.8,
                                 linewidth=3, 
                                 size=300,
                                 palette=c("#BFA6C9","#F5E0BA","#AED0DF","#9b5b33"),
                                 alpha=0.5,
                                 legend_loc="none",
                                 title="",
                                 add_margin=0,
                                 dpi=300)
dev.off()

pdf(paste0(mWAT_Ad_ASPC1_velocyto_prefix,"scvelo_PCA_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=2.5,
                                 arrow_size=1.6,
                                 density=1.8,
                                 linewidth=3, 
                                 size=500,
                                 palette=c("#BFA6C9","#F5E0BA","#AED0DF","#9b5b33"),
                                 alpha=0.5,
                                 legend_loc="none",
                                 title="",
                                 add_margin=0,
                                 dpi=300) ## other embedding
dev.off()

##save
adata$write(paste0(mWAT_Ad_ASPC1_velocyto_prefix,'scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(mWAT_Ad_ASPC1_velocyto_prefix,'scvelo_anndata.h5ad'))






########################CD HfiD & Ad ASPC-2###########################
mWAT_Ad_ASPC2_velocyto_prefix <- "Adipocyte_velocyto/mWAT_Ad_ASPC2_"

mWAT_Ad_ASPC_merged <- merge(mWAT_Ad_merged,y=mWAT_ASPC_merged,
                             add.cell.ids = NULL,
                             project = "CD_Ad_ASPC")
mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC <- mWAT_Ad_ASPC_merged@active.ident
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, Myannotation_Ad_ASPC %in% 
                                c("mWAT-ad1","mWAT-ad2","mWAT-ad3",
                                  "mWAT-ASPC2"))

mWAT_Ad_ASPC_merged <- FindVariableFeatures(mWAT_Ad_ASPC_merged,selection.method = "vst", nfeatures = 2000)

mWAT_Ad_ASPC_merged <- ScaleData(mWAT_Ad_ASPC_merged) %>% RunPCA(npcs = 20) %>% RunUMAP(dims = 1:10) %>% 
  RunTSNE(dims = 1:10) %>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))

Idents(mWAT_Ad_ASPC_merged) <- mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC
DimPlot(mWAT_Ad_ASPC_merged, reduction = "umap", label = F, repel = T, label.box = T, order = T, cols = hughie_color)
ggsave(paste0(mWAT_Ad_ASPC2_velocyto_prefix,"UMAP_Cluster_",Sys.Date(),".pdf"), width = 6, height = 5)

#scvelo
cell.chose <- as.data.frame(mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(mWAT_Ad_ASPC_merged@assays$RNA@meta.features)
feature.chose$feature <- row.names(feature.chose)

emat <- cbind(emat.CD.0w,emat.CD.4w,emat.CD.8w,
              emat.HfiD.0w,emat.HfiD.4w,emat.HfiD.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
emat <- emat[rownames(emat) %in% row.names(mWAT_Ad_ASPC_merged), ]

nmat <- cbind(nmat.CD.0w,nmat.CD.4w,nmat.CD.8w,
              nmat.HfiD.0w,nmat.HfiD.4w,nmat.HfiD.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% colnames(emat)]
nmat <- nmat[rownames(nmat) %in% row.names(mWAT_Ad_ASPC_merged), ]


cell.type <- cell.chose[colnames(emat),1]


#确保细胞数和emat，nmat一致
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, cells = colnames(emat))
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, features = rownames(emat))

#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings)))

#获取每个细胞在图上的坐标位置
emb <- mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings
pca <- mWAT_Ad_ASPC_merged@reductions$pca@cell.embeddings

Generate_scvelo_InputFiles(mWAT_Ad_ASPC_merged,OutPrefix = paste0(mWAT_Ad_ASPC2_velocyto_prefix,"scvelo_files_"))

counts_matrix <- GetAssayData(mWAT_Ad_ASPC_merged, assay='RNA', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(mWAT_Ad_ASPC_merged@meta.data)
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
pdf(paste0(mWAT_Ad_ASPC2_velocyto_prefix,"scvelo_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, 
                                 size=500,
                                 palette=c("#BFA6C9","#F5E0BA","#AED0DF","#f5cfe4"),
                                 alpha=0.5,
                                 legend_loc="none",
                                 title="",
                                 add_margin=0,
                                 dpi=300
                                 )
dev.off()

pdf(paste0(mWAT_Ad_ASPC2_velocyto_prefix,"scvelo_PCA_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=2.5,
                                 arrow_size=1.6,
                                 density=1.8,
                                 linewidth=3, 
                                 size=500,
                                 palette=c("#BFA6C9","#F5E0BA","#AED0DF","#f5cfe4"),
                                 alpha=0.5,
                                 legend_loc="none",
                                 title="",
                                 add_margin=0,
                                 dpi=300) ## other embedding
dev.off()

##save
adata$write(paste0(mWAT_Ad_ASPC2_velocyto_prefix,'scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(mWAT_Ad_ASPC2_velocyto_prefix,'scvelo_anndata.h5ad'))






########################CD HfiD & Ad ASPC-3###########################
mWAT_Ad_ASPC3_velocyto_prefix <- "Adipocyte_velocyto/mWAT_Ad_ASPC3_"

mWAT_Ad_ASPC_merged <- merge(mWAT_Ad_merged,y=mWAT_ASPC_merged,
                             add.cell.ids = NULL,
                             project = "CD_Ad_ASPC")
mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC <- mWAT_Ad_ASPC_merged@active.ident
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, Myannotation_Ad_ASPC %in% 
                                c("mWAT-ad1","mWAT-ad2","mWAT-ad3",
                                  "mWAT-ASPC3"))

mWAT_Ad_ASPC_merged <- FindVariableFeatures(mWAT_Ad_ASPC_merged,selection.method = "vst", nfeatures = 2000)

mWAT_Ad_ASPC_merged <- ScaleData(mWAT_Ad_ASPC_merged) %>% RunPCA(npcs = 20) %>% RunUMAP(dims = 1:10) %>% 
  RunTSNE(dims = 1:10) %>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))

Idents(mWAT_Ad_ASPC_merged) <- mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC
DimPlot(mWAT_Ad_ASPC_merged, reduction = "umap", label = F, repel = T, label.box = T, order = T, cols = hughie_color)
ggsave(paste0(mWAT_Ad_ASPC3_velocyto_prefix,"UMAP_Cluster_",Sys.Date(),".pdf"), width = 6, height = 5)

#scvelo
cell.chose <- as.data.frame(mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(mWAT_Ad_ASPC_merged@assays$RNA@meta.features)
feature.chose$feature <- row.names(feature.chose)

emat <- cbind(emat.CD.0w,emat.CD.4w,emat.CD.8w,
              emat.HfiD.0w,emat.HfiD.4w,emat.HfiD.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
emat <- emat[rownames(emat) %in% row.names(mWAT_Ad_ASPC_merged), ]

nmat <- cbind(nmat.CD.0w,nmat.CD.4w,nmat.CD.8w,
              nmat.HfiD.0w,nmat.HfiD.4w,nmat.HfiD.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% colnames(emat)]
nmat <- nmat[rownames(nmat) %in% row.names(mWAT_Ad_ASPC_merged), ]


cell.type <- cell.chose[colnames(emat),1]


#确保细胞数和emat，nmat一致
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, cells = colnames(emat))
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, features = rownames(emat))

#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings)))

#获取每个细胞在图上的坐标位置
emb <- mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings
pca <- mWAT_Ad_ASPC_merged@reductions$pca@cell.embeddings

Generate_scvelo_InputFiles(mWAT_Ad_ASPC_merged,OutPrefix = paste0(mWAT_Ad_ASPC3_velocyto_prefix,"scvelo_files_"))

counts_matrix <- GetAssayData(mWAT_Ad_ASPC_merged, assay='RNA', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(mWAT_Ad_ASPC_merged@meta.data)
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
pdf(paste0(mWAT_Ad_ASPC3_velocyto_prefix,"scvelo_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, 
                                 size=500,
                                 palette=c("#BFA6C9","#F5E0BA","#AED0DF","#8fa4ae"),
                                 alpha=0.5,
                                 legend_loc="none",
                                 title="",
                                 add_margin=0,
                                 dpi=300)
dev.off()

pdf(paste0(mWAT_Ad_ASPC3_velocyto_prefix,"scvelo_PCA_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=2.5,
                                 arrow_size=1.6,
                                 density=1.8,
                                 linewidth=3, 
                                 size=500,
                                 palette=c("#BFA6C9","#F5E0BA","#AED0DF","#8fa4ae"),
                                 alpha=0.5,
                                 legend_loc="none",
                                 title="",
                                 add_margin=0,
                                 dpi=300) ## other embedding
dev.off()

##save
adata$write(paste0(mWAT_Ad_ASPC3_velocyto_prefix,'scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(mWAT_Ad_ASPC3_velocyto_prefix,'scvelo_anndata.h5ad'))






########################CD & Ad ASPC############################
mWAT_Ad_ASPC2_CD_velocyto_prefix <- "Adipocyte_velocyto/mWAT_Ad_ASPC2_CD_"

mWAT_Ad_ASPC_merged <- merge(mWAT_Ad_merged,y=mWAT_ASPC_merged,
                             add.cell.ids = NULL,
                             project = "CD_Ad_ASPC")
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, group %in% "CD")
mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC <- mWAT_Ad_ASPC_merged@active.ident
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, Myannotation_Ad_ASPC %in% 
                                c("mWAT-ad1","mWAT-ad2","mWAT-ad3",
                                  "mWAT-ASPC2"))

mWAT_Ad_ASPC_merged <- FindVariableFeatures(mWAT_Ad_ASPC_merged,selection.method = "vst", nfeatures = 2000)

mWAT_Ad_ASPC_merged <- ScaleData(mWAT_Ad_ASPC_merged) %>% RunPCA(npcs = 20) %>% RunUMAP(dims = 1:10) %>% 
  RunTSNE(dims = 1:10) %>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))

Idents(mWAT_Ad_ASPC_merged) <- mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC

DimPlot(mWAT_Ad_ASPC_merged, reduction = "umap", label = F, repel = T, label.box = T, order = T, cols = hughie_color)
ggsave(paste0(mWAT_Ad_ASPC2_CD_velocyto_prefix,"UMAP_Cluster_",Sys.Date(),".pdf"), width = 6, height = 5)

#scvelo
cell.chose <- as.data.frame(mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(mWAT_Ad_ASPC_merged@assays$RNA@meta.features)
feature.chose$feature <- row.names(feature.chose)

emat <- cbind(emat.CD.0w,emat.CD.4w,emat.CD.8w,
              emat.HfiD.0w,emat.HfiD.4w,emat.HfiD.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
emat <- emat[rownames(emat) %in% row.names(mWAT_Ad_ASPC_merged), ]

nmat <- cbind(nmat.CD.0w,nmat.CD.4w,nmat.CD.8w,
              nmat.HfiD.0w,nmat.HfiD.4w,nmat.HfiD.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% colnames(emat)]
nmat <- nmat[rownames(nmat) %in% row.names(mWAT_Ad_ASPC_merged), ]


cell.type <- cell.chose[colnames(emat),1]


#确保细胞数和emat，nmat一致
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, cells = colnames(emat))
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, features = rownames(emat))

#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings)))

#获取每个细胞在图上的坐标位置
emb <- mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings
pca <- mWAT_Ad_ASPC_merged@reductions$pca@cell.embeddings

Generate_scvelo_InputFiles(mWAT_Ad_ASPC_merged,OutPrefix = paste0(mWAT_Ad_ASPC2_CD_velocyto_prefix,"scvelo_files_"))

counts_matrix <- GetAssayData(mWAT_Ad_ASPC_merged, assay='RNA', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(mWAT_Ad_ASPC_merged@meta.data)
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
pdf(paste0(mWAT_Ad_ASPC2_CD_velocyto_prefix,"scvelo_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300)
dev.off()

pdf(paste0(mWAT_Ad_ASPC2_CD_velocyto_prefix,"scvelo_PCA_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300) ## other embedding
dev.off()

##save
adata$write(paste0(mWAT_Ad_ASPC2_CD_velocyto_prefix,'scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(mWAT_Ad_ASPC_CD_velocyto_prefix,'scvelo_anndata.h5ad'))





########################HfiD & Ad ASPC############################
mWAT_Ad_ASPC2_HfiD_velocyto_prefix <- "Adipocyte_velocyto/mWAT_Ad_ASPC2_HfiD_"

mWAT_Ad_ASPC_merged <- merge(mWAT_Ad_merged,y=mWAT_ASPC_merged,
                             add.cell.ids = NULL,
                             project = "CD_Ad_ASPC")
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, group %in% "HfiD")
mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC <- mWAT_Ad_ASPC_merged@active.ident
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, Myannotation_Ad_ASPC %in% 
                                c("mWAT-ad1","mWAT-ad2","mWAT-ad3",
                                  "mWAT-ASPC2"))

mWAT_Ad_ASPC_merged <- FindVariableFeatures(mWAT_Ad_ASPC_merged,selection.method = "vst", nfeatures = 2000)

mWAT_Ad_ASPC_merged <- ScaleData(mWAT_Ad_ASPC_merged) %>% RunPCA(npcs = 20) %>% RunUMAP(dims = 1:10) %>% 
  RunTSNE(dims = 1:10) %>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))

Idents(mWAT_Ad_ASPC_merged) <- mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC

DimPlot(mWAT_Ad_ASPC_merged, reduction = "umap", label = F, repel = T, label.box = T, order = T, cols = hughie_color)
ggsave(paste0(mWAT_Ad_ASPC2_HfiD_velocyto_prefix,"UMAP_Cluster_",Sys.Date(),".pdf"), width = 6, height = 5)

#scvelo
cell.chose <- as.data.frame(mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(mWAT_Ad_ASPC_merged@assays$RNA@meta.features)
feature.chose$feature <- row.names(feature.chose)

emat <- cbind(emat.CD.0w,emat.CD.4w,emat.CD.8w,
              emat.HfiD.0w,emat.HfiD.4w,emat.HfiD.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
emat <- emat[rownames(emat) %in% row.names(mWAT_Ad_ASPC_merged), ]

nmat <- cbind(nmat.CD.0w,nmat.CD.4w,nmat.CD.8w,
              nmat.HfiD.0w,nmat.HfiD.4w,nmat.HfiD.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% colnames(emat)]
nmat <- nmat[rownames(nmat) %in% row.names(mWAT_Ad_ASPC_merged), ]


cell.type <- cell.chose[colnames(emat),1]


#确保细胞数和emat，nmat一致
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, cells = colnames(emat))
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, features = rownames(emat))

#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings)))

#获取每个细胞在图上的坐标位置
emb <- mWAT_Ad_ASPC_merged@reductions$umap@cell.embeddings
pca <- mWAT_Ad_ASPC_merged@reductions$pca@cell.embeddings

Generate_scvelo_InputFiles(mWAT_Ad_ASPC_merged,OutPrefix = paste0(mWAT_Ad_ASPC2_HfiD_velocyto_prefix,"scvelo_files_"))

counts_matrix <- GetAssayData(mWAT_Ad_ASPC_merged, assay='RNA', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(mWAT_Ad_ASPC_merged@meta.data)
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
pdf(paste0(mWAT_Ad_ASPC2_HfiD_velocyto_prefix,"scvelo_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300)
dev.off()

pdf(paste0(mWAT_Ad_ASPC2_HfiD_velocyto_prefix,"scvelo_PCA_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='Myannotation_Ad_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300) ## other embedding
dev.off()

##save
adata$write(paste0(mWAT_Ad_ASPC2_HfiD_velocyto_prefix,'scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(mWAT_Ad_ASPC_CD_velocyto_prefix,'scvelo_anndata.h5ad'))






########################CD & Ad############################
mWAT_Ad_CD_velocyto_prefix <- "Adipocyte_velocyto/mWAT_Ad_CD_"

mWAT_Ad_CD_merged <- subset(mWAT_Ad_merged, group %in% "CD")

#scvelo
cell.chose <- as.data.frame(mWAT_Ad_CD_merged$Myannotation_Ad)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(mWAT_Ad_CD_merged@assays$RNA@meta.features)
feature.chose$feature <- row.names(feature.chose)

emat <- cbind(emat.CD.0w,emat.CD.4w,emat.CD.8w,
              emat.HfiD.0w,emat.HfiD.4w,emat.HfiD.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
emat <- emat[rownames(emat) %in% row.names(mWAT_Ad_CD_merged), ]

nmat <- cbind(nmat.CD.0w,nmat.CD.4w,nmat.CD.8w,
              nmat.HfiD.0w,nmat.HfiD.4w,nmat.HfiD.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% colnames(emat)]
nmat <- nmat[rownames(nmat) %in% row.names(mWAT_Ad_CD_merged), ]


cell.type <- cell.chose[colnames(emat),1]


#确保细胞数和emat，nmat一致
mWAT_Ad_CD_merged <- subset(mWAT_Ad_CD_merged, cells = colnames(emat))
mWAT_Ad_CD_merged <- subset(mWAT_Ad_CD_merged, features = rownames(emat))

#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(mWAT_Ad_CD_merged@reductions$umap@cell.embeddings)))

#获取每个细胞在图上的坐标位置
emb <- mWAT_Ad_CD_merged@reductions$umap@cell.embeddings
pca <- mWAT_Ad_CD_merged@reductions$pca@cell.embeddings

Generate_scvelo_InputFiles(mWAT_Ad_CD_merged,OutPrefix = paste0(mWAT_Ad_CD_velocyto_prefix,"scvelo_files_"))

counts_matrix <- GetAssayData(mWAT_Ad_CD_merged, assay='RNA', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(mWAT_Ad_CD_merged@meta.data)
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
pdf(paste0(mWAT_Ad_CD_velocyto_prefix,"scvelo_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='Myannotation_Ad',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300)
dev.off()

pdf(paste0(mWAT_Ad_CD_velocyto_prefix,"scvelo_PCA_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='Myannotation_Ad',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300) ## other embedding
dev.off()

##save
adata$write(paste0(mWAT_Ad_CD_velocyto_prefix,'scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(mWAT_Ad_ASPC_CD_velocyto_prefix,'scvelo_anndata.h5ad'))






########################HfiD & Ad############################
mWAT_Ad_HfiD_velocyto_prefix <- "Adipocyte_velocyto/mWAT_Ad_HfiD_"

mWAT_Ad_HfiD_merged <- subset(mWAT_Ad_merged, group %in% "HfiD")

#scvelo
cell.chose <- as.data.frame(mWAT_Ad_HfiD_merged$Myannotation_Ad)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(mWAT_Ad_HfiD_merged@assays$RNA@meta.features)
feature.chose$feature <- row.names(feature.chose)

emat <- cbind(emat.CD.0w,emat.CD.4w,emat.CD.8w,
              emat.HfiD.0w,emat.HfiD.4w,emat.HfiD.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
emat <- emat[rownames(emat) %in% row.names(mWAT_Ad_HfiD_merged), ]

nmat <- cbind(nmat.CD.0w,nmat.CD.4w,nmat.CD.8w,
              nmat.HfiD.0w,nmat.HfiD.4w,nmat.HfiD.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% colnames(emat)]
nmat <- nmat[rownames(nmat) %in% row.names(mWAT_Ad_HfiD_merged), ]


cell.type <- cell.chose[colnames(emat),1]


#确保细胞数和emat，nmat一致
mWAT_Ad_HfiD_merged <- subset(mWAT_Ad_HfiD_merged, cells = colnames(emat))
mWAT_Ad_HfiD_merged <- subset(mWAT_Ad_HfiD_merged, features = rownames(emat))

#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(mWAT_Ad_HfiD_merged@reductions$umap@cell.embeddings)))

#获取每个细胞在图上的坐标位置
emb <- mWAT_Ad_HfiD_merged@reductions$umap@cell.embeddings
pca <- mWAT_Ad_HfiD_merged@reductions$pca@cell.embeddings

Generate_scvelo_InputFiles(mWAT_Ad_HfiD_merged,OutPrefix = paste0(mWAT_Ad_HfiD_velocyto_prefix,"scvelo_files_"))

counts_matrix <- GetAssayData(mWAT_Ad_HfiD_merged, assay='RNA', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(mWAT_Ad_HfiD_merged@meta.data)
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
pdf(paste0(mWAT_Ad_HfiD_velocyto_prefix,"scvelo_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='Myannotation_Ad',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300)
dev.off()

pdf(paste0(mWAT_Ad_HfiD_velocyto_prefix,"scvelo_PCA_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='Myannotation_Ad',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300) ## other embedding
dev.off()

##save
adata$write(paste0(mWAT_Ad_HfiD_velocyto_prefix,'scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(mWAT_Ad_HfiD_velocyto_prefix,'scvelo_anndata.h5ad'))






########################CD & ASPC############################
mWAT_ASPC_CD_velocyto_prefix <- "Adipocyte_velocyto/mWAT_ASPC_CD_"

mWAT_ASPC_CD_merged <- subset(mWAT_ASPC_merged, group %in% "CD")

#scvelo
cell.chose <- as.data.frame(mWAT_ASPC_CD_merged$Myannotation_ASPC)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(mWAT_ASPC_CD_merged@assays$RNA@meta.features)
feature.chose$feature <- row.names(feature.chose)

emat <- cbind(emat.CD.0w,emat.CD.4w,emat.CD.8w,
              emat.HfiD.0w,emat.HfiD.4w,emat.HfiD.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
emat <- emat[rownames(emat) %in% row.names(mWAT_ASPC_CD_merged), ]

nmat <- cbind(nmat.CD.0w,nmat.CD.4w,nmat.CD.8w,
              nmat.HfiD.0w,nmat.HfiD.4w,nmat.HfiD.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% colnames(emat)]
nmat <- nmat[rownames(nmat) %in% row.names(mWAT_ASPC_CD_merged), ]


cell.type <- cell.chose[colnames(emat),1]


#确保细胞数和emat，nmat一致
mWAT_ASPC_CD_merged <- subset(mWAT_ASPC_CD_merged, cells = colnames(emat))
mWAT_ASPC_CD_merged <- subset(mWAT_ASPC_CD_merged, features = rownames(emat))

#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(mWAT_ASPC_CD_merged@reductions$umap@cell.embeddings)))

#获取每个细胞在图上的坐标位置
emb <- mWAT_ASPC_CD_merged@reductions$umap@cell.embeddings
pca <- mWAT_ASPC_CD_merged@reductions$pca@cell.embeddings

Generate_scvelo_InputFiles(mWAT_ASPC_CD_merged,OutPrefix = paste0(mWAT_ASPC_CD_velocyto_prefix,"scvelo_files_"))

counts_matrix <- GetAssayData(mWAT_ASPC_CD_merged, assay='RNA', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(mWAT_ASPC_CD_merged@meta.data)
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
#clusters <- as.factor(adata$obs$Myannotation_ASPC_ASPC)
#rownames(emb) <- names(clusters) <- as.matrix(adata$obs_names$values)


## run scvelo dynamic model
scv$pp$filter_genes(adata) ## filter
scv$pp$moments(adata) ## normalize and compute moments
scv$tl$recover_dynamics(adata) ## model

## plot (creates pop up window)
scv$tl$velocity(adata, mode='dynamical')
scv$tl$velocity_graph(adata)
pdf(paste0(mWAT_ASPC_CD_velocyto_prefix,"scvelo_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='Myannotation_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300)
dev.off()

pdf(paste0(mWAT_ASPC_CD_velocyto_prefix,"scvelo_PCA_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='Myannotation_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300) ## other embedding
dev.off()

##save
adata$write(paste0(mWAT_ASPC_CD_velocyto_prefix,'scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(mWAT_ASPC_CD_velocyto_prefix,'scvelo_anndata.h5ad'))






########################HfiD & ASPC############################
mWAT_ASPC_HfiD_velocyto_prefix <- "Adipocyte_velocyto/mWAT_ASPC_HfiD_"

mWAT_ASPC_HfiD_merged <- subset(mWAT_ASPC_merged, group %in% "HfiD")

#scvelo
cell.chose <- as.data.frame(mWAT_ASPC_HfiD_merged$Myannotation_ASPC)
cell.chose$barcodes <- row.names(cell.chose)
feature.chose <- as.data.frame(mWAT_ASPC_HfiD_merged@assays$RNA@meta.features)
feature.chose$feature <- row.names(feature.chose)

emat <- cbind(emat.CD.0w,emat.CD.4w,emat.CD.8w,
              emat.HfiD.0w,emat.HfiD.4w,emat.HfiD.8w)
emat <- emat[unique(row.names(emat)),]

#emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
emat <- emat[rownames(emat) %in% row.names(mWAT_ASPC_HfiD_merged), ]

nmat <- cbind(nmat.CD.0w,nmat.CD.4w,nmat.CD.8w,
              nmat.HfiD.0w,nmat.HfiD.4w,nmat.HfiD.8w)
nmat <- nmat[unique(row.names(nmat)),]

nmat <- nmat[, colnames(nmat) %in% colnames(emat)]
nmat <- nmat[rownames(nmat) %in% row.names(mWAT_ASPC_HfiD_merged), ]


cell.type <- cell.chose[colnames(emat),1]


#确保细胞数和emat，nmat一致
mWAT_ASPC_HfiD_merged <- subset(mWAT_ASPC_HfiD_merged, cells = colnames(emat))
mWAT_ASPC_HfiD_merged <- subset(mWAT_ASPC_HfiD_merged, features = rownames(emat))

#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(mWAT_ASPC_HfiD_merged@reductions$umap@cell.embeddings)))

#获取每个细胞在图上的坐标位置
emb <- mWAT_ASPC_HfiD_merged@reductions$umap@cell.embeddings
pca <- mWAT_ASPC_HfiD_merged@reductions$pca@cell.embeddings

Generate_scvelo_InputFiles(mWAT_ASPC_HfiD_merged,OutPrefix = paste0(mWAT_ASPC_HfiD_velocyto_prefix,"scvelo_files_"))

counts_matrix <- GetAssayData(mWAT_ASPC_HfiD_merged, assay='RNA', slot='counts')

genes_attributes <- data.frame('gene'=rownames(counts_matrix))

ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(mWAT_ASPC_HfiD_merged@meta.data)
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
#clusters <- as.factor(adata$obs$Myannotation_ASPC_ASPC)
#rownames(emb) <- names(clusters) <- as.matrix(adata$obs_names$values)


## run scvelo dynamic model
scv$pp$filter_genes(adata) ## filter
scv$pp$moments(adata) ## normalize and compute moments
scv$tl$recover_dynamics(adata) ## model

## plot (creates pop up window)
scv$tl$velocity(adata, mode='dynamical')
scv$tl$velocity_graph(adata)
pdf(paste0(mWAT_ASPC_HfiD_velocyto_prefix,"scvelo_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color='Myannotation_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300)
dev.off()

pdf(paste0(mWAT_ASPC_HfiD_velocyto_prefix,"scvelo_PCA_",Sys.Date(),".pdf"), width = 9, height = 9)
scv$pl$velocity_embedding_stream(adata, basis='pca',
                                 color='Myannotation_ASPC',
                                 min_mass=3,
                                 arrow_size=1.3,
                                 density=1.8,
                                 linewidth=1.5, size=300) ## other embedding
dev.off()

##save
adata$write(paste0(mWAT_ASPC_HfiD_velocyto_prefix,'scvelo_anndata.h5ad'), compression='gzip')
#adata = scv$read(paste0(mWAT_ASPC_HfiD_velocyto_prefix,'scvelo_anndata.h5ad'))




