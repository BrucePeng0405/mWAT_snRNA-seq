
source("PZ_ScFunctions.R")

####################0 Basic stats#############################
GetCellComposition(mWAT_merged)
GetCellComposition(mWAT_Ad_merged)
GetCellComposition(mWAT_ASPC_merged)
GetCellComposition(mWAT_Immune_merged)

####################1 mWAT load & integration ##################
mWAT_prefix <- "All_Cells/mWAT_"
if(file.exists(paste0(mWAT_prefix,"Seurat.rds"))){
  mWAT_merged <- readRDS(paste0(mWAT_prefix,"Seurat.rds"))
} else {

  CD_HFD0w <- Read10X('raw/mWAT_CD/filtered_feature_bc_matrix/')
  CD_HFD0w_Obj <- CreateSeuratObject(counts = CD_HFD0w, project = "CD_HFD0w", min.cells = 10)
  CD_HFD4w <- Read10X('raw/mWAT_CD-4HFD/filtered_feature_bc_matrix/')
  CD_HFD4w_Obj <- CreateSeuratObject(counts = CD_HFD4w, project = "CD_HFD4w", min.cells = 10)
  CD_HFD8w <- Read10X('raw/mWAT_CD-8HFD/filtered_feature_bc_matrix/')
  CD_HFD8w_Obj <- CreateSeuratObject(counts = CD_HFD8w, project = "CD_HFD8w", min.cells = 10)
  HfiD_HFD0w <- Read10X('raw/mWAT_HFiD/filtered_feature_bc_matrix/')
  HfiD_HFD0w_Obj <- CreateSeuratObject(counts = HfiD_HFD0w, project = "HfiD_HFD0w", min.cells = 10)
  HfiD_HFD4w <- Read10X('raw/mWAT_HFiD-4HFD/filtered_feature_bc_matrix/')
  HfiD_HFD4w_Obj <- CreateSeuratObject(counts = HfiD_HFD4w, project = "HfiD_HFD4w", min.cells = 10)
  HfiD_HFD8w <- Read10X('raw/mWAT_HFiD-8HFD/filtered_feature_bc_matrix/')
  HfiD_HFD8w_Obj <- CreateSeuratObject(counts = HfiD_HFD8w, project = "HfiD_HFD8w", min.cells = 10)

  
  CD_HFD0w <- ReadCB_h5("raw/mWAT_CD/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  CD_HFD0w_Obj <-   CreateSeuratObject(counts = CD_HFD0w, project = "CD_HFD0w", min.cells = 10)
  CD_HFD4w <- ReadCB_h5("raw/mWAT_CD-4HFD/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  CD_HFD4w_Obj <-   CreateSeuratObject(counts = CD_HFD4w, project = "CD_HFD4w", min.cells = 10)
  CD_HFD8w <- ReadCB_h5("raw/mWAT_CD-8HFD/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  CD_HFD8w_Obj <-   CreateSeuratObject(counts = CD_HFD8w, project = "CD_HFD8w", min.cells = 10)
  HfiD_HFD0w <- ReadCB_h5("raw/mWAT_HFiD/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  HfiD_HFD0w_Obj <-   CreateSeuratObject(counts = HfiD_HFD0w, project = "HfiD_HFD0w", min.cells = 10)
  HfiD_HFD4w <- ReadCB_h5("raw/mWAT_HFiD-4HFD/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  HfiD_HFD4w_Obj <-   CreateSeuratObject(counts = HfiD_HFD4w, project = "HfiD_HFD4w", min.cells = 10)
  HfiD_HFD8w <- ReadCB_h5("raw/mWAT_HFiD-8HFD/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  HfiD_HFD8w_Obj <-   CreateSeuratObject(counts = HfiD_HFD8w, project = "HfiD_HFD8w", min.cells = 10)
  
  mWAT_merged <- RunSeurat(ObjList = c(CD_HFD0w_Obj, CD_HFD4w_Obj,CD_HFD8w_Obj,
                                       HfiD_HFD0w_Obj,HfiD_HFD4w_Obj,HfiD_HFD8w_Obj), 
                           OutPrefix = mWAT_prefix, NormalizationMethod = "VST",
                           removeDoublet = F,
                           Filter_nFeature_RNA = c(100, 8000), 
                           Filter_nCount_RNA = c(100, 50000), 
                           MitoPercent = 15)
  
  #load("raw/batchremoved_raw.Rdata")
  #mWAT_merged <- all
  #mWAT_merged <- RunPCA(mWAT_merged,npcs = 20) %>% RunUMAP(dims = 1:20) %>% 
  #  FindNeighbors(dims = 1:20) %>% FindClusters(resolution = 0.3)
  
  Idents(mWAT_merged) <- mWAT_merged$integrated_snn_res.0.25
  DimPlot(mWAT_merged, reduction = "umap", ncol = 2, label = F, split.by = "orig.ident", repel = T, label.box =T, order = T)
  ggsave(paste0(mWAT_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 12, height = 12)
  Plot_Cell_compoistion(Seurat_Obj = mWAT_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = mWAT_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = mWAT_merged, ColorUse = hughie_color, OutPrefix = mWAT_prefix)
}

PlotGenes <- c("Plin1","Acvr1c","Adipoq",#Adipocyte 7
               "Pdgfra","Nova1", # ASPC 11
               "Muc16","Msln", # Mesothelial Cell 10 8
               "Cyyr1","Cdh5",# Endothelial Cell 17
               "Prox1","Stab2","F8",#Lymphatic endo 14
               "Synpo2","Abcc9",#Pericyte 13 19
               "Ptprn2","Dmbt1", #Enteroendocrine cell 12 21
               "Ptprc","Dock2"# Pan immune others
)

(p1 <- VlnPlot(mWAT_merged, features = PlotGenes, fill.by = "ident", stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
#(p2 <- VlnPlot(mWAT_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mWAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1
ggsave(paste0(mWAT_prefix, "Before_annoataion_Selected_markers_VlnPlot_",Sys.Date(),".pdf"), width = 7, height = 8)

mWAT_merged <- RenameIdents(mWAT_merged,
                            `7` = "Adipocyte", 
                            `11` = "ASPC",
                            `10` = "Mesothelial cell",`8` = "Mesothelial cell",
                            `17` = "Endothelial cell", `14` = "Lymphatic endo",
                            `13` = "Pericyte",`19` = "Pericyte", 
                            `12`= "Enteroendocrine cell", `21` = "Enteroendocrine cell",
                            `0` = "Immune cell", `1`= "Immune cell",`2`= "Immune cell",
                            `3`= "Immune cell",`4`= "Immune cell",`5`= "Immune cell",
                            `6`= "Immune cell",`9`= "Immune cell",`15`= "Immune cell", 
                            `16`= "Immune cell", `18`= "Immune cell", `20`= "Immune cell"
                            )
mWAT_merged <- RenameIdents(mWAT_merged,
                            `11` = "Adipocyte", 
                            `9` = "ASPC",
                            `5` = "Mesothelial cell",
                            `13` = "Endothelial cell",
                            `21` = "Pericyte",
                            `6`= "Enteroendocrine cell", 
                            `0` = "Immune cell", `1`= "Immune cell",`2`= "Immune cell",
                            `3`= "Immune cell",`4`= "Immune cell",`7`= "Immune cell",
                            `8`= "Immune cell",`10`= "Immune cell",`12`= "Immune cell", 
                            `14`= "Immune cell", `15`= "Immune cell", `16`= "Immune cell",
                            `17` = "Immune cell", `18`= "Immune cell", `19`= "Immune cell",
                            `20` = "Immune cell"
                            )
mWAT_merged$Myannotation <- mWAT_merged@active.ident

mWAT_cellRank <- c("Adipocyte","ASPC","Mesothelial cell","Endothelial cell","Lymphatic endo",
                   "Pericyte","Enteroendocrine cell","Immune cell")
mWAT_colors <- c( "#AAD0AC","#EED0E0","#EFE2AA","#83B4EF", "#DBC9B3","#8ECFF8","#7EC0C3","#EBAEA9")
                  #,"#BFA6C9","#F5E0BA","#AED0DF","#89B780","#F5D8D0","#CB95BB","#808180FF")

mWAT_SAMcolors <- c("#A6CEE3", "#66A5CC","#267CB6",
                     "#A2D48E","#47A93A","#20854EFF")
names(mWAT_colors) <- mWAT_cellRank

Idents(mWAT_merged) <- factor(Idents(mWAT_merged), levels= mWAT_cellRank)
mWAT_merged$orig.ident <- factor(mWAT_merged$orig.ident, levels = c("CD_HFD0w","CD_HFD4w","CD_HFD8w",
                                                                    "HfiD_HFD0w","HfiD_HFD4w","HfiD_HFD8w"
                                                                    ))


PlotGenes <- c("Plin1","Acvr1c","Adipoq",#Adipocyte 7
               "Pdgfra","Nova1", # ASPC 11
               "Muc16","Msln", # Mesothelial Cell 10 8
               #"Cyyr1","Cdh5",# Endothelial Cell 17
               "Prox1","Stab2","F8",#Lymphatic endo 14
               "Synpo2","Abcc9",#Pericyte 13 19
               "Ptprn2","Dmbt1", #Enteroendocrine cell 12 21
               "Ptprc","Dock2"# Pan immune others
)

(p1 <- VlnPlot(mWAT_merged, features = PlotGenes, fill.by = "ident", cols = mWAT_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
#(p2 <- VlnPlot(mWAT_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mWAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1
ggsave(paste0(mWAT_prefix, "After_annotation_Selected_markers_VlnPlot_",Sys.Date(),".pdf"), width = 7, height = 8)

pdf(paste0(mWAT_prefix, "UMAP-ByClusters_",Sys.Date(),".pdf"), width = 14, height = 6)
p1 <- DimPlot(mWAT_merged, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T, cols = mWAT_colors, label.color = "white")
p2 <- DimPlot(mWAT_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = mWAT_SAMcolors)
p1 + p2
dev.off()
#ggsave(paste0(mWAT_prefix, "UMAP-ByClusters_",Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(mWAT_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = mWAT_colors, label.color = "white", order = T)
ggsave(paste0(mWAT_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 7, height = 9)

Plot_Cell_compoistion(Seurat_Obj = mWAT_merged, OutPrefix = mWAT_prefix, ColorUse1 = mWAT_colors, ColorUse2 = mWAT_SAMcolors)
Plot_Cell_compoistion(Seurat_Obj = subset(mWAT_merged, idents = c("Adipocyte","ASPC")), OutPrefix = paste0(mWAT_prefix,"Adipo_"), ColorUse1 = mWAT_colors, ColorUse2 = mWAT_SAMcolors)

TopMarkersGenes <- MarkersPlot(Seurat_Obj = mWAT_merged, ColorUse = mWAT_colors, OutPrefix = mWAT_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(mWAT_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4, cellRanks = mWAT_cellRank)
PlotObjMetrices(Seurat_Obj = mWAT_merged, OutPrefix= mWAT_prefix, ColorUse = mWAT_SAMcolors)


saveRDS(mWAT_merged, file = paste0(mWAT_prefix,"h5file_Seurat.rds"))

########################2 mWAT-Adipocytes############################
mWAT_Ad_prefix <- "Adipocyte/mWAT_Adipocytes_"
if(file.exists(paste0(mWAT_Ad_prefix,"Seurat.rds"))){
  mWAT_Ad_merged <- readRDS(paste0(mWAT_Ad_prefix,"Seurat.rds"))
} else {
  #houyu pepline
  mWAT_Ad_merged <- ScaleData(subset(mWAT_merged, idents = c("Adipocyte"))) %>% RunPCA(npcs = 20) %>% RunUMAP(dims = 1:20) %>% 
    FindNeighbors(dims = 1:20) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  #my pepline
  mWAT_Ad_merged <- FindNeighbors(subset(mWAT_merged, idents = c("Adipocyte")), dims = 1:20) %>%
    FindClusters(resolution = seq(from=0, by=0.05, length=10)) %>% RunUMAP(dims = 1:20)
  
  clustree(mWAT_Ad_merged)
  ggsave(paste0(mWAT_Ad_prefix, "clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
  
  
  
  Idents(mWAT_Ad_merged) <- mWAT_Ad_merged$RNA_snn_res.0.25
  DimPlot(mWAT_Ad_merged, reduction = "umap", label = F, repel = T, label.box = T, order = T, cols = hughie_color)
  ggsave(paste0(mWAT_Ad_prefix,"UMAP_",Sys.Date(),".pdf"), width = 6, height = 5)
  #Plot_Cell_compoistion(Seurat_Obj = mWAT_Ad_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = mWAT_Ad_prefix)
  #TopMarkersGenes <- MarkersPlot(Seurat_Obj = mWAT_Ad_merged, ColorUse = hughie_color, OutPrefix = mWAT_Ad_prefix, NMarkers = 40)
  
}

PlotGenes <- c("Adrb3","Lpl","Cd36","Adipoq",
               "Anxa2","Fasn",#lPRAT-ad1
               "Cmss1","Lars2","Camk1d","Cdk8",#lPRAT-ad4
               "Pparg","Acaca",#LGA
               "Nnat","Lrp3",#LSA
               "Cd36","Hif1a","Nedd9",#SLSA
               "Prune2",#Nature-mAd5
               "Mt2"#Nature-mAd6
               )
PlotGenes <- c("Hif1a","Nedd9","Rab7","Lep")

PlotGenes <- c("Lep","Palld","Hif1a","Adipoq", #SHA
               "Cmss1","Lars2","Camk1d","Cdk8", #Transient adipocytes
               "Acaca","Acly" #Lipogenensis aipocyte
               )
PlotGenes <- c("Plin1","Acvr1c","Adipoq","Pdgfra","Pdgfrb","Icam1","Nova1")
(p1 <- VlnPlot(mWAT_Ad_merged, features = PlotGenes, fill.by = "ident", cols = mWAT_Ad_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
p1
ggsave(paste0(mWAT_Ad_prefix, "Adipocyte_ASPC_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)


mWAT_Ad_colors <- c("#BFA6C9","#F5E0BA","#AED0DF")
(p1 <- VlnPlot(mWAT_Ad_merged, features = PlotGenes, fill.by = "ident", cols = mWAT_Ad_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(mWAT_Ad_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mWAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(mWAT_Ad_prefix, "Before_annotation_Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)

mWAT_Ad_merged <- RenameIdents(mWAT_Ad_merged, `0` = "mWAT-ad1", `1` = "mWAT-ad2", `2` = "mWAT-ad3")
mWAT_Ad_merged$Myannotation_Ad <- mWAT_Ad_merged@active.ident
mWAT_AdRanks <- c("mWAT-ad1","mWAT-ad2","mWAT-ad3")

#"#89B780,"#F5D8D0","#CB95BB","#808180FF")
names(mWAT_Ad_colors) <- mWAT_AdRanks
Idents(mWAT_Ad_merged) <- factor(Idents(mWAT_Ad_merged), levels = mWAT_AdRanks)

PlotGenes <- c("Lep","Adipoq", "Ebf2",#mAd1 Classic & Thermogenesis
               "Cmss1","Lars2","Camk1d","Cdk8", #mAd2 Oppose lipogenesis & control lipid handing
               "Acaca","Elovl6"#mAd3 Lipid biosynthesis
               )

(p1 <- VlnPlot(mWAT_Ad_merged, features = PlotGenes, fill.by = "ident", cols = mWAT_Ad_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(mWAT_Ad_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mWAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(mWAT_Ad_prefix, "SLSA_annotation_Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 14, height = 6)

(p1 <- DimPlot(mWAT_Ad_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = mWAT_Ad_colors, label.color = "white"))
(p2 <- DimPlot(mWAT_Ad_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = mWAT_SAMcolors))
p1 + p2
ggsave(paste0(mWAT_Ad_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(mWAT_Ad_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = mWAT_Ad_colors, label.color = "white", order = T)
ggsave(paste0(mWAT_Ad_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 7, height = 9)

Plot_Cell_compoistion(Seurat_Obj = mWAT_Ad_merged, OutPrefix = mWAT_Ad_prefix, ColorUse1 = mWAT_Ad_colors, ColorUse2 = mWAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = mWAT_Ad_merged, ColorUse = mWAT_Ad_colors, OutPrefix = mWAT_Ad_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(mWAT_Ad_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4,cellRanks = mWAT_AdRanks)

FeaturePlot(mWAT_Ad_merged, features =  c("Lpl","Acly","Anxa2","Cdk8"), min.cutoff = "q10", order = T, pt.size=0.1, ncol = 2)
ggsave(paste0(mWAT_Ad_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 6, height = 6)

saveRDS(mWAT_Ad_merged, file = paste0(mWAT_Ad_prefix,"Seurat.rds"))

###########markers
PlotGenes <- c( "Lrp5",  
               "Anxa1", "Adam9",    "Met", "Ptprk" ,
               "Lipe","Cdk8","Camk1d",
               "Acly", "Pdk1", "Elovl6", "Acat2"
               
)



mWAT_Ad_colors <- c("#BFA6C9","#F5E0BA","#AED0DF")
(p1 <- VlnPlot(mWAT_Ad_merged, features = PlotGenes, fill.by = "ident", cols = mWAT_Ad_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
#(p2 <- VlnPlot(mWAT_Ad_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mWAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1
ggsave(paste0(mWAT_Ad_prefix, "Selected_final_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)

PlotGenes <- c("Lep","Abcg1","Cd36","Hif1a","Anxa1","Lipe","Cdk8","Camk1d", "Aldh6a1","Apoe","Acly" ,"Dlat","Elovl6")
FeatureHeatmapPlot(Seurat_Obj = mWAT_Ad_merged,OutPrefix = mWAT_Ad_prefix,features = PlotGenes,cols = mWAT_Ad_colors, group.by = "Myannotation_Ad")

Selected_enrichment_Plot(GO_enrichment_file = "Adipocyte/mWAT_Adipocytes_Selectedmarkers_2024-08-15_GO-BP_list_2024-08-15.csv",
                         OutPrefix = mWAT_Ad_prefix)


################2.1.1 mWAT adipocytes Merge with Mandurp adipocytes  ################
AdMerged_prefix <- "Adipocyte_Merge/Mandrup2021_Ad_"
if(file.exists(paste0(Merged_prefix,"Seurat.rds"))){
  AdMerged_merged <- readRDS(paste0(AdMerged_prefix,"Seurat.rds"))
} else {
  MandrupAd <- readRDS("raw/Mandrup_Adipocytes_Seurat.rds")
  ObjList = c(mWAT_Ad_merged, MandrupAd)
  features <- SelectIntegrationFeatures(object.list = ObjList, nfeatures = 2000)
  anchors <-  FindIntegrationAnchors(object.list = ObjList, anchor.features = features)
  combined <- IntegrateData(anchorset = anchors)
  #houyu pepline
  AdMerged_merged <- ScaleData(combined) %>% RunPCA(npcs = 20) %>% RunUMAP(reduction = "pca", dims = 1:20) %>% 
    FindNeighbors(reduction = "pca", dims = 1:20) %>% FindClusters(resolution = seq(from=0, by=0.05, length=5))

  #my pepline
  AdMerged_merged <- ScaleData(combined) %>% RunPCA(npcs = 20) %>% FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = seq(from=0, by=0.05, length=10)) %>% RunUMAP(dims = 1:20)

  AdMerged_merged$orig.ident <- gsub("^.{8}$", "mWAT", AdMerged_merged$orig.ident)
  AdMerged_merged$orig.ident <- gsub("^.{10}$", "mWAT", AdMerged_merged$orig.ident)
  AdMerged_merged$orig.ident <- gsub("^.{6}$", "eWAT", AdMerged_merged$orig.ident)
  AdMerged_merged$orig.ident <- as.factor(AdMerged_merged$orig.ident)
  saveRDS(AdMerged_merged, file = paste0(AdMerged_prefix,"Seurat.rds"))
}

PlotGenes <- c("Plin1","Ucp1","Cidea","Adipoq","Cyp2e1","Aldh1a1",
               "Pparg","Acaca","Elovl6","Acly","Igf2r",
               "Nnat","Lrp3","Car3","Abcd2",
               "Hif1a","Nedd9","Gadd45g","Lep","Rab7")

Idents(AdMerged_merged) <- AdMerged_merged$integrated_snn_res.0.15
(p00 <- DimPlot(AdMerged_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = hughie_color, label.color = "white"))
(p01 <- DimPlot(AdMerged_merged, reduction = "umap", group.by = "orig.ident", label = F, repel = T, label.box = T, order = T, cols = c("#DC0000FF","#1F78B4")))
(p02 <- DimPlot(AdMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = hughie_color))
(pa <- VlnPlot(AdMerged_merged, features = PlotGenes, fill.by = "ident",cols = c("#DC0000FF","#1F78B4"), stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))
(p00 + p01)/p02
ggsave(paste0(AdMerged_prefix,"UMAP-mergedIdentity_",Sys.Date(),".pdf"), width = 12, height = 12)
Plot_Cell_compoistion(Seurat_Obj = AdMerged_merged, OutPrefix = AdMerged_prefix, ColorUse1 = rev(hughie_color[1:7]), ColorUse2 = c("#DC0000FF","#1F78B4"))

Idents(AdMerged_merged) <- AdMerged_merged$Subtype
Idents(AdMerged_merged) <- factor(Idents(AdMerged_merged), levels = c("LGA","LSA","SLSA"))
(p10 <- DimPlot(AdMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = c("#00B838","#F7766D","#629BFE")))
(pb <- VlnPlot(AdMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA")) 

Idents(AdMerged_merged) <- AdMerged_merged$Myannotation_Ad
Idents(AdMerged_merged) <- factor(Idents(AdMerged_merged), levels = mWAT_AdRanks)
(p20 <- DimPlot(AdMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = mWAT_Ad_colors))
(pc <- VlnPlot(AdMerged_merged, features = PlotGenes, fill.by = "ident",cols = mWAT_Ad_colors, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

(p00+p01)/(p10)/(p20)
ggsave(paste0(AdMerged_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 10, height = 15)

###################2.1.2 mWAT adipocytes Merge with eLife adipocytes #######################
AdMerged_prefix <- "Adipocyte_Merge/Houyu2024_Ad_"
if(file.exists(paste0(AdMerged_prefix,"Seurat.rds"))){
  AdMerged_merged <- readRDS(paste0(AdMerged_prefix,"Seurat.rds"))
} else {
  HouyuAd <- readRDS(paste0("raw/lPRAT_Adipocytes_Seurat.rds"))
  HouyuAd <- subset(HouyuAd, features = HouyuAd@assays[["integrated"]]@var.features)
  
  ObjList1 = c(mWAT_Ad_merged, HouyuAd)
  features1 <- SelectIntegrationFeatures(object.list = ObjList1, nfeatures = 2000)
  anchors1 <-  FindIntegrationAnchors(object.list = ObjList1, anchor.features = features1)
  combined <- IntegrateData(anchorset = anchors1)
  #houyu pepline
  AdMerged_merged <- ScaleData(combined) %>% RunPCA(npcs = 20) %>% RunUMAP(reduction = "pca", dims = 1:20) %>% 
    FindNeighbors(reduction = "pca", dims = 1:20) %>% FindClusters(resolution = seq(from=0, by=0.05, length=5))
  #my pepline
  AdMerged_merged <- ScaleData(combined) %>% RunPCA(npcs = 20) %>% FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = seq(from=0, by=0.05, length=10)) %>% RunUMAP(dims = 1:20)
  
  AdMerged_merged$orig.ident <- gsub("^.{8}$", "mWAT", AdMerged_merged$orig.ident)
  AdMerged_merged$orig.ident <- gsub("^.{10}$", "mWAT", AdMerged_merged$orig.ident)
  AdMerged_merged$orig.ident <- gsub("^.{9}$", "lPRAT", AdMerged_merged$orig.ident)
  
  AdMerged_merged$orig.ident <- as.factor(AdMerged_merged$orig.ident)
  saveRDS(AdMerged_merged, file = paste0(AdMerged_prefix,"Seurat.rds"))
}

PlotGenes <- c("Ebf2","Lep","Adipoq", #Lipolytic
               "Acaca","Elvl6",#Lipogenesis
               "Cmss1","Lars2","Camk1d","Cdk8"
)

Idents(AdMerged_merged) <- AdMerged_merged$integrated_snn_res.0.2
(p00 <- DimPlot(AdMerged_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = hughie_color, label.color = "white"))
(p01 <- DimPlot(AdMerged_merged, reduction = "umap", group.by = "orig.ident", label = F, repel = T, label.box = T, order = T, cols = c("#DC0000FF","#1F78B4")))
(p02 <- DimPlot(AdMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = hughie_color))
(pa <- VlnPlot(AdMerged_merged, features = PlotGenes, fill.by = "ident",cols = c("#DC0000FF","#1F78B4"), stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))
(p00 + p01)/p02
ggsave(paste0(AdMerged_prefix,"UMAP-mergedIdentity_",Sys.Date(),".pdf"), width = 12, height = 12)
Plot_Cell_compoistion(Seurat_Obj = AdMerged_merged, OutPrefix = AdMerged_prefix, ColorUse1 = rev(hughie_color[1:6]), ColorUse2 = c("#DC0000FF","#1F78B4"))

Idents(AdMerged_merged) <- AdMerged_merged$Subtype
Idents(AdMerged_merged) <- factor(Idents(AdMerged_merged), levels = c("lPRAT-ad1","lPRAT-ad2","lPRAT-ad3","lPRAT-ad4"))
(p10 <- DimPlot(AdMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = c("#729E4E","#FDBF6F","#A20056FF","#EE4C97FF")))
(pb <- VlnPlot(AdMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA")) 

Idents(AdMerged_merged) <- AdMerged_merged$Myannotation_Ad
Idents(AdMerged_merged) <- factor(Idents(AdMerged_merged), levels = mWAT_AdRanks)
(p20 <- DimPlot(AdMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = mWAT_Ad_colors))
(pc <- VlnPlot(AdMerged_merged, features = PlotGenes, fill.by = "ident",cols = mWAT_Ad_colors, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

(p00+p01)/(p10)/(p20)
ggsave(paste0(AdMerged_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 10, height = 15)

####################2.1.3 mWAT adipocytes Merge with Nature adipocytes ########################
AdMerged_prefix <- "Adipocyte_Merge/RosenNature2022_Ad_"
if(file.exists(paste0(AdMerged_prefix,"Seurat.rds"))){
  AdMerged_merged <- readRDS(paste0(AdMerged_prefix,"Seurat.rds"))
} else {
  RosenAd <- readRDS(paste0("raw/iWAT_Adipocytes_Seurat.rds"))
  #RosenAd <- subset(RosenAd, features = RosenAd@assays[["integrated"]]@var.features)
  
  ObjList1 = c(mWAT_Ad_merged, RosenAd)
  features1 <- SelectIntegrationFeatures(object.list = ObjList1, nfeatures = 2000)
  anchors1 <-  FindIntegrationAnchors(object.list = ObjList1, anchor.features = features1)
  combined <- IntegrateData(anchorset = anchors1)
  #houyu pepline
  AdMerged_merged <- ScaleData(combined) %>% RunPCA(npcs = 20) %>% RunUMAP(reduction = "pca", dims = 1:20) %>% 
    FindNeighbors(reduction = "pca", dims = 1:20) %>% FindClusters(resolution = seq(from=0, by=0.05, length=5))
  
  #my pepline
  AdMerged_merged <- ScaleData(combined) %>% RunPCA(npcs = 20) %>% FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = seq(from=0, by=0.05, length=10)) %>% RunUMAP(dims = 1:20)
  
  AdMerged_merged$orig.ident <- gsub("^.{8}$", "mWAT", AdMerged_merged$orig.ident)
  AdMerged_merged$orig.ident <- gsub("^.{10}$", "mWAT", AdMerged_merged$orig.ident)
  AdMerged_merged$orig.ident <- gsub("^.{2}$", "iWAT", AdMerged_merged$orig.ident)
  
  AdMerged_merged$orig.ident <- as.factor(AdMerged_merged$orig.ident)
  saveRDS(AdMerged_merged, file = paste0(AdMerged_prefix,"Seurat.rds"))
}

PlotGenes <- c("Ebf2","Lep","Adipoq", #Lipolytic
               "Acaca","Elvl6",#Lipogenesis
               "Cmss1","Lars2","Camk1d","Cdk8"
)

Idents(AdMerged_merged) <- AdMerged_merged$integrated_snn_res.0.05
(p00 <- DimPlot(AdMerged_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = hughie_color, label.color = "white"))
(p01 <- DimPlot(AdMerged_merged, reduction = "umap", group.by = "orig.ident", label = F, repel = T, label.box = T, order = T, cols = c("#DC0000FF","#1F78B4")))
(p02 <- DimPlot(AdMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = hughie_color))
(pa <- VlnPlot(AdMerged_merged, features = PlotGenes, fill.by = "ident",cols = c("#DC0000FF","#1F78B4"), stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))
(p00 + p01)/p02
ggsave(paste0(AdMerged_prefix,"UMAP-mergedIdentity_",Sys.Date(),".pdf"), width = 12, height = 12)
Plot_Cell_compoistion(Seurat_Obj = AdMerged_merged, OutPrefix = AdMerged_prefix, ColorUse1 = rev(hughie_color[1:6]), ColorUse2 = c("#DC0000FF","#1F78B4"))

Idents(AdMerged_merged) <- AdMerged_merged$Subtype
Idents(AdMerged_merged) <- factor(Idents(AdMerged_merged), levels = c("iWAT-mAd1","iWAT-mAd2","iWAT-mAd3","iWAT-mAd4","iWAT-mAd5","iWAT-mAd6"))
(p10 <- DimPlot(AdMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = c("#FFDF01","#FDAF05","#FC7E0D","#D1362D","#BB1B4B","#970667")))
(pb <- VlnPlot(AdMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA")) 

Idents(AdMerged_merged) <- AdMerged_merged$Myannotation_Ad
Idents(AdMerged_merged) <- factor(Idents(AdMerged_merged), levels = mWAT_AdRanks)
(p20 <- DimPlot(AdMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = mWAT_Ad_colors))
(pc <- VlnPlot(AdMerged_merged, features = PlotGenes, fill.by = "ident",cols = mWAT_Ad_colors, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

(p00+p01)/(p10)/(p20)
ggsave(paste0(AdMerged_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 10, height = 15)


#####################2.2 mWAT adipocytes monocle###################
AdMonocle_prefix <- "Adipocyte_monocle/mWAT_Adipocytes_"

#all smaple (CD & HfiD groups)
Idents(mWAT_Ad_merged) <- "Myannotation_Ad"
data <- as(as.matrix(mWAT_Ad_merged@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = mWAT_Ad_merged@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds)                        
mycds <- detectGenes(mycds, min_expr = 0.1)
mycds_expressed_genes <- row.names(mycds)
diff_test_res <- differentialGeneTest(mycds[mycds_expressed_genes,],
                                      fullModelFormulaStr = "~Myannotation_Ad")                      
diff_test_res[1:4,1:4]
write.csv(diff_test_res, file=paste0(AdMonocle_prefix,"diff_test_result.csv",sep=""))

ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
ordering_genes[1:5]
mycds <- setOrderingFilter(mycds, ordering_genes)
plot_ordering_genes(mycds)
ggsave(paste0(AdMonocle_prefix,"Ordering.genes_", Sys.Date(), ".pdf"), height = 7, width = 8)

mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree',verbose = T)
mycds <- orderCells(mycds)

saveRDS(mycds, file = paste0(AdMonocle_prefix, "mycds.rds"))


data <- GetAssayData(mWAT_Ad_merged, assay = 'RNA', slot = 'counts')
cell_metadata <- mWAT_Ad_merged@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
monocle3_cds <- new_cell_data_set(data,
                                  cell_metadata = cell_metadata,
                                  gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData ScaleData RunPCA
monocle3_cds <- preprocess_cds(monocle3_cds, num_dim = 100)
#umap降维
monocle3_cds <- reduce_dimension(monocle3_cds, preprocess_method = "PCA")
p1 <- plot_cells(monocle3_cds, reduction_method="UMAP", color_cells_by="Myannotation_Ad"
                 #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
)+ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
monocle3_cds.embed <- monocle3_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(mWAT_Ad_merged, reduction = "umap")
int.embed <- int.embed[rownames(monocle3_cds.embed),]
monocle3_cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(monocle3_cds, reduction_method="UMAP", color_cells_by="Myannotation_Ad"
                 #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
)+ggtitle('int.umap')
p1 + p2
ggsave(paste0(AdMonocle_prefix,"Monocle3_umap_embed_", Sys.Date(), ".pdf"), height = 7, width = 16)

## Monocle3聚类分区
monocle3_cds <- cluster_cells(monocle3_cds, resolution = 0.01)
p3 <- plot_cells(monocle3_cds, show_trajectory_graph = FALSE)+ggtitle("Label by clusterID")
p4 <- plot_cells(monocle3_cds, color_cells_by = "partition", show_trajectory_graph = FALSE)+ggtitle("Label by partitionID")
p3 + p4
ggsave(paste0(AdMonocle_prefix,"Monocle3_umap_cluster_", Sys.Date(), ".pdf"), height = 7, width = 16)


## 识别轨迹
monocle3_cds <- learn_graph(monocle3_cds)
p5 <- plot_cells(monocle3_cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                 label_branch_points = FALSE)
p5
ggsave(paste0(AdMonocle_prefix,"Monocle3_umap_trajectory_", Sys.Date(), ".pdf"), width = 5, height = 5)


##细胞按拟时排序
monocle3_cds <- order_cells(monocle3_cds)
p6 <- plot_cells(monocle3_cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                 label_leaves = F,  label_branch_points = FALSE,graph_label_size = 3,
                 cell_size = 1, alpha = 0.5, trajectory_graph_color = "black",trajectory_graph_segment_size = 1)
p6
ggsave(paste0(AdMonocle_prefix,"Monocle3_umap_trajectory_line_", Sys.Date(), ".pdf"), width = 5, height = 5)

#plot heatmap, all smaple (CD & HfiD groups)
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.000001))
nrow(subset(diff_test_res, qval < 0.000001))

bins <- monocle3_scale_to_100(mycds[sig_gene_names,])

pdf(paste0(AdMonocle_prefix,"monocel_heatmap_",Sys.Date(),".pdf"), width = 5, height = 15)

plot_pseudotime_heatmap(mycds[sig_gene_names,],
                        num_clusters = 3,
                        cores = 6,
                        show_rownames = F, return_heatmap = T,
                        add_annotation_col = bins,use_gene_short_name = T
)

dev.off()

df <- data.frame((mycds[sig_gene_names,]@phenoData@data))
df <- df[,c("Pseudotime", "Myannotation_Ad","timepoint")]

#"#A6CEE3", "#A2D48E",
#"#66A5CC","#47A93A",
#"#267CB6","#20854EFF"
color.celtype <- data.frame(celltype = c("mWAT-ad1","mWAT-ad2","mWAT-ad3"), 
                            color = c("#BFA6C9","#F5E0BA","#AED0DF"))  
mWAT_SAMcolors <- data.frame(timepoint = c("HFD0w","HFD4w","HFD8w"),
                             color = c("#A6CEE3","#66A5CC","#267CB6"))

df <- merge(df,color.celtype,by.x = "Myannotation_Ad", by.y = "celltype",drop = F)
df <- merge(df,mWAT_SAMcolors,by.x = "timepoint", by.y = "timepoint",drop = F)
df <- df[order(df$Pseudotime, decreasing = F),]

write.csv(df, file=paste0(AdMonocle_prefix,"ALL_colnames.csv",sep=""))

matrix <- as.matrix(df[,3])
pdf(paste0(AdMonocle_prefix,"ALL_monocel_heatmap_celltype_",Sys.Date(),".pdf"))
image(matrix, col = df$color.x)
dev.off()

pdf(paste0(AdMonocle_prefix,"ALL_monocel_heatmap_timepoint_",Sys.Date(),".pdf"))

image(matrix, col = df$color.y)
dev.off()


#plot heatmap CD group
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.000001))
nrow(subset(diff_test_res, qval < 0.000001))

bins <- monocle3_scale_to_100(mycds[sig_gene_names,grep("CD",mycds$group)])

pdf(paste0(AdMonocle_prefix,"ALL_CD_monocel_heatmap_",Sys.Date(),".pdf"), width = 5, height = 15)

plot_pseudotime_heatmap(mycds[sig_gene_names,grep("CD",mycds$group)],
                        num_clusters = 3,
                        cores = 6,
                        show_rownames = F, return_heatmap = T,
                        add_annotation_col = bins,use_gene_short_name = T
)

dev.off()


df <- data.frame((mycds[sig_gene_names,grep("CD",mycds$group)]@phenoData@data))
df <- df[,c("Pseudotime", "Myannotation_Ad","timepoint")]

#"#A6CEE3", "#A2D48E",
#"#66A5CC","#47A93A",
#"#267CB6","#20854EFF"
color.celtype <- data.frame(celltype = c("mWAT-ad1","mWAT-ad2","mWAT-ad3"), 
                            color = c("#BFA6C9","#F5E0BA","#AED0DF"))  
mWAT_SAMcolors <- data.frame(timepoint = c("HFD0w","HFD4w","HFD8w"),
                             color = c("#A6CEE3","#66A5CC","#267CB6"))

df <- merge(df,color.celtype,by.x = "Myannotation_Ad", by.y = "celltype",drop = F)
df <- merge(df,mWAT_SAMcolors,by.x = "timepoint", by.y = "timepoint",drop = F)
df <- df[order(df$Pseudotime, decreasing = F),]

write.csv(df, file=paste0(AdMonocle_prefix,"ALL_CD_group_colnames.csv",sep=""))

matrix <- as.matrix(df[,3])
pdf(paste0(AdMonocle_prefix,"ALL_CD_group_monocel_heatmap_celltype_",Sys.Date(),".pdf"))
image(matrix, col = df$color.x)
dev.off()

pdf(paste0(AdMonocle_prefix,"ALL_CD_group_monocel_heatmap_timepoint_",Sys.Date(),".pdf"))

image(matrix, col = df$color.y)
dev.off()

#plot heatmap HfiD group
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.000001))
nrow(subset(diff_test_res, qval < 0.000001))

bins <- monocle3_scale_to_100(mycds[sig_gene_names,grep("HfiD",mycds$group)])

pdf(paste0(AdMonocle_prefix,"ALL_HfiD_monocel_heatmap_",Sys.Date(),".pdf"), width = 5, height = 15)

plot_pseudotime_heatmap(mycds[sig_gene_names,grep("HfiD",mycds$group)],
                        num_clusters = 3,
                        cores = 6,
                        show_rownames = F, return_heatmap = T,
                        add_annotation_col = bins,use_gene_short_name = T
)

dev.off()
df <- data.frame((mycds[sig_gene_names,grep("HfiD",mycds$group)]@phenoData@data))
df <- df[,c("Pseudotime", "Myannotation_Ad","timepoint")]

#"#A6CEE3", "#A2D48E",
#"#66A5CC","#47A93A",
#"#267CB6","#20854EFF"
color.celtype <- data.frame(celltype = c("mWAT-ad1","mWAT-ad2","mWAT-ad3"), 
                            color = c("#BFA6C9","#F5E0BA","#AED0DF"))  
mWAT_SAMcolors <- data.frame(timepoint = c("HFD0w","HFD4w","HFD8w"),
                             color = c("#A2D48E","#47A93A","#20854EFF"))

df <- merge(df,color.celtype,by.x = "Myannotation_Ad", by.y = "celltype",drop = F)
df <- merge(df,mWAT_SAMcolors,by.x = "timepoint", by.y = "timepoint",drop = F)
df <- df[order(df$Pseudotime, decreasing = F),]

write.csv(df, file=paste0(AdMonocle_prefix,"ALL_HfiD_group_colnames.csv",sep=""))

matrix <- as.matrix(df[,3])
pdf(paste0(AdMonocle_prefix,"ALL_HfiD_group_monocel_heatmap_celltype_",Sys.Date(),".pdf"))
image(matrix, col = df$color.x)
dev.off()

pdf(paste0(AdMonocle_prefix,"ALL_HfiD_group_monocel_heatmap_timepoint_",Sys.Date(),".pdf"))

image(matrix, col = df$color.y)
dev.off()

#####################2.2.1 CD group mWAT adipocytes monocle###################

Idents(mWAT_Ad_merged) <- "group"
mWAT_Ad_merged_cd <- subset(mWAT_Ad_merged, idents = c("CD"))
Idents(mWAT_Ad_merged_cd) <- "Myannotation_Ad"

#Idents(mWAT_Ad_merged) <- "Myannotation_Ad"
data_cd <- as(as.matrix(mWAT_Ad_merged_cd@assays$RNA@counts), 'sparseMatrix')
pd_cd <- new('AnnotatedDataFrame', data = mWAT_Ad_merged_cd@meta.data)
fData_cd <- data.frame(gene_short_name = row.names(data_cd), row.names = row.names(data_cd))
fd_cd <- new('AnnotatedDataFrame', data = fData_cd)

mycds_cd <- newCellDataSet(data_cd,
                           phenoData = pd_cd,
                           featureData = fd_cd,
                           expressionFamily = negbinomial.size())

mycds_cd <- estimateSizeFactors(mycds_cd)
mycds_cd <- estimateDispersions(mycds_cd)                        
mycds_cd <- detectGenes(mycds_cd, min_expr = 0.1)
mycds_expressed_genes_cd <- row.names(mycds_cd)
diff_test_res_cd <- differentialGeneTest(mycds_cd[mycds_expressed_genes_cd,],
                                         fullModelFormulaStr = "~Myannotation_Ad")                      
diff_test_res_cd[1:4,1:4]
write.csv(diff_test_res_cd, file=paste0(AdMonocle_prefix,"CD_group_diff_test_result.csv",sep=""))
#diff_test_res_cd <- read.csv(file=paste0(AdMonocle_prefix,"CD_group_diff_test_result.csv",sep=""))

ordering_genes_cd <- row.names(subset(diff_test_res_cd, qval < 0.01))
ordering_genes_cd[1:5]
mycds_cd <- setOrderingFilter(mycds_cd, ordering_genes_cd)
plot_ordering_genes(mycds_cd)
ggsave(paste0(AdMonocle_prefix,"CD_group_Ordering.genes_", Sys.Date(), ".pdf"), height = 7, width = 8)

mycds_cd <- reduceDimension(mycds_cd, max_components = 2, method = 'DDRTree',verbose = T)
mycds_cd <- orderCells(mycds_cd)

saveRDS(mycds_cd, file = paste0(AdMonocle_prefix, "CD_group_mycds.rds"))
#mycds_cd <- readRDS("Adipocyte_monocle/mWAT_Adipocytes_CD_group_mycds.rds")

#plot heatmap of CD group
diff_test_res_cd <- diff_test_res_cd[order(diff_test_res_cd$qval, decreasing = F),]
sig_gene_names_cd <- row.names(diff_test_res_cd[1:600,])
#nrow(subset(diff_test_res_cd, qval < 0.0001))

bins <- monocle3_scale_to_100(mycds_cd[sig_gene_names_cd,])


pdf(paste0(AdMonocle_prefix,"CD_group_monocel_heatmap_",Sys.Date(),".pdf"), width = 5, height = 15)

plot_pseudotime_heatmap(mycds_cd[sig_gene_names_cd,],
                        num_clusters = 3,
                        cores = 6,
                        show_rownames = T, return_heatmap = T,
                        add_annotation_col = bins,use_gene_short_name = T
)

dev.off()

df <- data.frame((mycds_cd[sig_gene_names_cd,]@phenoData@data))
df <- df[,c("Pseudotime", "Myannotation_Ad","timepoint")]

#"#A6CEE3", "#A2D48E",
#"#66A5CC","#47A93A",
#"#267CB6","#20854EFF"
color.celtype <- data.frame(celltype = c("mWAT-ad1","mWAT-ad2","mWAT-ad3"), 
                            color = c("#BFA6C9","#F5E0BA","#AED0DF"))  
mWAT_SAMcolors <- data.frame(timepoint = c("HFD0w","HFD4w","HFD8w"),
                             color = c("#A6CEE3","#66A5CC","#267CB6"))

df <- merge(df,color.celtype,by.x = "Myannotation_Ad", by.y = "celltype",drop = F)
df <- merge(df,mWAT_SAMcolors,by.x = "timepoint", by.y = "timepoint",drop = F)
df <- df[order(df$Pseudotime, decreasing = F),]

write.csv(df, file=paste0(AdMonocle_prefix,"CD_group_colnames.csv",sep=""))

matrix <- as.matrix(df[,3])
pdf(paste0(AdMonocle_prefix,"CD_group_monocel_heatmap_celltype_",Sys.Date(),".pdf"))
image(matrix, col = df$color.x)
dev.off()

pdf(paste0(AdMonocle_prefix,"CD_group_monocel_heatmap_timepoint_",Sys.Date(),".pdf"))

image(matrix, col = df$color.y)
dev.off()


#####################2.2.2 HfiD group mWAT adipocytes monocle###################
Idents(mWAT_Ad_merged) <- "group"
mWAT_Ad_merged_HfiD <- subset(mWAT_Ad_merged, idents = c("HfiD"))
Idents(mWAT_Ad_merged_HfiD) <- "Myannotation_Ad"

#Idents(mWAT_Ad_merged) <- "Myannotation_Ad"
data_HfiD <- as(as.matrix(mWAT_Ad_merged_HfiD@assays$RNA@counts), 'sparseMatrix')
pd_HfiD <- new('AnnotatedDataFrame', data = mWAT_Ad_merged_HfiD@meta.data)
fData_HfiD <- data.frame(gene_short_name = row.names(data_HfiD), row.names = row.names(data_HfiD))
fd_HfiD <- new('AnnotatedDataFrame', data = fData_HfiD)

mycds_HfiD <- newCellDataSet(data_HfiD,
                             phenoData = pd_HfiD,
                             featureData = fd_HfiD,
                             expressionFamily = negbinomial.size())

mycds_HfiD <- estimateSizeFactors(mycds_HfiD)
mycds_HfiD <- estimateDispersions(mycds_HfiD)                        
mycds_HfiD <- detectGenes(mycds_HfiD, min_expr = 0.1)
mycds_expressed_genes_HfiD <- row.names(mycds_HfiD)
diff_test_res_HfiD <- differentialGeneTest(mycds_HfiD[mycds_expressed_genes_HfiD,],
                                           fullModelFormulaStr = "~Myannotation_Ad")                      
diff_test_res_HfiD[1:4,1:4]
write.csv(diff_test_res_HfiD, file=paste0(AdMonocle_prefix,"HfiD_group_diff_test_result.csv",sep=""))
#diff_test_res_HfiD <- read.csv(file=paste0(AdMonocle_prefix,"HfiD_group_diff_test_result.csv",sep=""))

ordering_genes_HfiD <- row.names(subset(diff_test_res_HfiD, qval < 0.01))
ordering_genes_HfiD[1:5]
mycds_HfiD <- setOrderingFilter(mycds_HfiD, ordering_genes_HfiD)
plot_ordering_genes(mycds_HfiD)
ggsave(paste0(AdMonocle_prefix,"HfiD_group_Ordering.genes_", Sys.Date(), ".pdf"), height = 7, width = 8)

mycds_HfiD <- reduceDimension(mycds_HfiD, max_components = 2, method = 'DDRTree',verbose = T)
mycds_HfiD <- orderCells(mycds_HfiD)

saveRDS(mycds_HfiD, file = paste0(AdMonocle_prefix, "HfiD_group_mycds.rds"))
#mycds_HfiD <- readRDS(file = paste0(AdMonocle_prefix, "HfiD_group_mycds.rds"))

#plot heatmap of HfiD group
diff_test_res_HfiD <- diff_test_res_HfiD[order(diff_test_res_HfiD$qval, decreasing = F),]
sig_gene_names_HfiD <- row.names(diff_test_res_HfiD[1:1000,])
nrow(subset(diff_test_res_HfiD, qval < 0.0001))

bins <- monocle3_scale_to_100(mycds_HfiD[sig_gene_names_HfiD,])


pdf(paste0(AdMonocle_prefix,"HfiD_group_monocel_heatmap_",Sys.Date(),".pdf"), width = 5, height = 15)

plot_pseudotime_heatmap(mycds_HfiD[sig_gene_names_HfiD,],
                        num_clusters = 3,
                        cores = 6,
                        show_rownames = T, return_heatmap = T,
                        add_annotation_col = bins,use_gene_short_name = T
)

dev.off()

df <- data.frame((mycds_HfiD[sig_gene_names_HfiD,]@phenoData@data))
df <- df[,c("Pseudotime", "Myannotation_Ad","timepoint")]

#"#A6CEE3", "#A2D48E",
#"#66A5CC","#47A93A",
#"#267CB6","#20854EFF"
color.celtype <- data.frame(celltype = c("mWAT-ad1","mWAT-ad2","mWAT-ad3"), 
                            color = c("#BFA6C9","#F5E0BA","#AED0DF"))  
mWAT_SAMcolors <- data.frame(timepoint = c("HFD0w","HFD4w","HFD8w"),
                             color = c("#A2D48E","#47A93A","#20854EFF"))

df <- merge(df,color.celtype,by.x = "Myannotation_Ad", by.y = "celltype",drop = F)
df <- merge(df,mWAT_SAMcolors,by.x = "timepoint", by.y = "timepoint",drop = F)
df <- df[order(df$Pseudotime, decreasing = F),]

write.csv(df, file=paste0(AdMonocle_prefix,"HfiD_group_colnames.csv",sep=""))

matrix <- as.matrix(df[,3])
pdf(paste0(AdMonocle_prefix,"HfiD_group_monocel_heatmap_celltype_",Sys.Date(),".pdf"))
image(matrix, col = df$color.x)
dev.off()

pdf(paste0(AdMonocle_prefix,"HfiD_group_monocel_heatmap_timepoint_",Sys.Date(),".pdf"))

image(matrix, col = df$color.y)
dev.off()



###################3 mWAT-ASPC ###################################
mWAT_ASPC_prefix <- "ASPC/mWAT_ASPC_"
if(file.exists(paste0(mWAT_ASPC_prefix,"Seurat.rds"))){
  mWAT_ASPC_merged <- readRDS(paste0(mWAT_ASPC_prefix,"Seurat.rds"))
} else {
  #my pepline
  mWAT_ASPC_merged <- FindNeighbors(subset(mWAT_merged, idents = c("ASPC")), dims = 1:20) %>%
    FindClusters(resolution = seq(from=0, by=0.05, length=10)) %>% RunUMAP(dims = 1:20)
  
  Idents(mWAT_ASPC_merged) <- mWAT_ASPC_merged$RNA_snn_res.0.2
  mWAT_ASPC_merged <- subset(mWAT_ASPC_merged, idents = c("0","1","2"))
  DimPlot(mWAT_ASPC_merged, reduction = "umap", label = F, repel = T, label.box = T, order = T, cols = hughie_color)
  ggsave(paste0(mWAT_ASPC_prefix,"UMAP_",Sys.Date(),".pdf"), width = 6, height = 5)
  
  #DimPlot(mWAT_ASPC_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = hughie_color)
  #Plot_Cell_compoistion(Seurat_Obj = mWAT_ASPC_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = mWAT_ASPC_prefix)
  #TopMarkersGenes <- MarkersPlot(Seurat_Obj = mWAT_ASPC_merged, ColorUse = hughie_color, OutPrefix = mWAT_ASPC_prefix)
}

mWAT_ASPC_merged <- RenameIdents(mWAT_ASPC_merged, `0` = "mWAT-ASPC1",`1` = "mWAT-ASPC2", `2` = "mWAT-ASPC3")
mWAT_ASPC_merged$Myannotation_ASPC <- mWAT_ASPC_merged@active.ident
mWAT_ASPCRanks <- c("mWAT-ASPC1","mWAT-ASPC2","mWAT-ASPC3")
mWAT_ASPC_colors <- c("#9b5b33","#f5cfe4","#AED0DF")
names(mWAT_ASPC_colors) <- mWAT_ASPCRanks
Idents(mWAT_ASPC_merged) <- factor(Idents(mWAT_ASPC_merged), levels = mWAT_ASPCRanks)

PlotGenes <- c("Pdgfra","Dpp4","Pi16","Anxa3",
               "Pdgfrb","Pparg","Cd36",
               "Hes1","Lox","Igf1","Lbp","Foxp2",
               "Lpl","Gata4","Fgf10",
               "Ebf1","Rock2","Zfp521","Bmp6","Ebf2",
               "Fbn1","Klf4","Klf2","Fn1","Loxl1")

PlotGenes <- c("Lpl","Pdgfra", #preadipocyte
               "Lars2","Camk1d","Cmss1","Cdk8",#ad2?
               "Dpp4","Pi16","Anxa3"# FAPs
               )

PlotGenes <- c("Plin1","Acvr1c","Adipoq","Pdgfra","Pdgfrb","Icam1","Nova1")
(p1 <- VlnPlot(mWAT_ASPC_merged, features = PlotGenes, fill.by = "ident", cols = mWAT_Ad_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
p1
ggsave(paste0(mWAT_ASPC_prefix, "Adipocyte_ASPC_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)

PlotGenes <- c("Pdgfra","Nova1", 
               "Lpl","Pdgfrb",   #preadipocyte
               "Cdk8","Camk1d","Cmss1",# aspc2
               "Dpp4","Pi16","Anxa3" #FAPs
)

(p1 <- VlnPlot(mWAT_ASPC_merged, features = PlotGenes, fill.by = "ident", cols = mWAT_ASPC_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(mWAT_ASPC_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mWAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(mWAT_ASPC_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)
FeatureHeatmapPlot(Seurat_Obj = mWAT_ASPC_merged,OutPrefix = mWAT_ASPC_prefix,features = PlotGenes,cols = mWAT_ASPC_colors, group.by = "Myannotation_ASPC")


(p1 <- DimPlot(mWAT_ASPC_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = mWAT_ASPC_colors, label.color = "white"))
(p2 <- DimPlot(mWAT_ASPC_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = mWAT_SAMcolors))
p1 + p2
ggsave(paste0(mWAT_ASPC_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(mWAT_ASPC_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = F, repel = T, label.box =T, cols = mWAT_ASPC_colors, label.color = "white", order = T)
ggsave(paste0(mWAT_ASPC_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 9, height = 5)

Plot_Cell_compoistion(Seurat_Obj = mWAT_ASPC_merged, OutPrefix = mWAT_ASPC_prefix, ColorUse1 = mWAT_ASPC_colors, ColorUse2 = mWAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = mWAT_ASPC_merged, ColorUse = mWAT_ASPC_colors, OutPrefix = mWAT_ASPC_prefix,NMarkers = 50)
DEG_enrichment(Seruat_DEG_file = paste0(mWAT_ASPC_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4, cellRanks = mWAT_ASPCRanks)

FeaturePlot(mWAT_ASPC_merged, features =  c("Ucp1","Cidea","Cyp2e1","Aldh1a1"), min.cutoff = "q10", order = T, pt.size=0.1, ncol = 2)
ggsave(paste0(mWAT_ASPC_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 6, height = 6)
saveRDS(mWAT_ASPC_merged, file = paste0(mWAT_ASPC_prefix,"Seurat.rds"))

Selected_enrichment_Plot(GO_enrichment_file = "ASPC/mWAT_ASPC_Selectedmarkers_2024-09-10_GO-BP_list_2024-09-10.csv",
                         OutPrefix = mWAT_ASPC_prefix)


###################3.1.1 ASPC Merge with Mandurp ASPC #########################
ASPCMerged_prefix <- "ASPC_Merge/Mandrup2021_ASPC_"
if(file.exists(paste0(ASPCMerged_prefix,"Seurat.rds"))){
  ASPCMerged_merged <- readRDS(paste0(ASPCMerged_prefix,"Seurat.rds"))
} else {
  MandrupASPC <- readRDS("raw/Mandrup_ASPC_Seurat.rds")
  #MandrupASPC <- subset(MandrupASPC, subset = Dataset == "LFD_R1")
  MandrupASPC$orig.ident <- "eWAT"
  #MandrupASPC$Subtype <- factor(MandrupASPC$Subtype, levels =c("FAP1","FAP2","FAP3","FAP4"))
  ObjList = c(mWAT_ASPC_merged, MandrupASPC)
  features <- SelectIntegrationFeatures(object.list = ObjList, nfeatures = 2000)
  anchors <-  FindIntegrationAnchors(object.list = ObjList, anchor.features = features)
  combined <- IntegrateData(anchorset = anchors)
  ASPCMerged_merged <- ScaleData(combined) %>% RunPCA(npcs = 30) %>%  RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>%   FindClusters(resolution = seq(from=0, by=0.05, length=5))

  ASPCMerged_merged$orig.ident <- gsub("^.{8}$", "mWAT", ASPCMerged_merged$orig.ident)
  ASPCMerged_merged$orig.ident <- gsub("^.{10}$", "mWAT", ASPCMerged_merged$orig.ident)
  
  ASPCMerged_merged$orig.ident <- factor(ASPCMerged_merged$orig.ident)
  saveRDS(ASPCMerged_merged, file = paste0(ASPCMerged_prefix,"Seurat.rds"))
}
ASPCMerged_merged$Myannotation_ASPC <- factor(ASPCMerged_merged$Myannotation_ASPC,level = mWAT_ASPCRanks)
PlotGenes <- c("Pdgfra","Dpp4","Pi16","Anxa3","Pdgfrb","Pparg","Cd36")
Idents(ASPCMerged_merged) <- ASPCMerged_merged$integrated_snn_res.0.2
(p00 <- DimPlot(ASPCMerged_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = hughie_color, label.color = "white"))
(p0 <- DimPlot(ASPCMerged_merged, reduction = "umap", group.by = "orig.ident", ncol = 1, label = F, repel = T, label.box = T, order = T, cols = c("#DC0000FF","#1F78B4")))
(p1 <- DimPlot(ASPCMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = c("#C67BFE","#FAA7A1","#7DAD02","#00BDC3")))
(pa <- VlnPlot(ASPCMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

Idents(ASPCMerged_merged) <- ASPCMerged_merged$Subtype
Idents(ASPCMerged_merged) <- factor(Idents(ASPCMerged_merged) , levels =c("FAP1","FAP2","FAP3","FAP4"))
(p2 <- DimPlot(ASPCMerged_merged, reduction = "umap",split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = c("#C67BFE","#FAA7A1","#7DAD02","#00BDC3")))
(pb <- VlnPlot(ASPCMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

Idents(ASPCMerged_merged) <- ASPCMerged_merged$Myannotation_ASPC
(p3 <- DimPlot(ASPCMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = mWAT_ASPC_colors))
(pc <- VlnPlot(ASPCMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

(p00+p0)/(p2)/(p3)
ggsave(paste0(ASPCMerged_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 10, height = 15)

#######################3.1.2 ASPC Merge with eLife ASPC ##############################
ASPCMerged_prefix <- "ASPC_Merge/Houyu2024_ASPC_"
if(file.exists(paste0(ASPCMerged_prefix,"Seurat.rds"))){
  ASPCMerged_merged <- readRDS(paste0(ASPCMerged_prefix,"Seurat.rds"))
} else {
  HouyuAd <- readRDS(paste0("raw/lPRAT_ASPC_Seurat.rds"))
  HouyuAd$orig.ident <- "lPRAT"
  #HouyuAd <- subset(HouyuAd, features = HouyuAd@assays[["integrated"]]@var.features)
  
  ObjList1 = c(mWAT_ASPC_merged, HouyuAd)
  features1 <- SelectIntegrationFeatures(object.list = ObjList1, nfeatures = 2000)
  anchors1 <-  FindIntegrationAnchors(object.list = ObjList1, anchor.features = features1)
  combined <- IntegrateData(anchorset = anchors1)
  #houyu pepline
  AdMerged_merged <- ScaleData(combined) %>% RunPCA(npcs = 20) %>% RunUMAP(reduction = "pca", dims = 1:20) %>% 
    FindNeighbors(reduction = "pca", dims = 1:20) %>% FindClusters(resolution = seq(from=0, by=0.05, length=5))
  
  ASPCMerged_merged <- ScaleData(combined) %>% RunPCA(npcs = 30) %>%  RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>%   FindClusters(resolution = seq(from=0, by=0.05, length=5))
  
  ASPCMerged_merged$orig.ident <- gsub("^.{8}$", "mWAT", ASPCMerged_merged$orig.ident)
  ASPCMerged_merged$orig.ident <- gsub("^.{10}$", "mWAT", ASPCMerged_merged$orig.ident)
  #ASPCMerged_merged$orig.ident <- gsub("^.{9}$", "lPRAT", ASPCMerged_merged$orig.ident)
  
  ASPCMerged_merged$orig.ident <- as.factor(ASPCMerged_merged$orig.ident)
  saveRDS(ASPCMerged_merged, file = paste0(ASPCMerged_prefix,"Seurat.rds"))
}

PlotGenes <- c("Ebf2","Lep","Adipoq", #Lipolytic
               "Acaca","Elvl6",#Lipogenesis
               "Cmss1","Lars2","Camk1d","Cdk8"
)

ASPCMerged_merged$Myannotation_ASPC <- factor(ASPCMerged_merged$Myannotation_ASPC,level = mWAT_ASPCRanks)
PlotGenes <- c("Pdgfra","Dpp4","Pi16","Anxa3","Pdgfrb","Pparg","Cd36")
Idents(ASPCMerged_merged) <- ASPCMerged_merged$integrated_snn_res.0.2
(p00 <- DimPlot(ASPCMerged_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = hughie_color, label.color = "white"))
(p0 <- DimPlot(ASPCMerged_merged, reduction = "umap", group.by = "orig.ident", ncol = 1, label = F, repel = T, label.box = T, order = T, cols = c("#DC0000FF","#1F78B4")))
(p1 <- DimPlot(ASPCMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = hughie_color))
(pa <- VlnPlot(ASPCMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

Idents(ASPCMerged_merged) <- ASPCMerged_merged$Subtype
Idents(ASPCMerged_merged) <- factor(Idents(ASPCMerged_merged) , levels =c("lPRAT-aspc1","lPRAT-aspc2"))
(p2 <- DimPlot(ASPCMerged_merged, reduction = "umap",split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = c("#2473AB","#D8823D")))
(pb <- VlnPlot(ASPCMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

Idents(ASPCMerged_merged) <- ASPCMerged_merged$Myannotation_ASPC
(p3 <- DimPlot(ASPCMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = mWAT_ASPC_colors))
(pc <- VlnPlot(ASPCMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

(p00+p0)/(p2)/(p3)
ggsave(paste0(ASPCMerged_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 10, height = 15)


######################3.2 mWAT ASPC monocle ######################
ASPCMonocle_prefix <- "ASPC_monocle/mWAT_ASPC_"

Idents(mWAT_ASPC_merged) <- "Myannotation_ASPC"
data <- as(as.matrix(mWAT_ASPC_merged@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = mWAT_ASPC_merged@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

mycds_ASPC <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

mycds_ASPC <- estimateSizeFactors(mycds_ASPC)
mycds_ASPC <- estimateDispersions(mycds_ASPC)                        
mycds_ASPC <- detectGenes(mycds_ASPC, min_expr = 0.1)
mycds_ASPC_expressed_genes <- row.names(mycds_ASPC)
diff_test_res_ASPC <- differentialGeneTest(mycds_ASPC[mycds_ASPC_expressed_genes,],
                                      fullModelFormulaStr = "~Myannotation_ASPC")                      
diff_test_res_ASPC[1:4,1:4]
write.csv(diff_test_res_ASPC, file=paste0(ASPCMonocle_prefix,"ALL_diff_test_result.csv",sep=""))

ordering_genes_ASPC <- row.names(subset(diff_test_res_ASPC, qval < 0.01))
ordering_genes_ASPC[1:5]
mycds_ASPC <- setOrderingFilter(mycds_ASPC, ordering_genes_ASPC)
plot_ordering_genes(mycds_ASPC)
ggsave(paste0(ASPCMonocle_prefix,"ALL_Ordering.genes_", Sys.Date(), ".pdf"), height = 7, width = 8)

mycds_ASPC <- reduceDimension(mycds_ASPC, max_components = 2, method = 'DDRTree',verbose = T)
mycds_ASPC <- orderCells(mycds_ASPC)

saveRDS(mycds_ASPC, file = paste0(ASPCMonocle_prefix, "ALL_mycds.rds"))


data <- GetAssayData(mWAT_ASPC_merged, assay = 'RNA', slot = 'counts')
cell_metadata_ASPC <- mWAT_ASPC_merged@meta.data
gene_annotation_ASPC <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation_ASPC) <- rownames(data)
monocle3_cds_ASPC <- new_cell_data_set(data,
                                  cell_metadata = cell_metadata_ASPC,
                                  gene_metadata = gene_annotation_ASPC)
#preprocess_cds函数相当于seurat中NormalizeData ScaleData RunPCA
monocle3_cds_ASPC <- preprocess_cds(monocle3_cds_ASPC, num_dim = 100)
#umap降维
monocle3_cds_ASPC <- reduce_dimension(monocle3_cds_ASPC, preprocess_method = "PCA")
p1 <- plot_cells(monocle3_cds_ASPC, reduction_method="UMAP", color_cells_by="Myannotation_ASPC"
                 #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
)+ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
monocle3_cds.embed_ASPC <- monocle3_cds_ASPC@int_colData$reducedDims$UMAP
int.embed_ASPC <- Embeddings(mWAT_ASPC_merged, reduction = "umap")
int.embed_ASPC <- int.embed_ASPC[rownames(monocle3_cds.embed_ASPC),]
monocle3_cds_ASPC@int_colData$reducedDims$UMAP <- int.embed_ASPC
p2 <- plot_cells(monocle3_cds_ASPC, reduction_method="UMAP", color_cells_by="Myannotation_ASPC"
                 #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
)+ggtitle('int.umap')
p1 + p2
ggsave(paste0(ASPCMonocle_prefix,"ALL_Monocle3_umap_embed_", Sys.Date(), ".pdf"), height = 7, width = 16)

## Monocle3聚类分区
monocle3_cds_ASPC <- cluster_cells(monocle3_cds_ASPC, resolution = 0.01)
p3 <- plot_cells(monocle3_cds_ASPC, show_trajectory_graph = FALSE)+ggtitle("Label by clusterID")
p4 <- plot_cells(monocle3_cds_ASPC, color_cells_by = "partition", show_trajectory_graph = FALSE)+ggtitle("Label by partitionID")
p3 + p4
ggsave(paste0(ASPCMonocle_prefix,"ALL_Monocle3_umap_cluster_", Sys.Date(), ".pdf"), height = 7, width = 16)


## 识别轨迹
monocle3_cds_ASPC <- learn_graph(monocle3_cds_ASPC)
p5 <- plot_cells(monocle3_cds_ASPC, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                 label_branch_points = FALSE)
p5
ggsave(paste0(ASPCMonocle_prefix,"ALL_Monocle3_umap_trajectory_", Sys.Date(), ".pdf"), width = 5, height = 5)


##细胞按拟时排序
# cds <- order_cells(cds) 存在bug，使用辅助线选择root细胞
#p   geom_vline(xintercept = seq(-7,-6,0.25))   geom_hline(yintercept = seq(0,1,0.25))
#embed <- data.frame(Embeddings(all, reduction = "umap"))
#embed <- subset(embed, umap_1 > -6.75 & umap_1 < -6.5 & umap_2 > 0.24 & umap_2 < 0.25)
#root.cell <- rownames(embed)
#monocle3_cds <- order_cells(monocle3_cds, root_cells = root.cell)
monocle3_cds_ASPC <- order_cells(monocle3_cds_ASPC)
p6 <- plot_cells(monocle3_cds_ASPC, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                 label_leaves = F,  label_branch_points = FALSE,graph_label_size = 3,
                 cell_size = 1, alpha = 0.5, trajectory_graph_color = "black",trajectory_graph_segment_size = 1)
p6
ggsave(paste0(ASPCMonocle_prefix,"ALL_Monocle3_umap_trajectory_line_", Sys.Date(), ".pdf"), width = 6, height = 5)

#plot heatmap 
sig_gene_names_ASPC <- row.names(subset(diff_test_res_ASPC, qval < 0.000001))
nrow(subset(diff_test_res_ASPC, qval < 0.000001))

bins <- monocle3_scale_to_100_ASPC(mycds_ASPC[sig_gene_names_ASPC,])

pdf(paste0(ASPCMonocle_prefix,"ALL_monocle_heatmap_",Sys.Date(),".pdf"), width = 5, height = 15)

plot_pseudotime_heatmap(mycds_ASPC[sig_gene_names_ASPC,],
                        num_clusters = 3,
                        cores = 6,
                        show_rownames = F, return_heatmap = T,
                        add_annotation_col = bins,use_gene_short_name = T
)

dev.off()

df <- data.frame((mycds_ASPC[sig_gene_names_ASPC,]@phenoData@data))
df <- df[,c("Pseudotime", "Myannotation_ASPC","timepoint")]

#"#A6CEE3", "#A2D48E",
#"#66A5CC","#47A93A",
#"#267CB6","#20854EFF"
color.celtype <- data.frame(celltype = c("mWAT-ASPC1","mWAT-ASPC2","mWAT-ASPC3"), 
                            color = c("#BFA6C9","#F5E0BA","#AED0DF"))  
mWAT_SAMcolors <- data.frame(timepoint = c("HFD0w","HFD4w","HFD8w"),
                             color = c("#A6CEE3","#66A5CC","#267CB6"))

df <- merge(df,color.celtype,by.x = "Myannotation_ASPC", by.y = "celltype",drop = F)
df <- merge(df,mWAT_SAMcolors,by.x = "timepoint", by.y = "timepoint",drop = F)
df <- df[order(df$Pseudotime, decreasing = F),]

write.csv(df, file=paste0(ASPCMonocle_prefix,"ALL_colnames.csv",sep=""))

matrix <- as.matrix(df[,3])
pdf(paste0(ASPCMonocle_prefix,"ALL_monocle_heatmap_celltype_",Sys.Date(),".pdf"))
image(matrix, col = df$color.x)
dev.off()

pdf(paste0(ASPCMonocle_prefix,"ALL_monocle_heatmap_timepoint_",Sys.Date(),".pdf"))

image(matrix, col = df$color.y)
dev.off()


#plot heatmap CD group
sig_gene_names_ASPC <- row.names(subset(diff_test_res_ASPC, qval < 0.000001))
nrow(subset(diff_test_res_ASPC, qval < 0.000001))

bins <- monocle3_scale_to_100_ASPC(mycds_ASPC[sig_gene_names_ASPC,grep("CD",mycds_ASPC$group)])

pdf(paste0(ASPCMonocle_prefix,"ALL_CD_monocle_heatmap_",Sys.Date(),".pdf"), width = 5, height = 15)

plot_pseudotime_heatmap(mycds_ASPC[sig_gene_names_ASPC,grep("CD",mycds_ASPC$group)],
                        num_clusters = 3,
                        cores = 6,
                        show_rownames = T, return_heatmap = T,
                        add_annotation_col = bins,use_gene_short_name = T
)

dev.off()


df <- data.frame((mycds_ASPC[sig_gene_names_ASPC,grep("CD",mycds_ASPC$group)]@phenoData@data))
df <- df[,c("Pseudotime", "Myannotation_ASPC","timepoint")]

#"#A6CEE3", "#A2D48E",
#"#66A5CC","#47A93A",
#"#267CB6","#20854EFF"
color.celtype <- data.frame(celltype = c("mWAT-ASPC1","mWAT-ASPC2","mWAT-ASPC3"), 
                            color = c("#BFA6C9","#F5E0BA","#AED0DF"))  
mWAT_SAMcolors <- data.frame(timepoint = c("HFD0w","HFD4w","HFD8w"),
                             color = c("#A6CEE3","#66A5CC","#267CB6"))

df <- merge(df,color.celtype,by.x = "Myannotation_ASPC", by.y = "celltype",drop = F)
df <- merge(df,mWAT_SAMcolors,by.x = "timepoint", by.y = "timepoint",drop = F)
df <- df[order(df$Pseudotime, decreasing = F),]

write.csv(df, file=paste0(ASPCMonocle_prefix,"ALL_CD_group_colnames.csv",sep=""))

matrix <- as.matrix(df[,3])
pdf(paste0(ASPCMonocle_prefix,"ALL_CD_group_monocle_heatmap_celltype_",Sys.Date(),".pdf"))
image(matrix, col = df$color.x)
dev.off()

pdf(paste0(ASPCMonocle_prefix,"ALL_CD_group_monocle_heatmap_timepoint_",Sys.Date(),".pdf"))

image(matrix, col = df$color.y)
dev.off()

#plot heatmap HfiD group
sig_gene_names_ASPC <- row.names(subset(diff_test_res_ASPC, qval < 0.000001))
nrow(subset(diff_test_res_ASPC, qval < 0.000001))

bins <- monocle3_scale_to_100_ASPC(mycds_ASPC[sig_gene_names_ASPC,grep("HfiD",mycds_ASPC$group)])

pdf(paste0(ASPCMonocle_prefix,"ALL_HfiD_monocle_heatmap_",Sys.Date(),".pdf"), width = 5, height = 15)

plot_pseudotime_heatmap(mycds_ASPC[sig_gene_names_ASPC,grep("HfiD",mycds_ASPC$group)],
                        num_clusters = 3,
                        cores = 6,
                        show_rownames = T, return_heatmap = T,
                        add_annotation_col = bins,use_gene_short_name = T
)

dev.off()
df <- data.frame((mycds_ASPC[sig_gene_names_ASPC,grep("HfiD",mycds_ASPC$group)]@phenoData@data))
df <- df[,c("Pseudotime", "Myannotation_ASPC","timepoint")]

#"#A6CEE3", "#A2D48E",
#"#66A5CC","#47A93A",
#"#267CB6","#20854EFF"
color.celtype <- data.frame(celltype = c("mWAT-ASPC1","mWAT-ASPC2","mWAT-ASPC3"), 
                            color = c("#BFA6C9","#F5E0BA","#AED0DF"))  
mWAT_SAMcolors <- data.frame(timepoint = c("HFD0w","HFD4w","HFD8w"),
                             color = c("#A2D48E","#47A93A","#20854EFF"))

df <- merge(df,color.celtype,by.x = "Myannotation_ASPC", by.y = "celltype",drop = F)
df <- merge(df,mWAT_SAMcolors,by.x = "timepoint", by.y = "timepoint",drop = F)
df <- df[order(df$Pseudotime, decreasing = F),]

write.csv(df, file=paste0(ASPCMonocle_prefix,"ALL_HfiD_group_colnames.csv",sep=""))

matrix <- as.matrix(df[,3])
pdf(paste0(ASPCMonocle_prefix,"ALL_HfiD_group_monocle_heatmap_celltype_",Sys.Date(),".pdf"))
image(matrix, col = df$color.x)
dev.off()

pdf(paste0(ASPCMonocle_prefix,"ALL_HfiD_group_monocle_heatmap_timepoint_",Sys.Date(),".pdf"))

image(matrix, col = df$color.y)
dev.off()

#####################3.2.1 CD group mWAT ASPCs monocle###################

Idents(mWAT_ASPC_merged) <- "group"
mWAT_ASPC_merged_cd <- subset(mWAT_ASPC_merged, idents = c("CD"))
Idents(mWAT_ASPC_merged_cd) <- "Myannotation_ASPC"

#Idents(mWAT_Ad_merged) <- "Myannotation_Ad"
data_cd_ASPC <- as(as.matrix(mWAT_ASPC_merged_cd@assays$RNA@counts), 'sparseMatrix')
pd_cd <- new('AnnotatedDataFrame', data = mWAT_ASPC_merged_cd@meta.data)
fData_cd <- data.frame(gene_short_name = row.names(data_cd_ASPC), row.names = row.names(data_cd_ASPC))
fd_cd <- new('AnnotatedDataFrame', data = fData_cd)

mycds_cd_ASPC <- newCellDataSet(data_cd_ASPC,
                           phenoData = pd_cd,
                           featureData = fd_cd,
                           expressionFamily = negbinomial.size())

mycds_cd_ASPC <- estimateSizeFactors(mycds_cd_ASPC)
mycds_cd_ASPC <- estimateDispersions(mycds_cd_ASPC)                        
mycds_cd_ASPC <- detectGenes(mycds_cd_ASPC, min_expr = 0.1)
mycds_expressed_genes_cd_ASPC <- row.names(mycds_cd_ASPC)
diff_test_res_cd_ASPC <- differentialGeneTest(mycds_cd_ASPC[mycds_expressed_genes_cd_ASPC,],
                                         fullModelFormulaStr = "~Myannotation_ASPC")                      
diff_test_res_cd_ASPC[1:4,1:4]
write.csv(diff_test_res_cd_ASPC, file=paste0(ASPCMonocle_prefix,"CD_group_diff_test_result.csv",sep=""))
#diff_test_res_cd <- read.csv(file=paste0(AdMonocle_prefix,"CD_group_diff_test_result.csv",sep=""))

ordering_genes_cd <- row.names(subset(diff_test_res_cd_ASPC, qval < 0.01))
ordering_genes_cd[1:5]
mycds_cd_ASPC <- setOrderingFilter(mycds_cd_ASPC, ordering_genes_cd)
plot_ordering_genes(mycds_cd_ASPC)
ggsave(paste0(ASPCMonocle_prefix,"CD_group_Ordering.genes_", Sys.Date(), ".pdf"), height = 7, width = 8)

mycds_cd_ASPC <- reduceDimension(mycds_cd_ASPC, max_components = 2, method = 'DDRTree',verbose = T)
mycds_cd_ASPC <- orderCells(mycds_cd_ASPC)

saveRDS(mycds_cd_ASPC, file = paste0(ASPCMonocle_prefix, "CD_group_mycds.rds"))
#mycds_cd <- readRDS("Adipocyte_monocle/mWAT_Adipocytes_CD_group_mycds.rds")

data <- GetAssayData(mWAT_ASPC_merged_cd, assay = 'RNA', slot = 'counts')
cell_metadata_ASPC <- mWAT_ASPC_merged_cd@meta.data
gene_annotation_ASPC <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation_ASPC) <- rownames(data)
monocle3_cds_ASPC <- new_cell_data_set(data,
                                       cell_metadata = cell_metadata_ASPC,
                                       gene_metadata = gene_annotation_ASPC)
#preprocess_cds函数相当于seurat中NormalizeData ScaleData RunPCA
monocle3_cds_ASPC <- preprocess_cds(monocle3_cds_ASPC, num_dim = 100)
#umap降维
monocle3_cds_ASPC <- reduce_dimension(monocle3_cds_ASPC, preprocess_method = "PCA")
p1 <- plot_cells(monocle3_cds_ASPC, reduction_method="UMAP", color_cells_by="Myannotation_ASPC"
                 #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
)+ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
monocle3_cds.embed_ASPC <- monocle3_cds_ASPC@int_colData$reducedDims$UMAP
int.embed_ASPC <- Embeddings(mWAT_ASPC_merged_cd, reduction = "umap")
int.embed_ASPC <- int.embed_ASPC[rownames(monocle3_cds.embed_ASPC),]
monocle3_cds_ASPC@int_colData$reducedDims$UMAP <- int.embed_ASPC
p2 <- plot_cells(monocle3_cds_ASPC, reduction_method="UMAP", color_cells_by="Myannotation_ASPC"
                 #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
)+ggtitle('int.umap')
p1 + p2
ggsave(paste0(ASPCMonocle_prefix,"CD_group_Monocle3_umap_embed_", Sys.Date(), ".pdf"), height = 7, width = 16)

## Monocle3聚类分区
monocle3_cds_ASPC <- cluster_cells(monocle3_cds_ASPC, resolution = 0.01)
p3 <- plot_cells(monocle3_cds_ASPC, show_trajectory_graph = FALSE)+ggtitle("Label by clusterID")
p4 <- plot_cells(monocle3_cds_ASPC, color_cells_by = "partition", show_trajectory_graph = FALSE)+ggtitle("Label by partitionID")
p3 + p4
ggsave(paste0(ASPCMonocle_prefix,"CD_group_Monocle3_umap_cluster_", Sys.Date(), ".pdf"), height = 7, width = 16)


## 识别轨迹
monocle3_cds_ASPC <- learn_graph(monocle3_cds_ASPC)
p5 <- plot_cells(monocle3_cds_ASPC, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                 label_branch_points = FALSE)
p5
ggsave(paste0(ASPCMonocle_prefix,"CD_group_Monocle3_umap_trajectory_", Sys.Date(), ".pdf"), width = 5, height = 5)


##细胞按拟时排序
# cds <- order_cells(cds) 存在bug，使用辅助线选择root细胞
#p   geom_vline(xintercept = seq(-7,-6,0.25))   geom_hline(yintercept = seq(0,1,0.25))
#embed <- data.frame(Embeddings(all, reduction = "umap"))
#embed <- subset(embed, umap_1 > -6.75 & umap_1 < -6.5 & umap_2 > 0.24 & umap_2 < 0.25)
#root.cell <- rownames(embed)
#monocle3_cds <- order_cells(monocle3_cds, root_cells = root.cell)
monocle3_cds_ASPC <- order_cells(monocle3_cds_ASPC)
p6 <- plot_cells(monocle3_cds_ASPC, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                 label_leaves = F,  label_branch_points = FALSE,graph_label_size = 3,
                 cell_size = 1, alpha = 0.5, trajectory_graph_color = "black",trajectory_graph_segment_size = 1)
p6
ggsave(paste0(ASPCMonocle_prefix,"CD_group_Monocle3_umap_trajectory_line_", Sys.Date(), ".pdf"), width = 5, height = 5)


#plot heatmap of CD group
diff_test_res_cd_ASPC <- diff_test_res_cd_ASPC[order(diff_test_res_cd_ASPC$qval, decreasing = F),]
sig_gene_names_cd_ASPC <- row.names(diff_test_res_cd_ASPC[1:600,])
#nrow(subset(diff_test_res_cd, qval < 0.0001))

bins <- monocle3_scale_to_100_ASPC(mycds_cd_ASPC[sig_gene_names_cd_ASPC,])


pdf(paste0(ASPCMonocle_prefix,"CD_group_monocle_heatmap_",Sys.Date(),".pdf"), width = 5, height = 15)

plot_pseudotime_heatmap(mycds_cd_ASPC[sig_gene_names_cd_ASPC,],
                        num_clusters = 3,
                        cores = 6,
                        show_rownames = T, return_heatmap = T,
                        add_annotation_col = bins,use_gene_short_name = T
)

dev.off()

df <- data.frame((mycds_cd_ASPC[sig_gene_names_cd_ASPC,]@phenoData@data))
df <- df[,c("Pseudotime", "Myannotation_ASPC","timepoint")]

#"#A6CEE3", "#A2D48E",
#"#66A5CC","#47A93A",
#"#267CB6","#20854EFF"
color.celtype <- data.frame(celltype = c("mWAT-ASPC1","mWAT-ASPC2","mWAT-ASPC3"), 
                            color = c("#BFA6C9","#F5E0BA","#AED0DF"))  
mWAT_SAMcolors <- data.frame(timepoint = c("HFD0w","HFD4w","HFD8w"),
                             color = c("#A6CEE3","#66A5CC","#267CB6"))

df <- merge(df,color.celtype,by.x = "Myannotation_ASPC", by.y = "celltype",drop = F)
df <- merge(df,mWAT_SAMcolors,by.x = "timepoint", by.y = "timepoint",drop = F)
df <- df[order(df$Pseudotime, decreasing = F),]

write.csv(df, file=paste0(ASPCMonocle_prefix,"CD_group_colnames.csv",sep=""))

matrix <- as.matrix(df[,3])
pdf(paste0(ASPCMonocle_prefix,"CD_group_monocle_heatmap_celltype_",Sys.Date(),".pdf"))
image(matrix, col = df$color.x)
dev.off()

pdf(paste0(ASPCMonocle_prefix,"CD_group_monocle_heatmap_timepoint_",Sys.Date(),".pdf"))

image(matrix, col = df$color.y)
dev.off()


#####################3.2.2 HfiD group mWAT ASPCs monocle###################

Idents(mWAT_ASPC_merged) <- "group"
mWAT_ASPC_merged_HfiD <- subset(mWAT_ASPC_merged, idents = c("HfiD"))
Idents(mWAT_ASPC_merged_HfiD) <- "Myannotation_ASPC"

#Idents(mWAT_Ad_merged) <- "Myannotation_Ad"
data_HfiD <- as(as.matrix(mWAT_ASPC_merged_HfiD@assays$RNA@counts), 'sparseMatrix')
pd_HfiD <- new('AnnotatedDataFrame', data = mWAT_ASPC_merged_HfiD@meta.data)
fData_HfiD <- data.frame(gene_short_name = row.names(data_HfiD), row.names = row.names(data_HfiD))
fd_HfiD <- new('AnnotatedDataFrame', data = fData_HfiD)

mycds_HfiD_ASPC <- newCellDataSet(data_HfiD,
                             phenoData = pd_HfiD,
                             featureData = fd_HfiD,
                             expressionFamily = negbinomial.size())

mycds_HfiD_ASPC <- estimateSizeFactors(mycds_HfiD_ASPC)
mycds_HfiD_ASPC <- estimateDispersions(mycds_HfiD_ASPC)                        
mycds_HfiD_ASPC <- detectGenes(mycds_HfiD_ASPC, min_expr = 0.1)
mycds_expressed_genes_HfiD_ASPC <- row.names(mycds_HfiD_ASPC)
diff_test_res_HfiD_ASPC <- differentialGeneTest(mycds_HfiD_ASPC[mycds_expressed_genes_HfiD_ASPC,],
                                           fullModelFormulaStr = "~Myannotation_ASPC")                      
diff_test_res_HfiD_ASPC[1:4,1:4]
write.csv(diff_test_res_HfiD_ASPC, file=paste0(ASPCMonocle_prefix,"HfiD_group_diff_test_result.csv",sep=""))
#diff_test_res_HfiD <- read.csv(file=paste0(AdMonocle_prefix,"HfiD_group_diff_test_result.csv",sep=""))

ordering_genes_HfiD <- row.names(subset(diff_test_res_HfiD_ASPC, qval < 0.01))
ordering_genes_HfiD[1:5]
mycds_HfiD_ASPC <- setOrderingFilter(mycds_HfiD_ASPC, ordering_genes_HfiD)
plot_ordering_genes(mycds_HfiD_ASPC)
ggsave(paste0(ASPCMonocle_prefix,"HfiD_group_Ordering.genes_", Sys.Date(), ".pdf"), height = 7, width = 8)

mycds_HfiD_ASPC <- reduceDimension(mycds_HfiD_ASPC, max_components = 2, method = 'DDRTree',verbose = T)
mycds_HfiD_ASPC <- orderCells(mycds_HfiD_ASPC)

saveRDS(mycds_HfiD_ASPC, file = paste0(ASPCMonocle_prefix, "HfiD_group_mycds.rds"))
#mycds_HfiD <- readRDS(file = paste0(AdMonocle_prefix, "HfiD_group_mycds.rds"))

data <- GetAssayData(mWAT_ASPC_merged_HfiD, assay = 'RNA', slot = 'counts')
cell_metadata_ASPC <- mWAT_ASPC_merged_HfiD@meta.data
gene_annotation_ASPC <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation_ASPC) <- rownames(data)
monocle3_HfiD_ASPC <- new_cell_data_set(data,
                                       cell_metadata = cell_metadata_ASPC,
                                       gene_metadata = gene_annotation_ASPC)
#preprocess_cds函数相当于seurat中NormalizeData ScaleData RunPCA
monocle3_HfiD_ASPC <- preprocess_cds(monocle3_HfiD_ASPC, num_dim = 100)
#umap降维
monocle3_HfiD_ASPC <- reduce_dimension(monocle3_HfiD_ASPC, preprocess_method = "PCA")
p1 <- plot_cells(monocle3_HfiD_ASPC, reduction_method="UMAP", color_cells_by="Myannotation_ASPC"
                 #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
)+ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
monocle3_cds.embed_ASPC <- monocle3_HfiD_ASPC@int_colData$reducedDims$UMAP
int.embed_ASPC <- Embeddings(mWAT_ASPC_merged_HfiD, reduction = "umap")
int.embed_ASPC <- int.embed_ASPC[rownames(monocle3_cds.embed_ASPC),]
monocle3_HfiD_ASPC@int_colData$reducedDims$UMAP <- int.embed_ASPC
p2 <- plot_cells(monocle3_HfiD_ASPC, reduction_method="UMAP", color_cells_by="Myannotation_ASPC"
                 #,trajectory_graph_color = c("#EFE2AA","#83B4EF", "#DBC9B3")
)+ggtitle('int.umap')
p1 + p2
ggsave(paste0(ASPCMonocle_prefix,"HfiD_group_Monocle3_umap_embed_", Sys.Date(), ".pdf"), height = 7, width = 16)

## Monocle3聚类分区
monocle3_HfiD_ASPC <- cluster_cells(monocle3_HfiD_ASPC, resolution = 0.01)
p3 <- plot_cells(monocle3_HfiD_ASPC, show_trajectory_graph = FALSE)+ggtitle("Label by clusterID")
p4 <- plot_cells(monocle3_HfiD_ASPC, color_cells_by = "partition", show_trajectory_graph = FALSE)+ggtitle("Label by partitionID")
p3 + p4
ggsave(paste0(ASPCMonocle_prefix,"HfiD_group_Monocle3_umap_cluster_", Sys.Date(), ".pdf"), height = 7, width = 16)


## 识别轨迹
monocle3_HfiD_ASPC <- learn_graph(monocle3_HfiD_ASPC)
p5 <- plot_cells(monocle3_HfiD_ASPC, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                 label_branch_points = FALSE)
p5
ggsave(paste0(ASPCMonocle_prefix,"HfiD_group_Monocle3_umap_trajectory_", Sys.Date(), ".pdf"), width = 5, height = 5)


##细胞按拟时排序
# cds <- order_cells(cds) 存在bug，使用辅助线选择root细胞
#p   geom_vline(xintercept = seq(-7,-6,0.25))   geom_hline(yintercept = seq(0,1,0.25))
#embed <- data.frame(Embeddings(all, reduction = "umap"))
#embed <- subset(embed, umap_1 > -6.75 & umap_1 < -6.5 & umap_2 > 0.24 & umap_2 < 0.25)
#root.cell <- rownames(embed)
#monocle3_cds <- order_cells(monocle3_cds, root_cells = root.cell)
monocle3_HfiD_ASPC <- order_cells(monocle3_HfiD_ASPC)
p6 <- plot_cells(monocle3_HfiD_ASPC, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                 label_leaves = F,  label_branch_points = FALSE,graph_label_size = 3,
                 cell_size = 1, alpha = 0.5, trajectory_graph_color = "black",trajectory_graph_segment_size = 1)
p6
ggsave(paste0(ASPCMonocle_prefix,"HfiD_group_Monocle3_umap_trajectory_line_", Sys.Date(), ".pdf"), width = 6, height = 5)


#plot heatmap of HfiD group
diff_test_res_HfiD_ASPC <- diff_test_res_HfiD_ASPC[order(diff_test_res_HfiD_ASPC$qval, decreasing = F),]
sig_gene_names_HfiD_ASPC <- row.names(diff_test_res_HfiD_ASPC[1:585,])
#nrow(subset(diff_test_res_HfiD_ASPC, qval < 0.0001))

bins <- monocle3_scale_to_100_ASPC(mycds_HfiD_ASPC[sig_gene_names_HfiD_ASPC,])

pdf(paste0(ASPCMonocle_prefix,"HfiD_group_monocel_heatmap_",Sys.Date(),".pdf"), width = 5, height = 15)

plot_pseudotime_heatmap(mycds_HfiD_ASPC[sig_gene_names_HfiD_ASPC,],
                        num_clusters = 3,
                        cores = 6,
                        show_rownames = T, 
                        return_heatmap = T,
                        add_annotation_col = bins,
                        use_gene_short_name = T
)

dev.off()

df <- data.frame((mycds_HfiD_ASPC[sig_gene_names_HfiD_ASPC,]@phenoData@data))
df <- df[,c("Pseudotime", "Myannotation_ASPC","timepoint")]

#"#A6CEE3", "#A2D48E",
#"#66A5CC","#47A93A",
#"#267CB6","#20854EFF"
color.celtype <- data.frame(celltype = c("mWAT-ASPC1","mWAT-ASPC2","mWAT-ASPC3"), 
                            color = c("#BFA6C9","#F5E0BA","#AED0DF"))  
mWAT_SAMcolors <- data.frame(timepoint = c("HFD0w","HFD4w","HFD8w"),
                             color = c("#A2D48E","#47A93A","#20854EFF"))

df <- merge(df,color.celtype,by.x = "Myannotation_ASPC", by.y = "celltype",drop = F)
df <- merge(df,mWAT_SAMcolors,by.x = "timepoint", by.y = "timepoint",drop = F)
df <- df[order(df$Pseudotime, decreasing = F),]

write.csv(df, file=paste0(ASPCMonocle_prefix,"HfiD_group_colnames.csv",sep=""))

matrix <- as.matrix(df[,3])
pdf(paste0(ASPCMonocle_prefix,"HfiD_group_monocel_heatmap_celltype_",Sys.Date(),".pdf"))
image(matrix, col = df$color.x)
dev.off()

pdf(paste0(ASPCMonocle_prefix,"HfiD_group_monocel_heatmap_timepoint_",Sys.Date(),".pdf"))

image(matrix, col = df$color.y)
dev.off()

#################4 mWAT-Ad2 ASPC2 ####################

mWAT_Ad_ASPC_merged <- merge(mWAT_Ad_merged,y=mWAT_ASPC_merged,
                             add.cell.ids = NULL,
                             project = "CD_Ad_ASPC")
mWAT_Ad_ASPC_merged$Myannotation_Ad_ASPC <- mWAT_Ad_ASPC_merged@active.ident
mWAT_Ad_ASPC_merged <- subset(mWAT_Ad_ASPC_merged, Myannotation_Ad_ASPC %in% 
                                c("mWAT-ad2","mWAT-ASPC2"))

mWAT_Ad2_ASPC2_prefix <- "ASPC/mWAT_Ad2_ASPC2_"

PlotGenes <- c("Plin1","Acvr1c","Pdgfra","Nova1",
               "Cdk8","Camk1d","Cmss1","Lars2")
(p1 <- VlnPlot(mWAT_Ad_ASPC_merged, features = PlotGenes, fill.by = "ident", cols = c("#F5E0BA","#f5cfe4"), stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
p1
ggsave(paste0(mWAT_Ad2_ASPC2_prefix, "Similar_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)



######################5 mWAT-Immune cell##############################

mWAT_immune_prefix <- "Immune/mWAT_immune_"
if(file.exists(paste0(mWAT_immune_prefix,"Seurat.rds"))){
  mWAT_Immune_merged <- readRDS(paste0(mWAT_immune_prefix,"Seurat.rds"))
} else {
  mWAT_Immune_merged <- ScaleData(subset(mWAT_merged, idents = c("Immune cell"))) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  clustree(mWAT_Immune_merged)
  ggsave(paste0(mWAT_immune_prefix, "clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
  
  Idents(mWAT_Immune_merged) <- mWAT_Immune_merged$RNA_snn_res.0.15
  DimPlot(mWAT_Immune_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = mWAT_Immune_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = mWAT_immune_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = mWAT_Immune_merged, ColorUse = hughie_color, OutPrefix = mWAT_immune_prefix)
}

PlotGenes <- c("Ptprc","Dock2", #all
               "Cd79b","Pax5",  #B cell
               "Cd3e", "Bcl11b",#T cell
               "Cxcr6",'Il7r',#CXCR6+ T cell
               "Tcf7", "Trbc2",#Treg cell 
               "Adgre1","Mrc1",#Macrophage
               "Csf3r",#Neutrophil
               "Plbd1",'Clec9a',   #Dendritic cell
               "Ccl5"#NK cell
               )
(p1 <- VlnPlot(mWAT_Immune_merged, features = PlotGenes, fill.by = "ident", cols = hughie_color, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
#(p2 <- VlnPlot(mWAT_Immune_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mWAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 
#+ p2
ggsave(paste0(mWAT_immune_prefix, "Before_annotation_Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)


mWAT_Immune_merged <- RenameIdents(mWAT_Immune_merged, 
                                   `0` = "B cell", `5` = "B cell", `10` = "B cell",`13` = "B cell",`14` = "B cell",
                                   `1` = "T cell", `2` = "T cell", `3` = "T cell", `12` = "T cell",
                                   `6` = "Macrophage",
                                   `7` = "Neutrophil",
                                   `8` = "Dendritic cell",`9` = "Dendritic cell",`11` = "Dendritic cell",
                                   `4` = "NK cell"
                                   )

mWAT_Immune_merged$Myannotation_Immune <- mWAT_Immune_merged@active.ident
mWAT_MacRanks <- c("B cell","T cell","Macrophage","Neutrophil","Dendritic cell","NK cell")
mWAT_Mac_colors <- hughie_color
names(mWAT_Mac_colors) <- mWAT_MacRanks
Idents(mWAT_Immune_merged) <- factor(Idents(mWAT_Immune_merged), levels = mWAT_MacRanks)

PlotGenes <- c("Ptprc","Dock2", #all
               "Pax5","Cd79b",  #B cell
               "Bcl11b","Tcf7", #T cell
               "Adgre1","Mrc1",#Macrophage
               "Csf3r",#Neutrophil
               "Plbd1",'Clec9a',   #Dendritic cell
               "Ccl5"#NK cell
)
mWAT_SAMcolors <- c("#A6CEE3", "#66A5CC","#267CB6",
                    "#A2D48E","#47A93A","#20854EFF")

(p1 <- VlnPlot(mWAT_Immune_merged, features = PlotGenes, fill.by = "ident", cols = mWAT_Mac_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(mWAT_Immune_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mWAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(mWAT_immune_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)

(p1 <- DimPlot(mWAT_Immune_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = hughie_color, label.color = "white"))
(p2 <- DimPlot(mWAT_Immune_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = mWAT_SAMcolors))
p1 + p2
ggsave(paste0(mWAT_immune_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(mWAT_Immune_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = F, repel = T, label.box =T, cols = mWAT_Mac_colors, label.color = "white", order = T)
ggsave(paste0(mWAT_immune_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 12, height = 5)

Plot_Cell_compoistion(Seurat_Obj = mWAT_Immune_merged, OutPrefix = mWAT_immune_prefix, ColorUse1 = mWAT_Mac_colors, ColorUse2 = mWAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = mWAT_Immune_merged, ColorUse = mWAT_Mac_colors, OutPrefix = mWAT_immune_prefix)
DEG_enrichment(Seruat_DEG_file = paste0(mWAT_immune_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4)

FeaturePlot(mWAT_Immune_merged, features =  c("Ucp1","Cidea","Cyp2e1","Aldh1a1"), min.cutoff = "q10", order = T, pt.size=0.1, ncol = 2)
ggsave(paste0(mWAT_immune_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 6, height = 6)
saveRDS(mWAT_Immune_merged, file = paste0(mWAT_immune_prefix,"Seurat.rds"))
