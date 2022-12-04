
setwd("/home/zhangjiayu/project/EC")

getwd()

library(Seurat)

library(harmony)

library(tidyverse)

library(patchwork)

library(dplyr)

library(future)

rm(list = ls())

plan()

scRNA = readRDS(file = "process/scRNAmerge.rds")

##数据标准化，寻找高变基因，数据中心化

scRNA <- NormalizeData(scRNA)

scRNA <- FindVariableFeatures(scRNA)

table(scRNA$tissuetype)

scRNA <- ScaleData(scRNA)

#PCA

scRNA <- RunPCA(scRNA, verbose = FALSE)

str(scRNA@assays)

str(scRNA@meta.data)

saveRDS(scRNA, file = "process/scRNApca.rds")

dir.create("process/pca")

plot1 = DimPlot(scRNA, reduction = "pca", group.by = "tissueorig")

ggsave("process/pca/pca1.pdf", plot = plot1, width = 6.5, height = 5)

plot2 = DimPlot(scRNA, reduction = "pca", group.by = "tissuetype")

ggsave("process/pca/pca2.pdf", plot = plot2, width = 6, height = 5)

plot3 = DimPlot(scRNA, reduction = "pca", group.by = "cancertype")

ggsave("process/pca/pca3.pdf", plot = plot3, width = 6, height = 5)

plot5 = DimPlot(scRNA, reduction = "pca", group.by = "Phase")

ggsave("process/pca/pca5.pdf", plot = plot5, width = 6, height = 5)

plot6 = DimPlot(scRNA, reduction = "pca")

ggsave("process/pca/pca6.pdf", plot = plot6, width = 6, height = 5)

plot6 = FeaturePlot(scRNA, features = "percent.mt")

ggsave("process/pca/pca7.pdf", plot = plot6, width = 6, height = 5)

plot6 = FeaturePlot(scRNA, features = "percent.rb")

ggsave("process/pca/pca8.pdf", plot = plot6, width = 6, height = 5)

plot6 = FeaturePlot(scRNA, features = "percent.HB")

ggsave("process/pca/pca9.pdf", plot = plot6, width = 6, height = 5)

plot6 = FeaturePlot(scRNA, features = "nFeature_RNA")

ggsave("process/pca/pca10.pdf", plot = plot6, width = 6, height = 5)

#Harmony

dir.create("process/harmony")

scRNA <- RunHarmony(scRNA, group.by.vars = "sampleid")

saveRDS(scRNA, file = "process/scRNAharmony.rds")

#UMAP, clustering

scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:50)

scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:50)

scRNA <- FindClusters(scRNA, resolution = 2)

saveRDS(scRNA, file = "process/scRNAharumap.rds")

# plot

dir.create("process/1.celltype")

dir.create("process/1.celltype/cluster")

scRNA = readRDS(file = "process/scRNAharumap.rds")

plot2 = DimPlot(scRNA, reduction = "umap", label = F, label.size = 3)

ggsave("process/1.celltype/cluster/UMAP.pdf", plot = plot2, width = 8, height = 5)

plot2 = DimPlot(scRNA, reduction = "umap", label = T, label.size = 2)

ggsave("process/1.celltype/cluster/UMAP0.pdf", plot = plot2, width = 8, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissueorig")

ggsave("process/1.celltype/cluster/UMAP1.pdf", plot = plot1, width = 6.5, height = 5)

plot3 = DimPlot(scRNA, reduction = "umap", group.by = "tissuetype")

ggsave("process/1.celltype/cluster/UMAP2.pdf", plot = plot3, width = 6, height = 5)

plot3 = DimPlot(scRNA, reduction = "umap", group.by = "cancertype")

ggsave("process/1.celltype/cluster/UMAP3.pdf", plot = plot3, width = 6, height = 5)

### plot-marker

select_genes <- c("EPCAM", "PTPRC", "PECAM1", "CSPG4", "COL1A1", "ACTA2")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 2, ncol = 3)

ggsave("process/1.celltype/cluster/1celltype_marker.pdf", p2, width = 20, height = 12)

select_genes <- c("EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "CDH1")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 2, ncol = 3)

ggsave("process/1.celltype/cluster/2epi.pdf", p2, width = 18, height = 12)

select_genes <- c("PECAM1", "CDH5", "VWF", "CLDN5")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 2, ncol = 2)

ggsave("process/1.celltype/cluster/3EC.pdf", p2, width = 12, height = 12)

select_genes1 <- c("CSPG4", "MCAM", "RGS5", "ACTA2", "TAGLN", "MYLK")

p1 <- FeaturePlot(scRNA, features = select_genes1, reduction = "umap", label = T, label.size = 2, ncol = 3)

ggsave("process/1.celltype/cluster/4pericyteSMC.pdf", p1, width = 18, height = 12)

select_genes2 <- c("ACTA2", "COL1A1", "FAP", "THY1", "PDPN", "S100A4", "CAV1", "VIM")

p1 <- FeaturePlot(scRNA, features = select_genes2, reduction = "umap", label = T, label.size = 2, ncol = 4)

ggsave("process/1.celltype/cluster/5CAF.pdf", p1, width = 24, height = 12)

select_genes2 <- c("CD248", "ACTA2", "COL1A1", "FAP", "THY1", "PDPN", "S100A4", "CAV1")

p1 <- FeaturePlot(scRNA, features = select_genes2, reduction = "umap", label = T, label.size = 2, ncol = 4)

ggsave("process/1.celltype/cluster/5CAF2.pdf", p1, width = 24, height = 12)

select_genes <- c("CD3D", "CD3G", "CD3E",
                  "CD4", "CD8A", "CD8B")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 2, ncol = 3)

ggsave("process/1.celltype/cluster/6T.pdf", p2, width = 18, height = 12)

select_genes <- c("CD3D", "CD3E", "CD4", "IL7R", "CD8A", "CD8B")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 2, ncol = 3)

ggsave("process/1.celltype/cluster/6T48.pdf", p2, width = 18, height = 12)

select_genes <- c("MS4A1", "CD79A", "NKG7", "GNLY", "CD14", "S100A8")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 2, ncol = 3)

ggsave("process/1.celltype/cluster/7BNM.pdf", p2, width = 18, height = 12)

select_genes <- c("FCGR3A", "CD68", "CSF1R", "CD163", "CD1C", "CD1A")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 2, ncol = 3)

ggsave("process/1.celltype/cluster/8MoMa.pdf", p2, width = 18, height = 12)

select_genes <- c("CPA3", "KIT", "FCGR3B")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 2, ncol = 3)

ggsave("process/1.celltype/cluster/9MastN.pdf", p2, width = 18, height = 6)

cellnumber = table(scRNA$seurat_clusters)

write.csv(cellnumber, "process/1.celltype/cluster/cellnumber.csv", row.names = F)

col.num <- length(levels(scRNA@active.ident))

violin <- VlnPlot(scRNA,
                  features = c("nFeature_RNA"),
                  cols = rainbow(col.num),
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 1) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("process/1.celltype/cluster/vlnplot_qc1.pdf", plot = violin, width = 15, height = 5)

violin <- VlnPlot(scRNA,
                  features = c("percent.mt", "percent.HB"),
                  cols = rainbow(col.num),
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 2) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("process/1.celltype/cluster/vlnplot_qc2.pdf", plot = violin, width = 20, height = 3)

# DEG

scRNA = readRDS(file = "process/scRNAharumap.rds")

plan("multicore", workers = 56)

options(future.globals.maxSize = 10 * 1024^3)

plan()

diff.wilcox = FindAllMarkers(scRNA)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

saveRDS(all.markers, file = "process/1.celltype/cluster/markers.rds")

#all.markers = readRDS(file = "process/1.celltype/cluster/markers.rds")

write.csv(all.markers, "process/1.celltype/cluster/diff_genes_wilcox.csv", row.names = F)

top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top20, "process/1.celltype/cluster/top20_diff_genes_wilcox.csv", row.names = F)

top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(top50, "process/1.celltype/cluster/top50_diff_genes_wilcox.csv", row.names = F)

table(scRNA$seurat_clusters)

n = length(unique(scRNA@meta.data$seurat_clusters)) #获取亚群个数

celltype = data.frame(ClusterID = 0:(n - 1),
                      celltype = "NA")  #构建数据框

## 判断亚群ID属于那类细胞

celltype[celltype$ClusterID %in% c(3, 71), 2] = "B_cells"

celltype[celltype$ClusterID %in% c(11, 31, 60), 2] = "B_Plasma_cells"

celltype[celltype$ClusterID %in% c(0, 9, 36, 62), 2] = "CD4+T_cells"

celltype[celltype$ClusterID %in% c(1, 2, 30, 64, 65), 2] = "CD8+T_cells"

celltype[celltype$ClusterID %in% c(16), 2] = "Dendritic_cells"

celltype[celltype$ClusterID %in% c(8, 22, 54, 59, 77), 2] = "Endothelial_cells"

celltype[celltype$ClusterID %in% c(7, 13, 15, 17, 20, 23, 24, 25, 33,
                                   37, 38, 39, 41, 43, 44, 45, 46, 47,
                                   48, 52, 53, 56, 57, 61, 67, 69, 73,
                                   76, 81), 2] = "Epithelial_cells"

celltype[celltype$ClusterID %in% c(10, 21, 26, 27, 28, 35, 40, 50, 63, 66, 68), 2] = "Fibroblasts"

celltype[celltype$ClusterID %in% c(4, 6, 14, 29, 32, 34, 49), 2] = "Macrophages"

celltype[celltype$ClusterID %in% c(19, 83), 2] = "Mast_cells"

celltype[celltype$ClusterID %in% c(5, 42, 55), 2] = "Monocytes"

celltype[celltype$ClusterID %in% c(18), 2] = "NK_cells"

celltype[celltype$ClusterID %in% c(12), 2] = "Pericytes"

## 重新赋值

scRNA@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]), "celltype"] <- celltype$celltype[i]}

table(scRNA@meta.data$celltype)

table(scRNA$seurat_clusters)

Cells.sub <- subset(scRNA@meta.data, celltype != "NA")

scRNAsub <- subset(scRNA, cells = row.names(Cells.sub))

table(scRNAsub$celltype)

scRNA = scRNAsub

table(Idents(scRNA))

table(scRNA$seurat_clusters)

Idents(scRNA) = scRNA$celltype

table(Idents(scRNA))

dir.create("keyobject")

saveRDS(scRNA, file = "keyobject/scRNActid.rds")

# DEG

plan("multicore", workers = 56)

options(future.globals.maxSize = 10 * 1024^3)

plan()

diff.wilcox = FindAllMarkers(scRNA)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

saveRDS(all.markers, file = "process/1.celltype/markers.rds")

write.csv(all.markers, "process/1.celltype/diff_genes_wilcox.csv", row.names = F)

top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top20, "process/1.celltype/top20_diff_genes_wilcox.csv", row.names = F)

top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(top50, "process/1.celltype/top50_diff_genes_wilcox.csv", row.names = F)

Cells.sub <- subset(scRNA@meta.data, celltype == "Endothelial_cells")

scRNAsub <- subset(scRNA, cells = row.names(Cells.sub))

table(scRNAsub$celltype)

saveRDS(scRNAsub, file = "keyobject/scRNAec.rds")


