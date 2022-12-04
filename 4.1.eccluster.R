
setwd("/home/zhangjiayu/project/EC")

getwd()

library(Seurat)

library(harmony)

library(tidyverse)

library(patchwork)

library(dplyr)

library(future)

rm(list = ls())

scRNA = readRDS(file = "keyobject/scRNAec.rds")

table(scRNA$celltype)

##数据标准化，寻找高变基因，数据中心化

scRNA <- NormalizeData(scRNA)

scRNA <- FindVariableFeatures(scRNA)

table(scRNA$tissuetype)

scRNA <- ScaleData(scRNA)

#PCA

scRNA <- RunPCA(scRNA, verbose = FALSE)

str(scRNA@assays)

str(scRNA@meta.data)

#Harmony

dir.create("process/2.ec")

scRNA <- RunHarmony(scRNA, group.by.vars = "sampleid")

#UMAP, clustering

scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:50)

scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:50)

scRNA <- FindClusters(scRNA, resolution = 0.2)

dir.create("process/2.ec/cluster")

plot2 = DimPlot(scRNA, reduction = "umap", label = F, label.size = 3)

ggsave("process/2.ec/cluster/UMAP.pdf", plot = plot2, width = 5.8, height = 5)

plot2 = DimPlot(scRNA, reduction = "umap", label = T, label.size = 2)

ggsave("process/2.ec/cluster/UMAP0.pdf", plot = plot2, width = 5.8, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissueorig")

ggsave("process/2.ec/cluster/UMAP1.pdf", plot = plot1, width = 6.5, height = 5)

plot3 = DimPlot(scRNA, reduction = "umap", group.by = "tissuetype")

ggsave("process/2.ec/cluster/UMAP2.pdf", plot = plot3, width = 6, height = 5)

plot3 = DimPlot(scRNA, reduction = "umap", group.by = "cancertype")

ggsave("process/2.ec/cluster/UMAP3.pdf", plot = plot3, width = 6, height = 5)

plot3 = DimPlot(scRNA, reduction = "umap", group.by = "Phase")

ggsave("process/2.ec/cluster/UMAPcc.pdf", plot = plot3, width = 6, height = 5)

###

select_genes <- c("EPCAM", "PTPRC", "PECAM1", "CSPG4", "COL1A1", "ACTA2")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = F, label.size = 3, ncol = 3)

ggsave("process/2.ec/cluster/1celltype_marker.pdf", p2, width = 20, height = 12)

select_genes <- c("EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "CDH1")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = F, label.size = 3, ncol = 3)

ggsave("process/2.ec/cluster/2epi.pdf", p2, width = 18, height = 12)

select_genes <- c("PECAM1", "CDH5", "VWF", "CLDN5")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = F, label.size = 3, ncol = 2)

ggsave("process/2.ec/cluster/3EC.pdf", p2, width = 12, height = 12)

select_genes1 <- c("CSPG4", "MCAM", "RGS5", "ACTA2", "TAGLN", "MYLK")

p1 <- FeaturePlot(scRNA, features = select_genes1, reduction = "umap", label = F, label.size = 3, ncol = 3)

ggsave("process/2.ec/cluster/4pericyteSMC.pdf", p1, width = 18, height = 12)

select_genes2 <- c("ACTA2", "COL1A1", "FAP", "THY1", "PDPN")

p1 <- FeaturePlot(scRNA, features = select_genes2, reduction = "umap", label = F, label.size = 3, ncol = 3)

ggsave("process/2.ec/cluster/5CAF.pdf", p1, width = 18, height = 12)

select_genes <- c("CD3D", "CD3G", "CD3E",
                  "CD68", "CD14", "CD79A")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = F, label.size = 3, ncol = 3)

ggsave("process/2.ec/cluster/6TM.pdf", p2, width = 18, height = 12)

table(scRNA$seurat_clusters)

saveRDS(scRNA, file = "keyobject/echarumap.rds")

scRNA = readRDS(file = "keyobject/echarumap.rds")

select_genes <- c("CXCR4", "APLNR", "GJA5", "ACKR1", "CA4", "PROX1")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 3, ncol = 3)

ggsave("process/2.ec/cluster2/ecsubmarker1.pdf", p2, width = 20, height = 12)

select_genes <- c("GJA5", "SOX17", "ACKR1", "VCAM1", "CA4", "PROX1")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 3, ncol = 3)

ggsave("process/2.ec/cluster2/ecsubmarker.pdf", p2, width = 20, height = 12)

select_genes <- c("CXCR4", "ESM1", "JAG1", "ENG", "PLVAP", "KDR")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 3, ncol = 3)

ggsave("process/2.ec/cluster2/ecsubtipimmarker.pdf", p2, width = 20, height = 12)

select_genes <- c("CA4", "EDNRB", "CAV1", "AKAP12", "TIMP3", "EDN1")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 3, ncol = 3)

ggsave("process/2.ec/cluster2/ecsubcapmarker.pdf", p2, width = 20, height = 12)

#DEG

plan()

plan("multicore", workers = 56)

options(future.globals.maxSize = 10 * 1024^3)

plan()

diff.wilcox = FindAllMarkers(scRNA)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

saveRDS(all.markers, file = "process/2.ec/cluster/markers.rds")

#all.markers = readRDS(file = "process/2.ec/cluster/markers.rds")

write.csv(all.markers, "process/2.ec/cluster/diff_genes_wilcox.csv", row.names = F)

top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top20, "process/2.ec/cluster/top20_deg.csv", row.names = F)

top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(top50, "process/2.ec/cluster/top50_deg_wilcox.csv", row.names = F)

table(scRNA$seurat_clusters)

features = c(top20$gene)

library(pheatmap) # 加载包

mk <- AverageExpression(
  scRNA,
  features = features,
  return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

p <- pheatmap(mfdf,
              main = "",
              scale = "row",
              angle_col = 270, # 设置显示角度
              #cellwidth = 20,cellheight = , # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 1, # 分别设置横向和纵向字体大小
              fontsize_col = 20)

ggsave("process/2.ec/cluster/top50mk.png", p, width = 10, height = 15)

col.num <- length(levels(scRNA@active.ident))

violin <- Seurat::VlnPlot(scRNA,
                          features = c("nFeature_RNA", "nCount_RNA"),
                          cols = rainbow(col.num),
                          pt.size = 0.01,
                          ncol = 2)

ggsave("process/2.ec/cluster/vlnplot_qca.png", plot = violin, width = 10, height = 6)

violin <- VlnPlot(scRNA,
                  features = c("percent.mt", "percent.HB", "percent.rb"),
                  cols = rainbow(col.num),
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 3) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("process/2.ec/cluster/vlnplot_qca1.png", plot = violin, width = 12, height = 5)

table(scRNA$seurat_clusters)

Cells.sub <- subset(scRNA@meta.data, seurat_clusters != 5 & seurat_clusters != 8 & seurat_clusters != 10 &
                    seurat_clusters != 11 & seurat_clusters != 12 & seurat_clusters != 13 & seurat_clusters != 14 & seurat_clusters != 15)

scRNAsub <- subset(scRNA, cells = row.names(Cells.sub))

table(scRNAsub$seurat_clusters)

scRNA = scRNAsub

saveRDS(scRNA, file = "keyobject/echuqc.rds")

dir.create("process/2.ec/cluster2")

##数据标准化，寻找高变基因，数据中心化

scRNA <- NormalizeData(scRNA)

scRNA <- FindVariableFeatures(scRNA)

table(scRNA$tissuetype)

scRNA <- ScaleData(scRNA)

#PCA

scRNA <- RunPCA(scRNA, verbose = FALSE)

#Harmony

scRNA <- RunHarmony(scRNA, group.by.vars = "sampleid")

saveRDS(scRNA, file = "keyobject/ecqcharmony.rds")

scRNA = readRDS(file = "keyobject/ecqcharmony.rds")

#UMAP, clustering

scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:50)

scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:50)

scRNA <- FindClusters(scRNA, resolution = 0.5)

dir.create("process/2.ec/cluster2")

plot2 = DimPlot(scRNA, reduction = "umap", label = F, label.size = 3)

ggsave("process/2.ec/cluster2/UMAP.png", plot = plot2, width = 5.8, height = 5)

plot2 = DimPlot(scRNA, reduction = "umap", label = T, label.size = 2)

ggsave("process/2.ec/cluster2/UMAP0.png", plot = plot2, width = 5.8, height = 5)

select_genes <- c("GJA5", "ACKR1", "CA4", "PROX1", "CXCR4", "ESM1")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 3, ncol = 3)

ggsave("process/2.ec/cluster2/ecsubmarkerid.pdf", p2, width = 20, height = 12)

select_genes <- c("CXCR4", "ESM1", "EDNRB", "EDN1", "AKAP12", "KDR")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 3, ncol = 3)

ggsave("process/2.ec/cluster2/ecsubmarkerid2.pdf", p2, width = 20, height = 12)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissueorig")

ggsave("process/2.ec/cluster3/UMAP1.png", plot = plot1, width = 6.5, height = 5)

plot3 = DimPlot(scRNA, reduction = "umap", group.by = "tissuetype")

ggsave("process/2.ec/cluster3/UMAP2.png", plot = plot3, width = 6, height = 5)

plot3 = DimPlot(scRNA, reduction = "umap", group.by = "cancertype")

ggsave("process/2.ec/cluster3/UMAP3.png", plot = plot3, width = 6, height = 5)

plot3 = DimPlot(scRNA, reduction = "umap", group.by = "Phase")

ggsave("process/2.ec/cluster3/UMAPcc.png", plot = plot3, width = 6, height = 5)

###

table(scRNA$tissueorig)

plan()

#deg_wilcox

diff.wilcox = FindAllMarkers(scRNA)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

saveRDS(all.markers, file = "process/2.ec/cluster3/markers.rds")

#all.markers = readRDS(file = "process/2.ec/cluster3/markers.rds")

write.csv(all.markers, "process/2.ec/cluster3/degs_wilcox.csv", row.names = F)

top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top20, "process/2.ec/cluster3/top20_degs.csv", row.names = F)

top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(top50, "process/2.ec/cluster3/top50_degs.csv", row.names = F)

table(scRNA$seurat_clusters)

features = c(top20$gene)

library(pheatmap) # 加载包

mk <- AverageExpression(
  scRNA,
  features = features,
  return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

p <- pheatmap(mfdf,
              main = "",
              scale = "row",
              angle_col = 270, # 设置显示角度
              #cellwidth = 20,cellheight = , # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 1, # 分别设置横向和纵向字体大小
              fontsize_col = 20)

ggsave("process/2.ec/cluster3/top20mk.png", p, width = 10, height = 15)

table(scRNA$seurat_clusters)

saveRDS(scRNA, file = "keyobject/ecqchumap.rds")





