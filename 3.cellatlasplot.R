
getwd()

library(Seurat)

library(tidyverse)

library(patchwork)

library(dplyr)

library(ggsci)

library(paletteer)

library(Scillus)

rm(list = ls())

scRNA = readRDS(file = "keyobject/scRNActid.rds")

table(scRNA$celltype)

table(scRNA$tissuetype)

dir.create("figure")

dir.create("figure/cellatlas")

dir.create("figure/cellatlas/celltype")

dir.create("figure/cellatlas/celltype/color2")

pal <- paletteer_d("ggsci::category20_d3")[c(4, 6, 12, 7, 9, 2,
                                             10, 3, 1, 8, 11, 13, 5)]

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "celltype", label = F, repel = T, cols = pal) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

plot1

ggsave("figure/cellatlas/celltype/color2/umapctc20.png", plot = plot1, width = 7, height = 5)

ggsave("figure/cellatlas/celltype/color2/umapctc20.pdf", plot = plot1, width = 7, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "celltype", label = F, repel = T) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

plot1

ggsave("figure/celltype/umapct.png", plot = plot1, width = 7, height = 5)

ggsave("figure/celltype/umapct.pdf", plot = plot1, width = 7, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

ggsave("figure/celltype/umapcluster.png", plot = plot1, width = 7.8, height = 5)

ggsave("figure/celltype/umapcluster.pdf", plot = plot1, width = 7.8, height = 5)

str(scRNA@meta.data)

pal <- paletteer_d("ggsci::category20_d3")[c(12, 3, 1, 9, 8, 7,
                                             6, 5, 4, 11, 2, 10)]

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissueorig", label = F, repel = T, cols = pal) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

plot1

ggsave("figure/cellatlas/celltype/color2/umaptissueorigc20.png", plot = plot1, width = 7, height = 5)

ggsave("figure/cellatlas/celltype/color2/umaptissueorigc20.pdf", plot = plot1, width = 7, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissueorig", label = F, repel = T) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

ggsave("figure/celltype/umaptissueorig.png", plot = plot1, width = 6.8, height = 5)

ggsave("figure/celltype/umaptissueorig.pdf", plot = plot1, width = 6.8, height = 5)

pal <- paletteer_d("ggsci::category20_d3")[c(1:20)]

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissuetype", label = F, repel = T, cols = pal) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

plot1

ggsave("figure/cellatlas/celltype/color2/umaptissuetypec20.png", plot = plot1, width = 6.3, height = 5)

ggsave("figure/cellatlas/celltype/color2/umaptissuetypec20.pdf", plot = plot1, width = 6.3, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissuetype", label = F, repel = T) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

ggsave("figure/celltype/umaptissuetype.png", plot = plot1, width = 6.3, height = 5)

ggsave("figure/celltype/umaptissuetype.pdf", plot = plot1, width = 6.3, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "cancertype", label = F, repel = T, cols = pal) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

ggsave("figure/celltype/umapcancertypec20.png", plot = plot1, width = 6.2, height = 5)

ggsave("figure/celltype/umapcancertypec20.pdf", plot = plot1, width = 6.2, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "cancertype", label = F, repel = T) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

ggsave("figure/celltype/umapcancertype.png", plot = plot1, width = 6.3, height = 5)

ggsave("figure/celltype/umapcancertype.pdf", plot = plot1, width = 6.3, height = 5)

table(scRNA$celltype)

markers = c("PTPRC", "CD79A", "MS4A1",
            "IGLC2", "IGHA1",
            "CD3D", "IL7R", "CD4",
            "CD8A", "CD8B",
            "CD1C",
            "PECAM1", "VWF",
            "EPCAM", "KRT8",
            "COL1A1", "DCN",
            "CD68", "MARCO",
            "CPA3", "KIT",
            "CD14", "S100A8",
            "NKG7", "GNLY",
            "RGS5", "ACTA2", "CSPG4")

plot1 = DotPlot(scRNA, features = markers, group.by = "celltype") + coord_flip() + theme_bw() + #去除背景，旋转图片
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=rel(1.1), face = "plain"),
        axis.text.y = element_text(size = rel(1.05), face = "plain")) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colors = c('#330066', '#336699', '#66CC66', '#FFCC33')) + #颜色渐变设置
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3))

ggsave("figure/celltype/cellmarker.png", plot1, width = 6, height = 7.5)

ggsave("figure/celltype/cellmarker.pdf", plot1, width = 6, height = 7.5)

select_genes <- c("EPCAM", "PTPRC", "PECAM1", "COL1A1", "CSPG4", "RGS5")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap",
                  label = F, label.size = 3, ncol = 3, cols = c("lightgrey", "red2"))

ggsave("figure/cellatlas/celltype/color2/celltype_marker.png", p2, width = 19, height = 12)

ggsave("figure/cellatlas/celltype/color2/celltype_marker.pdf", p2, width = 19, height = 12)

str(scRNA@meta.data)

scRNA$seurat_clusters = scRNA$celltype

plot1 = plot_stat(scRNA, plot_type = "prop_fill", group_by = "tissueorig", pal_setup = pal) +
  theme_bw() + coord_flip() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Cell type (Percentage)", x = "Tissue type") +
  guides(fill = guide_legend(title = "Cell types"))

ggsave("figure/celltype/celltypefreqc20.png", plot1, width = 7.2, height = 3.8)

ggsave("figure/celltype/celltypefreqc20.pdf", plot1, width = 7.2, height = 3.8)

pal <- paletteer_d("ggsci::category20_d3")[c(1:20)]

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissuetype",
                split.by = "tissuetype",
                label = F, repel = T, cols = pal) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

plot1

ggsave("figure/cellatlas/celltype/color2/umaptissuetypec20nt.png", plot = plot1, width = 10, height = 5)

ggsave("figure/cellatlas/celltype/color2/umaptissuetypec20nt.pdf", plot = plot1, width = 10, height = 5)

pal <- paletteer_d("ggsci::category20_d3")[c(4, 6, 12, 7, 9, 2,
                                             10, 3, 1, 8, 11, 13, 5)]

table(scRNA$celltype)

table(scRNA$seurat_clusters)

scRNA$seurat_clusters = scRNA$celltype

table(scRNA$seurat_clusters)

plot1 = plot_stat(scRNA, plot_type = "prop_fill", group_by = "tissueorig", pal_setup = pal) +
  theme_bw() + coord_flip() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Cell type (Percentage)", x = "Tissue type") +
  guides(fill = guide_legend(title = "Cell types"))

plot1

ggsave("figure/cellatlas/celltype/color2/celltypefreqc20.png", plot1, width = 6.5, height = 3.8)

ggsave("figure/cellatlas/celltype/color2/celltypefreqc20.pdf", plot1, width = 6.5, height = 3.8)

pal <- paletteer_d("ggsci::category20_d3")[c(4, 6, 12, 7, 9, 2,
                                             10, 3, 1, 8, 11, 13, 5)]

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "celltype",
                split.by = "tissuetype",
                label = F, repel = T, cols = pal) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

plot1

ggsave("figure/cellatlas/celltype/color2/umaptissuetypec20ct.png", plot = plot1, width = 10.5, height = 5)

ggsave("figure/cellatlas/celltype/color2/umaptissuetypec20ct.pdf", plot = plot1, width = 10, height = 5)

table(scRNA$celltype)

markers = c("PTPRC", "CD79A", "MS4A1",
            "IGLC2", "IGHA1",
            "CD3D", "IL7R", "CD4",
            "CD8A", "CD8B",
            "CD1C",
            "PECAM1", "VWF",
            "EPCAM", "KRT8",
            "COL1A1", "DCN",
            "CD68", "MARCO",
            "CPA3", "KIT",
            "CD14", "S100A8",
            "NKG7", "GNLY",
            "RGS5", "ACTA2", "CSPG4")

plot1 = DotPlot(scRNA, features = markers, group.by = "celltype") + theme_bw() + #去除背景，旋转图片
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=rel(1.1), face = "plain"),
        axis.text.y = element_text(size = rel(1.05), face = "plain")) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colors = c('#330066', '#336699', '#66CC66', '#FFCC33')) + #颜色渐变设置
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3))

plot1

ggsave("figure/celltype/cellmarkera.png", plot1, width = 9, height = 4.2)

ggsave("figure/celltype/cellmarkera.pdf", plot1, width = 9, height = 4.2)

select_genes <- c("CXCR4")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap",
                  label = F, cols = c("lightgrey", "red2"))

ggsave("figure/celltype/CXCR4umapcell.png", p2, width = 5.2, height = 5)

ggsave("figure/celltype/CXCR4umapcell.pdf", p2, width = 5.2, height = 5)

select_genes <- c("FOLH1")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap",
                  label = F, cols = c("lightgrey", "red2"))

ggsave("figure/celltype/FOLH1umapcell.png", p2, width = 5.2, height = 5)

ggsave("figure/celltype/FOLH1umapcell.pdf", p2, width = 5.2, height = 5)

surface = c("CA2", "INSR", "HTRA1", "MCAM", "CDH13", "PLVAP",
            "THY1", "TNFRSF4", "LAMA4", "PXDN", "FOLH1", "EDNRB")

p1 = VlnPlot(scRNA, features = surface, cols = pal, group.by = "celltype",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p1

ggsave(p1, file = "figure/celltype/surfacevio1allcell.png",
       width = 6.8, height = 5)

ggsave(p1, file = "figure/celltype/surfacevio1allcell.pdf",
       width = 6.8, height = 5)

p1 = VlnPlot(scRNA, features = "VEGFA", cols = pal, group.by = "celltype",
             stacked = T, pt.size = 0) #不显示坐标刻度

p1

ggsave("figure/celltype/vegfaexp.png", p1, width = 6.8, height = 3.9)

ggsave("figure/celltype/vegfaexp.pdf", p1, width = 6.8, height = 3.9)

p1 = VlnPlot(scRNA, features = "PGF", cols = pal, group.by = "celltype",
             stacked = T, pt.size = 0) #不显示坐标刻度

p1

ggsave("figure/celltype/pgfexp.png", p1, width = 6.8, height = 3.9)

ggsave("figure/celltype/pgfexp.pdf", p1, width = 6.8, height = 3.9)

select_genes <- c("PECAM1", "CLDN5", "CDH5", "PLVAP",
                  "VWF", "FLT1", "KDR", "FLT4")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap",
                  label = F, label.size = 3, ncol = 4,
                  cols = c("lightgrey", "red2"))

ggsave("figure/response/ecmkcellatlas.png", p2, width = 17, height = 8)

ggsave("figure/response/ecmkcellatlas.pdf", p2, width = 17, height = 8)

p1 = Seurat::VlnPlot(scRNA, features = "VEGFA", cols = pal, group.by = "celltype",
                     split.by = "tissueorig", pt.size = 0) #不显示坐标刻度

p1

ggsave("figure/response/vegfa-exp-cell-tissue.png", p1,
       width = 15, height = 3.5)

ggsave("figure/response/vegfa-exp-cell-tissue.pdf", p1,
       width = 15, height = 3.5)

p1 = Seurat::DotPlot(scRNA, features = "VEGFA",
                     group.by = "celltype",
                     split.by = "tissueorig") #不显示坐标刻度

p1

p1 = Seurat::VlnPlot(scRNA, features = "PGF", cols = pal, group.by = "celltype",
                     split.by = "tissueorig", pt.size = 0) #不显示坐标刻度

p1

ggsave("figure/response/pgf-exp-cell-tissue.png", p1,
       width = 15, height = 3.5)

ggsave("figure/response/pgf-exp-cell-tissue.pdf", p1,
       width = 15, height = 3.5)

pal <- paletteer_d("ggsci::category20_d3")[c(4, 6, 12, 7, 9, 2, 10,
                                             3, 1, 8, 11, 13, 5)]

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "celltype",
                split.by = "tissuetype",
                label = F, repel = T, cols = pal) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

plot1

ggsave("figure/cellatlas/umaptissuetypec20ct.png", plot = plot1, width = 10.5, height = 5)

ggsave("figure/cellatlas/umaptissuetypec20ct.pdf", plot = plot1, width = 10.5, height = 5)



