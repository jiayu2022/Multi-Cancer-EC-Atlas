
setwd("/home/zhangjiayu/project/EC")

getwd()

library(Seurat)

library(harmony)

library(tidyverse)

library(patchwork)

library(dplyr)

library(ggsci)

library(paletteer)

library(Scillus)

library(future)

rm(list = ls())

scRNA = readRDS(file = "keyobject/scRNActid.rds")

table(scRNA$celltype)

table(scRNA$seurat_clusters)

pal <- paletteer_d("ggsci::category20_d3")[c(1:20)]

dir.create("figure")

dir.create("figure/celltype")

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "celltype", label = F, repel = T, cols = pal) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

ggsave("figure/celltype/umapctc20.png", plot = plot1, width = 7, height = 5)

ggsave("figure/celltype/umapctc20.pdf", plot = plot1, width = 7, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "celltype", label = F, repel = T) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

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

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissueorig", label = F, repel = T, cols = pal) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

ggsave("figure/celltype/umaptissueorigc20.png", plot = plot1, width = 7, height = 5)

ggsave("figure/celltype/umaptissueorigc20.pdf", plot = plot1, width = 7, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissueorig", label = F, repel = T) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

ggsave("figure/celltype/umaptissueorig.png", plot = plot1, width = 6.8, height = 5)

ggsave("figure/celltype/umaptissueorig.pdf", plot = plot1, width = 6.8, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissuetype", label = F, repel = T, cols = pal) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

ggsave("figure/celltype/umaptissuetypec20.png", plot = plot1, width = 7, height = 5)

ggsave("figure/celltype/umaptissuetypec20.pdf", plot = plot1, width = 7, height = 5)

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

ggsave("figure/celltype/umapcancertypec20.png", plot = plot1, width = 7, height = 5)

ggsave("figure/celltype/umapcancertypec20.pdf", plot = plot1, width = 7, height = 5)

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
            "RGS5", "ACTA2")

plot1 = DotPlot(scRNA, features = markers, group.by = "celltype") + coord_flip() + theme_bw() + #去除背景，旋转图片
            theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=rel(1.1), face = "plain"),
            axis.text.y = element_text(size = rel(1.05), face = "plain")) +
            scale_color_gradientn(values = seq(0, 1, 0.2), colors = c('#330066', '#336699', '#66CC66', '#FFCC33')) + #颜色渐变设置
            labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3))

ggsave("figure/celltype/cellmarker.png", plot1, width = 6, height = 7.5)

ggsave("figure/celltype/cellmarker.pdf", plot1, width = 6, height = 7.5)

select_genes <- c("EPCAM", "PTPRC", "PECAM1", "COL1A1", "CSPG4", "RGS5")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = F, label.size = 3, ncol = 3)

ggsave("figure/celltype/celltype_marker.png", p2, width = 20, height = 12)

ggsave("figure/celltype/celltype_marker.pdf", p2, width = 20, height = 12)

str(scRNA@meta.data)

scRNA$seurat_clusters = scRNA$celltype

plot1 = plot_stat(scRNA, plot_type = "prop_fill", group_by = "tissueorig", pal_setup = pal) +
            theme_bw() + coord_flip() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            labs(y = "Cell type (Percentage)", x = "Tissue type") +
            guides(fill = guide_legend(title = "Cell types"))

ggsave("figure/celltype/celltypefreqc20.png", plot1, width = 6.2, height = 3.8)

ggsave("figure/celltype/celltypefreqc20.pdf", plot1, width = 6.2, height = 3.8)

col.num <- length(levels(scRNA@active.ident))

violin <- VlnPlot(scRNA,
                  features = c("nFeature_RNA", "nCount_RNA"),
                  group.by = "seurat_clusters",
                  pt.size = 0,
                  ncol = 2)

ggsave("process/1.celltype/vlnplot_qc.png", plot = violin, width = 20, height = 5)



