
setwd("/home/zhangjiayu/project/EC")

getwd()

library(Seurat)

library(tidyverse)

library(patchwork)

library(dplyr)

library(ggsci)

library(paletteer)

library(Scillus)

rm(list = ls())

pal <- paletteer_d("ggsci::category20_d3")[c(1:20)]

scRNA = readRDS(file = "keyobject/ecqchumap.rds")

table(scRNA$celltype)

scRNA@meta.data$ecsubtype = scRNA$seurat_clusters

table(scRNA$ecsubtype)

dir.create("figure/ecsubtype")

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "ecsubtype", label = F, repel = T, cols = pal) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

ggsave("figure/ecsubtype/umapctc20.png", plot = plot1, width = 6, height = 5)

ggsave("figure/ecsubtype/umapctc20.pdf", plot = plot1, width = 6, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "ecsubtype", label = T, repel = T) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

ggsave("figure/ecsubtype/umapcta.png", plot = plot1, width = 6, height = 5)

ggsave("figure/ecsubtype/umapcta.pdf", plot = plot1, width = 6, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissueorig", label = F, repel = T, cols = pal) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

ggsave("figure/ecsubtype/umaptissueorigc20.png", plot = plot1, width = 7, height = 5)

ggsave("figure/ecsubtype/umaptissueorigc20.pdf", plot = plot1, width = 7, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissueorig", label = F, repel = T) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

ggsave("figure/ecsubtype/umaptissueorig.png", plot = plot1, width = 6.8, height = 5)

ggsave("figure/ecsubtype/umaptissueorig.pdf", plot = plot1, width = 6.8, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissuetype", label = F, repel = T, cols = pal) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

ggsave("figure/ecsubtype/umaptissuetypec20.png", plot = plot1, width = 6.3, height = 5)

ggsave("figure/ecsubtype/umaptissuetypec20.pdf", plot = plot1, width = 6.3, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissuetype", label = F, repel = T) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

ggsave("figure/ecsubtype/umaptissuetype.png", plot = plot1, width = 6.3, height = 5)

ggsave("figure/ecsubtype/umaptissuetype.pdf", plot = plot1, width = 6.3, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "cancertype", label = F, repel = T, cols = pal) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

ggsave("figure/ecsubtype/umapcancertypec20.png", plot = plot1, width = 7, height = 5)

ggsave("figure/ecsubtype/umapcancertypec20.pdf", plot = plot1, width = 7, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "cancertype", label = F, repel = T) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

ggsave("figure/ecsubtype/umapcancertype.png", plot = plot1, width = 6.3, height = 5)

ggsave("figure/ecsubtype/umapcancertype.pdf", plot = plot1, width = 6.3, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "Phase", label = F, repel = T, cols = pal) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

ggsave("figure/ecsubtype/umapccphasec20.png", plot = plot1, width = 6.2, height = 5)

ggsave("figure/ecsubtype/umapccphasec20.pdf", plot = plot1, width = 6.2, height = 5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "Phase", label = F, repel = T) + #
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = ""))

ggsave("figure/ecsubtype/umapccphase.png", plot = plot1, width = 6.2, height = 5)

ggsave("figure/ecsubtype/umapccphase.pdf", plot = plot1, width = 6.2, height = 5)

select_genes <- c("EPCAM", "PTPRC", "PECAM1", "COL1A1", "CSPG4", "RGS5")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = F, label.size = 3, ncol = 3)

ggsave("figure/ecsubtype/celltype_marker.png", p2, width = 20, height = 12)

ggsave("figure/ecsubtype/celltype_marker.pdf", p2, width = 20, height = 12)

str(scRNA@meta.data)

plot1 = plot_stat(scRNA, plot_type = "prop_fill", group_by = "tissueorig", pal_setup = pal) +
            theme_bw() + coord_flip() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            labs(y = "EC subtype (Percentage)", x = "Tissue type") +
            guides(fill = guide_legend(title = "EC subtypes"))

ggsave("figure/ecsubtype/celltypefreqc20.png", plot1, width = 6.2, height = 3.8)

ggsave("figure/ecsubtype/celltypefreqc20.pdf", plot1, width = 6.2, height = 3.8)

all.markers = readRDS(file = "process/2.ec/cluster3/markers.rds")

top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)

dir.create("figure/ecsubmk")

write.csv(top15, "figure/ecsubmk/top15deg.csv", row.names = F)

features = c(top15$gene)

mk <- AverageExpression(scRNA,
                        group.by = "ecsubtype",
                        features = features,
                        return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

p <- pheatmap(mfdf,
              color = colorRampPalette(c("#0080FF", "white", "#FF8000"))(20),
              main = "",
              scale = "row",
              angle_col = 0, # 设置显示角度
              #cellwidth = 20,cellheight = , # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 3.5, # 分别设置横向和纵向字体大小
              fontsize_col = 6)

ggsave("figure/ecsubmk/top15mk.png", p, width = 3.5, height = 8)

ggsave("figure/ecsubmk/top15mk.pdf", p, width = 3.5, height = 8)

col.num <- length(levels(scRNA@active.ident))

violin <- Seurat::VlnPlot(scRNA,
                          features = c("nFeature_RNA", "nCount_RNA"),
                          cols = rainbow(col.num),
                          pt.size = 0.01,
                          ncol = 2)

ggsave("figure/ecsubtype/vlnplot_qca.png", plot = violin, width = 10, height = 6)

violin <- VlnPlot(scRNA,
                  features = c("percent.mt", "percent.HB", "percent.rb"),
                  cols = rainbow(col.num),
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 3) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("figure/ecsubtype/vlnplot_qca1.png", plot = violin, width = 12, height = 5)

select_genes <- c("CXCR4", "APLNR", "GJA5", "ACKR1", "CA4", "PROX1")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 3, ncol = 3)

ggsave("figure/ecsubtype/ecsubmarker.pdf", p2, width = 20, height = 12)

select_genes <- c("GJA5", "SOX17", "ACKR1", "VCAM1", "CA4", "PROX1")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 3, ncol = 3)

ggsave("figure/ecsubmk/ecsubmarker.pdf", p2, width = 20, height = 12)

select_genes <- c("CXCR4", "ESM1", "JAG1", "ENG", "PLVAP", "KDR")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 3, ncol = 3)

ggsave("figure/ecsubmk/ecsubtipimmarker.pdf", p2, width = 20, height = 12)

select_genes <- c("CA4", "EDNRB", "CAV1", "AKAP12", "TIMP3", "EDN1")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 3, ncol = 3)

ggsave("figure/ecsubmk/ecsubcapmarker.pdf", p2, width = 20, height = 12)

#ecsubid

table(scRNA$celltype)

table(scRNA$seurat_clusters)

n = length(unique(scRNA@meta.data$seurat_clusters)) #获取亚群个数

n

celltype = data.frame(ClusterID = 0:(n - 1),
                      celltype = "NA")  #构建数据框

## 判断亚群ID属于那类细胞

celltype[celltype$ClusterID %in% c(4, 7), 2] = "arterial_ECs"

celltype[celltype$ClusterID %in% c(2, 5), 2] = "capillary_ECs_1"

celltype[celltype$ClusterID %in% c(3, 9), 2] = "capillary_ECs_2"

celltype[celltype$ClusterID %in% c(6), 2] = "lymphatic_ECs"

celltype[celltype$ClusterID %in% c(0), 2] = "tip_like_ECs"

celltype[celltype$ClusterID %in% c(1, 8), 2] = "venous_ECs"

## 重新赋值

scRNA@meta.data$subtype = "NA"

for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]), "subtype"] <- celltype$celltype[i]}

table(scRNA@meta.data$subtype)

table(scRNA$seurat_clusters)

saveRDS(scRNA, file = "keyobject/ecsubid.rds")

table(Idents(scRNA))

Idents(scRNA) = scRNA$subtype

table(Idents(scRNA))

#deg_wilcox

plan()

plan("multicore", workers = 56)

options(future.globals.maxSize = 10 * 1024^3)

plan()

diff.wilcox = FindAllMarkers(scRNA)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

dir.create("process/2.ec/subid")

saveRDS(all.markers, file = "process/2.ec/subid/markers.rds")

#all.markers = readRDS(file = "process/2.ec/subid/markers.rds")

write.csv(all.markers, "process/2.ec/subid/degs_wilcox.csv", row.names = F)

top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)

write.csv(top15, "process/2.ec/subid/top15_degs.csv", row.names = F)

top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top20, "process/2.ec/subid/top20_degs.csv", row.names = F)

top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(top50, "process/2.ec/subid/top50_degs.csv", row.names = F)

select_genes <- c("CXCR4", "ESM1", "JAG1", "ENG", "PLVAP", "KDR")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 3, ncol = 3)

ggsave("process/2.ec/subid/ecsubtipimmarker.pdf", p2, width = 20, height = 12)

select_genes <- c("CA4", "EDNRB", "CAV1", "AKAP12", "TIMP3", "EDN1")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 3, ncol = 3)

ggsave("process/2.ec/subid/ecsubcapmarker.pdf", p2, width = 20, height = 12)

select_genes <- c("GJA5", "ACKR1", "CA4", "PROX1", "CXCR4", "ESM1")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 3, ncol = 3)

ggsave("process/2.ec/subid/ecsubmarkerid.pdf", p2, width = 20, height = 12)

select_genes <- c("CXCR4", "ESM1", "EDNRB", "EDN1", "AKAP12", "KDR")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", label = T, label.size = 3, ncol = 3)

ggsave("process/2.ec/subid/ecsubmarkerid2.pdf", p2, width = 20, height = 12)


