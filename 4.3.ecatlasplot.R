
getwd()

library(Seurat)

library(tidyverse)

library(patchwork)

library(dplyr)

library(ggsci)

library(paletteer)

library(Scillus)

library(pheatmap) # 加载包

library(dplyr)

library(MySeuratWrappers)

rm(list = ls())

scRNA = readRDS(file = "keyobject/ecsubid.rds")

table(scRNA$subtype)

table(scRNA$tissuetype)

dir.create("figure/ecsub")

dir.create("figure/ecsub/ecsubid")

dir.create("figure/ecsub/ecsubid/subid")

pal <- paletteer_d("ggsci::category20_d3")[c(5, 3, 4, 6, 2, 1)]

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "subtype", label = F, repel = T, cols = pal) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

plot1

ggsave("figure/ecsub/ecsubid/subid/umapctc20.png", plot = plot1, width = 6.8, height = 5)

ggsave("figure/ecsub/ecsubid/subid/umapctc20.pdf", plot = plot1, width = 6.8, height = 5)

scRNA$seurat_clusters = scRNA$subtype

plot1 = plot_stat(scRNA, plot_type = "prop_fill", group_by = "tissueorig", pal_setup = pal) +
  theme_bw() + coord_flip() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "EC subtype (Percentage)", x = "Tissue type") +
  guides(fill = guide_legend(title = "EC subtypes"))

plot1

ggsave("figure/ecsub/ecsubid/subid/celltypefreqc20.png", plot1, width = 6.8, height = 3.8)

ggsave("figure/ecsub/ecsubid/subid/celltypefreqc20.pdf", plot1, width = 6.8, height = 3.8)

top15 = read.csv(file = "figure/ecsub/ecsubid/top15a.csv")

features = c(top15$gene)

mk <- AverageExpression(scRNA,
                        group.by = "subtype",
                        features = features,
                        return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

p <- pheatmap(mfdf,
              color = colorRampPalette(c("#0080FF", "white", "#FF8000"))(20),
              main = "",
              scale = "row",
              angle_col = 90, # 设置显示角度
              #cellwidth = 20,cellheight = , # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 3.8, # 分别设置横向和纵向字体大小
              fontsize_col = 6)

ggsave("figure/ecsub/ecsubid/subid/top15mk.png", p, width = 3, height = 6.5)

ggsave("figure/ecsub/ecsubid/subid/top15mk.pdf", p, width = 3, height = 6.5)

str(scRNA@meta.data)

meta <- scRNA@meta.data[, c("tissueorig", "subtype")]

meta[1:5, 1:2]

stat = meta %>%
  group_by(tissueorig, subtype) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

write.csv(stat, "figure/ecsub/ecsubprop/ecsubfreq.csv", row.names = F)

stat = read.csv(file = "figure/ecsub/ecsubprop/ecsubfreqa.csv")

stat$tissue = stat$tissueorig

stat = separate(stat, col = tissue, into = c("cancertype", "tissuetype"), sep = "_")

str(stat)

stat$prop = stat$freq * 100

table(stat$subtype)

statsub = subset(stat, subtype == "venous_ECs")

plot1 = ggplot(statsub, aes(x = factor(cancertype), y = prop, fill = tissuetype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Venous EC / Total EC (Percentage)", x = "Cancer type") +
  theme_bw() + scale_fill_d3() +
  ggtitle("venous ECs") + theme(plot.title=element_text(hjust = 0.5))

plot1

ggsave("figure/ecsub/ecsubid/subid/vecfreq.png", plot1, width = 5, height = 3.5)

ggsave("figure/ecsub/ecsubid/subid/vecfreq.pdf", plot1, width = 5, height = 3.5)

statsub = subset(stat, subtype == "tip_like_ECs")

plot1 = ggplot(statsub, aes(x = factor(cancertype), y = prop, fill = tissuetype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Tip-like EC / Total EC (Percentage)", x = "Cancer type") +
  theme_bw() + scale_fill_d3() +
  ggtitle("tip-like ECs") + theme(plot.title=element_text(hjust = 0.5))

plot1

ggsave("figure/ecsub/ecsubprop/tipecfreq.png", plot1, width = 5, height = 3.5)

ggsave("figure/ecsub/ecsubprop/tipecfreq.pdf", plot1, width = 5, height = 3.5)

statsub = subset(stat, subtype == "lymphatic_ECs")

plot1 = ggplot(statsub, aes(x = factor(cancertype), y = prop, fill = tissuetype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Lymphatic EC / Total EC (Percentage)", x = "Cancer type") +
  theme_bw() + scale_fill_d3() +
  ggtitle("lymphatic ECs") + theme(plot.title=element_text(hjust = 0.5))

plot1

ggsave("figure/ecsub/ecsubid/subid/lecfreq.png", plot1, width = 5, height = 3.5)

ggsave("figure/ecsub/ecsubid/subid/lecfreq.pdf", plot1, width = 5, height = 3.5)

statsub = subset(stat, subtype == "capillary_ECs_1")

plot1 = ggplot(statsub, aes(x = factor(cancertype), y = prop, fill = tissuetype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Capillary EC-1 / Total EC (Percentage)", x = "Cancer type") +
  theme_bw() + scale_fill_d3() +
  ggtitle("capillary ECs-1") + theme(plot.title=element_text(hjust = 0.5))

plot1

ggsave("figure/ecsub/ecsubid/subid/cec1freq.png", plot1, width = 5, height = 3.5)

ggsave("figure/ecsub/ecsubid/subid/cec1freq.pdf", plot1, width = 5, height = 3.5)

statsub = subset(stat, subtype == "capillary_ECs_2")

plot1 = ggplot(statsub, aes(x = factor(cancertype), y = prop, fill = tissuetype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Capillary EC-2 / Total EC (Percentage)", x = "Cancer type") +
  theme_bw() + scale_fill_d3() +
  ggtitle("capillary ECs-2") + theme(plot.title=element_text(hjust = 0.5))

plot1

ggsave("figure/ecsub/ecsubid/subid/cec2freq.png", plot1, width = 5, height = 3.5)

ggsave("figure/ecsub/ecsubid/subid/cec2freq.pdf", plot1, width = 5, height = 3.5)

statsub = subset(stat, subtype == "arterial_ECs")

plot1 = ggplot(statsub, aes(x = factor(cancertype), y = prop, fill = tissuetype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Arterial EC / Total EC (Percentage)", x = "Cancer type") +
  theme_bw() + scale_fill_d3() +
  ggtitle("arterial ECs") + theme(plot.title=element_text(hjust = 0.5))

plot1

ggsave("figure/ecsub/ecsubid/subid/aecfreq.png", plot1, width = 5, height = 3.5)

ggsave("figure/ecsub/ecsubid/subid/aecfreq.pdf", plot1, width = 5, height = 3.5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "subtype",
                split.by = "tissuetype", label = F, repel = T, cols = pal) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

plot1

ggsave("figure/ecsub/ecsubid/subid/umapctc20nt.png", plot = plot1, width = 10.5, height = 5)

ggsave("figure/ecsub/ecsubid/subid/umapctc20nt.pdf", plot = plot1, width = 10.5, height = 5)

dir.create("figure/ecsub/ecatlas")

pal <- paletteer_d("ggsci::category20_d3")[c(1, 2)]

p1 = Seurat::VlnPlot(scRNA, features = "COL4A1", group.by = "cancertype",
                     split.by = "tissuetype", pt.size = 0, cols = pal)

p1

ggsave(p1, file = "figure/ecsub/ecatlas/col4a1c.png",
       width = 6, height = 3.5)

ggsave(p1, file = "figure/ecsub/ecatlas/col4a1c.pdf",
       width = 6, height = 3.5)

pal <- paletteer_d("ggsci::category20_d3")[c(5, 3, 4, 6, 2, 1)]

ECM = c("COL4A1", "COL4A2", "COL15A1", "COL18A1",
        "TIMP1", "TIMP2", "TIMP3")

p1 = VlnPlot(scRNA, features = ECM, cols = pal, group.by = "subtype",
             stacked = T, pt.size = 0,
             direction = "vertical", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  #不显示坐标刻度

p1

ggsave(p1, file = "figure/ecsub/ecatlas/ecmvio1.png",
       width = 6, height = 6)

ggsave(p1, file = "figure/ecsub/ecatlas/ecmvio1.pdf",
       width = 6, height = 6)

allgene = scRNA@assays$RNA@data@Dimnames[1]

str(allgene)

allgenes = allgene[[1]]

str(allgenes)

select_genes <- c("PECAM1", "CSPG4")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap",
                  blend = TRUE, label = F)

p2

ggsave("figure/ecsubtype/celltype_marker_pecam1_cspg4.png",
       p2, width = 15, height = 4)

ggsave("figure/ecsubtype/celltype_marker_pecam1_cspg4.pdf",
       p2, width = 15, height = 4)

select_genes <- c("PECAM1", "RGS5")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap",
                  blend = TRUE, label = F)

p2

ggsave("figure/ecsubtype/celltype_marker_pecam1_rgs5.png",
       p2, width = 15, height = 4)

ggsave("figure/ecsubtype/celltype_marker_pecam1_rgs5.pdf",
       p2, width = 15, height = 4)

select_genes <- c("PECAM1", "CSPG4", "RGS5")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap",
                  label = F, ncol = 3,
                  cols = c("lightgrey", "red2"))

ggsave("figure/ecsubtype/celltype_marker_pecam1_cspg4_rgs5.png",
       p2, width = 13, height = 4)

ggsave("figure/ecsubtype/celltype_marker_pecam1_cspg4_rgs5.pdf",
       p2, width = 13, height = 4)

p1 = VlnPlot(scRNA, features = "COL4A1", cols = pal,
             pt.size = 0, group.by = "subtype")

ggsave("figure/ecsubtype/ec_col4a1.png",
       p1, width = 8, height = 4.5)

ggsave("figure/ecsubtype/ec_col4a1.pdf",
       p1, width = 8, height = 4.5)

p1 = VlnPlot(scRNA, features = "COL4A2", cols = pal,
             pt.size = 0, group.by = "subtype")

ggsave("figure/ecsubtype/ec_col4a2.png",
       p1, width = 8, height = 4.5)

ggsave("figure/ecsubtype/ec_col4a2.pdf",
       p1, width = 8, height = 4.5)

select_genes <- c("PECAM1", "CSPG4", "VWF", "RGS5")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap",
                  label = F, ncol = 2,
                  cols = c("lightgrey", "red2"))

ggsave("figure/ecsubtype/celltype_marker_pecam1_vwf_cspg4_rgs5.png",
       p2, width = 9, height = 8)

ggsave("figure/ecsubtype/celltype_marker_pecam1_vwf_cspg4_rgs5.pdf",
       p2, width = 9, height = 8)

select_genes <- c("PECAM1", "VWF", "CSPG4", "RGS5")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap",
                  label = F, ncol = 4,
                  cols = c("lightgrey", "red2"))

ggsave("figure/ecsubtype/celltype_marker_pecam1_vwf_cspg4_rgs5_1.png",
       p2, width = 17, height = 4)

ggsave("figure/ecsubtype/celltype_marker_pecam1_vwf_cspg4_rgs5_1.pdf",
       p2, width = 17, height = 4)

allgene = scRNA@assays$RNA@data@Dimnames[1]

str(allgene)

allgenes = allgene[[1]]

str(allgenes)

HLA2 <- grep("^HLA-D", allgenes, value = T)

p1 = VlnPlot(scRNA, features = HLA2, cols = pal, group.by = "subtype",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
  theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
  theme(axis.ticks.x = element_line(size = 0.35))

p1

ggsave(p1, file = "figure/ecsubtype/ecsubhla2.png",
       width = 8, height = 5)

ggsave(p1, file = "figure/ecsubtype/ecsubhla2.pdf",
       width = 6, height = 3.5)

HLA2

hla2e = c("HLA-DRA", "HLA-DRB5", "HLA-DRB1", "HLA-DQA1",
          "HLA-DQB1", "HLA-DMB", "HLA-DMA",
          "HLA-DPA1", "HLA-DPB1")

p1 = VlnPlot(scRNA, features = hla2e, cols = pal, group.by = "subtype",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
  theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
  theme(axis.ticks.x = element_line(size = 0.35))

p1

ggsave(p1, file = "figure/ecsubtype/ecsubhla2e.png",
       width = 6.5, height = 5)

ggsave(p1, file = "figure/ecsubtype/ecsubhla2e.pdf",
       width = 6.5, height = 5)

mk <- AverageExpression(scRNA,
                        group.by = "subtype",
                        features = HLA2,
                        return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

p <- pheatmap(mfdf,
              color = colorRampPalette(c("#0080FF", "white", "#FF8000"))(20),
              main = "",
              scale = "row",
              angle_col = 90, # 设置显示角度
              #cellwidth = 20,cellheight = , # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 7.5, # 分别设置横向和纵向字体大小
              fontsize_col = 8)

ggsave("figure/ecsubtype/hla2hmp.png", p, width = 3.5, height = 5.5)

ggsave("figure/ecsubtype/hla2hmp.pdf", p, width = 3.5, height = 5.5)

pal <- paletteer_d("ggsci::category20_d3")[c(5, 3, 4, 6, 2, 1)]

p1 = VlnPlot(scRNA, features = "CD248", cols = pal, group.by = "subtype",
             stacked = T, pt.size = 0) + #横纵轴不标记任何东西
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  #不显示坐标刻度

p1

ggsave(p1, file = "figure/ecsub/ecatlas/ecmvio1.png",
       width = 6, height = 6)

ggsave(p1, file = "figure/ecsub/ecatlas/ecmvio1.pdf",
       width = 6, height = 6)

p1 = VlnPlot(scRNA, features = "FOLH1", cols = pal, group.by = "subtype",
             stacked = T, pt.size = 0) #不显示坐标刻度

p1

ggsave("figure/ecsubtype/psmaexp.png", p1, width = 6, height = 3.9)

ggsave("figure/ecsubtype/psmaexp.pdf", p1, width = 6, height = 3.9)

p1 = VlnPlot(scRNA, features = "CXCR4", cols = pal, group.by = "subtype",
             stacked = T, pt.size = 0) #不显示坐标刻度

p1

ggsave("figure/ecsubtype/cxcr4exp.png", p1, width = 6, height = 3.9)

ggsave("figure/ecsubtype/cxcr4exp.pdf", p1, width = 6, height = 3.9)

col.num <- length(levels(scRNA@active.ident))

violin <- Seurat::VlnPlot(scRNA,
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                                       "percent.HB", "percent.rb"),
                          cols = rainbow(col.num),
                          pt.size = 0.1, group.by = "subtype",
                          ncol = 5) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) 

ggsave("figure/ecsubtype/vlnplot_qc.png", plot = violin, width = 20, height = 6.5)

violin <- Seurat::VlnPlot(scRNA,
                          features = c("nFeature_RNA", "nCount_RNA"),
                          cols = rainbow(col.num),
                          pt.size = 0.01, group.by = "subtype",
                          ncol = 2)

ggsave("figure/ecsubtype/vlnplot_qca.png", plot = violin, width = 10, height = 6)

select_genes <- c("FOLH1")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap",
                  label = F, cols = c("lightgrey", "red2"))

ggsave("figure/ecsub/ecatlas/psmaumap.png", p2, width = 5.2, height = 5)

ggsave("figure/ecsub/ecatlas/psmaumap.pdf", p2, width = 5.2, height = 5)

pal <- paletteer_d("ggsci::category20_d3")[c(5, 3, 4, 6, 2, 1)]

surface = c("CA2", "INSR", "HTRA1", "MCAM", "CDH13", "PLVAP",
            "THY1", "TNFRSF4", "LAMA4", "PXDN", "FOLH1", "EDNRB")

p1 = VlnPlot(scRNA, features = surface, cols = pal, group.by = "subtype",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p1

ggsave(p1, file = "figure/ecsub/ecatlas/surfacevio1.png",
       width = 6.8, height = 3.8)

ggsave(p1, file = "figure/ecsub/ecatlas/surfacevio1.pdf",
       width = 6.8, height = 3.8)

violin <- Seurat::VlnPlot(scRNA,
                          features = "CXCR4",
                          cols = pal,
                          pt.size = 0, group.by = "subtype")

violin

ggsave("figure/ecsub/ecatlas/cxcr4vio.png", plot = violin, width = 8, height = 5)

ggsave("figure/ecsub/ecatlas/cxcr4vio.pdf", plot = violin, width = 8, height = 5)

select_genes <- c("CXCR4")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap",
                  label = F, cols = c("lightgrey", "red2"))

ggsave("figure/ecsub/ecatlas/CXCR4umap.png", p2, width = 5.2, height = 5)

ggsave("figure/ecsub/ecatlas/CXCR4umap.pdf", p2, width = 5.2, height = 5)

str(scRNA@meta.data)

pal <- paletteer_d("ggsci::category20_d3")[c(12:1)]

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "tissueorig", label = F, repel = T, cols = pal) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

plot1

ggsave("figure/ecsub/ecsubtype/umapcto20.png", plot = plot1, width = 6.6, height = 5)

ggsave("figure/ecsub/ecsubtype/umapcto20.pdf", plot = plot1, width = 6.6, height = 5)

surface = c("KCNE3", "NID2", "DLL4", "RAMP3", "EPLN")

p1 = VlnPlot(scRNA, features = surface, cols = pal, group.by = "subtype",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p1

ggsave(p1, file = "figure/ecsub/ecatlas/surfacevio1.png",
       width = 6.8, height = 3.8)

ggsave(p1, file = "figure/ecsub/ecatlas/surfacevio1.pdf",
       width = 6.8, height = 3.8)

select_genes <- c("GJA5", "ACKR1", "CA4", "PROX1",
                  "CXCR4", "ESM1", "ANGPT2")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap", ncol =4,
                  label = F, cols = c("lightgrey", "red2"))

ggsave("figure/ecsubtype/ecsubmk1.png", p2, width = 13, height = 6)

ggsave("figure/ecsubtype/ecsubmk1.pdf", p2, width = 13, height = 6)

dir.create("figure/response")

p1 = DotPlot(scRNA, features = "FOLH1", group.by = "tissuetype")

p1

ggsave("figure/response/psmatumorec.png", p1, width = 5, height = 5)

ggsave("figure/response/psmatumorec.pdf", p1, width = 5, height = 5)

p1 = DotPlot(scRNA, features = "FOLH1", group.by = "subtype")

p1

ggsave("figure/response/psmaecsub.png", p1, width = 6, height = 6)

ggsave("figure/response/psmaecsub.pdf", p1, width = 6, height = 6)

select_genes <- c("PECAM1", "CLDN5", "CDH5", "PLVAP",
                  "VWF", "FLT1", "KDR", "FLT4")

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "umap",
                  label = F, label.size = 3, ncol = 4,
                  cols = c("lightgrey", "red2"))

ggsave("figure/response/ecmkecatlas.png", p2, width = 17, height = 8)

ggsave("figure/response/ecmkecatlas.pdf", p2, width = 17, height = 8)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "subtype",
                split.by = "tissueorig", label = F, repel = T,
                ncol = 4, cols = pal) + #
  ggtitle(NULL) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
  theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
  guides(fill = guide_legend(title = ""))

plot1

ggsave("figure/response/umap-ecsub-tissueorig.png", plot = plot1,
       width = 12.5, height = 9)

ggsave("figure/response/umap-ecsub-tissueorig.pdf", plot = plot1,
       width = 12.5, height = 9)

#计数

str(scRNA@meta.data)

table(scRNA$tissueorig)

meta <- scRNA@meta.data[, c("tissueorig", "subtype")]

meta[1:5, 1:2]

stat = meta %>%
  group_by(tissueorig, subtype) %>%
  summarise(n = n())

dir.create("table")

write.csv(stat, "table/ecsub_tissueorig_count.csv", row.names = F)

table1 = table(scRNA$tissueorig)

write.csv(table1, "table/ec_tissueorig_count.csv", row.names = F)

top30 = read.csv(file = "figure/ecsubtype/tipec-top30.csv")

features = c(top30$gene)

mk <- AverageExpression(scRNA,
                        group.by = "subtype",
                        features = features,
                        return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

p <- pheatmap(mfdf,
              color = colorRampPalette(c("#0080FF", "white", "#FF8000"))(20),
              main = "",
              scale = "row",
              angle_col = 90, # 设置显示角度
              #cellwidth = 20,cellheight = , # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 5, # 分别设置横向和纵向字体大小
              fontsize_col = 6)

ggsave("figure/ecsubtype/tipec-top30.png", p, width = 3, height = 5)

ggsave("figure/ecsubtype/tipec-top30.pdf", p, width = 3, height = 5)


