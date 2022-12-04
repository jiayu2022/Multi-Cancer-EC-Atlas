
setwd("/home/zhangjiayu/project/EC")

getwd()

library(Seurat)

library(tidyverse)

library(GSVA)

library(msigdbr)

library(ggsci)

library(paletteer)

library(patchwork)

library(pheatmap) # 加载包

dir.create("figure/ecsub")

dir.create("figure/ecsub/gsva")

dir.create("figure/ecsub/gsva/tmp")

rm(list = ls())

###

scRNA = readRDS(file = "keyobject/ecsubid.rds")

table(Idents(scRNA))

library(future)

plan()

options(future.globals.maxSize = 10 * 1024^3)

plan("multicore", workers = 56)

#hallmarker

genesets <- msigdbr(species = "Homo sapiens", category = "H") %>% select("gs_name","gene_symbol") %>% as.data.frame()

genesets <- split(genesets$gene_symbol, genesets$gs_name)

str(genesets)

table(Idents(scRNA))

expr <- as.matrix(scRNA@assays$RNA@counts)

gsvahm = gsva(expr, genesets, method = "ssgsea", parallel.sz = 56)

saveRDS(gsvahm, file = "figure/ecsub/gsva/tmp/gsvahm.rds")

#gsvahm = readRDS(file = "figure/ecsub/gsva/tmp/gsvahm.rds")

str(gsvahm)

gsvahm[1:5, 1:10]

str(scRNA@meta.data)

table(scRNA$subtype)

meta <- scRNA@meta.data[, c("seurat_clusters", "subtype")]

meta[1:5, 1:2]

sch = CreateSeuratObject(gsvahm, project = "Hallmark", meta.data = meta)

table(Idents(sch))

Idents(sch) = sch@meta.data$subtype

table(Idents(sch))

saveRDS(sch, file = "figure/ecsub/gsva/tmp/sch.rds")

###

#sch = readRDS(file = "figure/ecsub/gsva/tmp/sch.rds")

str(sch)

str(sch@assays$RNA@data@Dimnames[1])

hm = sch@assays$RNA@data@Dimnames[1]

str(hm)

features = c(hm[[1]])

features

mk <- AverageExpression(sch,
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
              cluster_cols = FALSE, # 去掉横向、纵向聚类
              cluster_rows = TRUE,
              treeheight_col = 20, # 分别设置横、纵向聚类树高
              treeheight_row = 20,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 3.5, # 分别设置横向和纵向字体大小
              fontsize_col = 6)

ggsave("figure/ecsub/gsva/hallmark.png", p, width = 4.2, height = 6)

ggsave("figure/ecsub/gsva/hallmark.pdf", p, width = 4.2, height = 6)

table(scRNA$subtype)

#GO-BP

m_df = msigdbr(species = "Homo sapiens")

m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)

genesets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>% select("gs_name","gene_symbol") %>% as.data.frame()

genesets <- split(genesets$gene_symbol, genesets$gs_name)

str(genesets)

table(Idents(scRNA))

plan()

expr <- as.matrix(scRNA@assays$RNA@counts)

gsvago = gsva(expr, genesets, method = "ssgsea", parallel.sz = 56)

saveRDS(gsvago, file = "figure/ecsub/gsva/tmp/gsvago.rds")
 
gsvago = readRDS(file = "figure/ecsub/gsva/tmp/gsvago.rds")

str(gsvago)

gsvago[1:5, 1:10]

meta <- scRNA@meta.data[, c("seurat_clusters", "subtype")]

meta[1:5, 1:2]

scg = CreateSeuratObject(gsvago, project = "GO_BP", meta.data = meta)

table(Idents(scg))

Idents(scg) = scg@meta.data$subtype

table(Idents(scg))

saveRDS(scg, file = "figure/ecsub/gsva/tmp/scg.rds")

#scg = readRDS(file = "figure/ecsub/gsva/tmp/scg.rds")

plan()

table(Idents(scg))

diff.wilcox = FindAllMarkers(scg, logfc.threshold = 0.02)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

dir.create("figure/ecsub/gsva/go")

saveRDS(all.markers, file = "figure/ecsub/gsva/go/markers.rds")

#all.markers = readRDS(file = "figure/ecsub/gsva/go/markers.rds")

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(top10, "figure/ecsub/gsva/go/top10_degs.csv", row.names = F)

top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(top50, "figure/ecsub/gsva/go/top50_degs_wilcox.csv", row.names = F)

top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top20, "figure/ecsub/gsva/go/top20mk.csv", row.names = F)

goset = read.csv(file = "figure/ecsub/gsva/go/top10_degs.csv")

features = c(goset$gene)

features

mk <- AverageExpression(scg,
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
              cluster_cols = FALSE, # 去掉横向、纵向聚类
              cluster_rows = FALSE,
              treeheight_col = 20, # 分别设置横、纵向聚类树高
              treeheight_row = 20,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 3.5, # 分别设置横向和纵向字体大小
              fontsize_col = 6)

ggsave("figure/ecsub/gsva/go/gotop10mka.png", p, width = 6, height = 6)

ggsave("figure/ecsub/gsva/go/gotop10mka.pdf", p, width = 6, height = 6)
