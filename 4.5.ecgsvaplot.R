
getwd()

library(Seurat)

library(tidyverse)

library(GSVA)

library(msigdbr)

library(ggsci)

library(paletteer)

library(patchwork)

library(pheatmap) # 加载包

rm(list = ls())

dir.create("figure/ecsub/gsva")

###

sch = readRDS(file = "process/5.ec/gsva/hallmark/sch.rds")

table(sch@meta.data$subtype)

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
              fontsize_row = 4.5, # 分别设置横向和纵向字体大小
              fontsize_col = 6)

ggsave("figure/ecsub/gsva/hallmark/hallmark.png", p, width = 4.2, height = 6)

ggsave("figure/ecsub/gsva/hallmark/hallmark.pdf", p, width = 4.2, height = 6)

#GO-BP

scg = readRDS(file = "process/5.ec/gsva/go-bp/scg.rds")

table(scg@meta.data$subtype)

str(scg)

table(Idents(scg))

goset = read.csv(file = "process/5.ec/gsva/go-bp/top20mka.csv")

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
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              treeheight_col = 20, # 分别设置横、纵向聚类树高
              treeheight_row = 20,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 3, # 分别设置横向和纵向字体大小
              fontsize_col = 8)

ggsave("figure/ecsub/gsva/go/gotop20.png", p, width = 6, height = 8)

ggsave("figure/ecsub/gsva/go/gotop20.pdf", p, width = 6, height = 8)

goset = read.csv(file = "process/5.ec/gsva/go-bp/top10mk.csv")

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
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              treeheight_col = 20, # 分别设置横、纵向聚类树高
              treeheight_row = 20,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 5, # 分别设置横向和纵向字体大小
              fontsize_col = 8)

ggsave("figure/ecsub/gsva/go/gotop10.png", p, width = 8, height = 8)

ggsave("figure/ecsub/gsva/go/gotop10.pdf", p, width = 8, height = 8)

goset = read.csv(file = "process/5.ec/gsva/go-bp/top5select.csv")

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
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              treeheight_col = 20, # 分别设置横、纵向聚类树高
              treeheight_row = 20,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 6, # 分别设置横向和纵向字体大小
              fontsize_col = 8)

ggsave("figure/ecsub/gsva/go/goselect5.png", p, width = 6.8, height = 6)

ggsave("figure/ecsub/gsva/go/goselect5.pdf", p, width = 6.8, height = 6)



