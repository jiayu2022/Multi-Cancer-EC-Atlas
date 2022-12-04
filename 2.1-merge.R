
setwd("/home/zhangjiayu/project/EC")

getwd()

library(Seurat)

library(tidyverse)

library(patchwork)

library(dplyr)

library(future)

plan()

rm(list = ls())

scRNA1 = readRDS(file = "input/2020crcrc1.rds")

scRNA2 = readRDS(file = "input/2020crlc1.rds")

scRNA3 = readRDS(file = "input/2020crovc1.rds")

scRNA4 = readRDS(file = "input/2020pdacqc1.rds")

scRNA5 = readRDS(file = "input/2021gcqc1.rds")

scRNA6 = readRDS(file = "input/2021pnasrcc1.rds")

scRNA <- merge(scRNA1, y = c(scRNA2, scRNA3, scRNA4, scRNA5,
                             scRNA6))

Idents(scRNA)

str(scRNA@meta.data)

table(scRNA$tissueorig)

Idents(scRNA) = scRNA$tissuetype

table(Idents(scRNA))

str(scRNA@meta.data)

col.num <- length(levels(scRNA@active.ident))

violin <- VlnPlot(scRNA,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                               "percent.HB", "percent.rb"),
                  cols = rainbow(col.num),
                  pt.size = 0,
                  ncol = 5) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) 

dir.create("process")

dir.create("process/QC")

ggsave("process/QC/vlnplot_qc.pdf", plot = violin, width = 12, height = 5)

plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.ercc")

plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")

plot4 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.rb")

pearplot <- CombinePlots(plots = list(plot1, plot2, plot3, plot4), nrow = 1, legend = "none") 

ggsave("process/QC/vlnplot_qc1.pdf", plot = pearplot, width = 16, height = 5)

str(scRNA@meta.data)

#细胞周期评分

s.genes <- cc.genes$s.genes

g2m.genes <- cc.genes$g2m.genes

scRNA <- CellCycleScoring(scRNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#meta.data添加cancertype

table(scRNA$tissuetype)

table(scRNA$tissueorig)

table(scRNA$sampleid)

meta = scRNA@meta.data

cancertype = select(meta, tissueorig)

cancer = separate(cancertype, col = tissueorig, into = c("cancertype", "tissuetype"), sep = "_")

cancer = select(cancer, cancertype)

table(cancer$cancertype)

scRNA <- AddMetaData(scRNA, cancer)

table(scRNA$cancertype)

saveRDS(scRNA, file = "process/scRNAmerge.rds")


