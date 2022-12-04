
setwd("/home/zhangjiayu/project/EC/rawdata/2020NCpdac")

getwd()

library(Seurat)

library(tidyverse)

library(dplyr)

library(patchwork)

library(HGNChelper)

library(limma)

#library(hdf5r)

rm(list=ls())

scRNA = readRDS(file = "input/scRNAmerge.rds")

counts = scRNA@assays[["RNA"]]@counts

counts[1:6, 1:8]

counts = as.data.frame(counts)

counts[1:6, 1:8]

GeneUpdate = checkGeneSymbols(rownames(counts), unmapped.as.na = TRUE,
                              map = NULL, species = "human")

colnames(GeneUpdate)

table(GeneUpdate$Approved) #和HGNC不一致的基因数

sum(is.na(GeneUpdate$Suggested.Symbol)) #HGNC没有覆盖的基因数

GeneUpdate = checkGeneSymbols(rownames(counts), unmapped.as.na = FALSE,
                              map = NULL, species = "human")

sum(is.na(GeneUpdate$Suggested.Symbol))

Gene2up = select(GeneUpdate, rawgene = x, symbol = Suggested.Symbol)

Gene2up = column_to_rownames(Gene2up, var = "rawgene")

count2 = merge(Gene2up, counts, by = 0)

count2[1:6, 1:8]

count2 = count2[,-1]

count2[1:6, 1:8]

count2 = column_to_rownames(count2, var = "symbol")

count2[1:6, 1:8]

count2 = avereps(count2[, -1], ID = count2$symbol)

count2[1:6, 1:8]

class(count2)

meta = scRNA@meta.data

scRNA <- CreateSeuratObject(counts = count2, project = "2020_NC_PDAC", meta.data = meta,
                             min.cells = 3, min.features = 200)

table(scRNA$sampleid)

table(Idents(scRNA))

saveRDS(scRNA, file = "scRNAgeneup.rds")

# 质控

dir.create("QC")

scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

scRNA[["percent.ercc"]] <- PercentageFeatureSet(scRNA, pattern = "^ERCC-")

scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[SL]")

HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")

HB_m <- match(HB.genes, rownames(scRNA@assays$RNA))

HB.genes <- rownames(scRNA@assays$RNA)[HB_m]

HB.genes <- HB.genes[!is.na(HB.genes)]

scRNA[["percent.HB"]] <- PercentageFeatureSet(scRNA, features = HB.genes)

col.num <- length(levels(scRNA@active.ident))

violin <- VlnPlot(scRNA,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                               "percent.HB", "percent.rb"),
                  cols = rainbow(col.num),
                  pt.size = 0,
                  ncol = 5) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6)

plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.ercc")

plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")

plot4 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.rb")

pearplot <- CombinePlots(plots = list(plot1, plot2, plot3, plot4), nrow = 1,
                         legend = "none") 

ggsave("QC/pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5)

minGene = 500

maxGene = 6000

pctMT = 25

pctHB = 1

pctrb = 50

scRNA <- subset(scRNA, subset = nFeature_RNA > minGene
                & nFeature_RNA < maxGene
                & percent.mt < pctMT
                & percent.HB < pctHB
                & percent.rb < pctrb)

col.num <- length(levels(scRNA@active.ident))

violin <- VlnPlot(scRNA,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                               "percent.HB", "percent.rb"),
                  cols = rainbow(col.num),
                  pt.size = 0,
                  ncol = 5) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6)

saveRDS(scRNA, file = "/home/zhangjiayu/project/EC/input/2020pdacqc1.rds")


