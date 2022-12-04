
setwd("/home/zhangjiayu/project/EC")

getwd()

library(Seurat)

library(CellChat)

library(tidyverse)

library(patchwork)

library(dplyr)

library(ggalluvial)

library(future)

rm(list = ls())

options(stringsAsFactors = FALSE)

plan()

###

dir.create("process/7.tipeccellchat")

dir.create("process/7.tipeccellchat/input")

#ECsub

scRNA = readRDS(file = "keyobject/ecsubid.rds")

table(Idents(scRNA))

table(scRNA$subtype)

table(scRNA$tissueorig)

#crcn

scrnan = readRDS(file = "process/7.cellchat/input/crcn.rds")

table(scrnan$celltype)

scrnarmec = subset(scrnan, celltype != "Endothelial_cells")

table(scrnarmec$celltype)

ecsub = subset(scRNA, tissueorig == "CRC_Normal")

table(ecsub$subtype)

table(ecsub$tissueorig)

table(ecsub$celltype)

scrnarmec$subtype = scrnarmec$celltype

table(scrnarmec$subtype)

scrnansub = merge(scrnarmec, ecsub)

table(scrnansub$subtype)

saveRDS(scrnansub, file = "process/7.tipeccellchat/input/crcnsub.rds")

#crct

scrnan = readRDS(file = "process/7.cellchat/input/crct.rds")

table(scrnan$celltype)

scrnarmec = subset(scrnan, celltype != "Endothelial_cells")

table(scrnarmec$celltype)

ecsub = subset(scRNA, tissueorig == "CRC_Tumor")

table(ecsub$subtype)

table(ecsub$tissueorig)

table(ecsub$celltype)

scrnarmec$subtype = scrnarmec$celltype

table(scrnarmec$subtype)

scrnansub = merge(scrnarmec, ecsub)

table(scrnansub$subtype)

saveRDS(scrnansub, file = "process/7.tipeccellchat/input/crctsub.rds")

#gcn

scrnan = readRDS(file = "process/7.cellchat/input/gcn.rds")

table(scrnan$celltype)

scrnarmec = subset(scrnan, celltype != "Endothelial_cells")

table(scrnarmec$celltype)

ecsub = subset(scRNA, tissueorig == "GC_Normal")

table(ecsub$subtype)

table(ecsub$tissueorig)

table(ecsub$celltype)

scrnarmec$subtype = scrnarmec$celltype

table(scrnarmec$subtype)

scrnansub = merge(scrnarmec, ecsub)

table(scrnansub$subtype)

saveRDS(scrnansub, file = "process/7.tipeccellchat/input/gcnsub.rds")

#gct

scrnan = readRDS(file = "process/7.cellchat/input/gct.rds")

table(scrnan$celltype)

scrnarmec = subset(scrnan, celltype != "Endothelial_cells")

table(scrnarmec$celltype)

ecsub = subset(scRNA, tissueorig == "GC_Tumor")

table(ecsub$subtype)

table(ecsub$tissueorig)

table(ecsub$celltype)

scrnarmec$subtype = scrnarmec$celltype

table(scrnarmec$subtype)

scrnansub = merge(scrnarmec, ecsub)

table(scrnansub$subtype)

saveRDS(scrnansub, file = "process/7.tipeccellchat/input/gctsub.rds")

#lcn

scrnan = readRDS(file = "process/7.cellchat/input/lcn.rds")

table(scrnan$celltype)

scrnarmec = subset(scrnan, celltype != "Endothelial_cells")

table(scrnarmec$celltype)

ecsub = subset(scRNA, tissueorig == "LC_Normal")

table(ecsub$subtype)

table(ecsub$tissueorig)

table(ecsub$celltype)

scrnarmec$subtype = scrnarmec$celltype

table(scrnarmec$subtype)

scrnansub = merge(scrnarmec, ecsub)

table(scrnansub$subtype)

saveRDS(scrnansub, file = "process/7.tipeccellchat/input/lcnsub.rds")

#lct

scrnan = readRDS(file = "process/7.cellchat/input/lct.rds")

table(scrnan$celltype)

scrnarmec = subset(scrnan, celltype != "Endothelial_cells")

table(scrnarmec$celltype)

ecsub = subset(scRNA, tissueorig == "LC_Tumor")

table(ecsub$subtype)

table(ecsub$tissueorig)

table(ecsub$celltype)

scrnarmec$subtype = scrnarmec$celltype

table(scrnarmec$subtype)

scrnansub = merge(scrnarmec, ecsub)

table(scrnansub$subtype)

saveRDS(scrnansub, file = "process/7.tipeccellchat/input/lctsub.rds")

#ovcn

scrnan = readRDS(file = "process/7.cellchat/input/ovcn.rds")

table(scrnan$celltype)

scrnarmec = subset(scrnan, celltype != "Endothelial_cells")

table(scrnarmec$celltype)

ecsub = subset(scRNA, tissueorig == "OVC_Normal")

table(ecsub$subtype)

table(ecsub$tissueorig)

table(ecsub$celltype)

scrnarmec$subtype = scrnarmec$celltype

table(scrnarmec$subtype)

scrnansub = merge(scrnarmec, ecsub)

table(scrnansub$subtype)

saveRDS(scrnansub, file = "process/7.tipeccellchat/input/ovcnsub.rds")

#ovct

scrnan = readRDS(file = "process/7.cellchat/input/ovct.rds")

table(scrnan$celltype)

scrnarmec = subset(scrnan, celltype != "Endothelial_cells")

table(scrnarmec$celltype)

ecsub = subset(scRNA, tissueorig == "OVC_Tumor")

table(ecsub$subtype)

table(ecsub$tissueorig)

table(ecsub$celltype)

scrnarmec$subtype = scrnarmec$celltype

table(scrnarmec$subtype)

scrnansub = merge(scrnarmec, ecsub)

table(scrnansub$subtype)

saveRDS(scrnansub, file = "process/7.tipeccellchat/input/ovctsub.rds")

#pdacn

scrnan = readRDS(file = "process/7.cellchat/input/pdacn.rds")

table(scrnan$celltype)

scrnarmec = subset(scrnan, celltype != "Endothelial_cells")

table(scrnarmec$celltype)

ecsub = subset(scRNA, tissueorig == "PDAC_Normal")

table(ecsub$subtype)

table(ecsub$tissueorig)

table(ecsub$celltype)

scrnarmec$subtype = scrnarmec$celltype

table(scrnarmec$subtype)

scrnansub = merge(scrnarmec, ecsub)

table(scrnansub$subtype)

saveRDS(scrnansub, file = "process/7.tipeccellchat/input/pdacnsub.rds")

#pdact

scrnan = readRDS(file = "process/7.cellchat/input/pdact.rds")

table(scrnan$celltype)

scrnarmec = subset(scrnan, celltype != "Endothelial_cells")

table(scrnarmec$celltype)

ecsub = subset(scRNA, tissueorig == "PDAC_Tumor")

table(ecsub$subtype)

table(ecsub$tissueorig)

table(ecsub$celltype)

scrnarmec$subtype = scrnarmec$celltype

table(scrnarmec$subtype)

scrnansub = merge(scrnarmec, ecsub)

table(scrnansub$subtype)

saveRDS(scrnansub, file = "process/7.tipeccellchat/input/pdactsub.rds")

#rccn

scrnan = readRDS(file = "process/7.cellchat/input/rccn.rds")

table(scrnan$celltype)

scrnarmec = subset(scrnan, celltype != "Endothelial_cells")

table(scrnarmec$celltype)

ecsub = subset(scRNA, tissueorig == "RCC_Normal")

table(ecsub$subtype)

table(ecsub$tissueorig)

table(ecsub$celltype)

scrnarmec$subtype = scrnarmec$celltype

table(scrnarmec$subtype)

scrnansub = merge(scrnarmec, ecsub)

table(scrnansub$subtype)

saveRDS(scrnansub, file = "process/7.tipeccellchat/input/rccnsub.rds")

#rcct

scrnan = readRDS(file = "process/7.cellchat/input/rcct.rds")

table(scrnan$celltype)

scrnarmec = subset(scrnan, celltype != "Endothelial_cells")

table(scrnarmec$celltype)

ecsub = subset(scRNA, tissueorig == "RCC_Tumor")

table(ecsub$subtype)

table(ecsub$tissueorig)

table(ecsub$celltype)

scrnarmec$subtype = scrnarmec$celltype

table(scrnarmec$subtype)

scrnansub = merge(scrnarmec, ecsub)

table(scrnansub$subtype)

saveRDS(scrnansub, file = "process/7.tipeccellchat/input/rcctsub.rds")

