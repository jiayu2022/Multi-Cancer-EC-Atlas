
getwd()

library(Seurat)

library(tidyverse)

library(patchwork)

library(dplyr)

rm(list = ls())

# scRNA 取子集 ----

scRNA = readRDS(file = "keyobject/scRNActid.rds")

table(scRNA$celltype)

table(Idents(scRNA))

dir.create("process/6.cibersortx")

dir.create("process/6.cibersortx/sigmatrix")

dir.create("process/6.cibersortx/sigmatrix/input")

scsub = subset(scRNA, downsample = 500)

table(Idents(scsub))

dir.create("process/6.cibersortx/sigmatrix/tmp")

saveRDS(scsub, file = "process/6.cibersortx/sigmatrix/tmp/scsub500.rds")

#scsub = readRDS(file = "process/6.cibersortx/sigmatrix/tmp/scsub500.rds")

Cells.sub <- subset(scsub@meta.data, celltype != "Endothelial_cells")

sc <- subset(scsub, cells = row.names(Cells.sub))

table(Idents(sc))

saveRDS(sc, file = "process/6.cibersortx/sigmatrix/tmp/scsub500withoutec.rds")

#sc = readRDS(file = "process/6.cibersortx/sigmatrix/tmp/scsub500withoutec.rds")

ec = readRDS(file = "keyobject/ecsubid.rds")

table(Idents(ec))

table(ec$subtype)

Idents(ec) = ec$subtype

table(Idents(ec))

ecsub = subset(ec, downsample = 200)

table(Idents(ecsub))

saveRDS(ecsub, file = "process/6.cibersortx/sigmatrix/tmp/ecsub200.rds")

sccc = merge(sc, ecsub)

table(Idents(sccc))

saveRDS(sccc, file = "process/6.cibersortx/sigmatrix/input/scRNAsigmtx.rds")

scRNA2 = sccc

levels(scRNA2)

table(scRNA2$celltype)

table(Idents(scRNA2))

scRNA2$celltype = Idents(scRNA2)

table(scRNA2$celltype)

table(scRNA2$tissueorig)

saveRDS(scRNA2, file = "process/6.cibersortx/sigmatrix/input/pancancermtxscct.rds")

# scRNA 导出 data ----

rm(list = ls())

scRNA2 = readRDS(file = "process/6.cibersortx/sigmatrix/input/pancancermtxscct.rds")

sccounts = scRNA2@assays[["RNA"]]@counts

max(sccounts)

#counts <- GetAssayData(scRNA2, assay = "RNA", slot = "counts")

sccounts[1:5, 1:8]

exprdf = as.data.frame(sccounts)

exprdf[1:5, 1:8]

# 导出 ----

md = scRNA2@meta.data

md[1:6,1:10]

table(md$celltype)

ct = dplyr::select(md, celltype)

head(ct)

ct = t(ct)

ct = as.data.frame(ct)

ct[1:1, 1:5]

exprdf[1:5, 1:10]

rownames(ct) = "GeneSymbol"

exprdf = rbind(ct, exprdf)

exprdf[1:5, 1:8]

colnames(exprdf) = NULL

exprdf[1:5, 1:8]

write.table(exprdf,  # 要导出的数据
            file = "process/6.cibersortx/sigmatrix/input/pancasccounts.txt", # 指定导出数据的文件名称和格式
            sep = "\t", # 字段分隔符
            row.names = TRUE, # 逻辑词，是否将行名一起导出
            col.names = FALSE) # 逻辑词，是否将列名一起导出

refsample = read.table(file = "process/6.cibersortx/sigmatrix/input/pancasccounts.txt", # 指定导出数据的文件名称和格式
                       sep = "\t", head = TRUE, row.names = 1)

refsample[1:5, 1:8]

max(refsample)

# end ----



