
getwd()

library(tidyverse)

library(patchwork)

library(dplyr)

library(ggplot2)

library(ggthemes)

library(data.table)

library(scales)

rm(list = ls())

paad = read.csv(file = "tcga/paad/coxkm/sctcox.csv")

paad$sampletype = rep("PAAD", nrow(paad))

table(paad$sampletype)

ovc = read.csv(file = "tcga/ovc/coxkm/sctcox.csv")

ovc$sampletype = rep("OVC", nrow(ovc))

table(ovc$sampletype)

brca = read.csv(file = "tcga/brca/coxkm/sctcox.csv")

brca$sampletype = rep("BRCA",nrow(brca))

table(brca$sampletype)

lihc = read.csv(file = "tcga/lihc/coxkm/sctcox.csv")

lihc$sampletype = rep("LIHC", nrow(lihc))

table(lihc$sampletype)

coad = read.csv(file = "tcga/coad/coxkm/sctcox.csv")

coad$sampletype = rep("COAD", nrow(coad))

table(coad$sampletype)

read = read.csv(file = "tcga/read/coxkm/sctcox.csv")

read$sampletype = rep("READ", nrow(read))

table(read$sampletype)

prad = read.csv(file = "tcga/prad/coxkm/sctcox.csv")

prad$sampletype = rep("PRAD", nrow(prad))

table(prad$sampletype)

luad = read.csv(file = "tcga/luad/coxkm/sctcox.csv")

luad$sampletype = rep("LUAD", nrow(luad))

table(luad$sampletype)

lusc = read.csv(file = "tcga/lusc/coxkm/sctcox.csv")

lusc$sampletype = rep("LUSC", nrow(lusc))

table(lusc$sampletype)

blca = read.csv(file = "tcga/blca/coxkm/sctcox.csv")

blca$sampletype = rep("BLCA", nrow(blca))

table(blca$sampletype)

esca = read.csv(file = "tcga/esca/coxkm/sctcox.csv")

esca$sampletype = rep("ESCA", nrow(esca))

table(esca$sampletype)

gbm = read.csv(file = "tcga/gbm/coxkm/sctcox.csv")

gbm$sampletype = rep("GBM", nrow(gbm))

table(gbm$sampletype)

hnsc = read.csv(file = "tcga/hnsc/coxkm/sctcox.csv")

hnsc$sampletype = rep("HNSC", nrow(hnsc))

table(hnsc$sampletype)

thca = read.csv(file = "tcga/thca/coxkm/sctcox.csv")

thca$sampletype = rep("THCA", nrow(thca))

table(thca$sampletype)

stad = read.csv(file = "tcga/stad/coxkm/sctcox.csv")

stad$sampletype = rep("STAD", nrow(stad))

table(stad$sampletype)

cesc = read.csv(file = "tcga/cesc/coxkm/sctcox.csv")

cesc$sampletype = rep("CESC", nrow(cesc))

table(cesc$sampletype)

ucec = read.csv(file = "tcga/ucec/coxkm/sctcox.csv")

ucec$sampletype = rep("UCEC", nrow(ucec))

table(ucec$sampletype)

kirc = read.csv(file = "tcga/kirc/coxkm/sctcox.csv")

kirc$sampletype = rep("KIRC", nrow(kirc))

table(kirc$sampletype)

kich = read.csv(file = "tcga/kich/coxkm/sctcox.csv")

kich$sampletype = rep("KICH", nrow(kich))

table(kich$sampletype)

kirp = read.csv(file = "tcga/kirp/coxkm/sctcox.csv")

kirp$sampletype = rep("KIRP", nrow(kirp))

table(kirp$sampletype)

colnames(paad)

colnames(kirp)

l = list(blca, brca, cesc, coad, esca,
         gbm, hnsc, kich, kirc, kirp,
         lihc, luad, lusc, ovc, paad,
         prad, read, stad, thca, ucec)

panca = rbindlist(l, use.names=TRUE)

table(panca$sampletype)

dir.create("cox")

saveRDS(panca, file = "cox/pancact.rds")

#panca = readRDS(file = "cox/pancact.rds")

colnames(panca)

table(panca$sample)

panca$pstar <- ifelse(panca$pvalue < 0.05,
                      ifelse(panca$pvalue < 0.01,
                             ifelse(panca$pvalue < 0.001,"***", "**"), "*"), "")

write.csv(panca, file = "cox/pancact.csv")

panca = read.csv(file = "cox/pancacta.csv")

saveRDS(panca, file = "cox/pancact.rds")

panca = readRDS(file = "cox/pancact.rds")

p1 = ggplot(panca, aes(sampletype, id)) +
        geom_tile(aes(fill = HR), colour = "white", size = 0.1) +
        scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red", midpoint = 1) +
        geom_text(aes(label = pstar), col ="black", size = 5) +
        theme_minimal() +
        theme(axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 8))

p1

ggsave(p1, file = "cox/cox.png", width = 6.5, height = 5.5)

ggsave(p1, file = "cox/cox.pdf", width = 6.5, height = 5.5)

#panca = readRDS(file = "cox/pancact.rds")

table(panca$id)

pancaec = filter(panca, id == "arterialECs"|id == "capillaryECs1"|
                        id == "capillaryECs2"|id == "venousECs"|
                        id == "tipECs"|id == "lymphaticECs")

p1 = ggplot(pancaec, aes(sampletype, id)) +
  geom_tile(aes(fill = HR), colour = "white", size = 0.1) +
  scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red", midpoint = 1) +
  geom_text(aes(label = pstar), col ="black", size = 5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))

p1

#dir.create("tcga/coxec")

ggsave(p1, file = "coxec/cox.png", width = 6.5, height = 2.5)

ggsave(p1, file = "coxec/cox.pdf", width = 6.5, height = 2.5)



