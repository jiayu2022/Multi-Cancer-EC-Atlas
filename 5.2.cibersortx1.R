
getwd()

library(tidyverse)

library(patchwork)

library(dplyr)

library(limma)

library(ggplot2)

library(ggthemes)

library(ggpubr)

library(ggsci)

library(data.table)

rm(list = ls())

paad = read.csv(file = "tcga/paad/out/cellecsub.csv")

paad$sampletype = rep("PAAD", nrow(paad))

table(paad$sampletype)

ovc = read.csv(file = "tcga/ovc/out/cellecsub.csv")

ovc$sampletype = rep("OVC", nrow(ovc))

table(ovc$sampletype)

brca = read.csv(file = "tcga/brca/out/cellecsub.csv")

brca$sampletype = rep("BRCA",nrow(brca))

table(brca$sampletype)

lihc = read.csv(file = "tcga/lihc/out/cellecsub.csv")

lihc$sampletype = rep("LIHC", nrow(lihc))

table(lihc$sampletype)

coad = read.csv(file = "tcga/coad/out/cellecsub.csv")

coad$sampletype = rep("COAD", nrow(coad))

table(coad$sampletype)

read = read.csv(file = "tcga/read/out/cellecsub.csv")

read$sampletype = rep("READ", nrow(read))

table(read$sampletype)

prad = read.csv(file = "tcga/prad/out/cellecsub.csv")

prad$sampletype = rep("PRAD", nrow(prad))

table(prad$sampletype)

luad = read.csv(file = "tcga/luad/out/cellecsub.csv")

luad$sampletype = rep("LUAD", nrow(luad))

table(luad$sampletype)

lusc = read.csv(file = "tcga/lusc/out/cellecsub.csv")

lusc$sampletype = rep("LUSC", nrow(lusc))

table(lusc$sampletype)

blca = read.csv(file = "tcga/blca/out/cellecsub.csv")

blca$sampletype = rep("BLCA", nrow(blca))

table(blca$sampletype)

esca = read.csv(file = "tcga/esca/out/cellecsub.csv")

esca$sampletype = rep("ESCA", nrow(esca))

table(esca$sampletype)

gbm = read.csv(file = "tcga/gbm/out/cellecsub.csv")

gbm$sampletype = rep("GBM", nrow(gbm))

table(gbm$sampletype)

hnsc = read.csv(file = "tcga/hnsc/out/cellecsub.csv")

hnsc$sampletype = rep("HNSC", nrow(hnsc))

table(hnsc$sampletype)

thca = read.csv(file = "tcga/thca/out/cellecsub.csv")

thca$sampletype = rep("THCA", nrow(thca))

table(thca$sampletype)

stad = read.csv(file = "tcga/stad/out/cellecsub.csv")

stad$sampletype = rep("STAD", nrow(stad))

table(stad$sampletype)

cesc = read.csv(file = "tcga/cesc/out/cellecsub.csv")

cesc$sampletype = rep("CESC", nrow(cesc))

table(cesc$sampletype)

ucec = read.csv(file = "tcga/ucec/out/cellecsub.csv")

ucec$sampletype = rep("UCEC", nrow(ucec))

table(ucec$sampletype)

kirc = read.csv(file = "tcga/kirc/out/cellecsub.csv")

kirc$sampletype = rep("KIRC", nrow(kirc))

table(kirc$sampletype)

kich = read.csv(file = "tcga/kich/out/cellecsub.csv")

kich$sampletype = rep("KICH", nrow(kich))

table(kich$sampletype)

kirp = read.csv(file = "tcga/kirp/out/cellecsub.csv")

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

dir.create("input")

dir.create("out")

saveRDS(panca, file = "input/pancact.rds")

#panca = readRDS(file = "input/pancact.rds")

table(panca$sample)

colnames(panca)

plot.info_all <- panca[, c(26, 24, 25)]

p1 = ggboxplot(plot.info_all,
               x = "sampletype",
               y = "EC",
               color = "sample",
               fill = "white", #只需要修改这里，讲fill = "CellType"改为fill = "group"
               outlier.size = 0.5,
               xlab = "",
               ylab = "Fraction of  ECs") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5),
    axis.title.y = element_text(margin = margin(0, 0.35, 0, 0,'cm'))) +
  scale_color_d3() +
  stat_compare_means(aes(group = sample), label = "p.signif")

p1

ggsave("tcga/out/ecsuball.png", plot = p1, width = 8, height = 4)

ggsave("tcga/out/ecsuball.pdf", plot = p1, width = 8, height = 4)

colnames(panca)

plot.info_all <- panca[, c(26, 24, 5)]

p1 = ggboxplot(plot.info_all,
               x = "sampletype",
               y = "Pericytes",
               color = "sample",
               fill = "white", #只需要修改这里，讲fill = "CellType"改为fill = "group"
               outlier.size = 0.5,
               xlab = "",
               ylab = "Fraction of  Pericytes") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5),
    axis.title.y = element_text(margin = margin(0, 0.35, 0, 0,'cm'))) +
  scale_color_d3() +
  stat_compare_means(aes(group = sample), label = "p.signif")

p1

ggsave("tcga/out/perisuball.png", plot = p1, width = 8, height = 4)

ggsave("tcga/out/perisuball.pdf", plot = p1, width = 8, height = 4)

colnames(panca)

plot.info_all3 <- panca[, c(26, 24, 18)]

p3 = ggboxplot(plot.info_all3,
               x = "sampletype",
               y = "tip_like_ECs",
               color = "sample",
               fill = "white", #只需要修改这里，讲fill = "CellType"改为fill = "group"
               outlier.size = 0.5,
               xlab = "",
               ylab = "Fraction of  tip_ECs") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5),
    axis.title.y = element_text(margin = margin(0, 0.35, 0, 0,'cm'))) +
  scale_color_d3() +
  stat_compare_means(aes(group = sample), label = "p.signif")

p3

ggsave("tcga/out/tipec.png", plot = p3, width = 8, height = 4)

ggsave("tcga/out/tipec.pdf", plot = p3, width = 8, height = 4)

colnames(panca)

plot.info_all3 <- panca[, c(26, 24, 15)]

p3 = ggboxplot(plot.info_all3,
               x = "sampletype",
               y = "lymphatic_ECs",
               color = "sample",
               fill = "white", #只需要修改这里，讲fill = "CellType"改为fill = "group"
               outlier.size = 0.5,
               xlab = "",
               ylab = "Fraction of  lymphatic_ECs") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5),
    axis.title.y = element_text(margin = margin(0, 0.35, 0, 0,'cm'))) +
  scale_color_d3() +
  stat_compare_means(aes(group = sample), label = "p.signif")

p3

ggsave("tcga/out/lymphatic_ECs.png", plot = p3, width = 8, height = 4)

ggsave("tcga/out/lymphatic_ECs.pdf", plot = p3, width = 8, height = 4)

colnames(panca)

plot.info_all3 <- panca[, c(26, 24, 19)]

p3 = ggboxplot(plot.info_all3,
               x = "sampletype",
               y = "capillary_ECs_1",
               color = "sample",
               fill = "white", #只需要修改这里，讲fill = "CellType"改为fill = "group"
               outlier.size = 0.5,
               xlab = "",
               ylab = "Fraction of  capillary_ECs_1") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5),
    axis.title.y = element_text(margin = margin(0, 0.35, 0, 0,'cm'))) +
  scale_color_d3() +
  stat_compare_means(aes(group = sample), label = "p.signif")

p3

ggsave("tcga/out/capillary_ECs1.png", plot = p3, width = 8, height = 4)

ggsave("tcga/out/capillary_ECs1.pdf", plot = p3, width = 8, height = 4)

colnames(panca)

plot.info_all3 <- panca[, c(26, 24, 20)]

p3 = ggboxplot(plot.info_all3,
               x = "sampletype",
               y = "capillary_ECs_2",
               color = "sample",
               fill = "white", #只需要修改这里，讲fill = "CellType"改为fill = "group"
               outlier.size = 0.5,
               xlab = "",
               ylab = "Fraction of  capillary_ECs_2") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5),
    axis.title.y = element_text(margin = margin(0, 0.35, 0, 0,'cm'))) +
  scale_color_d3() +
  stat_compare_means(aes(group = sample), label = "p.signif")

p3

ggsave("tcga/out/capillary_ECs2.png", plot = p3, width = 8, height = 4)

ggsave("tcga/out/capillary_ECs2.pdf", plot = p3, width = 8, height = 4)

colnames(panca)

plot.info_all3 <- panca[, c(26, 24, 16)]

p3 = ggboxplot(plot.info_all3,
               x = "sampletype",
               y = "arterial_ECs",
               color = "sample",
               fill = "white", #只需要修改这里，讲fill = "CellType"改为fill = "group"
               outlier.size = 0.5,
               xlab = "",
               ylab = "Fraction of  arterial_ECs") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5),
    axis.title.y = element_text(margin = margin(0, 0.35, 0, 0,'cm'))) +
  scale_color_d3() +
  stat_compare_means(aes(group = sample), label = "p.signif")

p3

ggsave("tcga/out/arterial_ECs.png", plot = p3, width = 8, height = 4)

ggsave("tcga/out/arterial_ECs.pdf", plot = p3, width = 8, height = 4)

colnames(panca)

plot.info_all3 <- panca[, c(26, 24, 17)]

p3 = ggboxplot(plot.info_all3,
               x = "sampletype",
               y = "venous_ECs",
               color = "sample",
               fill = "white", #只需要修改这里，讲fill = "CellType"改为fill = "group"
               outlier.size = 0.5,
               xlab = "",
               ylab = "Fraction of  venous_ECs") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5),
    axis.title.y = element_text(margin = margin(0, 0.35, 0, 0,'cm'))) +
  scale_color_d3() +
  stat_compare_means(aes(group = sample), label = "p.signif")

p3

ggsave("tcga/out/venous_ECs.png", plot = p3, width = 8, height = 4)

ggsave("tcga/out/venous_ECs.pdf", plot = p3, width = 8, height = 4)

#panca = readRDS(file = "tcga/input/pancact.rds")

colnames(panca)

panca$periec = panca$Pericytes/panca$EC

colnames(panca)

plot.info_all5 <- panca[, c(26, 24, 27)]

write.csv(plot.info_all5, file = "tcga/out/peri_ECs.csv")

plot.info_all5 = read.csv(file = "tcga/out/peri_ECsa.csv")

p3 = ggboxplot(plot.info_all5,
               x = "sampletype",
               y = "periec",
               color = "sample",
               fill = "white", #只需要修改这里，讲fill = "CellType"改为fill = "group"
               outlier.size = 0.5,
               xlab = "",
               ylab = "Fraction of  Pericyte/EC") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5),
    axis.title.y = element_text(margin = margin(0, 0.35, 0, 0,'cm'))) +
  scale_color_d3() +
  stat_compare_means(aes(group = sample), label = "p.signif")

p3

ggsave("tcga/out/peri_ECs.png", plot = p3, width = 8, height = 4)

ggsave("tcga/out/peri_ECs.pdf", plot = p3, width = 8, height = 4)

colnames(panca)





