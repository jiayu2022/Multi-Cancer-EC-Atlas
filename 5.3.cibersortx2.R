
getwd()

library(ggplot2)

library(ggthemes)

library(dplyr)

library(data.table)

library(corrplot)

rm(list = ls())

paad = readRDS(file = "tcga/paad/immu/correlation.rds")

paad$cancertype = rep("PAAD", nrow(paad))

table(paad$cancertype)

ovc = readRDS(file = "tcga/ovc/immu/correlation.rds")

ovc$cancertype = rep("OVC", nrow(ovc))

table(ovc$cancertype)

brca = readRDS(file = "tcga/brca/immu/correlation.rds")

brca$cancertype = rep("BRCA", nrow(brca))

table(brca$cancertype)

lihc = readRDS(file = "tcga/lihc/immu/correlation.rds")

lihc$cancertype = rep("LIHC", nrow(lihc))

table(lihc$cancertype)

coad = readRDS(file = "tcga/coad/immu/correlation.rds")

coad$cancertype = rep("COAD", nrow(coad))

table(coad$cancertype)

read = readRDS(file = "tcga/read/immu/correlation.rds")

read$cancertype = rep("READ", nrow(read))

table(read$cancertype)

prad = readRDS(file = "tcga/prad/immu/correlation.rds")

prad$cancertype = rep("PRAD", nrow(prad))

table(prad$cancertype)

luad = readRDS(file = "tcga/luad/immu/correlation.rds")

luad$cancertype = rep("LUAD", nrow(luad))

table(luad$cancertype)

lusc = readRDS(file = "tcga/lusc/immu/correlation.rds")

lusc$cancertype = rep("LUSC", nrow(lusc))

table(lusc$cancertype)

blca = readRDS(file = "tcga/blca/immu/correlation.rds")

blca$cancertype = rep("BLCA", nrow(blca))

table(blca$cancertype)

esca = readRDS(file = "tcga/esca/immu/correlation.rds")

esca$cancertype = rep("ESCA", nrow(esca))

table(esca$cancertype)

gbm = readRDS(file = "tcga/gbm/immu/correlation.rds")

gbm$cancertype = rep("GBM", nrow(gbm))

table(gbm$cancertype)

hnsc = readRDS(file = "tcga/hnsc/immu/correlation.rds")

hnsc$cancertype = rep("HNSC", nrow(hnsc))

table(hnsc$cancertype)

thca = readRDS(file = "tcga/thca/immu/correlation.rds")

thca$cancertype = rep("THCA", nrow(thca))

table(thca$cancertype)

stad = readRDS(file = "tcga/stad/immu/correlation.rds")

stad$cancertype = rep("STAD", nrow(stad))

table(stad$cancertype)

cesc = readRDS(file = "tcga/cesc/immu/correlation.rds")

cesc$cancertype = rep("CESC", nrow(cesc))

table(cesc$cancertype)

ucec = readRDS(file = "tcga/ucec/immu/correlation.rds")

ucec$cancertype = rep("UCEC", nrow(ucec))

table(ucec$cancertype)

kirc = readRDS(file = "tcga/kirc/immu/correlation.rds")

kirc$cancertype = rep("KIRC", nrow(kirc))

table(kirc$cancertype)

kich = readRDS(file = "tcga/kich/immu/correlation.rds")

kich$cancertype = rep("KICH", nrow(kich))

table(kich$cancertype)

kirp = readRDS(file = "tcga/kirp/immu/correlation.rds")

kirp$cancertype = rep("KIRP", nrow(kirp))

table(kirp$cancertype)

colnames(paad)

colnames(kirp)

l = list(blca, brca, cesc, coad, esca,
         gbm, hnsc, kich, kirc, kirp,
         lihc, luad, lusc, ovc, paad,
         prad, read, stad, thca, ucec)

panca = rbindlist(l, use.names=TRUE)

table(panca$cancertype)

dir.create("cellfrac")

colnames(panca)

panca$pstar <- ifelse(panca$p.value < 0.05,
                     ifelse(panca$p.value < 0.01,
                            ifelse(panca$p.value < 0.001,"***", "**"), "*"), "")

write.csv(panca, file = "cellfrac/pancact.csv")

#panca = read.csv(file = "cellfrac/pancacta.csv")

table(panca$ieep_cells)

panca = filter(panca, ieep_cells != "Epithelial_cells")

table(panca$ieep_cells)

saveRDS(panca, file = "cellfrac/pancact.rds")

pancatb = tibble::as_tibble(panca)

table(panca$factori)

tip = filter(pancatb, factori == "tip_like_ECs")

p1 = ggplot(tip, aes(cancertype, ieep_cells)) +
  geom_tile(aes(fill = cor), colour = "white", size = 0.1) +
  scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red") +
  geom_text(aes(label = pstar), col ="black", size = 5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))+
  labs(title = "tip ECs", fill =paste0("Correlation")) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5))

p1

ggsave(p1, file = "cellfrac/tipeccor.png", width = 6.5, height = 3.2)

ggsave(p1, file = "cellfrac/tipeccor.pdf", width = 6.5, height = 3.2)

table(panca$factori)

v = filter(pancatb, factori == "venous_ECs")

p1 = ggplot(v, aes(cancertype, ieep_cells)) +
  geom_tile(aes(fill = cor), colour = "white", size = 0.1) +
  scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red") +
  geom_text(aes(label = pstar), col ="black", size = 5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))+
  labs(title = "venous_ECs", fill =paste0("Correlation")) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5))

p1

ggsave(p1, file = "tcga/cellfrac/veccor.png", width = 6.5, height = 3.45)

ggsave(p1, file = "tcga/cellfrac/veccor.pdf", width = 6.5, height = 3.45)

table(panca$factori)

a = filter(pancatb, factori == "arterial_ECs")

p1 = ggplot(a, aes(cancertype, ieep_cells)) +
  geom_tile(aes(fill = cor), colour = "white", size = 0.1) +
  scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red") +
  geom_text(aes(label = pstar), col ="black", size = 5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))+
  labs(title = "arterial_ECs", fill =paste0("Correlation")) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5))

p1

ggsave(p1, file = "tcga/cellfrac/accor.png", width = 6.5, height = 3.45)

ggsave(p1, file = "tcga/cellfrac/accor.pdf", width = 6.5, height = 3.45)

table(panca$factori)

c1 = filter(pancatb, factori == "capillary_ECs_1")

p1 = ggplot(c1, aes(cancertype, ieep_cells)) +
  geom_tile(aes(fill = cor), colour = "white", size = 0.1) +
  scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red") +
  geom_text(aes(label = pstar), col ="black", size = 5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))+
  labs(title = "capillary_ECs_1", fill =paste0("Correlation")) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5))

p1

ggsave(p1, file = "tcga/cellfrac/c1ccor.png", width = 6.5, height = 3.45)

ggsave(p1, file = "tcga/cellfrac/c1ccor.pdf", width = 6.5, height = 3.45)

table(panca$factori)

c2 = filter(pancatb, factori == "capillary_ECs_2")

p1 = ggplot(c2, aes(cancertype, ieep_cells)) +
  geom_tile(aes(fill = cor), colour = "white", size = 0.1) +
  scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red") +
  geom_text(aes(label = pstar), col ="black", size = 5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))+
  labs(title = "capillary_ECs_2", fill =paste0("Correlation")) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5))

p1

ggsave(p1, file = "tcga/cellfrac/c2ccor.png", width = 6.5, height = 3.45)

ggsave(p1, file = "tcga/cellfrac/c2ccor.pdf", width = 6.5, height = 3.45)

table(panca$factori)

l = filter(pancatb, factori == "lymphatic_ECs")

p1 = ggplot(l, aes(cancertype, ieep_cells)) +
  geom_tile(aes(fill = cor), colour = "white", size = 0.1) +
  scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red") +
  geom_text(aes(label = pstar), col ="black", size = 5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))+
  labs(title = "lymphatic_ECs", fill =paste0("Correlation")) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5))

p1

ggsave(p1, file = "tcga/cellfrac/lccor.png", width = 6.5, height = 3.45)

ggsave(p1, file = "tcga/cellfrac/lccor.pdf", width = 6.5, height = 3.45)

table(panca$factori)



