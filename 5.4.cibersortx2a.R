
getwd()

library(ggplot2)

library(ggthemes)

library(dplyr)

library(data.table)

library(corrplot)

rm(list = ls())

panca = read.csv(file = "cellfrac/pancact.csv")

table(panca$ieep_cells)

panca = filter(panca, ieep_cells != "Epithelial_cells" &
                      ieep_cells != "Fibroblasts" &
                      ieep_cells != "Pericytes")

table(panca$ieep_cells)

saveRDS(panca, file = "cellfrac/pancacti.rds")

dir.create("cellfrac2")

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

ggsave(p1, file = "tcga/cellfrac2/tipeccor.png", width = 6.5, height = 2.8)

ggsave(p1, file = "tcga/cellfrac2/tipeccor.pdf", width = 6.5, height = 2.8)

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

ggsave(p1, file = "tcga/cellfrac2/veccor.png", width = 6.5, height = 2.8)

ggsave(p1, file = "tcga/cellfrac2/veccor.pdf", width = 6.5, height = 2.8)

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

ggsave(p1, file = "tcga/cellfrac2/accor.png", width = 6.5, height = 2.8)

ggsave(p1, file = "tcga/cellfrac2/accor.pdf", width = 6.5, height = 2.8)

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

ggsave(p1, file = "tcga/cellfrac2/c1ccor.png", width = 6.5, height = 3.45)

ggsave(p1, file = "tcga/cellfrac2/c1ccor.pdf", width = 6.5, height = 3.45)

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

ggsave(p1, file = "tcga/cellfrac2/c2ccor.png", width = 6.5, height = 3.45)

ggsave(p1, file = "tcga/cellfrac2/c2ccor.pdf", width = 6.5, height = 3.45)

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

ggsave(p1, file = "tcga/cellfrac2/lccor.png", width = 6.5, height = 2.8)

ggsave(p1, file = "tcga/cellfrac2/lccor.pdf", width = 6.5, height = 2.8)

table(panca$factori)


