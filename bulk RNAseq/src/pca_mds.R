# Load libraries ----------------------------------------------------------

library(dplyr) 
library(stringr)
library(gplots)
library(ggplot2)
library(ggrepel)
library(biomaRt)
library(tximport)
library(limma)
library(edgeR)
library(Glimma)
library(ComplexHeatmap)
library(fgsea)
library(GSVA)
library(msigdbr)
library(WGCNA)

# Dimensionality reduction ------------------------------------------------

pc_data = as.data.frame(lcpm[, c(1:11, 22:26)])
pc_data = t(scale(t(pc_data)))
rna_pca = prcomp(x = t(pc_data))

varExp = (100*rna_pca$sdev^2)/sum(rna_pca$sdev^2)
varDF = data.frame(PCs = 1:length(varExp), varExp = varExp)

# extract the PCs and re-combine with labels
pcs = as.data.frame(rna_pca$x)
pcs_number = pcs %>% mutate(treatment = c(rep("Veh", 6), rep("27HC", 5), rep("GW3965", 5))) 
pcs_number$treatment = factor(pcs_number$treatment, 
                              levels = c("Veh","27HC","GW3965"))

# plot the data using the first 2 PCs
ggplot(data = pcs_number, aes(x = PC1, y = PC2, color = treatment)) + 
  geom_point(size = 2) + theme_classic() + 
  labs(x = paste0("PC1: ", round(varDF$varExp[1], 2), "% variance explained"), 
       y = paste0("PC2: ", round(varDF$varExp[2], 2), "% variance explained")) +
  scale_color_manual(values = c("navy", "darkred", "olivedrab3")) +
  theme_classic() +
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 10))  

rot = as.data.frame(rna_pca$rotation)
rot$ensembl_gene_id = rownames(rot)
rot = rot[rev(order(abs(rot$PC2))), ]
rot = rot[, c(1:2, 17)]
rot = rot %>% left_join(genes, by = "ensembl_gene_id")
pseudo = str_detect(rot$description, "pseudo", negate = TRUE)
pred = str_detect(rot$description, "predicted", negate = TRUE)

rot = rot[pseudo, ]
rot = rot[pred, ]

options(ggrepel.max.overlaps = Inf)

ggplot(rot[1:20, ], aes(x = PC1,y = PC2, label = external_gene_name)) + 
  geom_point() + theme_classic() + 
  geom_text_repel() + xlim(c(-0.065, 0.1)) +   
  labs(x = paste0("PC1: ",round(varExp[1],1),"%"),
       y = paste0("PC2: ",round(varExp[2],1),"%"), 
       title = "Loadings plot of PCA") + 
  theme(text = element_text(size = 10))

# mds
distance = dist(t(pc_data))
distance_mtx = as.matrix(distance)

mds = data.frame(cmdscale(distance_mtx))
mds = cbind(mds, c(rep("Veh", 6), rep("27HC", 5), rep("GW3965", 5)))
colnames(mds)[3] = "treatment"
relevel(as.factor(mds$treatment), ref = "Veh")
mds$treatment = factor(mds$treatment, levels = c("Veh",
                                                 "27HC",
                                                 "GW3965"))

ggplot(mds, aes(X1, X2, color = treatment)) + 
  geom_point(size = 2) +
  labs(x = "Dim 1", y = "Dim 2") +
  scale_color_manual(values = c("navy", "darkred", "orchid")) +
  theme_classic() +
  theme(axis.title = element_text(size = 11), axis.text = element_text(size = 10))
