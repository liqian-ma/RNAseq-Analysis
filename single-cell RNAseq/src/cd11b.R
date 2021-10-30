
# Load Libraries ----------------------------------------------------------


library(Seurat)
library(tidyverse)
library(biomaRt)


# Prep cd11b --------------------------------------------------------------

cd11b = readRDS("scaled.RDS")

# subset cd11b cells
cd11b = subset(cd11b, subset = Itgam > 1 & Cd3e <= 3)
cd11b2 = subset(cd11b, subset = Cd3e == 0)


cd11b2 = FindVariableFeatures(cd11b2)
cd11b2 = RunPCA(cd11b2, verbose = TRUE)

DimPlot(cd11b2, reduction = "pca", group.by = "treatment")
DimPlot(cd11b2, reduction = "pca", group.by = "sample_name")

ElbowPlot(cd11b2, 50)
cd11b2 = FindNeighbors(cd11b2, dims = 1:18, verbose = FALSE, reduction = "pca")
cd11b2 = FindClusters(cd11b2, verbose = FALSE, resolution = 0.3)
cd11b2 = RunUMAP(cd11b2, dims = 1:18, reduction = "pca")

Idents(cd11b2) = cd11b2@meta.data$RNA_snn_res.0.3
DimPlot(cd11b2, reduction = "umap", label = T, label.size = 8) + theme(legend.position = "")

Idents(cd11b2) = cd11b2@meta.data$RNA_snn_res.0.3
DimPlot(cd11b2, reduction = "umap", label = T, label.size = 8, split.by = "treatment") + 
  theme(legend.position = "")

FeaturePlot(cd11b2, features = c("Cd274", "Fcgr1", "Ly6c1", "Ly6g", "Klrk1", 
                                 "Cd19", "H2-Aa", "Prg2", "Adgre1", "Ccr3"
))

FeaturePlot(cd11b2, features = c("Cd274")) # PDL1+ cells
# cluster 4

FeaturePlot(cd11b2, features = c("Cd33", "Cd14")) # MHC-low cells
# 012

FeaturePlot(cd11b2, features = c("Klrk1")) # NK cells
# cluster 3

FeaturePlot(cd11b2, features = c("Ms4a1", "Cd19")) # Cd11b+ B cells
# cluster 7

FeaturePlot(cd11b2, features = c("Ly6g")) # neutrophils
# cluster 8

FeaturePlot(cd11b2, features = c("Cd86", "H2-Aa", "H2-Eb1")) # Ag presenting cells  #c7 + 5
FeaturePlot(cd11b2, features = c("Apoe", "Cd68", "H2-Aa")) # Ag-presenting macrophages
# cluster 5

FeaturePlot(cd11b2, features = c("Adgre1", "Irf8")) # Eosinophils
# https://www.nature.com/articles/s41586-019-1456-0?utm_source=other_website&utm_medium=display&utm_content=leaderboard&utm_campaign=JRCN_2_LW_X-moldailyfeed&error=cookies_not_supported&code=9cdf6076-c3c9-4c33-90fb-5851076e6cfd
# cluster 6

meta = cd11b2@meta.data %>% 
  mutate(barcodes = rownames(.),
         annotations = as.character(RNA_snn_res.0.3)) %>% 
  select(barcodes, annotations)

# annotate
pdl1 = meta %>% filter(annotations %in% c("4")) %>% pull(barcodes)
mmdsc0 = meta %>% filter(annotations %in% c("0")) %>% pull(barcodes)
mmdsc1 = meta %>% filter(annotations %in% c("1")) %>% pull(barcodes)
mmdsc2 = meta %>% filter(annotations %in% c("2")) %>% pull(barcodes)
nk = meta %>% filter(annotations %in% c("3")) %>% pull(barcodes)
bc = meta %>% filter(annotations %in% c("7")) %>% pull(barcodes)
gmdsc = meta %>% filter(annotations %in% c("8")) %>% pull(barcodes)
eosin = meta %>% filter(annotations %in% c("6")) %>% pull(barcodes)
agpmac = meta %>% filter(annotations %in% c("5")) %>% pull(barcodes)

## add to meta
meta = meta %>% mutate(manual_0.3_cd11b = case_when(
  barcodes %in% c(pdl1) ~"PDL1+ cells",
  barcodes %in% c(mmdsc0) ~"MHCIIlow monocytes 0",
  barcodes %in% c(mmdsc1) ~"MHCIIlow monocytes 1",
  barcodes %in% c(mmdsc2) ~"MHCIIlow monocytes 2",
  barcodes %in% c(nk) ~"NK cells",
  barcodes %in% c(bc) ~"B cells",
  barcodes %in% c(gmdsc) ~"Neutrophils",
  barcodes %in% c(eosin) ~"Eosinophils",
  barcodes %in% c(agpmac) ~"MHCIIhi macrophages",
  TRUE ~"unassigned cells"
))

# add to Seurat object
sample_vars = pull(meta, manual_0.3_cd11b)
names(sample_vars) = meta$barcodes

cd11b2 = AddMetaData(cd11b2, metadata = sample_vars, col.name = "manual_0.3_cd11b")

cd11b2$manual_0.3_cd11b = as.factor(cd11b2$manual_0.3_cd11b)
Idents(cd11b2) = "manual_0.3_cd11b"

DimPlot(cd11b2, reduction = "umap", split.by = "treatment") +
  scale_color_manual(name = "",
                     values = c("#A52A2A", "#D2B48C", "#6B8E23", "#32CD32",
                                "#4682B4", "#9370DB", "#DB7093", "#BC8F8F",
                                "#808080", "#FFB6C1"))

rm(list=setdiff(ls(), c("cd11b2")))

# Cell population difference ----------------------------------------------

meta = cd11b2@meta.data
total_cells = meta %>% group_by(sample_name) %>% summarise(total_n = n()) %>% ungroup()
subset_cells = meta %>% group_by(sample_name, manual_0.3_cd11b) %>% summarise(sub_n = n()) %>% ungroup()
subset_cells = left_join(subset_cells, total_cells, by = "sample_name")
subset_cells$percent = subset_cells$sub_n/subset_cells$total_n*100

subset_cells = arrange(subset_cells, manual_0.3_cd11b)

subset_cells2 = meta %>% group_by(treatment, manual_0.3_cd11b) %>% 
  summarise(mean_sub = n()/3) %>% ungroup()
total_cells2 = meta %>% group_by(treatment) %>% summarise(mean_total = n()/3) %>% ungroup()

subset_cells2 = left_join(subset_cells2, total_cells2, by = "treatment")
subset_cells2$percent = subset_cells2$mean_sub/subset_cells2$mean_total*100

ggplot(data = subset_cells2, aes(x = manual_0.3_cd11b, y = percent, fill = treatment)) +
  geom_col(position = "stack") + 
  scale_fill_brewer(palette = "Dark2") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

ggplot(data = subset_cells2, aes(x = treatment, y = percent, fill = manual_0.3_cd11b)) +
  geom_col(position = "stack") + 
  scale_fill_manual(name = "",
                    values = c("#A52A2A", "#D2B48C", "#6B8E23", "#32CD32",
                               "#4682B4", "#9370DB", "#DB7093", "#BC8F8F",
                               "#808080", "#FFB6C1")) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# heatmaps ----------------------------------------------------------------

cd11b2 = readRDS("cd11b2_annotated.RDS")

choles = c("Abca1", "Abcg1", "Srebf1", "Apoe")

# Monocyte0 population
Idents(cd11b2) = cd11b2@meta.data$sample_name
unique(cd11b2@meta.data$sample_name)

Idents(cd11b2) = cd11b2@meta.data$manual_0.3_cd11b
cd11b2_mdsc0 = subset(cd11b2, idents = c("MHClow monocytes 0"))

Idents(cd11b2_mdsc0) = cd11b2_mdsc0@meta.data$sample_name
avg = AverageExpression(cd11b2_mdsc0, return.seurat = TRUE)
avg_meta = avg@meta.data
avg_meta$sample_name = rownames(avg_meta)
avg_meta$treatment = str_extract(avg_meta$sample_name, ".*(?=\\.)")

# add to Seurat object
sample_vars = pull(avg_meta, treatment)
names(sample_vars) = avg_meta$sample_name

avg = AddMetaData(avg, metadata = sample_vars, col.name = "treatment")

Idents(avg) = "treatment"

DoHeatmap(avg, features = choles) + scale_fill_gradientn(colors = c("#4A6FE3", "#E2E2E2", "#D33F6A"), 
                                                         na.value = "white")

# Monocyte1 population

cd11b2_mdsc1 = subset(cd11b2, idents = c("MHClow monocytes 1"))

Idents(cd11b2_mdsc1) = cd11b2_mdsc1@meta.data$sample_name
avg = AverageExpression(cd11b2_mdsc1, return.seurat = TRUE)
avg_meta = avg@meta.data
avg_meta$sample_name = rownames(avg_meta)
avg_meta$treatment = str_extract(avg_meta$sample_name, ".*(?=\\.)")


# add to Seurat object
sample_vars = pull(avg_meta, treatment)
names(sample_vars) = avg_meta$sample_name

avg = AddMetaData(avg, metadata = sample_vars, col.name = "treatment")

Idents(avg) = "treatment"

DoHeatmap(avg, features = choles) + scale_fill_gradientn(colors = c("#4A6FE3", "#E2E2E2", "#D33F6A"), 
                                                         na.value = "white")

# Monocyte2 population

cd11b2_mdsc2 = subset(cd11b2, idents = c("MHClow monocytes 2"))

Idents(cd11b2_mdsc2) = cd11b2_mdsc2@meta.data$sample_name
avg = AverageExpression(cd11b2_mdsc2, return.seurat = TRUE)
avg_meta = avg@meta.data
avg_meta$sample_name = rownames(avg_meta)
avg_meta$treatment = str_extract(avg_meta$sample_name, ".*(?=\\.)")


# add to Seurat object
sample_vars = pull(avg_meta, treatment)
names(sample_vars) = avg_meta$sample_name

avg = AddMetaData(avg, metadata = sample_vars, col.name = "treatment")

Idents(avg) = "treatment"

DoHeatmap(avg, features = choles) + scale_fill_gradientn(colors = c("#4A6FE3", "#E2E2E2", "#D33F6A"), 
                                                         na.value = "white")

# Ag-pres macrophages population

Idents(cd11b2) = cd11b2@meta.data$sample_name
unique(cd11b2@meta.data$sample_name)

Idents(cd11b2) = cd11b2@meta.data$manual_0.3_cd11b
cd11b2_apm = subset(cd11b2, idents = c("Antigen-presenting macrophages"))

Idents(cd11b2_apm) = cd11b2_apm@meta.data$sample_name
avg = AverageExpression(cd11b2_apm, return.seurat = TRUE)
avg_meta = avg@meta.data
avg_meta$sample_name = rownames(avg_meta)
avg_meta$treatment = str_extract(avg_meta$sample_name, ".*(?=\\.)")


# add to Seurat object
sample_vars = pull(avg_meta, treatment)
names(sample_vars) = avg_meta$sample_name

avg = AddMetaData(avg, metadata = sample_vars, col.name = "treatment")

Idents(avg) = "treatment"

DoHeatmap(avg, features = choles) + scale_fill_gradientn(colors = c("#4A6FE3", "#E2E2E2", "#D33F6A"), 
                                                         na.value = "white")

# PDL1+ cells population

Idents(cd11b2) = cd11b2@meta.data$sample_name
unique(cd11b2@meta.data$sample_name)

Idents(cd11b2) = cd11b2@meta.data$manual_0.3_cd11b
cd11b2_pdl1 = subset(cd11b2, idents = c("PDL1+ cells"))

Idents(cd11b2_pdl1) = cd11b2_pdl1@meta.data$sample_name
avg = AverageExpression(cd11b2_pdl1, return.seurat = TRUE)
avg_meta = avg@meta.data
avg_meta$sample_name = rownames(avg_meta)
avg_meta$treatment = str_extract(avg_meta$sample_name, ".*(?=\\.)")


# add to Seurat object
sample_vars = pull(avg_meta, treatment)
names(sample_vars) = avg_meta$sample_name

avg = AddMetaData(avg, metadata = sample_vars, col.name = "treatment")

Idents(avg) = "treatment"

DoHeatmap(avg, features = choles) + scale_fill_gradientn(colors = c("#4A6FE3", "#E2E2E2", "#D33F6A"), 
                                                         na.value = "white")

# eosinophils

Idents(cd11b2) = cd11b2@meta.data$sample_name
unique(cd11b2@meta.data$sample_name)

Idents(cd11b2) = cd11b2@meta.data$manual_0.3_cd11b
cd11b2_tam = subset(cd11b2, idents = c("Eosinophils"))

Idents(cd11b2_tam) = cd11b2_tam@meta.data$sample_name
avg = AverageExpression(cd11b2_tam, return.seurat = TRUE)
avg_meta = avg@meta.data
avg_meta$sample_name = rownames(avg_meta)
avg_meta$treatment = str_extract(avg_meta$sample_name, ".*(?=\\.)")


# add to Seurat object
sample_vars = pull(avg_meta, treatment)
names(sample_vars) = avg_meta$sample_name

avg = AddMetaData(avg, metadata = sample_vars, col.name = "treatment")

Idents(avg) = "treatment"

DoHeatmap(avg, features = choles) + scale_fill_gradientn(colors = c("#4A6FE3", "#E2E2E2", "#D33F6A"), 
                                                         na.value = "white")

# DEG ---------------------------------------------------------------------

cd11b2$celltype_treat = paste(Idents(cd11b2), cd11b2$treatment, sep = "_")
ident_test = unique(cd11b2$celltype_treat)
vehicle_test = ident_test[grepl("_Vehicle$", ident_test)] %>% sort(.)
HC_test = ident_test[grepl("_27HC$", ident_test)] %>% sort(.)
gw3965_test = ident_test[grepl("_GW3965$", ident_test)] %>% sort(.)
cyp_test = ident_test[grepl("_GW273297X$", ident_test)] %>% sort(.)

Idents(cd11b2) = "celltype_treat"
marker_list = list()
for (i in 5:7) {
  marker = FindMarkers(cd11b2, ident.1 = HC_test[i], ident.2 = vehicle_test[i], 
                       verbose = FALSE, logfc.threshold = 0.1)
  marker = marker %>% filter(p_val_adj < 0.1)
  marker_list[[i]] = marker
  names(marker_list)[i] = paste0(HC_test[i], "vsVeh")
}

marker_list_gw = list()
for (i in 5:7) {
  marker = FindMarkers(cd11b2, ident.1 = gw3965_test[i], ident.2 = vehicle_test[i], 
                       verbose = FALSE, logfc.threshold = 0.1)
  marker = marker %>% filter(p_val_adj < 0.1)
  marker_list_gw[[i]] = marker
  names(marker_list_gw)[i] = paste0(gw3965_test[i], "vsVeh")
}

gw_mdsc0 = marker_list_gw[[5]]
hc_mdsc0 = marker_list[[5]]
gw_mdsc1 = marker_list_gw[[6]]
hc_mdsc1 = marker_list[[6]]
gw_mdsc2 = marker_list_gw[[7]]
hc_mdsc2 = marker_list[[7]]

length(intersect(rownames(gw_mdsc0), rownames(hc_mdsc0)))

# Venn Diagram ------------------------------------------------------------

siglist_27HC = list()
for (i in 1:3) {
  markers = marker_list_27hc[[i+4]]
  markers = markers %>% mutate(sig27HC = ifelse(p_val_adj < 0.1, 1, 0))
  markers$sig27HC[markers$p_val_adj < 0.1 & markers$avg_log2FC > 0] = 1
  markers$sig27HC[markers$p_val_adj < 0.1 & markers$avg_log2FC < 0] = -1
  siglist_27HC[[i]] = markers
  names(siglist_27HC)[i] = names(marker_list_27hc)[i+4]
}


siglist_GW = list()
for (i in 1:3) {
  markers = marker_list_gw[[i+4]]
  markers = markers %>% mutate(sigGW = ifelse(p_val_adj < 0.1, 1, 0))
  markers$sigGW[markers$p_val_adj < 0.1 & markers$avg_log2FC > 0] = 1
  markers$sigGW[markers$p_val_adj < 0.1 & markers$avg_log2FC < 0] = -1
  siglist_GW[[i]] = markers
  names(siglist_GW)[i] = names(marker_list_gw)[i+4]
}

df_mdsc0 = cbind(siglist_27HC[[1]]$sig27HC, siglist_GW[[1]]$sigGW)
colnames(df_mdsc0) = c("27HC", "GW3965")


df_mdsc1 = cbind(siglist_27HC[[2]]$sig27HC, siglist_GW[[2]]$sigGW)
colnames(df_mdsc1) = c("27HC", "GW3965")

df_mdsc2 = cbind(siglist_27HC[[3]]$sig27HC, siglist_GW[[3]]$sigGW)
colnames(df_mdsc2) = c("27HC", "GW3965")

library(limma)
vennDiagram(df_mdsc0[, c(1:2)], include = "down", 
            main = "downregulated genes", circle.col=c("darkred", "orchid"),
            names = c("27HC", "GW3965"), cex = 1)


vennDiagram(df_mdsc1[, c(1:2)], include = "down", 
            main = "downregulated genes", circle.col=c("darkred", "orchid"),
            names = c("27HC", "GW3965"), cex = 1)

vennDiagram(df_mdsc2[, c(1:2)], include = "down", 
            main = "downregulated genes", circle.col=c("darkred", "orchid"),
            names = c("27HC", "GW3965"), cex = 1)


# Enrichment Analysis -----------------------------------------------------

markers = load(file = "results/seurat/data/cd11b2_markers.RData")

gw_mdsc0 = marker_list_gw[[5]]
hc_mdsc0 = marker_list_27hc[[5]]
gw_mdsc1 = marker_list_gw[[6]]
hc_mdsc1 = marker_list_27hc[[6]]
gw_mdsc2 = marker_list_gw[[7]]
hc_mdsc2 = marker_list_27hc[[7]]

library(msigdbr)
# Retrieve mouse hallmark collection gene sets
m_df = msigdbr(species = "Mus musculus", category = "H")

# Use the gene sets data frame for fgsea
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
names(m_list) = str_remove_all(names(m_list), "HALLMARK_")

library(EGSEA)
# monocyte0 27HC
idx = buildCustomIdx(rownames(hc_mdsc0)[hc_mdsc0$p_val_adj < 0.05], 
                     gsets = m_list, anno = NULL, label = "custom",
                     name = "Hallmark", species = "Mouse", min.size = 1)

test = egsea.ora(rownames(hc_mdsc0)[hc_mdsc0$p_val_adj < 0.05], 
                 universe = rownames(cd11b2),
                 logFC = hc_mdsc0$avg_log2FC[hc_mdsc0$p_val_adj < 0.05], 
                 gs.annots = idx)

plot_test = test@results$custom$test.results$ExperimentalContrast
plot_test$pathway = rownames(plot_test)

plot_test$pathway = gsub("HALLMARK_", "", plot_test$pathway)
plot_test$pathway = gsub("_", " ", plot_test$pathway)
# plot_test$pathway = str_to_sentence(plot_test$pathway)

ggplot(plot_test, aes(reorder(pathway, avg.logfc.dir), avg.logfc.dir)) +
  geom_col(aes(fill = p.adj < 0.05)) + 
  scale_fill_manual(values = c("TRUE"="lightpink3", "FALSE"="lightblue4")) +
  coord_flip() +
  labs(x = "Pathway", y="Average logFC") +  theme_classic()


# monocyte0 gw3965
idx = buildCustomIdx(rownames(gw_mdsc0)[gw_mdsc0$p_val_adj < 0.05], 
                     gsets = m_list, anno = NULL, label = "custom",
                     name = "Hallmark", species = "Mouse", min.size = 1)

test = egsea.ora(rownames(gw_mdsc0)[gw_mdsc0$p_val_adj < 0.05], 
                 universe = rownames(cd11b2),
                 logFC = gw_mdsc0$avg_log2FC[gw_mdsc0$p_val_adj < 0.05], 
                 gs.annots = idx)

plot_test = test@results$custom$test.results$ExperimentalContrast
plot_test$pathway = rownames(plot_test)

plot_test$pathway = gsub("HALLMARK_", "", plot_test$pathway)
plot_test$pathway = gsub("_", " ", plot_test$pathway)

ggplot(plot_test, aes(reorder(pathway, avg.logfc.dir), avg.logfc.dir)) +
  geom_col(aes(fill = p.adj < 0.05)) + 
  scale_fill_manual(values = c("TRUE"="lightpink3", "FALSE"="lightblue4")) +
  coord_flip() +
  labs(x = "Pathway", y="Average logFC") +  theme_classic()

# Prepare data for RNA velocity -------------------------------------------

# meta data
meta = cd11b2@meta.data

customFun = function(DF, folder) {
  write.csv(DF, paste0("results/seurat/RNA velocity/", folder, 
                       "/", unique(DF$orig.ident), ".csv"), row.names = FALSE)
  return(DF)
}

meta %>% mutate(barcodes = str_extract(rownames(.), "(?<=_).*")) %>%
  mutate(barcodes = str_extract(barcodes, ".*(?=-)")) %>% 
  dplyr::select(orig.ident, barcodes) %>% group_by(orig.ident) %>% 
  do(customFun(., folder = "cells"))

Embeddings(cd11b, reduction = "umap") %>% as.data.frame() %>%
  mutate(barcodes = str_extract(rownames(.), "(?<=_).*")) %>% 
  mutate(barcodes = str_extract(barcodes, ".*(?=-)")) %>% 
  mutate(orig.ident = str_extract(rownames(.), "[0-9]+(?=_)")) %>%
  dplyr::select(barcodes, orig.ident, UMAP_1, UMAP_2) %>% group_by(orig.ident) %>% 
  do(customFun(.))

color = c("#A52A2A", "#D2B48C", "#6B8E23", "#32CD32",
          "#4682B4", "#9370DB", "#DB7093", "#BC8F8F",
          "#808080")
clusters = levels(cd11b2@meta.data$manual_0.3_cd11b)
cluster_col = data.frame(manual_0.3_cd11b = clusters, color = color)


meta %>% mutate(barcodes = str_extract(rownames(.), "(?<=_).*")) %>%
  mutate(barcodes = str_extract(barcodes, ".*(?=-)")) %>% 
  left_join(cluster_col, by = "manual_0.3_cd11b") %>%
  dplyr::select(orig.ident, barcodes, manual_0.3_cd11b, color) %>% group_by(orig.ident) %>% 
  do(customFun(., folder = "clusters"))
