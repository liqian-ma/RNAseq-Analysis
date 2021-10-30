# Load Libraries ----------------------------------------------------------

library(Seurat)
library(tidyverse)
library(Matrix)

# Read in cellranger matrices ---------------------------------------------

# read metadata into a table
sc_meta = read.table("~/scRNAseq/src/sc_meta.csv",
                     header = TRUE, sep = ",")

# read in 10x filtered matrices
files = file.path("~/scRNAseq/results/cellranger.compact/cellranger", 
                  sc_meta[, 1], "outs/filtered_feature_bc_matrix")
list10x = sapply(files,
                 function(x) Read10X(x))
names(list10x) = sc_meta[,1]

# filter lowly expressed genes
cellsPerGene = sapply(list10x, function(x) Matrix::rowSums(x > 0))
i.filter = which(rowSums(cellsPerGene) >= 15)
list10x_filt = sapply(list10x, function(x) x[i.filter,])

# make seurat object
seur_list = sapply(names(list10x_filt),
                   function(x) CreateSeuratObject(list10x_filt[[x]], project = x))
seur_all = merge(x = seur_list[[1]], y = seur_list[-1],
                 add.cell.ids = names(seur_list),
                 project = "4T1")

# number of cells per sample
(startNum <- table(seur_all$orig.ident))


# add mitochondrial gene percentage to metadata
seur_all[["percent.mt"]] = PercentageFeatureSet(seur_all, pattern = "^mt-")
head(seur_all@meta.data)

# add more experimental information to metadata
meta = seur_all@meta.data
colnames(sc_meta)[1] = "orig.ident"
sc_meta$orig.ident = as.character(sc_meta$orig.ident)
meta = left_join(meta, sc_meta, by = "orig.ident")
head(meta)
add_meta = meta[, 5:9]
rownames(add_meta) = rownames(seur_all@meta.data)
seur_all = AddMetaData(seur_all, add_meta)

head(seur_all@meta.data)


# Visualization prefiltering ----------------------------------------------

# visualize
seur_all$treatment = factor(seur_all$treatment, 
                            levels = c("Vehicle", "27HC", "GW3965", "GW273297X"))
seur_all$sample_name = factor(seur_all$sample_name, 
                              levels = c("Vehicle.1", "Vehicle.2", "Vehicle.3", 
                                         "27HC.1", "27HC.2", "27HC.3",
                                         "GW3965.1", "GW3965.2", "GW3965.3",
                                         "CYP.1", "CYP.2", "CYP.3"))
Idents(seur_all) = "sample_name"

# Filter samples ----------------------------------------------------------
umi = seur_all$nCount_RNA
UMI.threshold = 20000
mt.threshold = 7.5

seur_all = subset(seur_all,
                  subset = nCount_RNA < UMI.threshold & 
                    percent.mt < mt.threshold &
                    nFeature_RNA > 200)

(endNum <- table(seur_all$orig.ident))

# Data processing ---------------------------------------------------------

# normalize and scale data
seur_all = NormalizeData(seur_all)
# all.genes = rownames(seur_all)
seur_all = ScaleData(seur_all)

saveRDS(seur_all, "~/scRNAseq/results/scaled.RDS")


