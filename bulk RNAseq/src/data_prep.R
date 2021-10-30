
# DEG analysis using Limma-Voom --------------------------------------


# Load libraries ----------------------------------------------------------

# If do not have these libraries, uncomment and run the following code to download:
# install.packages(c("tidyverse", "ggplot2", "gplots"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c("biomaRt", "tximport", "limma", "edgeR", "Glimma", "ComplexHeatmap", "fgsea", "msigdbr"))

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

# Import data and pre-process ---------------------------------------------

# download mappings between transcript ID and gene ID
mm = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
tx2gene = getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = mm)
head(tx2gene)

# read in the sample information, relevel group factors
samples = read.table("Samples.txt", header = FALSE)
samples$V2 = relevel(as.factor(samples$V2), ref = "Vehicle")

# take the first column of table to build a file path to each quant/abundance file
files = file.path("kall_out_1", samples[, 1], "abundance.h5")
levels(samples$V2)[levels(samples$V2)=="27HC"] = "HC"
levels(samples$V2)[levels(samples$V2)=="GSK2033+27HC"] = "GSK_HC"
levels(samples$V2)[levels(samples$V2)=="GW3965+ICI"] = "GW_ICI"
levels(samples$V2)[levels(samples$V2)=="ICI+27HC"] = "ICI_HC"

# import Kallisto output and summarize to gene-level abundance
txi.kallisto = tximport(files, type = "kallisto", tx2gene = tx2gene, 
                        ignoreTxVersion = TRUE, 
                        countsFromAbundance = "lengthScaledTPM")
head(txi.kallisto)
names(files) = samples[ ,1]
head(txi.kallisto)


# Limma-voom analysis of DEG -----------------------------------------------------

# read in gene counts, define samples and groups
y = DGEList(txi.kallisto$counts, samples = samples, group = samples[, 2])

# TMM normalization
y = calcNormFactors(y)

# filter out the genes
keep = rowSums(cpm(y$counts) > 0) >= 1
y = y[keep, ]

# Processed data deposited to GEO
geo_processed = y$counts
samples_geo = samples
samples_geo$V2 = str_replace_all(samples_geo$V2, "HC", "27HC")
samples_geo$V2 = str_replace_all(samples_geo$V2, "\\_", " + ")
samples_geo$V2 = str_replace_all(samples_geo$V2, "GSK \\+ 27HC", "GSK2033 + 27HC")
samples_geo$V2 = str_replace_all(samples_geo$V2, "GW \\+ ICI", "GW3965 + 27HC")
samples_geo = samples_geo %>% group_by(V2) %>% 
  mutate(title = paste0(V2, ".", 1:length(V2)))

colnames(geo_processed) = samples_geo$title
geo_processed = as.data.frame(geo_processed)
geo_processed$ensembl_gene_id = rownames(geo_processed)
geo_processed = geo_processed %>% dplyr::select(ensembl_gene_id, everything())
write_delim(geo_processed, 
            "processed_data_bulkRNAseq.txt",
            delim = "\t")

# design matrix
design = model.matrix(~ 0 + samples[, 2])
colnames(design) = levels(y$samples$group)

# Transform count data to log2-counts per million (logCPM) 
# estimate the mean-variance relationship 
# use this to compute appropriate observation-level weights for linear modeling
v = voom(y, design)

# Fit linear model for each gene
fit = lmFit(v, design)

