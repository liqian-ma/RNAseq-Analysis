
# Load libraries ----------------------------------------------------------

library(WGCNA)
library(tidyverse)
library(biomaRt)
library(org.Mm.eg.db)
library(topGO) 

# WGCNA -------------------------------------------------------------------

table = table[,c(grep("logFC", names(table)), grep("FDR", names(table)))]
table = table %>% mutate(genes = rownames(table)) %>% filter_at(vars(4:5), any_vars(. < 0.05))
rownames(table) = table$genes

g.all = rownames(table)
datExpr = t(lcpm[g.all, c(1:11, 22:26)])
dim(datExpr)

sample_names = as.character(y$samples$group)[c(1:11, 22:26)]
sample_names[1:6] = paste0(sample_names[1:6], seq(1:6))
sample_names[7:11] = paste0("27HC", seq(1:5))
sample_names[12:16] = paste0(sample_names[12:16], seq(1:5))

rownames(datExpr) = sample_names

powers = c(c(1:10), seq(from = 12, to = 30, by = 2))


# calculate the power needed for a signed analysis
powers = c(c(1:10), seq(from = 12, to=30, by=2))

disableWGCNAThreads()
system.time(
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, 
                           corFnc = cor, networkType = "signed hybrid")
)

# adjacency matrix
adjacency_signed = adjacency(datExpr, power = 5, type = "signed")
TOM = TOMsimilarity(adjacency_signed) 
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
minModuleSize = 20

# module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

MEList = moduleEigengenes(datExpr, colors = dynamicColors) 
MEs = MEList$eigengenes

# calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
# plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")

# automatic merging
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

mergedColors = merge$colors
mergedMEs = merge$newMEs

# plot
plotDendroAndColors(geneTree, mergedColors,
                    "Merged modules", dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)

moduleColors = mergedColors

# numerical labels
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
table(moduleColors)

# plot heatmaps to see module patterns
hms = list()
for (i in seq_along(unique(moduleLabels))) {
  color = unique(moduleColors)
  label = unique(moduleLabels)
  col = color[i]
  lab = label[i]
  
  cols = moduleColors == col
  datM = datExpr[, cols]
  datM = t(scale(datM))
  hm = Heatmap(datM,
               col = color.scale,
               show_row_names = FALSE, show_column_names = FALSE, 
               cluster_rows = TRUE, cluster_columns = FALSE, show_column_dend = TRUE, 
               show_row_dend = TRUE, row_dend_reorder = TRUE, column_dend_reorder = TRUE, 
               column_title = paste0("Module ", i, " ", col, " : ", nrow(datM), " genes"),
               width = unit(100, "mm"), 
               top_annotation = colAnn, 
               heatmap_legend_param = list(title = "SD from mean", 
                                           at = c(-2, 0, 2), 
                                           direction = "vertical", 
                                           title_position = "topleft",
                                           legend_width = unit(2.2, "cm")))
  hms[[i]] = hm
}


module_genes = list()
for (i in seq_along(unique(moduleLabels))) {
  color = unique(moduleColors)
  label = unique(moduleLabels)
  col = color[i]
  lab = label[i]
  
  cols = moduleColors == col
  datM = datExpr[, cols]
  genesM = colnames(datM)
  module_genes[[i]] = genesM
  names(module_genes)[i] = col
}


# GO enrichment  ----------------------------------------------------------

# http://www.bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

# get all genes
all_genes = rownames(y)

# run GO for each module
list_go_bp = list()
all_genes = rownames(y)
for (i in seq_along(module_genes)) {
  gene_list = factor(as.integer(all_genes %in% module_genes[[i]]))
  names(gene_list) = all_genes
  
  GOdata =  new("topGOdata", ontology = "BP", 
                allGenes = gene_list, 
                geneSel = function(p) p < 1e-2, description = "Test", 
                annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "Ensembl")
  
  resultFisher = runTest(GOdata, algorithm = "classic", statistic = "fisher")
  wgcna = GenTable(GOdata, classicFisher = resultFisher, topNodes = 20)
  list_go_bp[[i]] = wgcna
}
