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

# Heatmap -----------------------------------------------------------------

my.contrasts = makeContrasts(
  HCvsVeh = HC - Vehicle,
  GWvsVeh = GW3965 - Vehicle,
  E2vsVeh = E2 - Vehicle,
  levels = design
)

# fit the model with contrast
fit.con = contrasts.fit(fit, my.contrasts)
fit.con = eBayes(fit.con)

con.coded = decideTests(fit.con, p.value = 0.05, n = Inf) 

g.all = rownames(con.coded[con.coded[, 1] != 0 | con.coded[, 2] != 0 ,])
heatgenes = lcpm[g.all, c(1:11, 22:26)]
heat.scaled = heatgenes %>% t() %>% scale() %>% t()

# set color scale of the heatmap
color.scale = colorpanel(1024, "#4A6FE3", "#E2E2E2", "#D33F6A")

# set column annotation by treatment groups
library(ComplexHeatmap)
ann = pcs_number$treatment  
ann = samples2$V2
# ann = ann %>% recode(., "HC" = "27HC") %>% data.frame(.)
ann = ann %>% recode(., "Veh" = "Vehicle") %>% data.frame(.)


colnames(ann) = "treatment"
colors = list("treatment"=c("Vehicle" = "navy", "27HC" = "darkred", "GW3965"="olivedrab3"))

colAnn = HeatmapAnnotation(df = ann, which = "col", col = colors, 
                           annotation_width=unit(c(1, 4), "cm"), gap = unit(1, "mm"))

# plot heatmap
p = Heatmap(heat.scaled,
            col = color.scale,
            show_row_names = FALSE, show_column_names = FALSE, 
            cluster_rows = TRUE, cluster_columns = FALSE, show_column_dend = TRUE, 
            show_row_dend = TRUE, row_dend_reorder = TRUE, column_dend_reorder = TRUE, 
            width = unit(100, "mm"), 
            top_annotation = colAnn, 
            heatmap_legend_param = list(title = "SD from mean", 
                                        at = c(-2, 0, 2), 
                                        direction = "vertical", 
                                        title_position = "topleft",
                                        legend_width = unit(2.2, "cm")))


# Volcano plot ------------------------------------------------------------

table = data.frame(a = 1:24671)
for (i in 1:3) {
  test = topTable(fit.con, sort.by = "none", adjust.method = "BH",
                  number = Inf, coef = i)
  test = test[c(1, 4)]
  test = test[sort(rownames(test)),]
  colnames(test)[1] = paste0("logFC.", colnames(fit.con$coefficients)[i])
  colnames(test)[2] = paste0("raw.Pvalue.", colnames(fit.con$coefficients)[i])
  table = cbind(table, test)
}
table = table[-1]

# separate p value FDR adjustment
pvals = grep("^raw", names(table))
p = table[, pvals]
for (j in 1:ncol(p)) {
  o = !is.na(p[, j])
  p[o, j] = p.adjust(p[o, j], method = "fdr")
}
colnames(p) = paste0("FDR.", colnames(fit.con$coefficients))
table = cbind(table, p)

# highlight genes of interest
nr_genes = c("Abca1", "Abcg1", "Srebf1", "Fasn", "Apoe", "Aacac", 
             "Scd1")

table.vol = table[,c(grep("logFC", names(table)), 
                     grep("FDR", names(table)))]
table.vol = table.vol %>% mutate(ensembl_gene_id = rownames(.)) %>% 
  left_join(., genes, by = "ensembl_gene_id")
table.vol$FDR.status.27HC = ifelse(table.vol$FDR.HCvsVeh < 0.01, "FDR < 0.05",
                                   "n.s")
table.vol$FDR.status.GW3965 = ifelse(table.vol$FDR.GWvsVeh < 0.01, "FDR < 0.05",
                                     "n.s")
table.vol$nrlabel[table.vol$external_gene_name %in% c(nr_genes)] = table.vol$external_gene_name[table.vol$external_gene_name %in% c(nr_genes)]


ggplot(table.vol, aes(logFC.HCvsVeh, -log10(FDR.HCvsVeh),  
                      label = nrlabel)) + 
  geom_point(aes(color = FDR.status.27HC)) + 
  geom_text_repel(box.padding = unit(0.45, "lines"), hjust = 1) + 
  labs(x = "logFC", y = "-log10(FDR)") +
  theme(legend.title = element_blank(),
        text = element_text(size = 20), axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  scale_color_manual(values=c("lightpink3", "lightblue4")) +
  theme_classic()


ggplot(table.vol, aes(logFC.GWvsVeh, -log10(FDR.GWvsVeh),  
                      label = nrlabel)) + 
  geom_point(aes(color = FDR.status.GW3965)) + 
  geom_text_repel(box.padding = unit(0.45, "lines"), hjust = 1) + 
  labs(x = "logFC", y = "-log10(FDR)") +
  theme(legend.title = element_blank(),
        text = element_text(size = 20), axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  scale_color_manual(values=c("lightpink3", "lightblue4")) +
  theme_classic()

# Venn Diagram ------------------------------------------------------------
con.coded = decideTests(fit.con, p.value = 0.05, n = Inf, method = "separate") 

# 27HC vs. GW
vennDiagram(con.coded[,c(1, 2)], include = "up", 
            main = "upregulated genes", circle.col=c("darkred", "olivedrab3"),
            names = c("27HC", "GW3965"), cex = 1)

vennDiagram(con.coded[,c(1, 2)], include = "down", 
            main = "downregulated genes", circle.col=c("darkred", "olivedrab3"),
            names = c("27HC", "GW3965"), cex = 1)
