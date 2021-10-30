library(Seurat)
library(tidyverse)

cd11b2 = readRDS("~/scRNAseq/results/cd11b2_annotated.RDS")

Idents(cd11b2) = cd11b2@meta.data$manual_0.3_cd11b
cd11b2$celltype_treat = paste(Idents(cd11b2), cd11b2$treatment, sep = "_")
ident_test = unique(cd11b2$celltype_treat)
vehicle_test = ident_test[grepl("_Vehicle$", ident_test)] %>% sort(.)
HC_test = ident_test[grepl("_27HC$", ident_test)] %>% sort(.)
gw3965_test = ident_test[grepl("_GW3965$", ident_test)] %>% sort(.)
cyp_test = ident_test[grepl("_GW273297X$", ident_test)] %>% sort(.)

Idents(cd11b2) = "celltype_treat"
marker_list_27hc = list()
for (i in 1:9) {
  marker = FindMarkers(cd11b2, ident.1 = HC_test[i], ident.2 = vehicle_test[i], 
                       verbose = FALSE, logfc.threshold = 0, min.pct = 0)
  marker_list_27hc[[i]] = marker
  names(marker_list_27hc)[i] = paste0(HC_test[i], "vsVeh")
}

marker_list_gw = list()
for (i in 1:9) {
  marker = FindMarkers(cd11b2, ident.1 = gw3965_test[i], ident.2 = vehicle_test[i], 
                       verbose = FALSE, logfc.threshold = 0, min.pct = 0)
  marker_list_gw[[i]] = marker
  names(marker_list_gw)[i] = paste0(gw3965_test[i], "vsVeh")
}

marker_list_cyp = list()
for (i in 1:9) {
  marker = FindMarkers(cd11b2, ident.1 = cyp_test[i], ident.2 = vehicle_test[i], 
                       verbose = FALSE, logfc.threshold = 0, min.pct = 0)
  marker_list_cyp[[i]] = marker
  names(marker_list_cyp)[i] = paste0(marker_list_cyp[i], "vsVeh")
}

rm(cd11b2)

save.image("~/scRNAseq/results/cd11b2_markers.RData")
