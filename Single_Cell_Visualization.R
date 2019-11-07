
cell_type <- dataForTSNEPlot$Encode_Blueprint_Subtype
cell_type

levels(cell_type) <- seq(nlevels(cell_type)) #make cluster names short as numeric

cell_type_label <- paste(cell_type, dataForTSNEPlot$Encode_Blueprint_Subtype, sep = ':')

cell_type_label <- factor(cell_type_label, levels=str_sort(unique(cell_type_label),numeric = TRUE))


col.Cluster <- list(`Cell type`=colorspace::qualitative_hcl(nlevels(cell_type), palette="Dark 3"))
names(col.Cluster$`Cell type`) <- levels(cell_type)
#col.Signatures <- list(Signatures=colorspace::qualitative_hcl(nlevels(sig.ord$CellType), palette="Dark 3"))
#names(col.Signatures$Signatures) <- levels(sig.ord$CellType)
col.Cluster

sig.exprs <- seurat.integrated[['RNA']]@data[markers.to.plot,dataForTSNEPlot$barcode]
sig.exprs



#annoColors <- list(
#  `Sample Type`=c(Normal='darkolivegreen',
#                  Dysplasia='lightcoral'))


#topAnnotation = HeatmapAnnotation(`Sample Type`=phenoData[,'SampleType'],
#                                  `Histology Grade`=phenoData[,'HistologyGrade'],
#                                  col=annoColors,
#                                  simple_anno_size_adjust = TRUE,
#                                  #annotation_height = c(1,1),
#                                  height = unit(8, "mm"),
#                                  #summary = anno_summary(height = unit(4, "cm")),
#                                  show_legend = c("bar" = TRUE),
#                                  show_annotation_name = F)


ht <- Heatmap(as.matrix(sig.exprs), name='Log2(TP10K)', show_row_names = TRUE, show_column_names = FALSE, 
              top_annotation = HeatmapAnnotation(`Cell type` = cell_type, 
                                                 annotation_legend_param = 
                                                   list(`Cell type` = list(at = levels(cell_type),
                                                                           labels = levels(cell_type_label))),
                                                 col=col.Cluster,
                                                 show_legend=TRUE, 
                                                 show_annotation_name=FALSE),
              #top_annotation = HeatmapAnnotation(Cluster = cell_type, col=col.Cluster, show_legend=FALSE, show_annotation_name=FALSE),
              #left_annotation = rowAnnotation(Signatures = sig.ord$CellType, col=col.Signatures, show_legend=FALSE, show_annotation_name=FALSE),
              #left_annotation = rowAnnotation(Signatures = sig.ord$CellType, show_legend=FALSE, show_annotation_name=FALSE),
              row_names_side = "left", row_names_gp = gpar(fontsize = 10), 
              row_title_rot = 0, row_title_gp = gpar(fontsize = 9, fontface = 'bold'),
              column_title_gp = gpar(fontsize = 8),
              #row_split = sig.ord$CellType,
              column_split = cell_type,
              cluster_rows = FALSE, cluster_columns = FALSE, 
              show_column_dend = FALSE, show_row_dend = FALSE,
              cluster_column_slices = FALSE,
              show_heatmap_legend = TRUE,
              row_order = markers.to.plot,
              #heatmap_legend_param = list(legend_direction='horizontal'),
              col = colorRamp2(c(0, floor(max(sig.exprs))), c("grey95", "red")))
ht
