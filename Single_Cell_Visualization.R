
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



### Violin plot for marker genes

dataForViolinPlot <- data.frame(expr=as.numeric(t(sig.exprs)),
                                cell.type=rep(cell_type_label, length(markers.to.plot)),
                                gene=rep(markers.to.plot, each=length(cell_type_label)),
                                stringsAsFactors = F)

dataForViolinPlot$gene <- factor(dataForViolinPlot$gene, levels=markers.to.plot)

ggplot(data=dataForViolinPlot, aes(x=cell.type, y=expr)) +
  geom_violin(aes(fill=cell.type, color=cell.type), lwd=0.1, #size=0.5,
              outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
              outlier.fill = NA, width=0.8) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  scale_x_discrete(position = 'bottom') +
  #geom_jitter(size=0.01, width = 0.1) +
  facet_wrap(~gene, nrow=1, strip.position = 'bottom') +
  labs(x='', y=expression('Log'[2]*'(TP10K+1)')) +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = 'none') +
  #theme_set(theme_minimal()) #
  theme(axis.title=element_blank(),
        
        axis.text = element_text(color='black', size=14),
        axis.text.x = element_blank(),
        strip.text = element_text(angle = 45, size=12, 
                                  hjust=0.3),
        strip.background = element_blank()) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(panel.spacing = unit(0, "lines")) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))



########### Cell composition

dataForBarPlot$Encode_Blueprint_Subtype <- as.character(dataForBarPlot$Encode_Blueprint_Subtype)

cells.for.plot <- c('CD4+ Naive T cells','CD4+ Tcm','CD4+ Tem',
           'CD8+ Naive T cells','CD8+ Tcm','CD8+ Tem',
           'Tregs','NK cells','Naive B-cells','Memory B-cells',
           'Class-switched memory B-cells')

dataForBarPlot$Encode_Blueprint_Subtype


idx <- which(! dataForBarPlot$Encode_Blueprint_Subtype %in% cells.for.plot)

dataForBarPlot$Encode_Blueprint_Subtype[idx] <- 'Others'


idx <- which(dataForBarPlot$Encode_Blueprint_Subtype %in% cells.for.plot)
idx

unique(dataForBarPlot$Encode_Blueprint_Subtype)

dataForBarPlot$Encode_Blueprint_Subtype <- factor(dataForBarPlot$Encode_Blueprint_Subtype,
                                                  levels=rev(c(cells.for.plot,'Others')))

cellColor <- colorPanel22[1:12]
cellColor

#### stacked
ggplot(data=dataForBarPlot, aes(x=sample, y=proportion, 
                                fill=Encode_Blueprint_Subtype)) +
  geom_bar(stat='identity', width=0.8) + #coord_flip()
  #scale_y_discrete(limits=rev(unique(dataForBarPlot$Encode_Blueprint_Subtype))) +
  #ylim(0,50) +
  #geom_errorbar(aes(ymin=expr, ymax=expr+sd), width=.2, size=0.5, #expr-sd
  #              position=position_dodge(.9)) +
  labs(x='', y=expression('Proportion of Cells (%)')) +
  #scale_y_continuous(trans = 'sqrt',
  #                   breaks = c(0,2.5,50,250,750),
  #                   labels = c(0,2.5,50,250,750)) +
  #scale_y_sqrt() +
  #scale_y_continuous(trans='log2') +
  #scale_fill_manual(values = rep('black',nrow(dataForBarPlot))) +
  scale_fill_manual(values = cellColor) +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = 'right') +
  theme(axis.title=element_text(size=16),
        axis.text = element_text(color='black', size=14),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank()) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 1, unit = "cm"))

# 800*500



dataForBubblePlot <- data.frame(dataForBubblePlot, stringsAsFactors = F)
dataForBubblePlot

colnames(dataForBubblePlot) <- c('ifng', 'total', 'ratio', 'cell', 'sample', 'treatment', 'mean.all', 'mean.expressed')

dataForBubblePlot$total <- as.numeric(as.character(dataForBubblePlot$total))
dataForBubblePlot$ifng <- as.numeric(as.character(dataForBubblePlot$ifng))

dataForBubblePlot$mean.all <- round(as.numeric(as.character(dataForBubblePlot$mean.all)),2)
dataForBubblePlot$mean.expressed <- round(as.numeric(as.character(dataForBubblePlot$mean.expressed)),2)


dataForBubblePlot$ratio <- round(as.numeric(as.character(dataForBubblePlot$ratio))*100,2)
dataForBubblePlot

dataForBubblePlot$label <- paste(paste0(dataForBubblePlot$ratio, '%'), ' (', paste(dataForBubblePlot$ifng, dataForBubblePlot$total, sep=' / '), ')', sep='')
dataForBubblePlot$label

exprMax <- max(dataForBubblePlot$mean.all)
exprMax

dataForBubblePlot

cols = brewer.pal(4, "Reds")
cols = colorRampPalette(cols)(10)
col_fun = colorRampPalette(rev(c(cols[10],cols[1])), space = "Lab")(2)

ggplot(dataForBubblePlot, mapping=aes(x=cell, y=treatment, #y=-log10(Benjamini), #y=Fold.Enrichment
                                      color=mean.all,size=ratio)) +
  geom_point()+ #coord_flip() +
  scale_x_discrete(limits=unique(dataForBubblePlot$cell)) +
  scale_y_discrete(limits=rev(unique(dataForBubblePlot$treatment))) +
  #scale_x_discrete(limits=Order)+
  scale_colour_gradientn(limits=c(0,exprMax+0.07),
                         colors= c(col_fun[1],col_fun[2])) + #
  #facet_wrap(~Comparison) +
  #facet_grid(Regulation~Comparison) + # scales=free
  xlab('')+ylab('') +
  guides(size = guide_legend(order=2, title='Percent\nExpressed'),
         colour = guide_colourbar(order=1, title = 'Average\nExpression')) + #'P Value\n(Benjamini)'))
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.minor = element_blank(),
                   #panel.border = element_rect(colour='black'),
                   panel.border = element_blank(),
                   panel.background = element_blank()) +
  ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size=20)) +
  theme(axis.text=element_text(size=14, color='black'),
        axis.text.x =element_text(size=14, color='black', angle=90, hjust=1),
        axis.title=element_text(size=15)) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  theme(strip.text = element_text(size = 14),
        legend.key.size = unit(0.8,'cm'))
