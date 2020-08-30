
dataForBarPlot <- annotation %>% group_by(celltype) %>%
  summarise(freq=length(celltype))

o <- order(dataForBarPlot$freq, decreasing = T)

dataForBarPlot$celltype <- factor(dataForBarPlot$celltype, levels=dataForBarPlot$celltype[o])


ggplot(data=dataForBarPlot, aes(x=celltype, y=freq, fill=celltype, color=celltype)) +
  geom_bar(stat='identity', width=.6) + #coord_flip()
  #geom_errorbar(aes(ymin=expr, ymax=expr+sd), width=.2, size=0.5, #expr-sd
  #              position=position_dodge(.9)) +
  labs(x='', y=expression('Number of Cells')) +
  #scale_y_continuous(trans = 'sqrt',
  #                   breaks = c(0,2.5,50,250,750),
  #                   labels = c(0,2.5,50,250,750)) +
  #scale_y_sqrt() +
  #scale_y_continuous(trans='log2') +
  #scale_fill_manual(values = rep('black',nrow(dataForBarPlot))) +
  #scale_color_manual(values = rep('black',nrow(dataForBarPlot))) +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = 'none') +
  theme(axis.title=element_text(size=16),
        axis.text = element_text(color='black', size=14),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank()) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))




cell.type <- factor(annotation$celltype, levels=unique(annotation$celltype))
cell.type

markers.to.plot <- idx


dataForViolinPlot <- data.frame(expr=as.numeric(t(sig.exprs)),
                                cell.type=rep(cell.type, length(markers.to.plot)), #cell_type_label
                                gene=rep(markers.to.plot, each=length(cell.type)),
                                stringsAsFactors = F)

dataForViolinPlot$gene <- factor(dataForViolinPlot$gene, levels=markers.to.plot) 


# ggplot(data=dataForViolinPlot, aes(x=cell.type, y=expr)) +
#   geom_violin(aes(fill=cell.type, color=cell.type), lwd=0.1, #size=0.5,
#               outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
#               outlier.fill = NA, width=0.8) +
#   coord_flip() +
#   scale_y_continuous(position = "right") +
#   scale_x_discrete(position = 'bottom') +
#   #geom_jitter(size=0.01, width = 0.1) +
#   facet_wrap(~gene, nrow=1, strip.position = 'bottom') +
#   labs(x='', y=expression('Log'[2]*'(TP10K+1)')) +
#   theme_bw()+
#   theme(legend.title = element_blank(),
#         legend.position = 'none') +
#   #theme_set(theme_minimal()) #
#   theme(axis.title=element_blank(),
#         
#         axis.text = element_text(color='black', size=14),
#         axis.text.x = element_blank(),
#         strip.text = element_text(angle = 45, size=12, 
#                                   hjust=0.3),
#         strip.background = element_blank()) +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank()) +
#   theme(panel.spacing = unit(0, "lines")) +
#   theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))



dataForBubblePlot <- dataForViolinPlot %>% group_by(gene, cell.type) %>% 
  summarise(mean.all=mean(expr), 
            mean.expressed=ifelse(sum(expr>0)>0, sum(expr)/sum(expr>0), 0),
            percent.expressed=sum(expr>0)/length(expr)*100)


exprMax <- max(dataForBubblePlot$mean.all)
exprMax

dataForBubblePlot$mean.expressed

dataForBubblePlot

library(RColorBrewer)

cols = brewer.pal(4, "Reds")
cols = colorRampPalette(cols)(10)
col_fun = colorRampPalette(rev(c(cols[10],cols[1])), space = "Lab")(2)

#col_fun = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
#col_fun

dataForBubblePlot

ggplot(dataForBubblePlot, mapping=aes(x=cell.type, y=gene, #y=-log10(Benjamini), #y=Fold.Enrichment
                                      color=mean.all,size=percent.expressed)) +
  geom_point()+ #coord_flip() +
  scale_x_discrete(limits=unique(dataForBubblePlot$cell.type)) +
  scale_y_discrete(limits=rev(unique(dataForBubblePlot$gene))) +
  #scale_x_discrete(limits=Order)+
  scale_colour_gradientn(limits=c(0,exprMax+0.07),
                         colors= c(col_fun[1],col_fun[2])) + #col_fun
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
        axis.text.x =element_text(size=14, color='black', angle=45, hjust = 1),
        axis.title=element_text(size=15)) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  theme(strip.text = element_text(size = 14),
        legend.key.size = unit(0.8,'cm'))
