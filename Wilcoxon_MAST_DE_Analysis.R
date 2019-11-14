
####################################################################
####### STEP 7. Differential Gene Expression Analysis

seurat.integrated@meta.data$sample_cell <- paste(seurat.integrated@meta.data$sample,
                                                 seurat.integrated@meta.data$EncodeBlueprint_Subtype,
                                                 sep='_')
seurat.integrated@meta.data$sample_cell

feature.names <- readRDS('data/rData/feature_names.rds')
feature.names


######## Seurat

cells <- c("CD4+ Naive T cells","CD4+ Tem","CD4+ Tcm",
           "CD8+ Naive T cells", "CD8+ Tem","CD8+ Tcm",
           "Tregs","NK cells")

Idents(seurat.integrated) <- 'EncodeBlueprint_Subtype'
Idents(seurat.integrated)


#Idents(seurat.integrated) <- 'sample_cell'
#Idents(seurat.integrated)


######## wilcoxon

stimSample <- 'S2_GEX'
ctrlSample <- 'S1_GEX'

for (cell in cells) {
  print (cell)
  
  ident.1 <- paste(stimSample, cell, sep='_')
  ident.2 <- paste(ctrlSample, cell, sep='_')
  
  seurat.integrated.subset <- subset(seurat.integrated,
                                     subset = Encode_Blueprint_Subtype == cell)
  
  
  expr.data <- as.matrix(seurat.integrated.subset[['RNA']]@data)
  meta.data <- seurat.integrated.subset@meta.data
  
  cell1 <- rownames(meta.data)[meta.data$sample_cell==ident.1]
  cell2 <- rownames(meta.data)[meta.data$sample_cell==ident.2]
  
  
  #if (length(cell1) < 3 | length(cell2) < 3) {
  #  next
  #}
  
  #idx <- which(rowSums(expr.data>1) > ncol(expr.data)*0.05)
  
  idx <- which(rowSums(expr.data>0.5) > ncol(expr.data)*0.1)
  genes <- rownames(expr.data)[idx]
  
  #genes <- rownames(expr.data)
  
  
  mito <- grep(pattern = "^MT-", x = genes, value = FALSE)
  ribo <- grep(pattern = "^RPS|^RPL|^MRPL", x = genes, value = FALSE)
  
  genes <-genes[-c(mito, ribo)]
  
  group <- factor(meta.data$sample_cell, levels=c(ident.1, ident.2))
  
  degList <- pbmclapply(
    X = genes,
    FUN = function(gene) {wilcoxFun(expr.data, group, gene, cell1, cell2)},
    mc.cores = 4
  )
  
  degTable <- Reduce(rbind, degList)
  
  colnames(degTable) <- c(paste0('n_', c(stimSample, ctrlSample)),
                          paste0('mean_', c(stimSample, ctrlSample)),
                          't', 'pVal')
  
  degTable <- data.frame(degTable)
  
  rownames(degTable) <- genes#rownames(seurat.integrated@data)
  
  idx <- which(rowSums(is.na(degTable))>0 | rowSums(abs(degTable)==Inf)>0)
  
  if (length(idx) > 0) {
    degTable <- degTable[-idx,]
  }
  
  degTable$FDR <- p.adjust(degTable$pVal, method = 'BH')
  degTable$Bonferroni <- p.adjust(degTable$pVal, method = 'bonferroni')
  
  degTable$logFC_TP10K <- degTable[,3]-degTable[,4]
  
  degTable <- degTable[order(degTable$FDR, decreasing = F),]
  
  #cell <- gsub('/', '-', cell, fixed=T)
  
  degTable$Ensembl <- feature.names[rownames(degTable),]$V1
  
  fl = paste(stimSample, ctrlSample, cell, 'Wilcox_filter_TP10K0.5_Cell10_noMito_noRibo', sep='_')
  #fl = paste(stimSample, ctrlSample, cell, 'Wilcox_test_noFilter_noMito_noRibo', sep='_')
  write.table(degTable, file=paste0('report/', fl, '.txt'),
              sep='\t', quote=F)
  
}

View(degTable)


#idx <- which(rowSums(expr.data>0.5) > ncol(expr.data)*0.1)
#idx

#ribo <- grep(pattern = "^RPS|^RPL|^MRPL", x = genes, value = FALSE)
#ribo





### MAST

stimSample <- 'S2_GEX'
ctrlSample <- 'S1_GEX'

for (cell in cells) {
      print (cell)
      
      ident.1 <- paste(stimSample, cell, sep='_')
      ident.2 <- paste(ctrlSample, cell, sep='_')
      
      seurat.integrated.subset <- subset(seurat.integrated,
                                         subset = Encode_Blueprint_Subtype == cell)
      
      expr.data <- as.matrix(seurat.integrated.subset[['RNA']]@data)
      raw.data <- as.matrix(seurat.integrated.subset[['RNA']]@counts)
      meta.data <- seurat.integrated.subset@meta.data
      
      cell1 <- rownames(meta.data)[meta.data$sample_cell==ident.1]
      cell2 <- rownames(meta.data)[meta.data$sample_cell==ident.2]
      
      
      #if (length(cell1) < 3 | length(cell2) < 3) {
      #  next
      #}
      
      #idx <- which(rowSums(expr.data>1) > ncol(expr.data)*0.05)
      
      idx <- which(rowSums(expr.data>0.5) > ncol(expr.data)*0.1)
      
      genes <- rownames(expr.data)[idx]
      
      mito <- grep(pattern = "^MT-", x = genes, value = FALSE)
      ribo <- grep(pattern = "^RPS|^RPL|^MRPL", x = genes, value = FALSE)
      
      genes <-genes[-c(mito, ribo)]
      
      expr.data <- expr.data[genes,]
      raw.data <- raw.data[genes,]
      
      group <- factor(meta.data$sample, levels=c(stimSample, ctrlSample))
      names(group) <- rownames(meta.data)
      
      ###
      idx <- which(meta.data$sample %in% c(stimSample, ctrlSample))
      
      L <- list()
      
      L$count <- raw.data[,idx]
      L$tpm <- expr.data[,idx]
      L$condt <- group[idx]
      
      grp <- L$condt
      cdr <- scale(colMeans(L$count > 0))
      
      dge <- DGEList(counts=L$count)
      dge <- edgeR::calcNormFactors(dge)
      cpms <- edgeR::cpm(dge)
      
      sca <- FromMatrix(exprsArray = log2(cpms + 1), 
                        cData = data.frame(wellKey = names(grp), 
                                           grp = grp, cdr = cdr))
      zlmdata <- zlm(~cdr + grp, sca) # zlm.SingleCellAssay
      mast <- lrTest(zlmdata, "grp")
      
      degTable = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"],
                            row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
      
      genes <- rownames(degTable)
      
      exprSumList <- pbmclapply(
        X = genes,
        FUN = function(gene) {exprSummaryFun(expr.data, gene, cell1, cell2)},
        mc.cores = 4
      )
      
      
      exprSumTable <- Reduce(rbind, exprSumList)
      degTable <- cbind(exprSumTable, degTable)
      degTable <- cbind(length(cell2), degTable)
      degTable <- cbind(length(cell1), degTable)
      colnames(degTable)[1:4] <- c(paste0('n_', c(stimSample, ctrlSample)),
                                   paste0('mean_', c(stimSample, ctrlSample)))
      
      rownames(degTable) <- genes
      
      degTable$logFC_TP10K <- degTable[,3]-degTable[,4]
      
      degTable$Ensembl <- feature.names[rownames(degTable),]$V1
      
      degTable$FDR <- p.adjust(degTable$pval, method='BH')
      degTable$Bonferroni <- p.adjust(degTable$pval, method='bonferroni')
      
      
      tSumList <- pbmclapply(
        X = genes,
        FUN = function(gene) {ttestFun(expr.data, group, gene)},
        mc.cores = 4
      )
      
      tSumTable <- Reduce(rbind, tSumList)
      
      degTable$t <- tSumTable[,1]
      degTable$ttest.pVal <- tSumTable[,2]
      degTable$ttest.FDR <- p.adjust(degTable$ttest.pVal, method='BH')
      
      
      ### order
      degTable <- degTable[order(degTable$FDR, decreasing = F),]
      
      fl = paste(stimSample, ctrlSample, cell, 'MASTcpmDetRate_filter_TP10K0.5_Cell5_noMito_noRibo', sep='_')
      
      write.table(degTable, file=paste0('report/MAST/', fl, '.txt'),
                  sep='\t', quote=F)
      
}

#lrtGroup <- paste0('grp', stimSample)
#summaryCond <- summary(zlmdata, doLRT=lrtGroup)
#summaryDt <- summaryCond$datatable

#fcHurdle <- merge(summaryDt[contrast==lrtGroup & component=='logFC', .(primerid, coef)], #, ci.hi, ci.lo #logFC coefficients
#                  summaryDt[contrast==lrtGroup & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
#                  by='primerid')

#fcHurdle[,FDR:=p.adjust(`Pr(>Chisq)`, 'BH')]

#fcHurdle <- data.frame(fcHurdle, row.names = 1)

#colnames(fcHurdle) <- c('logFC','pval')
#fcHurdle$FDR <- p.adjust(fcHurdle$pval, method = 'BH')


## barplot

degForPlot <- c()

for (cell in cells) {
  
  print (cell)
  
  fl = paste(stimSample, ctrlSample, cell, 'Wilcox_filter_TP10K0.5_Cell10_noMito_noRibo', sep='_')
  fl <- file.path('report', paste0(fl, '.txt'))
  
  degTable <- read.table(fl, sep='\t', stringsAsFactors = F, header = T)
  
  upDEGs <- rownames(degTable[degTable$FDR<=0.05 & degTable$logFC>log2(1.2),])
  downDEGs <- rownames(degTable[degTable$FDR<=0.05 & degTable$logFC<log2(1.2)*-1,])
  
  degForPlot <- rbind(degForPlot, c(length(upDEGs), 'Up-regulated', cell))
  degForPlot <- rbind(degForPlot, c(length(downDEGs), 'Down-regulated', cell))
  
}

degForPlot <- data.frame(degForPlot)
colnames(degForPlot) <- c('Num','Regulation','Cell')
degForPlot$Num <- as.numeric(degForPlot$Num)
degForPlot$Regulation <- factor(degForPlot$Regulation, 
                                levels=c('Up-regulated','Down-regulated'))

degForPlot$Cell <- factor(degForPlot$Cell, cells)

ggplot(data=degForPlot, aes(x=Cell, y=Num, fill=Regulation)) + 
  geom_bar(stat='identity', width=.8, position=position_dodge()) + #coord_flip()
  #geom_text(aes(label=Num), hjust=0.5, vjust=-0.5, size=4, position=position_dodge(width=1)) +
  #geom_errorbar(aes(ymin=expr, ymax=expr+sd), width=.2, size=0.5, #expr-sd
  #              position=position_dodge(.9)) +
  labs(x='', y='Number of DEGs') + 
  ylim(0,30) +
  theme_bw()+
  #scale_y_continuous(trans = 'sqrt',
  #                   breaks = c(0,2.5,50,250,750),
  #                   labels = c(0,2.5,50,250,750)) +
  #scale_y_sqrt() +
  #scale_y_continuous(trans='log2') +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.spacing.x = unit(0.2, "cm"),
        legend.position = 'top') +
  theme(axis.title=element_text(size=16), 
        axis.text = element_text(color='black', size=14),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank()) #+
#theme(plot.margin =  unit(c(0.5,0.5,0.5,3), "cm"))



###
stimSample <- 'S2_GEX'
ctrlSample <- 'S1_GEX'


for (cell in cells) {
  
  print (cell)
  
  fl = paste(stimSample, ctrlSample, cell, 'Wilcox_filter_TP10K0.5_Cell10_noMito_noRibo', sep='_')
  fl <- file.path('report', paste0(fl, '.txt'))
  
  degTable <- read.table(fl, sep='\t', stringsAsFactors = F, header = T)
  
  upDEGs <- degTable[degTable$FDR<=0.05 & degTable$logFC>log2(1.2),]
  downDEGs <- degTable[degTable$FDR<=0.05 & degTable$logFC<log2(1.2)*-1,]
  
  fl = paste(stimSample, ctrlSample, 'Wilcox_filter_TP10K0.5_Cell10_noMito_noRibo', sep='_')
  
  degTable <- rbind(upDEGs, downDEGs)
  
  if (nrow(degTable)==0) {
    next
  } else {
    write.xlsx(degTable, file = file.path('report', paste0(fl, '.xlsx')), sheetName = cell, append = TRUE)
  }
  
}


### Visualization

### scatter plot

Idents(seurat.integrated) <- 'Encode_Blueprint_Subtype'
Idents(seurat.integrated)

stimSample <- 'S4_GEX'
ctrlSample <- 'S3_GEX'


for (cell in cells) {
  print (cell)
  
  seurat.integrated.subset <- subset(seurat.integrated,
                                     subset = sample %in% c(ctrlSample, stimSample) & 
                                       Encode_Blueprint_Subtype == cell)
  
  Idents(seurat.integrated.subset) <- 'sample'
  
  fl = paste(stimSample, ctrlSample, cell, 'Wilcox_filter_TP10K0.5_Cell10_noMito_noRibo', sep='_')
  fl <- file.path('report', paste0(fl, '.txt'))
  
  degTable <- read.table(fl, sep='\t', stringsAsFactors = F, header = T)
  
  avg <- AverageExpression(seurat.integrated.subset, show.progress = FALSE)
  
  avg.t.cells <- log2(avg$RNA+1)
  avg.t.cells
  
  #avg.t.cells <- log1p(AverageExpression(seurat.integrated.subset, verbose = FALSE)$RNA)
  
  avg.t.cells <- avg.t.cells[rownames(degTable),]
  colnames(avg.t.cells) <- c('ctrlSample', 'stimSample')
  
  avg.t.cells$gene <- rownames(avg.t.cells)
  avg.t.cells$gene <- gsub('[dup]','',avg.t.cells$gene, fixed=T)
  
  genes.to.label1 <- rownames(degTable[degTable$FDR<=0.05 & degTable$logFC_TP10K>log2(1.2),])
  genes.to.label2 <- rownames(degTable[degTable$FDR<=0.05 & degTable$logFC_TP10K<log2(1.2)*-1,])
  
  print (length(genes.to.label1))
  print (length(genes.to.label2))
  print (genes.to.label2)
  #genes.to.label1 <- rownames(deg[deg$p_val<=0.01 & deg$avg_logFC>0,])
  #genes.to.label2 <- rownames(deg[deg$p_val<=0.01 & deg$avg_logFC<0,])
  
  avg.t.cells$sig <- 'NS'
  if (length(genes.to.label1)>0) {
    avg.t.cells[genes.to.label1,]$sig <- 'UP'
  }
  
  if (length(genes.to.label2)>0) {
    avg.t.cells[genes.to.label2,]$sig <- 'DOWN'
  }
  
  avg.t.cells$sig <- factor(avg.t.cells$sig)
  
  nx <- degTable[1,2]
  ny <- degTable[1,1]
  
  if (length(genes.to.label1)>0 & length(genes.to.label2)>0) {
    val <- c("green3", "black", "red")
  } else if (length(genes.to.label1)>0 & length(genes.to.label2)==0) {
    val <- c("black", "red")
  } else if (length(genes.to.label1)==0 & length(genes.to.label2)>0) {
    val <- c("green3", "black")
  } else if (length(genes.to.label1)==0 & length(genes.to.label2)==0) {
    val <- c("black")
  }
  
  
  p1 <- ggplot(avg.t.cells, aes(ctrlSample, stimSample)) + geom_point(aes(color=sig), size=0.5) + ggtitle(cell) +
    geom_point(data = subset(avg.t.cells, sig == 'UP'),
               aes(ctrlSample, stimSample, color=sig), size=0.5) +
    geom_point(data = subset(avg.t.cells, sig == 'DOWN'),
               aes(ctrlSample, stimSample, color=sig), size=0.5) +
    
    scale_color_manual(values = val) +
    geom_abline(intercept = 0, slope=1, color='blue', linetype='dashed')+
    labs(x=paste0(ctrlSample, ' (N=', nx, ')'),y=paste0(stimSample, ' (N=', ny, ')')) +
    
    geom_text_repel(data=subset(avg.t.cells, avg.t.cells$gene %in% genes.to.label1), 
                    aes(label=gene), segment.alpha = 0.4,size = 3, color='red',
                    nudge_y = 1) +
    geom_text_repel(data=subset(avg.t.cells, avg.t.cells$gene %in% genes.to.label2), 
                    aes(label=gene), segment.alpha = 0.4,size = 3, color='darkgreen',
                    nudge_x = 1) +
    theme(legend.position = 'none') +
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.text = element_text(size=14),
          legend.position = 'none') +
    theme(axis.title=element_text(size=16),
          axis.text = element_text(color='black', size=14)) +
    theme(plot.title = element_text(hjust = 0.5, size=16, face='bold')) +
    theme(axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major = element_blank()) #+
  #theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 1, unit = "cm"))
  
  
  png(filename = gsub('.txt','.repel.png', fl), width = 500, height = 500, res = 100)
  print (p1)
  dev.off()
  
  
}


library(ggrepel)




## volcano

logFcThreshold <- log2(1.2)
logFcThreshold <- log2(1.2)

for (cell in cells) {
  
  print (cell)
  
  fl = paste(stimSample, ctrlSample, cell, 'Wilcox_test_filter_TP10K0.5_Cell10_noMito_noRibo', sep='_')
  fl <- file.path('report', paste0(fl, '.txt'))
  
  degTable <- read.table(fl, sep='\t', stringsAsFactors = F, header = T)
  
  rownames(degTable) <- gsub('[dup]',' ',rownames(degTable), fixed=T)
  
  genes.to.label1 <- rownames(degTable[degTable$FDR<=0.05 & degTable$logFC>log2(1.2),])
  genes.to.label2 <- rownames(degTable[degTable$FDR<=0.05 & degTable$logFC<log2(1.2)*-1,])
  
  avg.t.cells <- degTable
  
  avg.t.cells$sig <- 'NS'
  if (length(genes.to.label1)>0) {
    avg.t.cells[genes.to.label1,]$sig <- 'UP'
  }
  
  if (length(genes.to.label2)>0) {
    avg.t.cells[genes.to.label2,]$sig <- 'DOWN'
  }
  
  avg.t.cells$sig <- factor(avg.t.cells$sig)
  
  avg.t.cells$gene <- rownames(avg.t.cells)
  
  nx <- degTable[1,2]
  ny <- degTable[1,1]
  
  
  if (length(genes.to.label1)>0 & length(genes.to.label2)>0) {
    val <- c("green3", "black", "red")
  } else if (length(genes.to.label1)>0 & length(genes.to.label2)==0) {
    val <- c("black", "red")
  } else if (length(genes.to.label1)==0 & length(genes.to.label2)>0) {
    val <- c("green3", "black")
  } else if (length(genes.to.label1)==0 & length(genes.to.label2)==0) {
    val <- c("black")
  }
  
  
  p1 <- ggplot(avg.t.cells, aes(x = logFC, y = -log10(FDR))) +
    geom_point(aes(color=sig), alpha=1, size=0.6) +
    geom_vline(xintercept = c(-logFcThreshold, logFcThreshold), 
               color='darkgreen', linetype='dashed') +
    geom_hline(yintercept = -log10(0.05), color='darkgreen',linetype='dashed')+
    scale_color_manual(values = val) +
    xlab(expression('log'[2]*'(Fold Change)')) + ylab(expression('-log'[10]*'(Adjusted P Value)')) +
    ylim(0,7) +
    
    geom_text_repel(data=subset(avg.t.cells, avg.t.cells$gene %in% genes.to.label1), 
                    aes(label=gene), segment.alpha = 0.4,size = 3, color='red') +
    geom_text_repel(data=subset(avg.t.cells, avg.t.cells$gene %in% genes.to.label2), 
                    aes(label=gene), segment.alpha = 0.4,size = 3, color='darkgreen') +
    
    
    #geom_text_repel(data=subset(pairwiseComparisonStatistics.forPlot, pairwiseComparisonStatistics.forPlot$Gene %in% ovlp), 
    #                aes(label=Gene), segment.alpha = 0.4,size = 3, color='blue') +
    #geom_text_repel(data = subset(comparisonData, abs(ddCt) >= log2threshold & -log10(Pvalue) >= -log10(pvalThreshold)), segment.alpha = 0.4,
    #                aes(label = Gene), size = annotSize) +
    #geom_hline(yintercept = -log10(adjPvalThreshold), color = 'red', linetype = 'dashed') +
    #geom_vline(xintercept = -logFcThreshold, color = 'red', linetype = 'dashed') +
    #geom_vline(xintercept = logFcThreshold, color = 'red', linetype = 'dashed') +
    #theme(panel.background = element_rect(color = 'black', size = 2)) +
    
  #theme_bw() +
  #theme(axis.line = element_line(colour = "black"),
  #      panel.grid.major = element_blank(),
  #      panel.grid.minor = element_blank(),
  #      panel.border = element_rect(colour='black'),
  #      panel.background = element_blank()) +
  theme(legend.position="none") +
    theme(axis.text=element_text(size=14), 
          axis.title=element_text(size=16),
          strip.text = element_text(size=16))
  
  
  png(filename = gsub('.txt','.Volcano.png', fl), width = 500, height = 500, res = 100)
  print (p1)
  dev.off()
}



########## bubble plot

stimSample <- 'S2_GEX'
ctrlSample <- 'S1_GEX'

cells <- c('CD4+ Naive T cells','CD4+ Tcm','CD4+ Tem',
           'CD8+ Naive T cells','CD8+ Tcm','CD8+ Tem',
           'Tregs','NK cells')



deg <- c()
for (cell in cells) {
  
  print (cell)
  
  fl = paste(stimSample, ctrlSample, cell, 'Wilcox_filter_TP10K0.5_Cell10_noMito_noRibo', sep='_')
  fl <- file.path('report', paste0(fl, '.txt'))
  
  degTable <- read.table(fl, sep='\t', stringsAsFactors = F, header = T)
  
  degTable$Gene <- rownames(degTable)
  degTable$Cell <- cell
  
  deg <- rbind(deg, degTable)
  
}

deg$Cell <- factor(deg$Cell, levels=cells)

#deg <- deg[with(deg, order(Cell, logFC_TP10K, FDR)),]
#deg <- deg[with(deg, order( decreasing = T)),]
deg


degForPlot <- deg

degForPlot <- degForPlot[degForPlot$FDR<=0.05 & degForPlot$logFC_TP10K<=-log2(1.2),]

#degForPlot <- degForPlot[degForPlot$FDR<=0.05 & degForPlot$logFC_TP10K>=log2(1.2),]
degForPlot

#degForPlot <- degForPlot[degForPlot$p_val_adj<=0.05,]
#degForPlot[order(degForPlot$avg_logFC, decreasing = T),]

dim(degForPlot)

pval <- 0.05

p = ggplot(data=degForPlot, mapping=aes(x=Gene, y=Cell, #y=-log10(Benjamini), #y=Fold.Enrichment
                                        color=FDR,size=abs(logFC_TP10K)))
p + geom_point() + coord_flip() + 
  scale_x_discrete(limits=rev(unique(degForPlot$Gene))) +  # REV
  #scale_x_discrete(limits=unique(degForPlot$gene)) +
  scale_colour_gradientn(limits=c(0,pval),
                         colors= c("red","yellow","green")) + #
  #facet_wrap(~Comparison) +
  #facet_grid(Regulation~Comparison) + # scales=free
  xlab('')+ylab('') +
  guides(shape = guide_legend(order=1, title = 'Absolute\nlogFC'),
         colour = guide_colourbar(order=1, title = 'FDR'),
         size = guide_legend(order=2, title = 'Absolute\nlogFC')) + 
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour='black'),
                   panel.background = element_blank()) +
  ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size=20)) +
  theme(axis.text=element_text(size=14, color='black'), 
        axis.text.x =element_text(size=16, color='black', angle = 45, hjust=1),
        axis.title=element_text(size=15)) + 
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14), 
        legend.key.size = unit(0.8,'cm'))



#####

genes <- unique(degForPlot$Gene)
genes


celltypes <- unique(degForPlot$Cell)

fcMatrix <- matrix(rep(NA,length(genes)*length(celltypes)), nrow=length(genes), ncol=length(celltypes))
fcMatrix

rownames(fcMatrix) <- genes
colnames(fcMatrix) <- celltypes

sigMatrix <- matrix(rep(NA,length(genes)*length(celltypes)), nrow=length(genes), ncol=length(celltypes))
sigMatrix

rownames(sigMatrix) <- genes
colnames(sigMatrix) <- celltypes


genes

deg <- c()
for (cell in celltypes) {
  
  print (cell)
  
  fl = paste(stimSample, ctrlSample, cell, 'Wilcox_filter_TP10K0.5_Cell10_noMito_noRibo', sep='_')
  fl <- file.path('report', paste0(fl, '.txt'))
  
  degTable <- read.table(fl, sep='\t', stringsAsFactors = F, header = T)
  
  fcMatrix[genes, cell] <- degTable[genes,'logFC_TP10K']
  sigMatrix[genes, cell] <- ifelse(degTable[genes,'FDR'] <= 0.05, '*', '')
  
}
fcMatrix
sigMatrix





library(circlize)
library(ComplexHeatmap)
col_fun = colorRampPalette(rev(c("red",'white','blue')), space = "Lab")(100)

col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

ht <- Heatmap(fcMatrix,
              #name = 'Expression',
              
              # COLOR
              #col = colorRampPalette(rev(c("red",'white','blue')), space = "Lab")(100),
              col=col_fun,
              na_col = 'grey',
              rect_gp = gpar(col = "grey", lwd = 0.5),
              
              # MAIN PANEL
              column_title = NULL,
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              column_names_rot = 90,
              column_names_gp = gpar(fontsize = 10),
              row_names_gp = gpar(fontsize = 10),
              #column_names_max_height = unit(3, 'cm'),
              #column_split = factor(phenoData$Day,
              #                      levels=str_sort(unique(phenoData$Day), numeric = T)),
              
              #column_order = rownames(phenoData),
              
              # ANNOTATION
              #top_annotation = topAnnotation,
              
              # ADD TEXT
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(sigMatrix[i, j], x, y, gp = gpar(fontsize = 10, col = "black", fill = 'black'))
              },
              
              # LEGEND
              #heatmap_legend_param = list(
              #  at = c(-1.5, -1, 0, 1, 1.5),
              #labels = c("Negative",'Positive'),
              #  title = "",
              #  title_position = 'leftcenter-rot',
              #  legend_height = unit(3, "cm"),
              #  adjust = c("right", "top")
              #)
)

draw(ht,annotation_legend_side = "right",row_dend_side = "left")
