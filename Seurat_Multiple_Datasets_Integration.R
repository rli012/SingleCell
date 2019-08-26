setwd('C:\\Users/rli3/Desktop/share/BioinfoHub/')

#library(devtools)
#devtools::install_github(repo = "satijalab/seurat", ref = "develop",
#                         INSTALL_opts = c('--no-lock'))

library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(reshape2)
source('Single_Cell_Functions.R')


#################################################################
######## STEP 1. Read Raw Data

fds <- list.dirs('data/CellRanger')
fds <- fds[grep('filtered_feature_bc_matrix', fds)]
fds

num.samples <- length(fds)

metadata <- c()
mat.list <- list()

for (i in 1:num.samples) {
  matrixDir <- fds[i]
  sam = strsplit(matrixDir, '/')[[1]][3]
  print (sam)
  
  mat.list[[sam]] <- Read10XRL(data.dir = matrixDir)
  colnames(mat.list[[sam]]) <- gsub('(\\w+)-1', paste0(i, '_\\1'), colnames(mat.list[[sam]]))
  
  metadata <- rbind(metadata, data.frame(barcode=colnames(mat.list[[sam]]), sample=sam))

}

metadata[1:5,]



##################################################################
####### STEP 2. Create Seurat Object, Cell QC, Data Normalization

min.cells <- 5
seurat.list <- list()

for (i in 1:length(mat.list)) {
  seurat.list[[i]] <- CreateSeuratObject(counts = mat.list[[i]], min.cells = min.cells)
  seurat.list[[i]]@meta.data$sample <- names(mat.list)[i]
  seurat.list[[i]]$percent_mt <- PercentageFeatureSet(seurat.list[[i]], pattern = "^MT-")
  
  p <- VlnPlot(seurat.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3)
  
  png(filename = paste0('report/Cell_QC_', names(mat.list)[i], '.png'), width = 700, height = 500)
  print (p)
  dev.off()
  
  seurat.list[[i]] <- subset(seurat.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mt < 5)
  
  message(paste0(ncol(seurat.list[[i]][['RNA']]@counts), '/', ncol(mat.list[[i]]), ' of cells are kept'))
  
  seurat.list[[i]] <- NormalizeData(seurat.list[[i]], normalization.method = "LogNormalize",
                                    scale.factor = 10000,verbose = TRUE)
  seurat.list[[i]] <- FindVariableFeatures(seurat.list[[i]], selection.method = "vst", 
                                           nfeatures = 2000, verbose = TRUE)
  
  seurat.list[[i]] <- SCTransform(seurat.list[[i]], vars.to.regress = "percent_mt", 
                                  verbose = TRUE)
  
}
names(seurat.list) <- names(mat.list)
seurat.list[[i]]@meta.data[1:5,]


dataForBarPlot <- c()
for (i in 1:num.samples) {
  dataForBarPlot <- rbind(dataForBarPlot, 
                          c(ncol(mat.list[[i]]), ncol(seurat.list[[i]][['RNA']]@counts)))
}
dataForBarPlot


dataForBarPlot <- data.frame(freq=unlist(c(dataForBarPlot)),
                             group=rep(c('Before QC','After QC'), each=num.samples),
                             sample=rep(names(mat.list), num.samples))

dataForBarPlot$group <- factor(dataForBarPlot$group, levels=c('Before QC','After QC'))


ggplot(data=dataForBarPlot, aes(x=sample, y=freq, 
                                fill=group)) +
  geom_bar(stat='identity', width=0.8, position='dodge') + #coord_flip()
  ylim(0,6000) +
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
        legend.text = element_text(size=12),
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


####
names(seurat.list)
DefaultAssay(seurat.list[[1]])

reference.list <- seurat.list[c('S1_GEX','S2_GEX')]
seurat.features <- SelectIntegrationFeatures(object.list = reference.list, nfeatures = 3000)
reference.list <- PrepSCTIntegration(object.list = reference.list, anchor.features = seurat.features, 
                                  verbose = TRUE)

seurat.anchors <- FindIntegrationAnchors(object.list = reference.list, normalization.method = "SCT", 
                                       anchor.features = seurat.features)
seurat.integrated <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT")



#================================================================================
###
reference.list <- seurat.list[c('S1_GEX','S2_GEX')]
seurat.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
seurat.integrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:30)
seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)

#?IntegrateData()
?MergeSeurat()

### No Batch Correction
#seurat.integrated <- merge(x=seurat.list[[1]], 
#                           y=c(seurat.list[[2]],seurat.list[[3]],seurat.list[[4]],
#                               seurat.list[[5]],seurat.list[[6]],seurat.list[[7]],
#                               seurat.list[[8]],seurat.list[[9]],seurat.list[[10]]),
#                           merge.data=FALSE)

#seurat.integrated <- NormalizeData(seurat.integrated, verbose = TRUE)
#seurat.integrated

#================================================================================



# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
names(seurat.integrated)

dim(seurat.integrated[['RNA']]@data)
dim(seurat.integrated[['integrated']]@scale.data)

seurat.integrated[['RNA']]@data[1:10,1:10]
seurat.integrated[['RNA']]@data[1:10,1:10] == seurat.list[[1]][['RNA']]@data[1:10,1:10]

DefaultAssay(seurat.integrated)
if (DefaultAssay(seurat.integrated) != 'integrated') {
  DefaultAssay(seurat.integrated) <- "integrated"
}


# Run the standard workflow for visualization and clustering
seurat.integrated <- RunPCA(seurat.integrated, npcs = 30, verbose = FALSE)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat.integrated), 5)
seurat.integrated@meta.data

seurat.integrated <- RunTSNE(seurat.integrated, reduction = "pca", dims = 1:30, seed.use = '777')
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:30, seed.use = '777')

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.5)

seurat.integrated@meta.data[1:5,]

p1 <- DimPlot(seurat.integrated, reduction = "tsne", group.by = "sample")
p2 <- DimPlot(seurat.integrated, reduction = "tsne", group.by = "seurat_clusters")
#p2 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
#              repel = TRUE) + NoLegend()
plot_grid(p1, p2)

DimPlot(seurat.integrated, reduction = "tsne", split.by = "sample")

saveRDS(seurat.integrated, file='data/rData/Seurat_Integrated.rds')



####################################################################
####### STEP 5. Find Cluster Biomarkers

DefaultAssay(seurat.integrated) <- "RNA"

cell.markers <- FindConservedMarkers(seurat.integrated, ident.1 = 6, 
                                     grouping.var = "sample", 
                                     verbose = TRUE)
cell.markers[1:20,]

idx <- which(!grepl('^RPL|^RPS|^MT',rownames(cell.markers)))
head(cell.markers[idx,])

markers.to.plot <- c('CD3E','CD4','CD8B','CCR7',
                     'SELL','IL7R','FOXP3','GNLY',
                     'NKG7','TRAC', 'TRDC', 'TRGC1')

FeaturePlot(object = seurat.integrated,
            features = markers.to.plot,
            min.cutoff = "q9", 
            cols = c("lightgrey", "blue"),
            ncol = 4,
            pt.size = 0.1)

seurat.integrated@meta.data
#Idents(seurat.integrated) <- seurat.integrated@meta.data$seurat_clusters

DotPlot(seurat.integrated, 
        features = rev(markers.to.plot), 
        cols = c('blue', 'red'), 
        dot.scale = 8)+ #, 
  #split.by = "sample") + 
  RotatedAxis()



#seurat.integrated <- RenameIdents(seurat.integrated, 
#                                `0` = "CD14 Mono", 
#                                `1` = "CD4 Naive T", 
#                                `2` = "CD4 Memory T", 
#                                `3` = "CD16 Mono", 
#                                `4` = "B", 
#                                `5` = "CD8 T", 
#                                `6` = "T activated", 
#                                `7` = "NK", 
#                                `8` = "DC", 
#                                `9` = "B Activated", 
#                                `10` = "Mk", 
#                                `11` = "pDC")

###
new.cluster.ids <- c("CD14 Mono", "CD4 Naive T", "CD4 Memory T",
                     "T activated", "CD8 T", 'Tregs', 'NK',  
                     "CD16 Mono", "B",  
                     "NK", "DC", "B Activated")

levels(seurat.integrated)

names(new.cluster.ids) <- levels(seurat.integrated)
seurat.integrated <- RenameIdents(seurat.integrated, new.cluster.ids)


seurat.integrated@meta.data$cell_type <- seurat.integrated@active.ident #Idents(seurat.integrated)
seurat.integrated@meta.data

DimPlot(seurat.integrated, label = TRUE, label.size = 5, pt.size = 1) + NoLegend() +
  theme(axis.title = element_text(size=16, face='bold'))



###
cells <- c("CD14 Mono", "CD4 Naive T", "CD4 Memory T", 
           "CD16 Mono", "B", "CD8 T", "T activated", 
           "NK", "DC", "B Activated", "Mk", "pDC")

Idents(seurat.integrated) <- factor(Idents(seurat.integrated), 
                                    levels = cells)

DotPlot(seurat.integrated, 
        features = rev(markers.to.plot), 
        cols = c('blue', 'red'), 
        dot.scale = 8)+ #, 
        #split.by = "sample") + 
  RotatedAxis()



####################################################################
####### STEP 6. SingleR Annotation

devtools::install_github('dviraran/SingleR')

### Single cell
singlerCellAnnotation = CreateSinglerObject(seurat.integrated[['RNA']]@counts, annot = NULL, project.name='PBMC', 
                                        min.genes = 0, technology = "10X", species = "Human", citation = "",
                                        ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                                        fine.tune = T, do.signatures = T, clusters = NULL, do.main.types = T, 
                                        reduce.file.size = T, numCores = 4)

#saveRDS(singlerAnnotation, file='data/rData/singlerAnnotation.rds')

enblueAnnotation <- data.frame(subtype=singlerCellAnnotation$singler[[2]]$SingleR.single$labels,
                               maintype=singlerAnnotation$singler[[2]]$SingleR.single.main$labels)


hpcaAnnotation <- data.frame(subtype=singlerAnnotation$singler[[1]]$SingleR.single$labels,
                             maintype=singlerAnnotation$singler[[1]]$SingleR.single.main$labels)

rownames(hpcaAnnotation) == rownames(enblueAnnotation)

barcode <- rownames(seurat.integrated@meta.data)

seurat.integrated@meta.data[, 'EncodeBlueprint_Subtype'] <- capitalizeRL(as.character(enblueAnnotation[barcode,'subtype']))
seurat.integrated@meta.data[, 'EncodeBlueprint_Maintype'] <- capitalizeRL(as.character(enblueAnnotation[barcode,'maintype']))

seurat.integrated@meta.data[, 'HPCA_Subtype'] <- capitalizeRL(as.character(hpcaAnnotation[barcode,'subtype']))
seurat.integrated@meta.data[, 'HPCA_Maintype'] <- capitalizeRL(as.character(hpcaAnnotation[barcode,'maintype']))

## TO DO
## Define cells


original.cell <- c('CD4+ T-cells', 'CD8+ T-cells', 'naive B-cells')
target.cell <- c('CD4+ Naive T cells', 'CD8+ Naive T cells', 'Naive B cells')

seurat.integrated@meta.data$EncodeBlueprint_Subtype <- RenameCellRL(seurat.integrated@meta.data$EncodeBlueprint_Subtype, 
                                                                    original.cell, target.cell)

seurat.integrated@meta.data$HPCA_Subtype <- RenameCellRL(seurat.integrated@meta.data$HPCA_Subtype, 
                                                         original.cell, target.cell)


### Cluster
clusters <- seurat.integrated@meta.data$seurat_clusters
singlerClusterAnnotation = CreateSinglerObject(seurat.integrated[['RNA']]@counts, annot = NULL, project.name='PBMC', 
                                               min.genes = 0, technology = "10X", species = "Human", citation = "",
                                               ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                                               fine.tune = T, do.signatures = T, clusters = clusters, 
                                               do.main.types = T, reduce.file.size = T, numCores = 4)

# BLUEPRINT+ENCODE
original.cell <- 1:length(singlerClusterAnnotation$singler[[2]]$SingleR.clusters$labels)-1
original.cell

target.cell <- singlerClusterAnnotation$singler[[2]]$SingleR.clusters$labels[,1]
target.cell

seurat.integrated@meta.data$EncodeBlueprint_Subtype_Cluster <- RenameCellRL(as.numeric(as.character(clusters)), 
                                                                            original.cell, target.cell)

# HPCA
original.cell <- 1:length(singlerClusterAnnotation$singler[[1]]$SingleR.clusters$labels)-1
original.cell

target.cell <- singlerClusterAnnotation$singler[[1]]$SingleR.clusters$labels[,1]
target.cell

seurat.integrated@meta.data$HPCA_Subtype_Cluster <- RenameCellRL(as.numeric(as.character(clusters)), 
                                                                  original.cell, target.cell)

# 
original.cell <- c('CD4+ T-cells', 'CD8+ T-cells', 'naive B-cells')
target.cell <- c('CD4+ Naive T cells', 'CD8+ Naive T cells', 'Naive B cells')

seurat.integrated@meta.data$EncodeBlueprint_Subtype_Cluster <- RenameCellRL(seurat.integrated@meta.data$EncodeBlueprint_Subtype_Cluster, 
                                                                            original.cell, target.cell)

seurat.integrated@meta.data$HPCA_Subtype_Cluster <- RenameCellRL(seurat.integrated@meta.data$HPCA_Subtype_Cluster, 
                                                                 original.cell, target.cell)



p1 <- DimPlot(seurat.integrated, reduction = "tsne", 
              group.by = "EncodeBlueprint_Subtype_Cluster", 
              label = TRUE)

p2 <- DimPlot(seurat.integrated, reduction = "tsne", 
              group.by = "seurat_clusters",
              label = TRUE)

plot_grid(p1, p2)

p3 <- DimPlot(seurat.integrated, reduction = "tsne", 
              group.by = "HPCA_Subtype_Cluster", 
              label = TRUE)

p4 <- DimPlot(seurat.integrated, reduction = "tsne", 
              group.by = "seurat_clusters",
              label = TRUE)

plot_grid(p3, p4)


##### Visualization

## Cell type composition
annotation <- data.frame(seurat.integrated@meta.data, 
                         seurat.integrated@reductions$tsne@cell.embeddings,
                         stringsAsFactors = F)
annotation$barcode <- rownames(annotation)
rownames(seurat.integrated@meta.data) == rownames(seurat.integrated@reductions$tsne@cell.embeddings)

total <- nrow(annotation)
total

dataForBarPlot <- annotation %>% group_by(EncodeBlueprint_Subtype) %>%
  summarise(freq=length(EncodeBlueprint_Subtype), total=total, proportion=freq/total*100)

annotation

o <- order(dataForBarPlot$freq, decreasing = T)

dataForBarPlot$EncodeBlueprint_Subtype <- factor(dataForBarPlot$EncodeBlueprint_Subtype, 
                                                  levels=dataForBarPlot$EncodeBlueprint_Subtype[o])

# Proportion of cells
ggplot(data=dataForBarPlot, aes(x=EncodeBlueprint_Subtype, y=proportion, 
                                fill=EncodeBlueprint_Subtype, color=EncodeBlueprint_Subtype)) +
  geom_bar(stat='identity', width=.6) + #coord_flip()
  ylim(0,40) +
  #geom_errorbar(aes(ymin=expr, ymax=expr+sd), width=.2, size=0.5, #expr-sd
  #              position=position_dodge(.9)) +
  labs(x='', y=expression('Proportion of Cells (%)')) +
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
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 1, unit = "cm"))


# Number of cells
ggplot(data=dataForBarPlot, aes(x=EncodeBlueprint_Subtype, y=freq, 
                                fill=EncodeBlueprint_Subtype, color=EncodeBlueprint_Subtype)) +
  geom_bar(stat='identity', width=.6) + #coord_flip()
  ylim(0,6000) +
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
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 1, unit = "cm"))


###
cells <- unique(annotation$EncodeBlueprint_Subtype)
samples <- unique(annotation$sample)
samples

dataForBarPlot <- data.frame(sample=rep(samples,each=length(cells)),
                             EncodeBlueprint_Subtype=rep(cells, length(samples)),
                             freq=0,
                             proportion=0,
                             stringsAsFactors = F)
dataForBarPlot

rownames(dataForBarPlot) <- paste0(dataForBarPlot$sample, '_', dataForBarPlot$EncodeBlueprint_Subtype)
dataForBarPlot


total <- annotation %>% group_by(sample) %>%
  summarise(freq=length(EncodeBlueprint_Subtype), sam=length(unique(EncodeBlueprint_Subtype)))
total

freq <- annotation %>% group_by(sample, EncodeBlueprint_Subtype) %>%
  summarise(freq=length(EncodeBlueprint_Subtype))

freq$proportion <- freq$freq/rep(total$freq, total$sam)*100
freq

freq <- data.frame(freq, stringsAsFactors = F)

rownames(freq) <- paste0(freq$sample, '_', freq$EncodeBlueprint_Subtype)
freq


dataForBarPlot[rownames(freq),]$freq <- freq$freq
dataForBarPlot[rownames(freq),]$proportion <- freq$proportion

dataForBarPlot


o <- order(dataForBarPlot$freq[dataForBarPlot$sample=='S1_GEX'], decreasing = T)

dataForBarPlot$EncodeBlueprint_Subtype <- factor(dataForBarPlot$EncodeBlueprint_Subtype, 
                                                  levels=dataForBarPlot$EncodeBlueprint_Subtype[dataForBarPlot$sample=='S1_GEX'][o])


ggplot(data=dataForBarPlot, aes(x=EncodeBlueprint_Subtype, y=freq, 
                                fill=sample)) +
  geom_bar(stat='identity', width=0.8, position='dodge') + #coord_flip()
  ylim(0,2000) +
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
        legend.text = element_text(size=12),
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



####
ggplot(data=dataForBarPlot, aes(x=EncodeBlueprint_Subtype, y=proportion, 
                                fill=sample)) +
  geom_bar(stat='identity', width=0.8, position='dodge') + #coord_flip()
  ylim(0,50) +
  #geom_errorbar(aes(ymin=expr, ymax=expr+sd), width=.2, size=0.5, #expr-sd
  #              position=position_dodge(.9)) +
  labs(x='', y=expression('Proportion of Cells (%)')) +
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



################


dataForTSNEPlot <- annotation
colnames(dataForTSNEPlot)

#========================================
# A subset of cells
cells <- c('CD4+ Naive T cells','CD4+ Tcm','CD4+ Tem',
           'CD8+ Naive T cells','CD8+ Tcm','CD8+ Tem',
           'Tregs','Naive B-cells','Memory B-cells',
           'Class-switched memory B-cells','NK cells')

dataForTSNEPlot <- annotation[annotation$EncodeBlueprint_Subtype %in% cells,]
colnames(dataForTSNEPlot)
#=========================================

dataForTSNEPlot$EncodeBlueprint_Subtype <- factor(dataForTSNEPlot$EncodeBlueprint_Subtype, 
                                                   levels=dataForBarPlot$EncodeBlueprint_Subtype[o])

colorPanel22 <- c('#F8766D','#CD9600','#FFDF00','#7CAE00','#00A9FF','#00BFC4','#C77CFF','#FF61CC',
                  'brown1','orange','limegreen','royalblue','aquamarine','purple','magenta', # gold
                  '#AE3121','#875700','#3C6D00','#006CC2','#007F85','#8D26C8','#BB008D') # '#FFB90F'

cellNum <- length(unique(dataForTSNEPlot$EncodeBlueprint_Subtype))
cellNum

cellColor <- colorPanel22[1:cellNum]

### 
p2 <- ggplot(dataForTSNEPlot) + geom_point(aes(x=tSNE_1, y=tSNE_2, 
                                               color=EncodeBlueprint_Subtype), size=0.3) + 
  #scale_fill_manual(values = colors) +
  scale_color_manual(values=cellColor) +
  #scale_color_hue() +
  labs(x='tSNE_1', y='tSNE_2') +
  guides(colour = guide_legend(override.aes = list(size=2),
                               ncol=1)) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour='white'),
                   panel.background = element_blank()) +
  theme(axis.text=element_text(size=14, color='black'),
        axis.title=element_text(size=16)) +
  theme(legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = 'right')

p2




###

markers.to.plot <- c('CD3E','CD4','CD8B','CCR7','SELL','IL7R',
                     'FOXP3','GNLY','NKG7','TRAC', 'TRDC', 'TRGC1')

### Manual FeaturePlot

exprDa <- as.matrix(seurat.integrated[['RNA']]@data[markers.to.plot,annotation$barcode])
exprDa[1:3,1:10]
dim(exprDa)

rownames(seurat.integrated[['tsne']]@cell.embeddings)==annotation$barcode

df = data.frame(seurat.integrated[['tsne']]@cell.embeddings,
                celltype=annotation$EncodeBlueprint_Subtype,
                t(exprDa))

df[1:5,]


df = melt(df,id.vars = c('tSNE_1','tSNE_2','celltype'))
df[1:5,]

ggplot(df,aes(x=tSNE_1,y=tSNE_2,color=value)) +
  geom_point(size=0.1)+
  scale_color_gradient(low="lightgrey", high="blue") +
  facet_wrap(~variable,ncol=4) +
  labs(color = 'Normalized\nExpression', x='tSNE_1', y='tSNE_2') +
  theme_classic()+
  theme(strip.background = element_blank(), 
        strip.text = element_text(size=14, face='bold'),
        axis.text=element_text(size=14, color='black'),
        axis.title=element_text(size=16),
        legend.title = element_text(size=8))


### Expression of Individual Genes




dfDa <- df %>% group_by(variable, celltype) %>%
  summarise(num=sum(value!=0), total=length(value), proportion=num/total*100)



dataForBarPlot <- dfDa[dfDa$variable=='GNLY',]


o <- order(dataForBarPlot$proportion, decreasing = T)

dataForBarPlot$celltype <- factor(dataForBarPlot$celltype, 
                                  levels = dataForBarPlot$celltype[o])

dataForBarPlot$proportion


###
ggplot(data=dataForBarPlot, aes(x=celltype, y=proportion, fill='#4285F4')) +
  geom_bar(stat='identity', width=.8) + #coord_flip()
  #scale_x_discrete(limits=rev(dataForBarPlot$celltype)) +
  #geom_errorbar(aes(ymin=expr, ymax=expr+sd), width=.2, size=0.5, #expr-sd
  #              position=position_dodge(.9)) +
  labs(x='', y=expression('Proportion of Cells (%)')) +
  #ylim(0,60) +
  #facet_wrap(~variable, nrow=1, scales='free_x') +
  #scale_y_continuous(trans = 'sqrt',
  #                   breaks = c(0,2.5,50,250,750),
  #                   labels = c(0,2.5,50,250,750)) +
  #scale_y_sqrt() +
  #scale_y_continuous(trans='log2') +
  #scale_y_continuous(breaks=seq(0,60,10),
  #                   labels=seq(0,60,10),
  #                   limits = c(0,60)) +
  scale_fill_manual(values = rep('#4285F4',nrow(dataForBarPlot))) +
  #scale_color_manual(values = rep('black',nrow(dataForBarPlot))) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour='black'),
                   panel.background = element_blank()) +
  theme(axis.text=element_text(size=14, color='black'),
        axis.text.x =element_text(size=14, angle = 45, color='black', hjust = 1),
        axis.title.x =element_blank(),
        axis.title.y =element_text(size=16),
        strip.text = element_text(face = 'bold', size=12)) +
  theme(legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = 'none') + 
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 2.75, unit = "cm"))



dataForBoxPlot <- df[df$variable=='GNLY',]

colnames(dataForBoxPlot) <- c('x','y','group','gene','expr')


exprMean <- dataForBoxPlot %>% group_by(group) %>% summarise(meanExpr=mean(expr))

pro <- dataForBoxPlot %>% group_by(group) %>%
  summarise(num=sum(expr!=0), total=length(expr), proportion=round(num/total*100,3))

exprMean$group==pro$group
pro

dataForBoxPlot$group <- paste0(dataForBoxPlot$group, ' (', pro$proportion[match(dataForBoxPlot$group, pro$group)], '%)')

exprMean$group <- paste0(exprMean$group,' (', pro$proportion, '%)')





#dataForBoxPlot$group <- factor(dataForBoxPlot$group,
#                               levels=exprMean$group[order(dataForBoxPlot[match(exprMean$group,dataForBoxPlot$group),]$dataset,exprMean$meanExpr)])

dataForBoxPlot$group <- factor(dataForBoxPlot$group,
                               levels=exprMean$group[order(exprMean$meanExpr, decreasing = T)])


#dataForBoxPlot$group <- factor(dataForBoxPlot$group,
#                               levels=pro$group[order(pro$proportion, decreasing = T)])


ggplot(data=dataForBoxPlot, aes(x=group, y=expr)) +
  geom_violin(aes(fill=group),#size=0.5,
              outlier.shape = NA, outlier.size = NA,#outlier.colour = 'black',
              outlier.fill = NA, width=0.8) +
  #coord_flip() +
  geom_jitter(size=0.01, width = 0.1) +
  #facet_wrap(~study, nrow=3, scales='free') +
  labs(x='', y=expression('Log'[2]*'(Normalized UMI Count)')) +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = 'none') +
  #theme_set(theme_minimal()) #
  theme(axis.title=element_text(size=18),
        axis.text = element_text(color='black', size=14),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.text = element_text(size=14)) +
  theme(axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 4.25, unit = "cm"))



###

expr.data <- seurat.integrated[['RNA']]@data
meta.data <- seurat.integrated@meta.data
meta.data


colnames(expr.data) == rownames(meta.data)

nkt <- c('CD4+ Naive T cells','CD4+ Tcm','CD4+ Tem',
         'CD8+ Naive T cells','CD8+ Tcm','CD8+ Tem',
         'Tregs','NK cells')

d <- c()
for (sam in samples) {
  
  idx <- which(meta.data$sample==sam)
  
  ifng <- expr.data['IFNG',idx]
  ifng <- ifelse(ifng==0,0,1)
  
  cell.type <- meta.data$EncodeBlueprint_Subtype[idx]
  #clonotype <- tcr.contig[[sam]][ovlp, 'label']
  
  for (cell in nkt) {
    print (cell)
    cell.num <- sum(cell.type==cell)
    ifng.num <- sum(cell.type==cell & ifng==1)
    ratio <- ifng.num/cell.num
    
    d <- rbind(d, c(ifng.num, cell.num, ratio, cell, sam))
    
  }
  
}

d <- data.frame(d)

colnames(d) <- c('ifng', 'total', 'ratio', 'cell', 'sample')

d$ratio <- round(as.numeric(as.character(d$ratio))*100,2)
d

d$label <- paste(paste0(d$ratio, '%'), ' (', paste(d$ifng, d$total, sep=' / '), ')', sep='')
d$label


d$ratio

ggplot(data=d, aes(x=cell, y=ratio, fill=sample)) + 
  geom_bar(stat='identity', width=.8, position=position_dodge2(reverse=T)) + coord_flip() +
  scale_x_discrete(limits=rev(unique(d$cell))) +
  #geom_errorbar(aes(ymin=expr, ymax=expr+sd), width=.2, size=0.5, #expr-sd
  #              position=position_dodge(.9)) +
  labs(x='', y='Proportion of Cells Expressing IFNG (%)') + 
  geom_text(aes(label = label), hjust=-0.05, size=3, 
            position=position_dodge2(width=0.8, reverse=T)) +
  #ylim(0,10) +
  theme_bw()+
  #scale_y_continuous(trans = 'sqrt',
  #                   breaks = c(0,2.5,50,250,750),
  #                   labels = c(0,2.5,50,250,750)) +
  #scale_y_sqrt() +
  #scale_y_continuous(trans='log2') +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.position = c(0.85,0.85)) +
  theme(axis.title=element_text(size=16), 
        axis.text = element_text(color='black', size=14),
        axis.text.x = element_text(angle = 0, hjust=0.5)) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank()) +
  theme(plot.margin =  margin(t = 0.25, r = 2, b = 0.25, l = 0.25, unit = "cm"))



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
  
  L <- list()
  
  L$count <- raw.data
  L$tpm <- expr.data
  L$condt <- group
  
  grp <- L$condt
  cdr <- scale(colMeans(L$tpm > 0))
  sca <- FromMatrix(exprsArray = log2(L$tpm + 1), 
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
  degTable <- degTable[order(degTable$FDR, decreasing = F),]
  
  fl = paste(stimSample, ctrlSample, cell, 'MASTtpmDetRate_filter_TP10K0.5_Cell10_noMito_noRibo', sep='_')
  
  write.table(degTable, file=paste0('report/', fl, '.txt'),
              sep='\t', quote=F)
  
}


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



###
FeaturePlot(seurat.integrated, 
            features = c("CD3D", "GNLY", "IFI6"), 
            split.by = "sample", 
            max.cutoff = 3, 
            cols = c("grey", "red"))


plots <- VlnPlot(seurat.integrated, 
                 features = c("IFI6", "ISG15"), 
                 split.by = "sample", 
                 group.by = "EncodeBlueprint_Subtype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)




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


####################################################################
####### STEP 8. Metascape GO Enrichment Analysis

regulation <- 'Down'

fls <- file.path('report/Enrichment', list.files('report/Enrichment/'))
#idx <- which(grepl('S2_GEX_S1_GEX', fls) & grepl('Up', fls))
idx <- grepl(regulation, fls)
idx

fls <- fls[idx]

commonTerms <- c()

for (fl in fls) {
  
  print (fl)
  
  kegg <- read_excel(fl, sheet='Enrichment')
  idx <- grep('Summary', kegg$GroupID)+1
  kegg <- kegg[idx,]
  
  commonTerms <- c(commonTerms, kegg$Term)
  
}

unique(commonTerms)


###

stimSample <- 'S2_GEX'
ctrlSample <- 'S1_GEX'


dataForBubblePlot <- c()

for (cell in cells) {
  
  print (cell)
  
  fl = paste('Metascape', stimSample, ctrlSample, cell, regulation, sep='_')
  fl <- file.path('report/Enrichment', paste0(fl, '.xlsx'))
  
  if (!file.exists(fl)) {
    next
  }
  
  kegg <- read_excel(fl, sheet='Enrichment')
  kegg$Comparison <- cell
  
  idx <- which(grepl('Member', kegg$GroupID) & kegg$Term %in% commonTerms)
  
  kegg <- kegg[idx,]
  
  dataForBubblePlot <- rbind(dataForBubblePlot, kegg)
  
}

dataForBubblePlot

unique(dataForBubblePlot$Term)





colnames(dataForBubblePlot)[6] <- 'FDR'

dataForBubblePlot$FDR <- 10^dataForBubblePlot$FDR

dataForBubblePlot <- dataForBubblePlot[dataForBubblePlot$FDR<=0.05,]
dataForBubblePlot$Count <- unlist(lapply(strsplit(dataForBubblePlot$Symbols, ','),length))
dataForBubblePlot$Term <- paste0(dataForBubblePlot$Term, '~', dataForBubblePlot$Description)
dataForBubblePlot$Comparison <- factor(dataForBubblePlot$Comparison, levels=cells)
dataForBubblePlot

dataForBubblePlot <- dataForBubblePlot %>% group_by(Comparison) %>% dplyr::arrange(FDR, .by_group = TRUE)

write.xlsx(data.frame(dataForBubblePlot), file='dataForBubblePlot.xlsx')

pval <- 0.05

ggplot(dataForBubblePlot, mapping=aes(x=Term, y=Comparison, #y=-log10(Benjamini), #y=Fold.Enrichment
                                      color=FDR,size=Count)) +
  geom_point()+ coord_flip() +
  scale_x_discrete(limits=rev(unique(dataForBubblePlot$Term))) +
  #scale_x_discrete(limits=Order)+
  scale_colour_gradientn(limits=c(0,pval),
                         colors= c("red","yellow","green")) + #
  #facet_wrap(~Comparison) +
  #facet_grid(Regulation~Comparison) + # scales=free
  xlab('')+ylab('') +
  guides(shape = guide_legend(order=1),
         colour = guide_colourbar(order=2, title = 'FDR')) + #'P Value\n(Benjamini)'))
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour='black'),
                   panel.background = element_blank()) +
  ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size=20)) +
  theme(axis.text=element_text(size=14, color='black'),
        axis.text.x =element_text(size=14, color='black', angle=45, hjust=1),
        axis.title=element_text(size=15)) +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14),
        legend.key.size = unit(0.8,'cm'))



####################################################################
####### STEP 9. Gene Set Enrichment Analysis (GSEA)




