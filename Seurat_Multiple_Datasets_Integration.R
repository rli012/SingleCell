setwd('C:\\Users/rli3/Desktop/share/BioinfoHub/')

#library(devtools)
#devtools::install_github(repo = "satijalab/seurat", ref = "develop",
#                         INSTALL_opts = c('--no-lock'))

library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
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

#saveRDS(seurat.integrated, file='data/rData/Seurat_Integrated.rds')
seurat.integrated@meta.data[1:5,]

p1 <- DimPlot(seurat.integrated, reduction = "tsne", group.by = "sample")
p2 <- DimPlot(seurat.integrated, reduction = "tsne", group.by = "seurat_clusters")
#p2 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
#              repel = TRUE) + NoLegend()
plot_grid(p1, p2)


