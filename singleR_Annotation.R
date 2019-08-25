####################################################################
####### STEP 6. SingleR Annotation

devtools::install_github('dviraran/SingleR')
library(singleR)

seurat.integrated <- readRDS('data/rData/seurat_integrated.rds')

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



########################################
### Create Reference

name = 'FANTOM5'
expr = as.matrix(expr) # the expression matrix
types = as.character(types) # a character list of the types. Samples from the same type should have the same name.
main_types = as.character(main_types) # a character list of the main types. 
fantom5 = list(name=name,
               data = expr, 
               types=types, 
               main_types=main_types)

# if using the de method, we can predefine the variable genes
fantom5$de.genes = CreateVariableGeneSet(expr,types,200)
fantom5$de.genes.main = CreateVariableGeneSet(expr,main_types,300)

# if using the sd method, we need to define an sd threshold
sd = rowsSd(expr)
sd.thres = sort(sd, decreasing = T)[4000] # or any other threshold
fantom5$sd.thres = sd.thres

saveRDS(fantom5,file='SingleR_FANTOM5_Reference.rds') 
# it is best to name the object and the file with the same name.


###########
fantom5 <- readRDS(file='SingleR_FANTOM5_Reference.rds') 

### Single cell
singlerCellAnnotation = CreateSinglerObject(seurat.integrated[['RNA']]@counts, annot = NULL, project.name='PBMC', 
                                            min.genes = 0, technology = "10X", species = "Human", citation = "",
                                            ref.list = list(fantom5), normalize.gene.length = F, variable.genes = "de",
                                            fine.tune = T, do.signatures = T, clusters = NULL, do.main.types = T, 
                                            reduce.file.size = T, numCores = 4)
