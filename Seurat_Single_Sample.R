
################################################

######## STEP 1. Read Raw Data

setwd('~/Projects/Infrastructure_20190109/BioinfoHub/pages/BioinfoHub/')

Read10XRL <- function(data.dir=NULL) {
  
  #barcode.path <- file.path(data.dir, "barcodes.tsv.gz")
  #features.path <- file.path(data.dir, "features.tsv.gz")
  #matrix.path <- file.path(data.dir, "matrix.mtx.gz")
  
  barcode.path <- file.path(data.dir, "barcodes.tsv")
  features.path <- file.path(data.dir, "genes.tsv")
  matrix.path <- file.path(data.dir, "matrix.mtx")
  
  mat <- readMM(file = matrix.path)
  
  feature.names = read.delim(features.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V2
  
  dupIdx <- which(duplicated(rownames(mat)))
  
  rownames(mat)[dupIdx] <- paste0(rownames(mat)[dupIdx], '[dup]')
  
  return (mat)
  
}


matrixDir <- paste0('data/from10X/pbmc4k_filtered_gene_bc_matrices/GRCh38/')

mat <- Read10XRL(data.dir = matrixDir)
mat[1:5,1:5]
dim(mat)



####### STEP 2. Create Seurat Object and Cell QC

###### Control sample

# Initialize the Seurat object with the raw (non-normalized) data
# Keep all genes expressed in >= 5 cells (~0.1% of the data)
ctrlSample <- '4kPBMC'

ctrl <- CreateSeuratObject(counts = mat, project = ctrlSample, min.cells = 5)
ctrl@meta.data$sample <- ctrlSample
ctrl@meta.data

ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^MT-")
names(ctrl@meta.data)


VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



# Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


# Number of genes expressed by a cell [200, 2500]
# Percentage of total mitochondrial RNA in a cell < 5%

ctrl <- subset(ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

message(paste0(ncol(ctrl[['RNA']]@counts), '/', ncol(mat), ' of cells are kept'))



# Normalization
ctrl <- NormalizeData(object = ctrl, normalization.method = "LogNormalize", 
                      scale.factor = 10000)


# Find variable genes
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ctrl), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ctrl)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))



# === TODO ===
# What is vst ?
# Understand Variable Genes


# Scaling the data
all.genes <- rownames(ctrl)
ctrl <- ScaleData(ctrl, features = all.genes)
ctrl


####### STEP 3. PCA & UMAP & tSNE

### Perform linear dimensional reduction

ctrl <- RunPCA(ctrl, features = VariableFeatures(object = ctrl))

# Examine and visualize PCA results a few different ways
print(ctrl[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(ctrl, dims = 1:2, reduction = "pca")

DimPlot(ctrl, reduction = "pca")

DimHeatmap(ctrl, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(ctrl, dims = 1:15, cells = 500, balanced = TRUE)



######################################################################
### Determine the ‘dimensionality’ of the dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc)

ctrl <- FindNeighbors(ctrl, dims = 1:10)
ctrl <- FindClusters(ctrl, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(ctrl), 5)



##################################################################
### Run non-linear dimensional reduction (UMAP/tSNE)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
ctrl <- RunUMAP(ctrl, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")


saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")



#######################################################################
### t-SNE and Clustering
ctrl <- RunTSNE(ctrl, reduction = "pca", dims = 1:20)

DimPlot(ctrl, reduction = "tsne")

#saveRDS(ctrl, file='data/rData/PBMC68K.rds')

ctrl@reductions$tsne@cell.embeddings


ctrl <- FindNeighbors(ctrl, reduction = "pca", dims = 1:20)
ctrl <- FindClusters(ctrl, #reduction.type = "cca.aligned", 
                     resolution = 0.5, dims.use = 1:20, force.recalc = TRUE)


# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(ctrl, reduction = "tsne")




