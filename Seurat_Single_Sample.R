library(Seurat)

#################################################################
######## STEP 1. Read Raw Data

setwd('C:\\Users/rli3/Desktop/share/BioinfoHub/')

source('Single_Cell_Functions.R')

matrixDir <- paste0('data/pbmc3k_filtered_gene_bc_matrices/hg19/')
mat <- Read10XRL(data.dir = matrixDir)
mat[1:5,1:5]
dim(mat)

### ONE LINE COMMAND
ctrl <- CreateSeuratObject(mat) %>% 
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent_mt") %>% 
  SCTransform(vars.to.regress = "percent.mt") %>% 
  RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  RunTSNE(dims = 1:30) %>% FindClusters()


##################################################################
####### STEP 2. Create Seurat Object and Cell QC

# Initialize the Seurat object with the raw (non-normalized) data
# Keep genes expressed in >= 3 cells (~0.1% of the data)

project <- 'pbmc3k'
ctrlSample <- 'CTRL'

ctrl <- CreateSeuratObject(counts = mat, project = project, min.cells = 3)
ctrl@meta.data$sample <- ctrlSample
head(ctrl@meta.data)
ctrl[['RNA']]@counts[1:5,1:5]

ctrl[["percent_mt"]] <- PercentageFeatureSet(ctrl, pattern = "^MT-")
ctrl$percent_mt <- PercentageFeatureSet(ctrl, pattern = "^MT-")
names(ctrl@meta.data)


VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3)



# Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "percent_mt")
plot2 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# The number of unique genes detected in each cell [200, 2500]
  # Low-quality cells or empty droplets will often have very few genes
  # Cell doublets or multiplets may exhibit an aberrantly high gene count
# Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
# # Percentage of total mitochondrial RNA in a cell < 5%
  # Low-quality / dying cells often exhibit extensive mitochondrial contamination


ctrl <- subset(ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mt < 5)
VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3)

message(paste0(ncol(ctrl[['RNA']]@counts), '/', ncol(mat), ' of cells are kept'))


##########################################################################
####### STEP 3. Normalization, Feature Selection

# === TODO ===
# vst ?
# CombinePlots doesn't work

###
### Log Normalization (For DE Analysis)
ctrl <- NormalizeData(object = ctrl, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
ctrl[['RNA']]@data[1:5,1:5]


# Find variable genes
# Features that exhibit high cell-to-cell variation in the dataset
# Focusing on these genes in downstream analysis helps to 
# highlight biological signal in single-cell datasets.
# The variable genes will be used in downstream analysis, like PCA

ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ctrl), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ctrl)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, 
                     xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))


###
# Scaling the data
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# Scaling is an essential step in the Seurat workflow, but only on genes that will be used as input to PCA
# Seurat heatmaps (produced as shown below with DoHeatmap) require genes in the heatmap to be scaled

all.genes <- rownames(ctrl)
ctrl <- ScaleData(ctrl, features = all.genes)
ctrl[['RNA']]@scale.data[1:5,1:5]


# remove unwanted sources of variation
# ‘regress out’ heterogeneity associated with (for example) 
# cell cycle stage, or mitochondrial contamination

ctrl <- ScaleData(ctrl, vars.to.regress = "percent_mt")
###


###
### Run sctransform (For PCA, Clustering, etc.)
# Replaces NormalizeData, ScaleData, and FindVariableFeatures
ctrl <- SCTransform(ctrl, vars.to.regress = "percent_mt", verbose = FALSE)
ctrl[['SCT']]@scale.data[1:5,1:5]



####################################################################
####### STEP 4. PCA & UMAP & tSNE

DefaultAssay(ctrl) <- "SCT"

###
### Perform linear dimensional reduction
ctrl <- RunPCA(ctrl, features = VariableFeatures(object = ctrl))

# Examine and visualize PCA results a few different ways
print(ctrl[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(ctrl, dims = 1:2, reduction = "pca")

DimPlot(ctrl, reduction = "pca")

DimHeatmap(ctrl, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(ctrl, dims = 1:15, cells = 500, balanced = TRUE)



### Determine the ‘dimensionality’ of the dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
ctrl <- JackStraw(ctrl, num.replicate = 100)
ctrl <- ScoreJackStraw(ctrl, dims = 1:20)

# sharp drop-off in significance after the first 10-12 PCs.
JackStrawPlot(ctrl, dims = 1:15)

# we can observe an ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.
ElbowPlot(ctrl)

# Cluster the cells
ctrl <- FindNeighbors(ctrl, dims = 1:10)
#ctrl <- FindClusters(ctrl, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(ctrl), 5)



###
### Run non-linear dimensional reduction (UMAP/tSNE)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
ctrl <- RunUMAP(ctrl, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(ctrl, reduction = "umap")


###
### t-SNE and Clustering
ctrl <- RunTSNE(ctrl, reduction = "pca", dims = 1:10)

DimPlot(ctrl, reduction = "tsne")

#saveRDS(ctrl, file='data/rData/ctrl3k.rds')

ctrl@reductions$tsne@cell.embeddings
ctrl@meta.data


ctrl <- FindNeighbors(ctrl, reduction = "pca", dims = 1:10)
ctrl <- FindClusters(ctrl, #reduction.type = "cca.aligned", 
                     resolution = 0.5, dims.use = 1:10, force.recalc = TRUE)


# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(ctrl, reduction = "tsne", label = TRUE)



####################################################################
####### STEP 5. Find Cluster Biomarkers

# find all markers of cluster 1
cluster1.markers <- FindMarkers(ctrl, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(ctrl, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
ctrl.markers <- FindAllMarkers(ctrl, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ctrl.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(ctrl, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(ctrl, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(ctrl, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(ctrl, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))

top10 <- ctrl.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(ctrl, features = top10$gene) + NoLegend()

###
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(ctrl)
ctrl <- RenameIdents(ctrl, new.cluster.ids)
DimPlot(ctrl, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(ctrl, file = "data/rData/pbmc3k_final.rds")



####################################################################
####### STEP 6. SingleR Annotation






####################################################################
####### STEP 7. Summary of Visualization


