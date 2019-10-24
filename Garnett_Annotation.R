
############################################################################################

BiocManager::install('batchelor')

BiocManager::install('IRanges')

devtools::install_github('cole-trapnell-lab/monocle3')

# Install a few Garnett dependencies:
BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db'))

# Install Garnett
devtools::install_github("cole-trapnell-lab/garnett", ref="monocle3")


library(garnett)

pbmc_classifier <- readRDS('script/hsPBMC_20191017.RDS')
pbmc_classifier

pbmc_cds <- new_cell_data_set(as(seurat.integrated[['RNA']]@counts, "dgCMatrix"),
                              cell_metadata = seurat.integrated@meta.data,
                              gene_metadata = data.frame(gene_short_name=rownames(seurat.integrated[['RNA']]@counts),
                                                         Symbol=rownames(seurat.integrated[['RNA']]@counts),
                                                         row.names=rownames(seurat.integrated[['RNA']]@counts),
                                                         stringsAsFactors = F))


library(org.Hs.eg.db)
pbmc_cds <- classify_cells(pbmc_cds, pbmc_classifier,
                           db = org.Hs.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "SYMBOL")


View(data.frame(pData(pbmc_cds)))

table(pData(pbmc_cds)$cell_type)
table(pData(pbmc_cds)$cluster_ext_type)
