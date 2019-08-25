

# Read Expression Matrix

Read10XRL <- function(data.dir=NULL) {
  
  files <- list.files(data.dir)
  
  barcode <- files[grep('barcodes', files)]
  mtx <- files[grep('matrix', files)]
  feature <- files[grep('gene|feature', files)]
  #feature <- setdiff(files, c(barcode, mtx))
  
  barcode.path <- file.path(data.dir, barcode)
  feature.path <- file.path(data.dir, feature)
  matrix.path <- file.path(data.dir, mtx)
  
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