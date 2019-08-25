

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
  
  feature.names = read.delim(feature.path, 
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


# Google colors
google.red <- '#EA4335'
google.yellow <- '#FBBC05'
google.green <- '#34A853'
google.blue <- '#4285F4'



#########
## The first letter of the first world
capitalizeRL <- function(x) {
  substr(x,1,1) <- toupper(substr(x,1,1))
  return(x)
}

RenameCellRL <- function(x, original.cell, target.cell) {
  n <- length(original.cell)
  
  for (i in 1:n) {
    x[x==original.cell[i]] <- target.cell[i]
  }
  
  return (x)
}


#############
fantom5 <- readRDS(file='SingleR_FANTOM5_Reference.rds') 



