# Google colors
google.colors <- c('#EA4335', '#FBBC05', '#34A853', '#4285F4')

google.red <- '#EA4335'
google.yellow <- '#FBBC05'
google.green <- '#34A853'
google.blue <- '#4285F4'


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
  
   if (length(dupIdx)!=0) {
    rownames(mat)[dupIdx] <- paste0(rownames(mat)[dupIdx], '[dup]')
  }
  
  return (mat)
  
}


### DEG functions

wilcoxFun <- function(expr.data, group, gene, cell1, cell2) {
  
  meanX <- log2(mean(expm1(expr.data[gene,cell1])))
  meanY <- log2(mean(expm1(expr.data[gene,cell2])))
  
  test <- wilcox.test(expr.data[gene,] ~ group)
  
  pVal <- test$p.value
  #est <- as.numeric(format(tTest$estimate, digits = 3))
  t <- as.numeric(format(test$statistic, digits = 3))
  
  #degTable <- rbind(degTable, c(length(cell1), length(cell2), meanX, meanY, t, pVal))
  return (c(length(cell1), length(cell2), meanX, meanY, t, pVal))
  
}



exprSummaryFun <- function(expr.data, gene, cell1, cell2) {
  
  meanX <- log2(mean(expm1(expr.data[gene,cell1])))
  meanY <- log2(mean(expm1(expr.data[gene,cell2])))
  
  return (c(meanX, meanY))
  
}



ttestFun <- function(expr.data, group, gene) {
  test <- t.test(expr.data[gene,] ~ group)
  ttest.pVal <- test$p.value
  t <- as.numeric(format(test$statistic, digits = 3))
  return (c(t, ttest.pVal))
}


#########
## The first letter of the first world
capitalizeFun <- function(x, to.lower=TRUE) {
  if (!is.character(x)) {
    x <- as.character(x)
  }
  
  if (to.lower==TRUE) {
    x <- tolower(x)
  }
  
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


#################
fantom5 <- readRDS(file='SingleR_FANTOM5_Reference.rds') 


#################
#Identify differential expressed genes across conditions

LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0, 
                       adj.y.s = 0, text.size = 2.5, segment.size = 0.1) {
  for (i in genes) {
    x1 <- exp.mat[i, 1]
    y1 <- exp.mat[i, 2]
    plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t, 
                            label = i, size = text.size)
    plot <- plot + annotate("segment", x = x1 + adj.x.s, xend = x1, y = y1 + 
                              adj.y.s, yend = y1, size = segment.size)
  }
  return(plot)
}

LabelUR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05, 
                    adj.r.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = adj.r.t, 
                    adj.y.s = adj.u.s, adj.x.s = adj.r.s, ...))
}

LabelUL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05, 
                    adj.l.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = -adj.l.t, 
                    adj.y.s = adj.u.s, adj.x.s = -adj.l.s, ...))
}


