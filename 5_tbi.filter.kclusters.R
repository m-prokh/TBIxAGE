# this script is used to take results from k-means clustering and
# perform filtering and merging of the results to produce a final network

library(Seurat)
library(ggplot2)
library(ggpubr)
library(factoextra)
source('alra.R')
DefaultAssay(obj) = 'RNA'
library(WGCNA)

 #allowWGCNAThreads() # req in R studio
options(stringsAsFactors = FALSE)

# minimum size of a gene co-expression cluster
minModSize = 22
# minimum correlation between a gene and the cluster centroid
kAEtoStay = 0.15 #0.15, 0.25! 0.3, 0.35
# threshold for merging clusters with strongly correlated patterns
# the "height" corresponds to 1 - Pearson's r
cutHeight = 0.25 #22, 37 0.25!, 0.3, 0.37, 40, 42

# expression seurat obj
objs = readRDS('tbi.14celltypes.seurat.alra.rds')
types = levels(objs)

#load kmeans
kmns.files = Sys.glob('*kmeans.rds')
kmns = list()
types = gsub( '\\.(.*)','',kmns.files)
n = length(types)
for( i in 1:n ) {
  kmns[[i]] = readRDS( kmns.files[i] )
}
names(kmns) = types

for( i in 1:n ) {
  obj = subset(x = objs, idents = types[i])
  impute1 = NormalizeData(obj)
  datExpr0 = as.matrix( GetAssayData(impute1) )
  # Z score
  datExpr=sweep(datExpr0,1,apply(datExpr0,1,mean),"-")
  indx.sd=(apply(datExpr0,1,sd))==0 # these will produce NAs
  datExpr=sweep(datExpr,1,apply(datExpr,1,sd),"/")
  datExpr[indx.sd,]=0
  if(sum(is.na(datExpr))!=0){print("NAs in exprsMTX.Z Zscores!")}
  kAE = cor( t(km$centers) , t(datExpr) )
  colors = km$cluster
  for( i in 1:length(colors) ) {
    r.i = kAE[ colors[i] , i ]
    if( r.i < kAEtoStay ) colors[i] = 0
  }
  size = table( colors )
  too.small = as.numeric(names(which(size<minModSize)))
  colors[ colors %in% too.small ] = 0

  centers = sapply( sort(unique(colors)) , function(i)
    colMeans(datExpr[ colors == i , ]) )

  colnames(centers) = paste( 'AE' , sort(unique(colors)) , sep = '' )
  r = cor( centers )
  d = as.dist( 1 - r )
  hc = hclust( d , method = 'average' ) #'average', 'ward.D2', 'mcquitty'
  # cophenetic distance between two observations that have been clustered:
  #   intergroup dissimilarity at which the two observations are first combined
  #   into a single cluster
  sort(unique(cophenetic(hc)))

  cl = cutree( hclust( d , method = 'average' ) , h = cutHeight )

  mergeColors = rep(NA,length(colors))
  for( i in 1:max(cl) ) {
    idx = as.numeric( gsub( 'AE','', names(cl)[ cl == i ] ))
    mergeColors[ colors %in% idx ] = i
  }
  mergeColors = mergeColors - 1
  names(mergeColors) = names(colors)
  print(table(mergeColors))
  norm = t(datExpr0)
  norm = norm[ , rownames(impute1) ]
  res = alra( norm )
  alra.out = res[[3]]
  rownames(alra.out) = colnames(obj)

  datExpr2 = alra.out
  clgenes = names(km$cluster)
  datExpr2 = datExpr2[,clgenes]
  reordr = match(colnames(datExpr2), clgenes)

  MEs = moduleEigengenes( datExpr2 , mergeColors )

  saveRDS( mergeColors , paste( types[i] , 'k50.mm22r15ch25.merged.clusters.rds' , sep='.'))
  saveRDS( MEs , paste( types[i] , 'k50.merged.MEs.rds' , sep='.'))
  saveRDS( datExpr2 , paste( 'tbi', types[i] , 'ALRA.rds', sep='.'))
  meta = obj@meta.data
  meta = cbind( meta , MEs$eigengenes , MEs$averageExpr )
  obj@meta.data = meta
  saveRDS( obj , file = paste('tbi', types[i] , 'wMEs.rds', sep='.'))
}
