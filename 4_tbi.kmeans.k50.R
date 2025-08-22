# calculate kmeans clusters

library( Seurat )
# the network reconstruction is then performed using the smoothed data (after knn-smoothing)
objs = readRDS('/tbi/tbi.14celltypes.seurat.alra.rds')
types = levels(objs)

for (i in 1:length(types)){
  obj = subset(x = objs, idents = types[i])
  obj = NormalizeData(obj)
  datExpr0 = as.matrix( GetAssayData(obj) )
  # Z score
  datExpr=sweep(datExpr0,1,apply(datExpr0,1,mean),"-")
  indx.sd=(apply(datExpr0,1,sd))==0 # these will produce NAs
  datExpr=sweep(datExpr,1,apply(datExpr,1,sd),"/")
  datExpr[indx.sd,]=0
  if(sum(is.na(datExpr))!=0){print("NAs in exprsMTX.Z Zscores!")}
  #kmeans
  k = 50
  cl=kmeans(x=datExpr,centers=k,iter.max=10000,nstart=10,alg="Lloyd")
  saveRDS( cl ,  paste( types[i] , 'kmeans.k50.rds' , sep='.'))
}
