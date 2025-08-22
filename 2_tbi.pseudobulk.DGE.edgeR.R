library( edgeR )

meta = read.csv('tbi/tbi_snRNA_meta.csv')
rownames(meta) = meta$sample_name
Group <- factor(paste(meta$group,meta$age,sep=".")) 
meta <- cbind(meta,Group=Group)

expr.list = readRDS('/tbi/MoTBI.AllSamples.Integrated.Cleaned.Annotated.CellRanger.pseudobulk_counts.rds')

dge.list = list()

types = names(expr.list)
n = length(types)

# which cell types have enough data

lib.sizes = list()
for( i in 1:n ) {
  lib.sizes[[i]] = colSums( expr.list[[i]] )
}
lib.sizes = do.call( cbind , lib.sizes )
colnames(lib.sizes) = types
min.size = apply( lib.sizes , 2 , min )
nsmall = colSums( lib.sizes < 10000 ) 

drop.types = which( nsmall > 0 )

for( i in 1:n ) {
  if( i %in% drop.types ) {
    dge.list[[i]] = NULL
    next
  }
  cat( types[i] , '\n' )
  expr = expr.list[[i]]
  meta = meta[ colnames(expr) , ]
  keepCols = which( colSums(expr) >= 10000 )
  m = meta[ keepCols, ] 
  e = expr[ , keepCols]
  keepRows = filterByExpr( e , group = m$group , min.count = 2 )
  #keepRows = which( rowSums(expr) >= 2 )
  e = e[ keepRows , ]
  logcpm = cpm( e , log = T )
  centered = t(scale(t(logcpm)))
  pca = prcomp( centered )
  pdf( paste( types[i] , 'mdsplot.pdf' , sep = '.' ) )
  plotMDS( logcpm )
  plot( x = pca$rotation[,1] , y = pca$rotation[,2] , 
        xlab = 'PC1' , ylab = 'PC2' , type = 'n' )
  text( x = pca$rotation[,1] , y = pca$rotation[,2] ,
        labels = rownames(pca$rotation) )
  plot( x = pca$rotation[,3] , y = pca$rotation[,4] , 
        xlab = 'PC3' , ylab = 'PC4' , type = 'n' )
  text( x = pca$rotation[,3] , y = pca$rotation[,4] ,
        labels = rownames(pca$rotation) )
  dev.off()
  z = sapply( 1:ncol(pca$rotation) , function(i) 
	      pca$rotation[,i] / sd(pca$rotation[,1]) )
  maxz = apply( abs(z[,1:5]) , 1 , max )
  outlier = maxz > 4 #3
  e = e[,outlier==F]
  m = m[outlier==F,]
  mt = colSums( e[ grep('^MT-',rownames(e)) , ] )
  m$percent.mt = mt / colSums(e)
  m$mt2 = m$percent.mt ^ 2
  y = DGEList( counts = e , group = m$Group )
  y <- calcNormFactors(y)
  design = model.matrix( ~ 0 + Group, data = m )
  colnames(design) = levels(Group)#gsub('group','',colnames(design)[1:3] )
  y = estimateDisp( y , design , robust = T )
  fit = glmFit(y,design)
  sham.3movs18mo = glmLRT( fit , design , 
	contrast = c(1,-1,0,0) )
  TBI.3movs18mo = glmLRT( fit , design , 
	contrast = c(0,0,1,-1) )
  TBIvssham.3mo = glmLRT( fit , design , 
	contrast = c(0,-1,0,1) )
  TBIvssham.18mo = glmLRT( fit , design , 
	contrast = c(-1,0,1,0) )
  TBIvssham.bothages = glmLRT( fit , design , 
	contrast = c(-0.5,-0.5,0.5,0.5) )
  bothgroups.3movs18mo = glmLRT( fit , design , 
	contrast = c(0.5,-0.5,0.5,-0.5) )
  res1 = as.data.frame(topTags( sham.3movs18mo , n = Inf ))
  res2 = as.data.frame(topTags( TBI.3movs18mo , n = Inf ))
  res3 = as.data.frame(topTags( TBIvssham.3mo , n = Inf ))
  res4 = as.data.frame(topTags( TBIvssham.18mo , n = Inf ))
  res5 = as.data.frame(topTags( TBIvssham.bothages , n = Inf ))
  res6 = as.data.frame(topTags( bothgroups.3movs18mo , n = Inf ))
  res1$Contrast = 'sham.3movs18mo'
  res2$Contrast = 'TBI.3movs18mo'
  res3$Contrast = 'TBIvssham.3mo'
  res4$Contrast = 'TBIvssham.18mo'
  res5$Contrast = 'TBIvssham.bothages'
  res6$Contrast = 'bothgroups.3movs18mo'
  res6$CellType = types[i]
  res6$Gene = rownames(res6)
  res5$CellType = types[i]
  res5$Gene = rownames(res5)
  res4$CellType = types[i]
  res4$Gene = rownames(res4)
  res1$CellType = types[i]
  res1$Gene = rownames(res1)
  res2$CellType = types[i]
  res2$Gene = rownames(res2)
  res3$CellType = types[i]
  res3$Gene = rownames(res3)
  res = rbind( res1 , res2 , res3 , res4, res5, res6 )
  res = res[ order( res$PValue ) , ]
  dge.list[[i]] = res

}

dge = do.call( rbind , dge.list )
dge = dge[ order( dge$PValue ) , ]
write.table( dge , 'tbi/tbi.pseudobulk.edgeR.rmOutliers.dge.2024.01.30.txt' )


