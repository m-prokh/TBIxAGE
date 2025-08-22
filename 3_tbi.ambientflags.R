dge = read.delim('tbi.pseudobulk.edgeR.rmOutliers.dge.2024.01.30.txt')

dge$direction = NA
for (i in 1:length(dge$logFC)){
  if (dge$logFC[i] < 0 ){
    dge$direction[i] <- "down"
  } else if (dge$logFC[i] > 0 ){
    dge$direction[i] <- "up"
  } 
}

colnames(dge)[8] = 'gene'

g = unique(dge$Gene )
n = length(g)
logcpm.max = sapply( 1:n , function(i)
  max( dge$logCPM[ dge$gene == g[i] ] ) )
logcpm.max = data.frame(
  gene = g , logcpm.max )

dge2 = merge( dge , logcpm.max , by = 'gene' )
dge2$logcpm.diff = dge2$logCPM - dge2$logcpm.max
dge2$flag_ambient = dge2$logcpm.diff < -3.32

dge2 = merge( dge , dge2[,c(1,2,3,8,10,11,12)] , by = c('gene','CellType') )
dge2 = dge2[ order( dge2$PValue ) , ]

write.table( dge2 , file = 'dge.with_ambient_flag.txt' )

###########

dge2 = read.table('dge.with_ambient_flag.txt' , header=T)

dge2$logP = -log10( dge2$PValue)

dge2$logP.signed = dge2$logP * sign( dge2$logFC.x)

clean = dge2[ dge2$flag_ambient == F , ]

write.table( clean , file = 'dge.no_ambient.txt' )
