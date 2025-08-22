library( Seurat )

# test for overlap between modules and DEGs in each cell type

# load modules
setwd("/Users/mariya/tbi_local/colorfiles")
color.files = Sys.glob('*clusters.rds')
color.list = list()
types = gsub( '\\.(.*)','',color.files) #the first part of the filename before<.>
n = length(types)
for( i in 1:n ) {
  color.list[[i]] = readRDS( color.files[i] )
}
names(color.list) = types

# load DEGs
dge = read.table('tbi/tbi.pseudobulk.edgeR.rmOutliers.dge.2024.01.30.txt',
                 header=TRUE)

types = names(color.list) #celltypes
n = length(types)
df.all = list()

for( j in 1:6 ) {
 contrast = c('bothgroups.3movs18mo', 'sham.3movs18mo', 'TBI.3movs18mo', 
              'TBIvssham.bothages', 'TBIvssham.3mo', 'TBIvssham.18mo')[j]
 df = list()
 for( i in 1:n ) {
  dge.i = dge[ dge$Contrast == contrast & dge$CellType == types[i] , ]
  colors = color.list[[i]]
  m = max(colors)
  u = intersect( names(colors) , dge.i$Gene )
  up = dge.i$Gene[ dge.i$PValue < 0.05 & dge.i$logFC > 0 ]
  dn = dge.i$Gene[ dge.i$PValue < 0.05 & dge.i$logFC < 0 ]
  p.dn = rep(NA,m)
  p.up = rep(NA,m)
  or.dn = rep(NA,m)
  or.up = rep(NA,m)
  o.dn = rep(NA,m)
  o.up = rep(NA,m)
  nSetGenes = rep(NA,m)
  nUp = rep(length(up),m)
  nDown = rep(length(dn),m)
  for( k in 1:m ) {
      set = intersect( names(colors)[colors==k] , u )
      if( length(set) == 0 ) next
      nSetGenes[k] = length(set)
      t.up = table( u %in% up , u %in% set )
      test.up = fisher.test( t.up , alternative = 'greater' )
      p.up[k] = test.up$p.value
      or.up[k] = test.up$estimate
      o.up[k] = t.up[2,2]
      t.dn = table( u %in% dn , u %in% set )
      test.dn = fisher.test( t.dn , alternative = 'greater' )
      p.dn[k] = test.dn$p.value
      or.dn[k] = test.dn$estimate
      o.dn[k] = t.dn[2,2]
  }
      df.i.up = data.frame(
        celltype = types[i] ,
        set = paste('M',1:m,sep='') ,
	contrast = contrast ,
        dir = 'up' ,
        nDEGs = nUp ,
        nSetGenes ,
        nOverlap = o.up ,
        OR = or.up ,
        pval = p.up )
      df.i.dn = data.frame(
        celltype = types[i] ,
        set = paste('M',1:m,sep='') ,
	contrast = contrast ,
        dir = 'down' ,
        nDEGs = nDown ,
        nSetGenes ,
        nOverlap = o.dn ,
        OR = or.dn ,
        pval = p.dn )
  df.i = rbind( df.i.up , df.i.dn )
  df[[i]] = df.i
 }
 df = do.call( rbind , df )
 df = df[ order( df$pval ) , ]
 df$fdr = p.adjust( df$pval , method = 'fdr' )
 df$bonf = p.adjust( df$pval , method = 'bonf' )
 df.all[[j]] = df
}

df.all = do.call( rbind , df.all )
df.all = df.all[ order( df.all$pval ) , ]

write.table( df.all , quote=F , sep='\t' , row.names=F , 
             file = 'tbi.module.deg.enrichment.fisher.txt' )