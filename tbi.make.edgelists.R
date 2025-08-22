# create list of pairwise correlations among the top 50 genes in each module

library( Seurat )
obj = readRDS("tbi/tbi.microglia.w26MEs.rds")
colors = readRDS('mic.k50.merged.clusters.rds')
MM = readRDS('mic_module_membership.kAE.rds')

obj = NormalizeData(obj)
datExpr = t(as.matrix(GetAssayData(obj)))

n = max(colors)
nGenes = 75
rthresh = 0.5
edgelists = list()
for( i in 1:n ) {
  mod = names(colors)[ colors == i ]
  MM.mod = MM[ mod , ]
  ranked = MM.mod[ order( MM.mod[,paste('M',i,sep='')] ,
			  decreasing = T ) , ]
  top = head( rownames(ranked) , nGenes )
  r = cor( datExpr[ , top ] )
  diag(r) = 0
  r[ lower.tri(r) ] = 0
  e = which( r > rthresh , arr.ind = T )
  edges = data.frame(
    		gene1 = colnames(r)[e[,1]] ,
		gene2 = colnames(r)[e[,2]] ,
		weight = r[ e ] )
  write.table( edges , sep = '\t' , quote=F , row.names=F ,
	       file = paste( 'mic.net.edgelist.top75genes.r0.5.M' , i , '.txt' , 
			    sep = '' ) )
  edgelists[[i]] = edges
}
names(edgelists) = paste('M',1:n,sep='')

saveRDS( edgelists , file = 'mic.net.edgelists.top75genes.r0.5.rds' )


