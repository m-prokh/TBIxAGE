# conda activate /conda_envs/seurat5

library( Seurat )
library( edgeR )
library( Signac )

obj = readRDS('~/TBI/CellRanger5/MoTBI.AllSamples.Integrated.Cleaned.Annotated.CellRanger.rds')

DefaultAssay(obj) = 'RNA'

types = levels(obj)

elist = list()

samples = unique( obj$sample_name )

for( i in 1:length(types) ) {
  cat( types[i] , '\n' )
  expr = matrix( 0 , ncol = length(donors) , nrow = nrow(obj) )
  sub = subset( obj , idents = types[i] )
  for( j in 1:length(samples) ) {
    idx = which( sub$sample_name == samples[j] )
    if(length(idx)==0) next
    sub.j = subset( sub , cells = idx )
    expr[,j] = rowSums( GetAssayData( sub.j , slot = 'counts' ) )
    rownames(expr) = rownames(obj)
    colnames(expr) = samples
  }
  elist[[i]] = expr
}

names(elist) = types

saveRDS( elist , file = 'MoTBI.AllSamples.Integrated.Cleaned.Annotated.CellRanger.pseudobulk_counts.rds')


