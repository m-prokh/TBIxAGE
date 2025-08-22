# test for overlap between modules and DEGs in each cell type
library( limma )

# load modules
dgColors <- readRDS("~/tbi/k50.merged.clusters.DG.mm22r15ch32.rds")
# load DEGs
  
dge_DG = read.delim('pseudobulk.mdgeR.rmOutliers.dge.txt')
dge = dge_DG
#temporarily
color.list[[1]] = mergeColors 
types = cellt
names(color.list) = types
n = length(types)

df.all = list()
for( j in 1:6 ) {
 contrast = c('bothgroups.3movs18mo', 'sham.3movs18mo', 'TBI.3movs18mo', 'TBIvssham.bothages', 'TBIvssham.3mo', 'TBIvssham.18mo')[j]
 df = list()
 dge.i = dge[ dge$Contrast == contrast, ]
 dge.i$score = sign(dge.i$logFC) * -log10(dge.i$PValue)
 colors = dgColors
 m = max(colors)
 p.dn = rep(NA,m)
 p.up = rep(NA,m)
 nGenes = rep(NA,m)
 for( k in 1:m ) {
   set = names(colors)[colors==k]
   idx = dge.i$Gene %in% set
   nGenes[k] = length(set)
   p.up[k] = geneSetTest(
                statistics = dge.i$score ,
                index = idx ,
                alternative = 'up' )
   p.dn[k] = geneSetTest(
                statistics = dge.i$score ,
                index = idx ,
                alternative = 'down' )
   df.i.up = data.frame(
        celltype = "EN-DG" ,
        set = paste('M',1:m,sep='') ,
	contrast = contrast ,
        dir = 'up' ,
        nGenes ,
        pval = p.up )
      df.i.dn = data.frame(
        celltype = "EN-DG",
        set =  paste('M',1:m,sep=''),
	contrast = contrast ,
        dir = 'down' ,
        nGenes ,
        pval = p.dn )
      df.i = rbind( df.i.up , df.i.dn )
  df[[i]] = df.i
 }
 df = do.call( rbind , df )
 df$fdr = p.adjust( df$pval , method = 'fdr' )
 df$bonf = p.adjust( df$pval , method = 'bonf' )
 df = df[ order( df$pval ) , ]
 df.all[[j]] = df
}

df.all = do.call( rbind , df.all )
df.all = df.all[ order( df.all$pval ) , ]
 
write.table( df.all , quote=F , sep='\t' , row.names=F , 
	     file = 'tbi.DG_module_deg_enrichment.geneSetTest.txt' )




