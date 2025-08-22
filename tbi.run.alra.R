library( Seurat )
library(rsvd)
library( edgeR )
#import single cell pseudobulk counts as objs
# #objs = readRDS('MoTBI.AllSamples.Integrated.Cleaned.Annotated.CellRanger.rds')
objs = readRDS('/tbi/MoTBI.AllSamples.Integrated.Cleaned.Annotated.CellRanger.pseudobulk_counts.rds')

### ALRA (Kluger method)
(.packages())
source('/ALRA/alra.R')

#types = unique( objs$CellType )
types = unique(Idents(objs))
types = names(objs)
n = length(types)

alra.list = list()
for( i in 1:length(types) ) {
  cat( i , '\n' )
  tmp = subset( obj , idents = types[i] )
  norm = t(as.matrix(LayerData(tmp,layer='data')))
  alra.out = alra( norm )
  alra.list[[i]] = alra.out[[3]]
}

meta = readRDS('meta_for_alra.rds')
meta = obj@meta.data
mm = meta[ colnames(obj) , ]
DefaultAssay(objs) = 'RNA'
for( i in 1:n ) {
  cat( types[i] , '\n' )
  obj = subset( objs, idents = types[i] ) #obj = objs[[i]]
  if( ncol(obj) > 25000 ) {obj = subset( obj , downsample = 25000 )}
  obj$percent.mt = PercentageFeatureSet( obj , pattern = '^MT-' )
  norm = t(as.matrix(GetAssayData(obj, assay = "RNA")))
  pct = colSums( norm > 0 ) / nrow(norm)
  norm = norm[ , pct > 0.01 ]
  res = alra( norm )
  alra.out = res[[3]]
  rownames(alra.out) = colnames(obj)
  meta = obj@meta.data
  alra.obj = CreateSeuratObject( t(alra.out) , meta = meta)
  alra.obj = NormalizeData( alra.obj )
  Idents(alra.obj) = types[i]
  alra.obj = FindVariableFeatures(alra.obj)
  alra.obj = ScaleData( alra.obj)
  alra.obj = RunPCA(alra.obj)
  alra.obj = RunUMAP( alra.obj , dims = 1:10 )
  alra.obj = FindNeighbors( alra.obj , dims = 1:10 )
  alra.obj = FindClusters( alra.obj)
  alra.list[[i]] = alra.obj
}
#meta=obj@meta.data
#norm=t(as.matrix(LayerData(obj, layer='data')))
names(alra.list) = types

#########
norm = t(as.matrix(GetAssayData(obj)))
pct = colSums(norm>0)/nrow(norm)
norm = norm[ , pct > 0.03 ]

res = alra( norm )
alra.out = res[[3]]
rownames(alra.out) = colnames(obj)

m = obj@meta.data
alra = CreateSeuratObject( counts = t(alra.out) , meta.data = m )
alra = NormalizeData( alra )
Idents(alra) = Idents(obj)

alra = FindVariableFeatures(alra)
alra = ScaleData( alra )
alra = RunPCA(alra)
alra = RunUMAP( alra , dims = 1:10 )#, umap.method = "umap-learn",
alra = FindNeighbors( alra , dims = 1:10 )
alra = FindClusters( alra )


saveRDS( alra , file = 'tbi.rna.pseudobulk.allcells.ALRA.seurat.rds' )

#### find markers and plot features

# Define cell types and associated features for plotting
cell_types_features <- list(
  "Microglia" = c('Siglech','Nav3', 'C1qb', 'Trem2', 'Apoe', 'Cacna1a', 'Mrc1', 'Dab2', 'H2-Eb1', 'Cd74', 'Grn'),
  "Neuron_Dentate_C1ql2" = c('Marc2','Rbfox3', 'Dcx', 'Prox1', 'Calb1')
)

# Define cell types for marker finding and output file names
cell_types_files <- list(
  "Microglia" = "tbi.Microglia.ALRA.markers.txt",
  "Neuron_Dentate_C1ql2" = "tbi.Dentate_C1ql2.ALRA.markers.txt",
  "Neuron_CajalRhetzius_Lhx1" = "tbi.Neuron_CajalRhetzius_Lhx1.ALRA.markers.txt",
  "Neuron_Subiculum_Slc17a6" = "tbi.Neuron_Subiculum_Slc17a6.ALRA.markers.txt",
  "Endothelial-Mural-Fibroblast" = "tbi.Endothelial-Mural-Fibroblast.ALRA.markers.txt",
  "Neuron_Subiculum_Entorhinal_Nxph3" = "tbi.Neuron_Subiculum_Entorhinal_Nxph3.ALRA.markers.txt",
  "Neuron_CA1_Subiculum_Postsubiculum_Entorhinal" = "tbi.Neuron_CA1_Subiculum_Postsubiculum_Entorhinal.ALRA.markers.txt",
  "Neuron_CA2CA3_Pvrl3-Rgs15-Calb2" = "tbi.Neuron_CA2CA3_Pvrl3-Rgs15-Calb2.ALRA.markers.txt",
  "Astrocytes" = "tbi.Astrocytes.ALRA.markers.txt",
  "Oligodendrocytes" = "tbi.Oligodendrocytes.ALRA.markers.txt",
  "Polydendrocytes" = "tbi.Polydendrocytes.ALRA.markers.txt",
  "Interneuron_Gad2-Lamp5" = "tbi.Interneuron_Gad2-Lamp5.ALRA.markers.txt",
  "Interneuron_Gad2-Sst" = "tbi.PInterneuron_Gad2-Sst.ALRA.markers.txt",
  "Interneuron_Gad2-Vip" = "tbi.Interneuron_Gad2-Vip.ALRA.markers.txt"
)

# Plot and marker finding for cell types with features
for (ct in names(cell_types_features)) {
  alra <- alra.list[[ct]]
  markers <- FindAllMarkers(alra, only.pos = TRUE, max.cells.per.ident = 1000)
  features <- cell_types_features[[ct]]
  pdf(paste0("tbi.ALRA.UMAP.", ct, "features.pdf"))
  print(DimPlot(alra, label = TRUE) + NoLegend())
  print(FeaturePlot(alra, label = TRUE, combine = FALSE, features = features))
  # If Microglia or DG, do blend plot for specific features
  if (ct == "Microglia") {
    print(FeaturePlot(alra, label = TRUE, combine = FALSE, features = c('Apoe', 'Grn'), blend = TRUE))
  }
  if (ct == "Neuron_Dentate_C1ql2") {
    print(FeaturePlot(alra, label = TRUE, combine = FALSE, features = c('Dcx', 'Calb1'), blend = TRUE))
  }
  dev.off()
  write.table(markers, quote = FALSE, sep = '\t', row.names = FALSE, file = cell_types_files[[ct]])
}

# Marker finding for remaining cell types
for (ct in setdiff(names(cell_types_files), names(cell_types_features))) {
  alra <- alra.list[[ct]]
  markers <- FindAllMarkers(alra, only.pos = TRUE, max.cells.per.ident = 1000)
  write.table(markers, quote = FALSE, sep = '\t', row.names = FALSE, file = cell_types_files[[ct]])
}
