library(monocle3)
library(dplyr)
DefaultAssay(obj) <- "RNA"
datExpr = obj@assays$RNA@counts
annot = as.factor(rownames(datExpr), )
annot = data.frame(gene_short_name = as.factor(rownames(datExpr)), row.names=rownames(datExpr))
meta$sample_id = meta$sample_name #meta can't have sample_name as variable
meta$cell_name = rownames(meta)
meta = meta[,-4]
dgT <- new_cell_data_set(datExpr,
                         cell_metadata = meta,
                         gene_metadata = annot)
mrkrlist0 = c("Marc2", "Prox1", "Sox2", "Rbfox1", "Rbfox3", "Vangl2", "Gfap", 
              "Eomes", "Sox4", "Dcx", "Calb2", "Calb1", "Bend6", "Plppr4", "Fgf13",
              "Grm7", "Robo2", "Slc17a7", "Kcnq5", "Kcnj6", "Gabra1", "Grin2a",
              "Grin1","Grin2b", "Gabbr1", "Gabbr2", "Gria2", "Snap25", "Tuba1a",
              "Cobl", "Cryab", "Lpar1", "Ncam1", "Igfbp7", "Prrx1", "Slc1a2", 
              "Slc1a3", "Cd44", "Aqp4", "Sox5", "Sox6")
mrkrlist1 = c("Calb1", "Bcl11b", "Kcnq5",
             "Kcnj6", "Gabra1", "Gria2", "Grin2a", "Grin1", "Prox1", "Rbfox3", "Rbfox1","Slc1a3",
             "Dgkh", "Ncam1", "Cntnap5a", "C1ql2", "Grm7", "Robo2", "Penk",
             "Arc", "Ptgs2", "Gabbr2", "Gabbr1", "Dcx", "Fgf13", "Plppr4", "Slc1a2", 
             "Slc17a7", "Grin2b", "Cobl", "Snap25", "Bend6") #final
mrkrlist2 = c("Bcl11b", "Kcnq5", "Kcnj6", "Prox1", "Rbfox1", "Slc1a3", "Ncam1",
              "Cntnap5a", "C1ql2", "Grm7", "Robo2", "Penk", "Gabbr1", "Grin2b", 
              "Snap25", "Bend6")
mrkrlist3 = c("Bcl11b", "Kcnj6", "Prox1", "Rbfox1", "Slc1a3", "Dcx", "Ascl1",
              "Cntnap5a", "C1ql2", "Grm7", "Robo2", "Penk", "Gabbr1", "Bend6")
## Step 1: Normalize and pre-process the data
dgT <- preprocess_cds(dgT, num_dim = 100)
plot_pc_variance_explained(dgT)

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(dgT, preprocess_method = 'PCA', reduction_method = 'UMAP')
plot_cells(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "identz")

## Step 4: Cluster the cells
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by="identz")
plot_cells(cds, genes=mrkrlist2)
## Step 5: Learn a graph
cds <- learn_graph(cds)

## Step 6: Order cells
cds <- order_cells(cds)

plot_cells(cds)

######## with LSI preprocessing
cds <- preprocess_cds(dgT, method = "LSI", num_dim = 100)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds, preprocess_method = 'LSI', reduction_method = 'UMAP')
plot_cells(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "identz")
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by="identz")
plot_cells(cds, genes=mrkrlist2)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds)

######## with louvain clustering
cds2 <- reduce_dimension(dgT, preprocess_method = 'PCA', reduction_method = 'UMAP')
cds2 <- cluster_cells(cds2, reduction_method = "UMAP", 
                      cluster_method = "louvain", 
                      k=50, 
                      num_iter = 10, verbose=T)

plot_cells(cds2, label_groups_by_cluster=FALSE,  color_cells_by = "identz")
cds2 <- cluster_cells(cds2)
cds2 <- learn_graph(cds2)
plot_cells(cds2, color_cells_by="identz")
plot_cells(cds2, genes=mrkrlist2)
cds2 <- order_cells(cds2)
plot_cells(cds2)
saveRDS(cds2, "pseudoDG_pca100_umap_louvain.rds")

######## with leiden clustering at diff resolutions
cds <- reduce_dimension(dgT, preprocess_method = 'PCA', reduction_method = 'UMAP')
cds2 <- cluster_cells(cds, reduction_method = "UMAP", 
                      k=30, 
                      num_iter = 10,
                      resolution = 1e-3,
                      verbose=T)
cds2 <- learn_graph(cds2)
plot_cells(cds2, color_cells_by="identz")
cds2 <- order_cells(cds2)
plot_cells(cds2, color_cells_by="identz")
saveRDS(cds2, "pseudoDG_pca100_umap_leivenk30res001.rds")

########
########
#########
detach("package:monocle3", unload=TRUE)
library(monocle)
library(dplyr)
cds4 <- newCellDataSet(datExpr,
                       phenoData = dg_meta,
                       featureData = annot)



cds4 <- reduce_dimension(dgT, preprocess_method = 'PCA', reduction_method = 'UMAP')

## Step 4: Cluster the cells
cds4 <- cluster_cells(cds4)

## Step 5: Learn a graph


cds4 <- partitionCells(cds4)
cds4 <- learnGraph(cds4,
                  max_components = 3,
                  RGE_method = 'SimplePPT',
                  partition_component = F,
                  verbose = F)

cds4 <- learnGraph(cds,  RGE_method = 'DDRTree', do_partition=FALSE)
cds4 <- order_cells(cds4)


marker_test_res <- top_markers(cds2, group_cells_by="partition", 
                               reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds2,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)
plt2 = plot_cells(cds2,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 1) + 
  viridis::scale_color_viridis()
plt3 = plot_cells(cds2, color_cells_by="identz", 
           label_branch_points = FALSE,
           cell_size = 1,
           label_cell_groups = F) + 
  scale_color_manual(values = dg_cols2)
plt3 = plot_cells(cds2, reduction_method = "Aligned",
                  color_cells_by="identz", 
                  label_branch_points = FALSE,
                  cell_size = 1,
                  label_cell_groups = F) + 
  scale_color_manual(values = dg_cols2)

plt4 = plot_cells(cds2, color_cells_by="group_age", 
           label_branch_points = FALSE,
           cell_size = 1,
           label_cell_groups = F) + 
  scale_color_manual(values = c("#403872FF", "#60CEACFF","#BD1655FF","#F6B893FF"))
fill3 = c("#49C1ADFF","#3E356BFF", "#F4875EFF", "#961C5BFF")
library(ragg)
agg_png("dentatepsed2.png", width = 8.4, height = 6.8, units = "in", res = 300)
plt3
dev.off()
png("dentate_monocleumap_pseudotime.png", width = 8.4, height = 6.8, units = "in", res = 300)
plt2
dev.off()
png("dentate_monocleumap_groupage.png", width = 8.4, height = 6.8, units = "in", res = 300)
plt4
dev.off()
################################################################################
################################################################################
p3 <- plot_cell_trajectory(cds2, color_by = "subtypes") + #, theta = -45
  scale_color_manual(values = dg_cols2)
png("mic_monoclepseudotime.png", width = 6, height = 4, units = 'in', res = 300)
p3
dev.off()

plot_cell_trajectory(mic, color_by = "State")

p3 <- plot_cell_trajectory(mic, color_by = "Pseudotime")
png("mic_monoclepseudotime_bypseudotime.png", width = 6, height = 4, units = 'in', res = 300)
p3
dev.off()

fill2 = c("#F4875EFF", "#49C1ADFF", "#961C5BFF", "#3E356BFF")
p3 <- plot_cell_trajectory(mic, color_by = "grp_age") + 
  scale_color_manual(values = fill2)
png("mic_monoclepseudotime_byGrpAge.png", width = 6, height = 4, units = 'in', res = 300)
p3
dev.off()


dg_cols2 = c("dg6" = "#51367A", 
             "dg5" = "#A04295", 
             "dg4" = "#E2B5D5", 
             "dg3" = "#6699CC",
             "dg2" = "#FF88FF",
             "dg1" = "#A1E2C7")

plt5 = plot_genes_in_pseudotime(cds2[mrkrlist3,], 
                         color_cells_by="identz",
                         ncol = 7,
                         panel_order = mrkrlist3)  + 
  scale_color_manual(values = dg_cols2)

png("dg_pseudotimegenetrajectories2.png", width = 12, height = 3.2, units = 'in', res = 300)
plt5
dev.off()



###########
cds3 <- reduce_dimension(dgT, preprocess_method = 'PCA', 
                        max_components = 3,
                        reduction_method = 'UMAP')
cds3 <- cluster_cells(cds3, reduction_method = "UMAP", 
                      k=30, 
                      num_iter = 10,
                      resolution = 1e-3,
                      verbose=T)
cds3 <- learn_graph(cds3)
plot_cells_3d(cds3, color_cells_by="identz")
plot_cells(cds3, color_cells_by="identz", cell_size = 0.5) + 
  scale_color_manual(values = dg_cols2)
cds3 <- order_cells(cds3)


dg_cols2 = c("dg6" = "#51367A", 
             "dg5" = "#A04295", 
             "dg4" = "#E2B5D5", 
             "dg3" = "#6699CC",
             "dg2" = "#FF88FF",
             "dg1" = "#A1E2C7")
plot_cells_3d(
  cds3,
  dims = c(1, 2, 3),
  color_cells_by = "identz",
  cell_size = 50,
  color_palette = dg_cols2,
  reduction_method = c("UMAP"))
