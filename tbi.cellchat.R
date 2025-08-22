#devtools::install_github("jinworks/CellChat")
library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
# reticulate::use_python("/Users/suoqinjin/anaconda3/bin/python", required=T)
library(NMF)
library(circlize)
library(ComplexHeatmap)
library(ggalluvial)

objs = readRDS('allcellsSansMixed_seurat.rds')
data.in <- objs[["RNA"]]@data # normalized data matrix
meta = objs@meta.data# a dataframe with rownames containing cell meta data
meta$labels = Idents(objs)

my_cols2 = c("#045C4A", "#999933", "#DDCC77",
            "#882255", "#6C37C7",
            "#AA4499", "#6699CC", "#88CCEE",
            "#44AA99", "#CC6677", "#AA4466",
            "#332288", "#117733", "#47910E",
            "#661100", "#51367A", "#A04295",
            "#E2B5D5", "#6699CC", "#FF88FF",
            "#A1E2C7")
names(my_cols2) = c("IN-Lamp5", "IN-Sst", "IN-Vip", "EN-CajalRhetzius", "EN-Nxph3", "EN-L3EC",
                    "EN-CA1", "EN-CA2CA3", "Astro", "Ependyma", "CP", "EndoMural",
                    "OPC", "Oligo", "Micro", "dg6", "dg5", "dg4", "dg3", "dg2", "dg1")
my_cols2 = my_cols2[c(9,4,7,8,11,21,20,19,18,17,16,12,10,6,1,15,5,14,13,2,3)]

##############################################################################
##############################################################################
for (grp_name in unique(meta$group_age)) {
  cell.use = rownames(meta)[meta$group_age == grp_name] # extract the cell names from disease data
  # Subset the input data for CellChat analysis
  data.input = data.in[, cell.use]
  meta2 = meta[cell.use, ]
  meta2$labels = meta2$identzz
  meta2$samples = meta2$sample_name #cellchat will look for a column in meta called samples

  # create a new CellChat object from Seurat obj
  cellchat <- createCellChat(object = data.input, meta = meta2, group.by = "labels")

  # subset the expression data of signaling genes for saving computation cost
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multisession", workers = 6) # do parallel
  options(future.globals.maxSize= 1572864000) #891289600 850 MiB
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  levels(cellchat@idents)
  groupSize <- as.numeric(table(cellchat@idents))

  ptm = Sys.time()
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  # triMean is used for calculating the average gene expression per cell group

  cellchat <- filterCommunication(cellchat, min.cells = 10)

  df.net <- subsetCommunication(cellchat)
  # all the inferred cell-cell communications at the level of ligands/receptor

  #compute the communication probability on signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)

  cellchat <- aggregateNet(cellchat)
  ptm = Sys.time()
  execution.time = Sys.time() - ptm

  groupSize <- as.numeric(table(cellchat@idents))
  pdf(paste0("cellchat_circleplots_", grp_name, ".pdf"), width = 12, height = 6)
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                   weight.scale = T, label.edge= F,
                   color.use = my_cols2,
                   title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                   color.use = my_cols2,
                   weight.scale = T, label.edge= F,
                   title.name = "Interaction weights/strength")
  dev.off()

  mat <- cellchat@net$weight

  png(paste0("tbi_cellchatpercell_", grp_name, "1.png"), width = 13, height = 13, units = 'in', res = 300)
  par(mfrow = c(4,4), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, color.use = my_cols2, vertex.weight = groupSize,
                     weight.scale = T, edge.weight.max = max(mat),
                     title.name = rownames(mat)[i])
  }
  dev.off()

  # show all the significant interactions (L-R pairs) associated with certain signaling pathways
  # dg to 6 targets
  png(paste0("chordgene_NRG_", grp_name, "_DGto6sources.png"), width = 11, height = 12, units = 'in', res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(6:11),
                       targets.use = c(1,15,18:21),
                       signaling = "NRG",
                       legend.pos.x = 8,
                       color.use = my_cols2)
  dev.off()

  png(paste0("chordgene_SEMA3_", grp_name, "_6sources.png"), width = 11, height = 12, units = 'in', res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(1,15,18:21),
                       targets.use = c(6:11),
                       signaling = "SEMA3",
                       color.use = my_cols2,
                       lab.cex = 0.5,
                       show.legend = F)
  dev.off()

  # show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
  png(paste0("cellchatDGchordgene_", grp_name, "_6sourcestoDG.png"), width = 11, height = 12, units = 'in', res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(1,15,18:21), targets.use = c(6:11),
                       slot.name = "netP", legend.pos.x = 10, color.use = my_cols2)
  dev.off()

  png(paste0("cellchatDGchordgene_", grp_name, "_4sourcestoDG.png"), width = 11, height = 12, units = 'in', res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(15,19:21), targets.use = c(6:11),
                       slot.name = "netP", legend.pos.x = 10, color.use = my_cols2)
  dev.off()

  pairLR.use <- extractEnrichedLR(cellchat, signaling = c("GABA-A", "GABA-B"))
  png(paste0("cellchatDGchordgene_", grp_name, "_6sourcestoDG_gaba2.png"), width = 11, height = 12, units = 'in', res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(1,15,18:21), targets.use = c(6:11),
                       pairLR.use = pairLR.use, slot.name = "netP", legend.pos.x = 10, color.use = my_cols2)
  dev.off()

  png(paste0("cellchatDGchordgene_", grp_name, "_sourcestoDG_gaba3_test2.png"), width = 11, height = 12, units = 'in', res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(1:21), targets.use = c(1:21),
                       signaling = c("GABA-A", "GABA-B"),
                       slot.name = "netP", show.legend = F, color.use = my_cols2)
  dev.off()
  # The following line references cellchat3, which may not exist in this loop context. If needed, adapt accordingly.
  # png(paste0("cellchatDGchordgene_sham3mo_sourcestoDG_gaba2_test.png"), width = 11, height = 12, units = 'in', res = 300)
  # netVisual_chord_gene(cellchat3, sources.use = c(1:21), targets.use = c(1:21),
  #                      signaling = c("GABA-A", "GABA-B"), slot.name = "netP", show.legend = F, color.use = my_cols2)
  # dev.off()

  png(paste0("cellchatDGchordgene_", grp_name, "_sourcestoDG_gabaglutNRG.png"), width = 11, height = 12, units = 'in', res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(1,2,9:21), targets.use = c(3:8),
                       signaling = c("Glutamate", "GABA-A", "GABA-B", "NRG"),
                       slot.name = "netP", show.legend = F, color.use = my_cols2)
  dev.off()

  png(paste0("cellchatDGchordgene_", grp_name, "_sourcestoDG_glut.png"), width = 11, height = 12, units = 'in', res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(1,2,9:21), targets.use = c(3:8),
                       signaling = c("Glutamate"), slot.name = "netP", show.legend = F, color.use = my_cols2)
  dev.off()

  png(paste0("cellchatDGchordgene_", grp_name, "_DGto6sources.png"), width = 11, height = 12, units = 'in', res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(3:8), targets.use = c(1,2,9:21),
                       slot.name = "netP", legend.pos.x = 10, color.use = my_cols2)
  dev.off()
  png(paste0("cellchatDGchordgene_", grp_name, "_DGto6sources_gaba.png"), width = 11, height = 12, units = 'in', res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(3:8), targets.use = c(1,2,9:21),
                       signaling = c("GABA-A", "GABA-B"), slot.name = "netP", legend.pos.x = 10, color.use = my_cols2)
  dev.off()
  png(paste0("cellchatDGchordgene_", grp_name, "_DGtosources_glut.png"), width = 11, height = 12, units = 'in', res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(3:8), targets.use = c(1,2,9:21),
                       signaling = c("Glutamate"), slot.name = "netP", show.legend = F, color.use = my_cols2)
  dev.off()

  png(paste0("cellchatDGchordgene_", grp_name, "_DGto6sources_gaba.png"), width = 11, height = 12, units = 'in', res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(6:11), targets.use = c(1,15,18:21),
                       pairLR.use = pairLR.use, slot.name = "netP", legend.pos.x = 10, color.use = my_cols2)
  dev.off()

  png(paste0("cellchatDGchordgene_", grp_name, "_DGto4sources.png"), width = 11, height = 12, units = 'in', res = 300)
  netVisual_chord_gene(cellchat, sources.use = c(6:11), targets.use = c(15,19:21),
                       slot.name = "netP", legend.pos.x = 10, color.use = my_cols2)
  dev.off()

  pdf(paste0(grp_name, "_GeneExpression_PTN.pdf"), width = 8, height = 5)
  plotGeneExpression(cellchat, signaling = "PTN", enriched.only = TRUE, type = "violin", color.use = my_cols2)
  dev.off()
  pdf(paste0(grp_name, "_GeneExpression_SLITRK.pdf"), width = 8, height = 5)
  plotGeneExpression(cellchat, signaling = "SLITRK", enriched.only = TRUE, type = "violin", color.use = my_cols2)
  dev.off()
  pdf(paste0(grp_name, "_GeneExpression_SLIT.pdf"), width = 8, height = 5)
  plotGeneExpression(cellchat, signaling = "SLIT", enriched.only = TRUE, type = "violin", color.use = my_cols2)
  dev.off()
  pdf(paste0(grp_name, "_GeneExpression_SEMA3.pdf"), width = 8, height = 5)
  plotGeneExpression(cellchat, signaling = "SEMA3", enriched.only = TRUE, type = "violin", color.use = my_cols2)
  dev.off()
  pdf(paste0(grp_name, "_GeneExpression_CypA.pdf"), width = 8, height = 5)
  plotGeneExpression(cellchat, signaling = "CypA", enriched.only = TRUE, type = "violin", color.use = my_cols2)
  dev.off()
  pdf(paste0(grp_name, "_GeneExpression_FGF.pdf"), width = 8, height = 5)
  plotGeneExpression(cellchat, signaling = "FGF", enriched.only = TRUE, type = "violin", color.use = my_cols2)
  dev.off()
  pdf(paste0(grp_name, "_GeneExpression_IGF.pdf"), width = 8, height = 5)
  plotGeneExpression(cellchat, signaling = "IGF", enriched.only = TRUE, type = "violin", color.use = my_cols2)
  dev.off()

  ptm = Sys.time()
  # Compute the network centrality scores
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
  cellchat@netP$pathways
  pathways.show = c("CypA", "FGF", "IGF", "NRG", "PTN", "SEMA3", "SLIT", "SLITRK")
  pdf(paste0(grp_name, "_SignalingRoleNetwork_pathways.show.pdf"), width = 9, height = 2.5)
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show,
                                    font.size = 10,
                                    color.use = my_cols2, color.heatmap = "PuBuGn")
  dev.off()

  pathways.show2 = c("WNT", "TGFb", "COMPLEMENT")
  pdf(paste0(grp_name, "_SignalingRoleNetwork_pathways.show2_WNT_TGFb_COMPLEMENT.pdf"), width = 9, height = 2.5)
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show2,
                                    font.size = 10,
                                    color.use = my_cols2, color.heatmap = "PuBuGn")
  dev.off()

  pathways.show2 = c("GRN", "ncWNT", "ACTIVIN", "OPIOID", "ANGPTL", "EGF")
  pdf(paste0(grp_name, "_SignalingRoleNetwork_pathways.show2_GRN_ncWNT_ACTIVIN_OPIOID_ANGPTL_EGF.pdf"), width = 9, height = 2.5)
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show2,
                                    font.size = 10,
                                    color.use = my_cols2, color.heatmap = "PuBuGn")
  dev.off()

  selectK(cellchat, pattern = "outgoing")
  nPatterns = 3
  dev.off()
  cellchat <- identifyCommunicationPatterns(cellchat,
                                            pattern = "outgoing",
                                            k = nPatterns,
                                            color.use = my_cols2,
                                            height = 9)
  # river plot outgoing comm
  pdf(paste0("cellchat_riverplot_", grp_name, "_outgoing.pdf"), width = 9, height = 6)
  netAnalysis_river(cellchat, pattern = "outgoing", color.use = c(my_cols2))
  dev.off()

  selectK(cellchat, pattern = "incoming")
  nPatterns = 4
  dev.off()
  cellchat <- identifyCommunicationPatterns(cellchat,
                                            pattern = "incoming",
                                            k = nPatterns,
                                            color.use = my_cols2,
                                            height = 9)
  # river plot incoming comm
  pdf(paste0("cellchat_riverplot_", grp_name, "_incoming.pdf"), width = 9, height = 6)
  netAnalysis_river(cellchat, pattern = "incoming", color.use = c(my_cols2))
  dev.off()

  # save cellchat object
  write.csv(df.net, paste0('cellchat_allLRcomms_', grp_name, '_subDG.csv'))
  saveRDS(cellchat, paste0('cellchat_', grp_name, '_dgsubs.rds'))
}
