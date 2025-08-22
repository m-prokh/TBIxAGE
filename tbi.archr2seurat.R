### archr2seurat
packages <- c("ArchR", "ArchRtoSignac", "Seurat", "Signac", "stringr")
loadinglibrary(packages)
lapply(packages, library, character.only=T)
tbi <- loadArchRProject(path = "/TBI/snATAC/Aggr_ArchR/Proj3-ArchRProj", force = FALSE, showLogo = FALSE)
pkm <- getPeakMatrix(tbi)
## Extract Ensembl gene annotation and convert to UCSC style if needed
library(EnsDb.Mmusculus.v79) # Ensembl database to convert to human hg38. Install what is appropriate for your analysis
annotations <- getAnnotation(reference = EnsDb.Mmusculus.v79, refversion = "mm10")
fragments_dirs <- list(
  "/local/projects-t3/idea/parthbhatia/TBI/snATAC/CellRangerV2/TBI_18Mo_12/",
  "/local/projects-t3/idea/parthbhatia/TBI/snATAC/CellRangerV2/TBI_3Mo_23/",
  "/local/projects-t3/idea/parthbhatia/TBI/snATAC/CellRangerV2/TBI_18Mo_8/",
  "/local/projects-t3/idea/parthbhatia/TBI/snATAC/CellRangerV2/Sham_18Mo_18/",
  "/local/projects-t3/idea/parthbhatia/TBI/snATAC/CellRangerV2/Sham_3Mo_36/",
  "/local/projects-t3/idea/parthbhatia/TBI/snATAC/CellRangerV2/TBI_3Mo_227/",
  "/local/projects-t3/idea/parthbhatia/TBI/snATAC/CellRangerV2/Sham_3Mo_31/",
  "/local/projects-t3/idea/parthbhatia/TBI/snATAC/CellRangerV2/Sham_18Mo_14/"
)
# created ArchR2Signac2 function by modifying code for ArchR2Signac
seurat_atac <- ArchR2Signac2(
  ArchRProject = tbi,
  refversion = "mm10",
  fragments_dir = fragments_dirs,
  pm = pkm, 
  fragments_fromcellranger = "Yes",
  fragments_file_extension = NULL,
  annotation = annotations
)
saveRDS(seurat_atac, file = "tbi_seurat_from_archr.rds")

## Transfer ArchRProject gene score matrix to Signac SeuratObject
gsm <- getGeneScoreMatrix(ArchRProject = tbi, SeuratObject = seurat_atac)
seurat_atac[['RNA']] <- CreateAssayObject(counts = gsm)
saveRDS(gsm, file = "gsm_from_archr.rds")

## Transfer ArchRProject dimension reduction ("IterativeLSI" and "Harmony") and UMAP to Signac SeuratObject
seurat_atac <- addDimRed(
  ArchRProject = tbi,
  SeuratObject = seurat_atac,
  addUMAPs = "UMAP",
  reducedDims = "IterativeLSI"
) # default is "IterativeLSI"

#add both 'Harmony' and ‘IterativeLSI’:
seurat_atac <- addTwoDimRed(
  ArchRProject = tbi,
  SeuratObject = seurat_atac,
  addUMAPs = "UMAP",
  reducedDims1 = "IterativeLSI",
  reducedDims2 = "Harmony" # IterativeLSI2 or Harmony
)
saveRDS(seurat_atac, file="tbi_seurat_from_archr_wLSIwHarmony.rds")
save.image()