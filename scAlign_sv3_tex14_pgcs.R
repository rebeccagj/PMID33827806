library(scAlign)
library(Seurat)
library(ggplot2)
library(DittoSeq)
library(SingleCellExperiment)
library(gridExtra)
library(cowplot)
library(dplyr)

# This file works with first set of tex14 fastqs 
# (eg, more shallowly sequenced one, as no difference was 
# found in Seurat analysis of the more deeply sequenced fastqs)

# create seurat objects ####
read_create_meta = function(dirvar,projvar,cellname){
  mat = Read10X(dirvar)
  seuratobj = CreateSeuratObject(counts = mat, project = projvar)
  seuratobj = RenameCells(seuratobj, add.cell.id = cellname)
  seuratobj[["percent.blood"]] = PercentageFeatureSet(object = seuratobj, pattern = "^Hb")
  seuratobj[["percent.mito"]] = PercentageFeatureSet(object = seuratobj, pattern = "^mt-")
  seuratobj[["percent.ribo"]] = PercentageFeatureSet(object = seuratobj, pattern = c("^Rps"))
  seuratobj
}

tex14wt_bh = read_create_meta("tex14wt_bh/outs/filtered_gene_bc_matrices/mm10","tex14wt_bh","1")
tex14ko_bh = read_create_meta("tex14ko_bh/outs/filtered_gene_bc_matrices/mm10","tex14ko_bh","2")

tex14wt_bh
# An object of class Seurat 
# 27998 features across 12585 samples within 1 assay 
# Active assay: RNA (27998 features)
tex14ko_bh
# An object of class Seurat 
# 27998 features across 12585 samples within 1 assay 
# Active assay: RNA (27998 features)

tex14_bh_gcs = readRDS("tex14_bh_gcs.rds") 
# this file has been provided and was generated as previously described in 
# tex14_wt_seuratv3.R and tex14_ko_seuratv3.R, lines 110-118

tex14wt_bh_pgcs = subset(tex14wt_bh, cells = tex14_gcs)
tex14ko_bh_pgcs = subset(tex14ko_bh, cells = tex14_gcs)

tex14wt_bh_pgcs = NormalizeData(tex14wt_bh_pgcs)
tex14ko_bh_pgcs = NormalizeData(tex14ko_bh_pgcs)

tex14wt_bh_pgcs = FindVariableFeatures(tex14wt_bh_pgcs, selection.method = "vst", nfeatures = 3000)
tex14ko_bh_pgcs = FindVariableFeatures(tex14ko_bh_pgcs, selection.method = "vst", nfeatures = 3000)

tex14wt_bh_pgcs@assays$RNA@var.features = append(tex14wt_bh_pgcs@assays$RNA@var.features, c('Dazl'))
tex14ko_bh_pgcs@assays$RNA@var.features = append(tex14ko_bh_pgcs@assays$RNA@var.features, c('Dazl'))

tex14wt_bh_pgcs <- ScaleData(tex14wt_bh_pgcs, features = rownames(tex14wt_bh_pgcs@assays$RNA@counts))
tex14ko_bh_pgcs <- ScaleData(tex14ko_bh_pgcs, features = rownames(tex14wt_bh_pgcs@assays$RNA@counts))

genes.use = Reduce(intersect, 
                   list(VariableFeatures(tex14wt_bh_pgcs),
                        VariableFeatures(tex14ko_bh_pgcs),
                        rownames(tex14wt_bh_pgcs),
                        rownames(tex14ko_bh_pgcs)))
grep("Dazl", genes.use)

# port out to scAlign ####
seuratv3_to_scefxn = function(seurobj,genelist) {
  sce = SingleCellExperiment(
    assays = list(counts = seurobj@assays$RNA@counts[genelist,colnames(seurobj@assays$RNA@scale.data)],
                  logcounts = seurobj@assays$RNA@data[genelist,colnames(seurobj@assays$RNA@scale.data)],
                  scale.data = seurobj@assays$RNA@scale.data[genelist,colnames(seurobj@assays$RNA@scale.data)]),
    colData = seurobj@meta.data[,1:3]
  )
  sce
}

library(SingleCellExperiment)
tex14wt_bh_pgcs_sce = seuratv3_to_scefxn(tex14wt_bh_pgcs,genes.use)
tex14ko_bh_pgcs_sce = seuratv3_to_scefxn(tex14ko_bh_pgcs,genes.use)

scAlign_tex14_bh_pgcs = scAlignCreateObject(sce.objects = list("wt"=tex14wt_bh_pgcs_sce, "ko"=tex14ko_bh_pgcs_sce),
                                            labels = list("wt", "ko"),
                                            data.use="scale.data",
                                            pca.reduce = TRUE,
                                            pcs.compute = 50,
                                            cca.reduce = TRUE,
                                            ccs.compute = 15,
                                            project.name = "scAlign_tex14_bh_pgcs")
saveRDS(scAlign_tex14_bh_pgcs, "scAlign_tex14_bh_pgcs.rds")

results.dir = "."
#dir.create(results.dir)

## Run scAlign with high_var_genes as input to the encoder (alignment) and logcounts with the decoder (projections).
scAlign_tex14_bh_pgcs = scAlign(scAlign_tex14_bh_pgcs,
                                options=scAlignOptions(steps=5000, log.every=5000, norm=TRUE, batch.norm.layer=TRUE, early.stop=FALSE, architecture="small"),
                                encoder.data="scale.data",
                                decoder.data="logcounts",
                                supervised='none',
                                run.encoder=TRUE,
                                run.decoder=TRUE,
                                log.dir=file.path(results.dir, 'models','gene_input'),
                                device="CPU")

## Additional run of scAlign with PCA, the early.stopping heuristic terminates the training procedure too early with PCs as input so it is disabled.
scAlign_tex14_bh_pgcs = scAlign(scAlign_tex14_bh_pgcs,
                                options=scAlignOptions(steps=15000, log.every=1000, norm=TRUE, batch.norm.layer=TRUE, early.stop=FALSE, architecture="small"),
                                encoder.data="PCA",
                                supervised='none',
                                run.encoder=TRUE,
                                run.decoder=FALSE,
                                log.dir=file.path(results.dir, 'models','pca_input'),
                                device="CPU")

## Additional run of scAlign with CCA
scAlign_tex14_bh_pgcs = scAlign(scAlign_tex14_bh_pgcs,
                                options=scAlignOptions(steps=10000, log.every=1000, norm=TRUE, batch.norm.layer=TRUE, early.stop=FALSE, architecture="small"),
                                encoder.data="CCA",
                                supervised='none',
                                run.encoder=TRUE,
                                run.decoder=FALSE,
                                log.dir=file.path(results.dir, 'models','cca_input'),
                                device="CPU")

saveRDS(scAlign_tex14_bh_pgcs, "scAlign_tex14_bh_pgcs.rds")

## Plot aligned data in tSNE space, when the data was processed in three different ways: 1) either using the original gene inputs, 2) after PCA dimensionality reduction for preprocessing, or 3) after CCA dimensionality reduction for preprocessing. Cells here are colored by dataset.
set.seed(1)
gene_plot = PlotTSNE(scAlign_tex14_bh_pgcs, "ALIGNED-GENE", title="scAlign-Gene", perplexity=30)
pca_plot = PlotTSNE(scAlign_tex14_bh_pgcs, "ALIGNED-PCA", title="scAlign-PCA", perplexity=30)
cca_plot = PlotTSNE(scAlign_tex14_bh_pgcs, "ALIGNED-CCA", title="scAlign-CCA", perplexity=30)
wt2ko = PlotTSNE(scAlign_tex14_bh_pgcs, "wt2ko", title = "scAlign-wt2ko", perplexity=30)
ko2wt = PlotTSNE(scAlign_tex14_bh_pgcs, "ko2wt", title = "scAlign-ko2wt", perplexity=30)
legend = get_legend(PlotTSNE(scAlign_tex14_bh_pgcs, "ALIGNED-GENE", title="scAlign-Gene", legend="right", max_iter=1))
combined_plot = grid.arrange(wt2ko, ko2wt, nrow = 2)

# Take these integrations back to Seurat for further analysis ####
scAlign_tex14_bh_pgcs = readRDS("scAlign_tex14_bh_pgcs.rds")
scAlign_tex14_bh_pgcs_seurat = as.Seurat(scAlign_tex14_bh_pgcs, counts = "counts", scale.data = "scale.data")
scAlign_tex14_bh_pgcs_seurat@reductions
saveRDS(scAlign_tex14_bh_pgcs_seurat, "scAlign_tex14_bh_pgcs_seurat.rds")

scAlign_tex14_bh_pgcs_seurat = FindNeighbors(scAlign_tex14_bh_pgcs_seurat, reduction = "ko2wt", dims = 1:2210)
scAlign_tex14_bh_pgcs_seurat = FindClusters(scAlign_tex14_bh_pgcs_seurat, reduction = "ko2wt", resolution = 0.8)
scAlign_tex14_bh_pgcs_seurat = RunTSNE(scAlign_tex14_bh_pgcs_seurat, reduction = "ko2wt", dims = 1:2210)

cexSize = 1
x = DimPlot(scAlign_tex14_bh_pgcs_seurat, 
            group.by = "orig.ident", 
            reduction = "tsne", 
            order = c('tex14ko','tex14wt'), 
            cols = c('royalblue1','black'), 
            pt.size = cexSize) + ggtitle('Genotype')
y = DimPlot(scAlign_tex14_bh_pgcs_seurat, 
            pt.size = cexSize, label = TRUE) + ggtitle('Cluster Resolution 0.8')
z = FeaturePlot(scAlign_tex14_bh_pgcs_seurat, 
                features = "Dazl", 
                pt.size = cexSize)
CombinePlots(list(z,x,y), ncol = 3)
saveRDS(scAlign_tex14_bh_pgcs_seurat, "scAlign_tex14_bh_pgcs_seurat.rds")


# Manually Selecting Cells ####
pluripotentSubset = CellSelector(plot = x) #manually downloaded gate and highlight image from plot > export
meioticSubset_wt = CellSelector(plot = x) #manually downloaded gate and highlight image from plot > export
meioticSubset_ko = CellSelector(plot = x) #manually downloaded gate and highlight image from plot > export

# RDS files from this code are provided 
saveRDS(pluripotentSubset, "pluripotentSubset.rds")
saveRDS(meioticSubset_wt, "meioticSubset_wt.rds")
saveRDS(meioticSubset_ko, "meioticSubset_ko.rds")

# Marker Comparison Tests ####
#wt vs ko meiotic subsets
#pluri vs mei in wt
#mei wt vs everything
#mei ko vs everything
#pluri wt everything

markers_rn2col = function(scobj, id1, id2, test) {
  markers = FindMarkers(scobj, ident.1 = id1, ident.2 = id2, test.use = test, logfc.threshold = 0.1)
  markers = tibble::rownames_to_column(markers)
  head(markers)
  markers = markers %>% dplyr::arrange(desc(avg_logFC))
  return(markers)
}

wt_vs_ko_mitoticSubsetMarkers = markers_rn2col(scAlign_tex14_bh_pgcs_seurat, id1 = meioticSubset_wt, id2 = meioticSubset_ko, "bimod")
pluri_vs_mei_WTSubsetMarkers = markers_rn2col(scAlign_tex14_bh_pgcs_seurat, id1 = pluripotentSubset, id2 = meioticSubset_wt, "bimod")
mei_WT_vs_allSubsetMarkers = markers_rn2col(scAlign_tex14_bh_pgcs_seurat, id1 = meioticSubset_wt, id2 = NULL, "bimod")
mei_KO_vs_allSubsetMarkers = markers_rn2col(scAlign_tex14_bh_pgcs_seurat, id1 = meioticSubset_ko, id2 = NULL, "bimod")
pluri_WT_vs_allSubsetMarkers = markers_rn2col(scAlign_tex14_bh_pgcs_seurat, id1 = pluripotentSubset, id2 = NULL, "bimod")

write.csv(wt_vs_ko_mitoticSubsetMarkers,
          file = paste0(results.dir,
                        "mei-pluri-subset-subset/",
                        Sys.Date(),"-",
                        deparse(substitute(scAlign_tex14_bh_pgcs_seurat)),"-",
                        deparse(substitute(wt_vs_ko_mitoticSubsetMarkers)),".csv"), 
          quote = F)

write.csv(pluri_vs_mei_WTSubsetMarkers,
          file = paste0(results.dir,
          "mei-pluri-subset-subset/",
          Sys.Date(),"-",
          deparse(substitute(scAlign_tex14_bh_pgcs_seurat)),"-",
          deparse(substitute(pluri_vs_mei_WTSubsetMarkers)),".csv"), 
          quote = F)

write.csv(mei_WT_vs_allSubsetMarkers, 
          file = paste0(results.dir,
                        "mei-pluri-subset-subset/",
                        Sys.Date(),"-",
                        deparse(substitute(scAlign_tex14_bh_pgcs_seurat)),"-",
                        deparse(substitute(mei_WT_vs_allSubsetMarkers)),".csv"),
          quote = F)

write.csv(mei_KO_vs_allSubsetMarkers, 
          file = paste0(results.dir,
                        "mei-pluri-subset-subset/",
                        Sys.Date(),"-",
                        deparse(substitute(scAlign_tex14_bh_pgcs_seurat)),"-",
                        deparse(substitute(mei_KO_vs_allSubsetMarkers)),".csv"),
          quote = F)

write.csv(pluri_WT_vs_allSubsetMarkers, 
          file = paste0(results.dir,
                        "mei-pluri-subset-subset/",
                        Sys.Date(),"-",
                        deparse(substitute(scAlign_tex14_bh_pgcs_seurat)),"-",
                        deparse(substitute(pluri_WT_vs_allSubsetMarkers)),".csv"),
          quote = F)