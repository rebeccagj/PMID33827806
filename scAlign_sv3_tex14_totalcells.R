library(scAlign)
library(Seurat)
library(ggplot2)
library(DittoSeq)
library(SingleCellExperiment)
library(gridExtra)
library(cowplot)

# This file works with the same libraries, but a second more deeply sequenced run

# create seurat objects ####
read_create_meta = function(dirvar,projvar,cellname){
  mat = Read10X(dirvar)
  seuratobj = CreateSeuratObject(counts = mat, project = projvar)
  seuratobj = RenameCells(seuratobj, add.cell.id = cellname)
  seuratobj[['percent.blood']] = PercentageFeatureSet(object = seuratobj, pattern = '^Hb')
  seuratobj[['percent.mito']] = PercentageFeatureSet(object = seuratobj, pattern = '^mt-')
  seuratobj[['percent.ribo']] = PercentageFeatureSet(object = seuratobj, pattern = c('^Rps'))
  seuratobj
}

tex14wt = read_create_meta('tex14wt_biohub/outs/filtered_gene_bc_matrices/mm10/','tex14wt','1')
tex14ko = read_create_meta('tex14ko_biohub/outs/filtered_gene_bc_matrices/mm10/','tex14ko','2')

tex14wt
# An object of class Seurat 
# 31053 features across 14311 samples within 1 assay 
# Active assay: RNA (31053 features)
tex14ko
# An object of class Seurat 
# 31053 features across 18287 samples within 1 assay 
# Active assay: RNA (31053 features)

# Subsetting Values Determined in tex14_wt_seuratv3.R and tex14_ko_seuratv3.R
tex14wt = subset(tex14wt, subset = nCount_RNA < 30000 & nFeature_RNA > 1000 & percent.mito < 5 & percent.blood < .25)
tex14ko = subset(tex14ko, subset = nCount_RNA < 30000 & nFeature_RNA > 1000 & percent.mito < 5 & percent.blood < .25)

tex14wt = NormalizeData(tex14wt)
tex14ko = NormalizeData(tex14ko)

tex14wt = FindVariableFeatures(tex14wt, selection.method = 'vst', nfeatures = 3000)
tex14ko = FindVariableFeatures(tex14ko, selection.method = 'vst', nfeatures = 3000)

tex14wt <- ScaleData(tex14wt)
tex14ko <- ScaleData(tex14ko)

genes.use = Reduce(intersect, list(VariableFeatures(tex14wt),
                                   VariableFeatures(tex14ko),
                                   rownames(tex14wt),
                                   rownames(tex14ko)))

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

tex14_wt_sce = seuratv3_to_scefxn(tex14wt,genes.use)
tex14_ko_sce = seuratv3_to_scefxn(tex14ko,genes.use)

scAlign_tex14 = scAlignCreateObject(sce.objects = list('wt'=tex14_wt_sce, 'ko'=tex14_ko_sce),
                                 labels = list('wt', 'ko'),
                                 data.use='scale.data',
                                 pca.reduce = TRUE,
                                 pcs.compute = 50,
                                 cca.reduce = TRUE,
                                 ccs.compute = 15,
                                 project.name = 'scAlign_tex14')

results.dir = paste0('scAlign-',Sys.Date(),'/')
dir.create(results.dir)

## Run scAlign with high_var_genes as input to the encoder (alignment) and logcounts with the decoder (projections).
scAlign_tex14 = scAlign(scAlign_tex14,
                     options=scAlignOptions(steps=5000, log.every=5000, norm=TRUE, batch.norm.layer=TRUE, early.stop=FALSE, architecture='small'),
                     encoder.data='scale.data',
                     decoder.data='logcounts',
                     supervised='none',
                     run.encoder=TRUE,
                     run.decoder=TRUE,
                     log.dir=file.path(results.dir, 'models','gene_input'),
                     device='CPU')

## Additional run of scAlign with PCA, the early.stopping heuristic terminates the training procedure too early with PCs as input so it is disabled.
scAlign_tex14 = scAlign(scAlign_tex14,
                     options=scAlignOptions(steps=15000, log.every=1000, norm=TRUE, batch.norm.layer=TRUE, early.stop=FALSE, architecture='small'),
                     encoder.data='PCA',
                     supervised='none',
                     run.encoder=TRUE,
                     run.decoder=FALSE,
                     log.dir=file.path(results.dir, 'models','pca_input'),
                     device='CPU')

## Additional run of scAlign with CCA
scAlign_tex14 = scAlign(scAlign_tex14,
                     options=scAlignOptions(steps=10000, log.every=1000, norm=TRUE, batch.norm.layer=TRUE, early.stop=FALSE, architecture='small'),
                     encoder.data='CCA',
                     supervised='none',
                     run.encoder=TRUE,
                     run.decoder=FALSE,
                     log.dir=file.path(results.dir, 'models','cca_input'),
                     device='CPU')

saveRDS(scAlign_tex14, 'scAlign_tex14.rds')

## Plot aligned data in tSNE space, when the data was processed in three different ways: 1) either using the original gene inputs, 2) after PCA dimensionality reduction for preprocessing, or 3) after CCA dimensionality reduction for preprocessing. Cells here are colored by dataset.
set.seed(5678)
gene_plot = PlotTSNE(scAlign_tex14_bh_pgcs, "ALIGNED-GENE", title="scAlign-Gene", perplexity=30)
pca_plot = PlotTSNE(scAlign_tex14_bh_pgcs, "ALIGNED-PCA", title="scAlign-PCA", perplexity=30)
cca_plot = PlotTSNE(scAlign_tex14_bh_pgcs, "ALIGNED-CCA", title="scAlign-CCA", perplexity=30)
wt2ko = PlotTSNE(scAlign_tex14_bh_pgcs, "wt2ko", title = "scAlign-wt2ko", perplexity=30)
ko2wt = PlotTSNE(scAlign_tex14_bh_pgcs, "ko2wt", title = "scAlign-ko2wt", perplexity=30)
legend = get_legend(PlotTSNE(scAlign_tex14_bh_pgcs, "ALIGNED-GENE", title="scAlign-Gene", legend="right", max_iter=1))
combined_plot = grid.arrange(wt2ko, ko2wt, nrow = 2)

# Take these integrations back to Seurat for further analysis ####
scAlign_tex14_seuratObj = as.Seurat(scAlign_tex14, 
                                    counts = 'counts', 
                                    scale.data = 'scale.data')
scAlign_tex14_seuratObj@reductions
saveRDS(scAlign_tex14_seuratObjObj, 'scAlign_tex14_seuratObj.rds')

# Total Cells Clustering and Dimension Reduction ####
scAlign_tex14_seuratObj = FindNeighbors(scAlign_tex14_seuratObj, reduction = 'ko2wt', dims = 1:2390)
scAlign_tex14_seuratObj = FindClusters(scAlign_tex14_seuratObj, reduction = 'ko2wt')
scAlign_tex14_seuratObj = RunTSNE(scAlign_tex14_seuratObj, reduction = 'ko2wt', dims = 1:2390)

colswtko = c('black','#007FFF')
origidplot = DimPlot(scAlign_tex14_seuratObj,
                     group.by = 'orig.ident', 
                     reduction = 'tsne', 
                     order = c('tex14ko','tex14wt'), 
                     cols = colswtko) + ggtitle('Genotype')
clusterplot = DimPlot(scAlign_tex14_seuratObj,
                      label = TRUE) + ggtitle('Cluster Resolution 0.8')
pgcmarker = FeaturePlot(scAlign_tex14_seuratObj, 
                      features = 'Dazl',
                      pt.size = cexSize)
CombinePlots(list(origidplot,clusterplot,pgcmarker), ncol = 3)
saveRDS(scAlign_tex14_seuratObj, 'scAlign_tex14_seuratObj.rds')