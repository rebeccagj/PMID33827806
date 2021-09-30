library(Seurat)
library(dplyr)
library(rlang)
library(ggplot2)

tex14_wt_data = Read10X(data.dir = 'tex14wt/outs/filtered_gene_bc_matrices/mm10/')
tex14_wt = CreateSeuratObject(counts = tex14_wt_data, min.cells = 3, project = 'tex14wt', assay = 'RNA')
tex14_ko_data = Read10X(data.dir = 'tex14ko/outs/filtered_gene_bc_matrices/mm10/')
tex14_ko = CreateSeuratObject(counts = tex14_ko_data, min.cells = 3, project = 'tex14ko', assay = 'RNA')

rm(list=ls(pattern='_data'))

# Inital QC ####
#create metadata columns
tex14_wt[['percent.blood']] = PercentageFeatureSet(object = tex14_wt, pattern = '^Hb')
tex14_wt[['percent.mito']] = PercentageFeatureSet(object = tex14_wt, pattern = '^mt-')
tex14_wt[['percent.ribo']] = PercentageFeatureSet(object = tex14_wt, pattern = '^Rps')

tex14_ko[['percent.blood']] = PercentageFeatureSet(object = tex14_ko, pattern = '^Hb')
tex14_ko[['percent.mito']] = PercentageFeatureSet(object = tex14_ko, pattern = '^mt-')
tex14_ko[['percent.ribo']] = PercentageFeatureSet(object = tex14_ko, pattern = '^Rps')

#pre filtered scatter and violin plots
metafeatures = colnames(tex14_wt@meta.data)
VlnPlot(tex14_wt, features = metafeatures[2:6], group.by = "orig.ident", ncol = 5)
VlnPlot(tex14_ko, features = metafeatures[2:6], group.by = "orig.ident", ncol = 5)

preplot1a = FeatureScatter(object = tex14_wt, feature1 = metafeatures[2], feature2 = metafeatures[3])
preplot1b = FeatureScatter(object = tex14_wt, feature1 = metafeatures[2], feature2 = metafeatures[4])
preplot1c = FeatureScatter(object = tex14_wt, feature1 = metafeatures[2], feature2 = metafeatures[5])

preplot2a = FeatureScatter(object = tex14_ko, feature1 = metafeatures[2], feature2 = metafeatures[3])
preplot2b = FeatureScatter(object = tex14_ko, feature1 = metafeatures[2], feature2 = metafeatures[4])
preplot2c = FeatureScatter(object = tex14_ko, feature1 = metafeatures[2], feature2 = metafeatures[5])

CombinePlots(list(preplot1a, preplot1b, preplot1c,
                  preplot2a, preplot2b, preplot2c),
             ncol = 3)

# Filter Cells ####
tex14_wt = subset(x = tex14_wt, subset = nCount_RNA < 25000 & nFeature_RNA > 200 & percent.mito < 7.0 & percent.blood < 0.5)
tex14_ko = subset(x = tex14_ko, subset = nCount_RNA < 25000 & nFeature_RNA > 200 & percent.mito < 7.0 & percent.blood < 0.5)

#post plots
VlnPlot(tex14_wt, features = metafeatures[2:6], group.by = "orig.ident", ncol = 5)
VlnPlot(tex14_ko, features = metafeatures[2:6], group.by = "orig.ident", ncol = 5)

postplot1a = FeatureScatter(object = tex14_wt, feature1 = metafeatures[2], feature2 = metafeatures[3])
postplot1b = FeatureScatter(object = tex14_wt, feature1 = metafeatures[2], feature2 = metafeatures[4])
postplot1c = FeatureScatter(object = tex14_wt, feature1 = metafeatures[2], feature2 = metafeatures[5])

postplot2a = FeatureScatter(object = tex14_ko, feature1 = metafeatures[2], feature2 = metafeatures[3])
postplot2b = FeatureScatter(object = tex14_ko, feature1 = metafeatures[2], feature2 = metafeatures[4])
postplot2c = FeatureScatter(object = tex14_ko, feature1 = metafeatures[2], feature2 = metafeatures[5])

CombinePlots(list(postplot1a, postplot1b, postplot1c,
                  postplot2a, postplot2b, postplot2c),
             ncol = 3)

# Normalize Data and Find Variable Genes ####
tex14_wt = NormalizeData(object = tex14_wt, normalization.method = 'LogNormalize', scale.factor = 10000)
tex14_ko = NormalizeData(object = tex14_ko, normalization.method = 'LogNormalize', scale.factor = 10000)
#Normalized values are stored in  pbmc[['RNA']]@data.

tex14_wt = FindVariableFeatures(object = tex14_wt, selection.method = 'vst', nfeatures = 2000)
tex14_ko = FindVariableFeatures(object = tex14_ko, selection.method = 'vst', nfeatures = 2000)

#identify the 10 most highly variable genes
top10_wt = head(VariableFeatures(object = tex14_wt), 10)
top10_ko = head(VariableFeatures(object = tex14_ko), 10)

#plot variable features with and without labels
plot1d = VariableFeaturePlot(object = tex14_wt)
plot1e = LabelPoints(plot = plot1d, points = top10_wt) #removal of repel = true fixes 'Viewport has zero dimension(s)' error
CombinePlots(list(plot1d, plot1e))

plot2d = VariableFeaturePlot(object = tex14_ko)
plot2e = LabelPoints(plot = plot2d, points = top10_ko) #removal of repel = true fixes 'Viewport has zero dimension(s)' error
CombinePlots(list(plot2d, plot2e))

plot1e = LabelPoints(plot = plot1d, points = top10_wt) + ggtitle('WT')
plot2e = LabelPoints(plot = plot2d, points = top10_ko) + ggtitle('KO')
CombinePlots(list(plot1e, plot2e))

# Integrate Data ####
tex14.anchors = FindIntegrationAnchors(list(tex14_wt, tex14_ko), dims = 1:30)
tex14.integrated = IntegrateData(tex14.anchors, dims = 1:30)
tex14.integrated = ScaleData(object = tex14.integrated, features = tex14.integrated@assays$RNA@data@Dimnames[[1]])
tex14.integrated = RunPCA(object = tex14.integrated, npcs = 30)

#investigate PCs
tex14.integrated = JackStraw(object = tex14.integrated, num.replicate = 100)
tex14.integrated = ScoreJackStraw(tex14.integrated, dims = 1:20)
JackStrawPlot(tex14.integrated, dims = 1:20)
ElbowPlot(tex14.integrated) + ggtitle('All Cells from tex14 wt and ko')

saveRDS(tex14.integrated, 'tex14.integrated.rds')
#tex14.integrated = readRDS('tex14.integrated.rds')

# Total Cells Clustering ####
dims.used = 1:16
tex14.integrated = FindNeighbors(object = tex14.integrated, dims = dims.used)
tex14.integrated = FindClusters(tex14.integrated, resolution = seq(0.4,1,0.2))
#tex14.integrated@meta.data$integrated_snn_res.* 
#here is where different res's are stored now

metafeatures = colnames(tex14.integrated@meta.data)
metafeatures
# [1] "orig.ident"            
# [2] "nCount_RNA"            
# [3] "nFeature_RNA"          
# [4] "percent.blood"         
# [5] "percent.mito"          
# [6] "percent.ribo"          
# [7] "integrated_snn_res.0.4"
# [8] "integrated_snn_res.0.6"
# [9] "integrated_snn_res.0.8"
# [10] "integrated_snn_res.1"  
# [11] "seurat_clusters"

tex14.integrated = RunTSNE(object = tex14.integrated, reduction = 'pca', dims = dims.used)
tex14.integrated = RunUMAP(object = tex14.integrated, reduction = 'pca', dims = dims.used)
saveRDS(tex14.integrated, 'tex14.integrated.rds')
tex14.integrated = readRDS('tex14.integrated.rds')

# Generate TSNE/UMAP for 'orig ident' and clustering of total cells ####
cols30 = c('#339966','#800000','#808080','#000080','#99CCFF',
           '#9bc39b','#008000','#9999FF','#3366FF','#FF99CC',
           '#00FF00','#FF9900','#993366','#FF6600','#CC99FF',
           '#0000FF','#808000','#660066','#008080','#FFCC99',
           '#FF00FF','#666699','#FF8080','#c3c374','#3366FF',
           '#00FFFF','#008080','#0066CC','#000000','#FF0000')
colswtko = c('black','#007FFF')

metafeatures[10]
# [1] "integrated_snn_res.1"
metafeatures[1]
# [1] "orig.ident"

DimPlot(tex14.integrated, 
        reduction = 'tsne', 
        group.by = metafeatures[10], 
        label = T,
        cols = cols30) + ggtitle('tex14 wt and ko germ cells and somatic niche') + NoAxes()

DimPlot(tex14.integrated,
        reduction = 'tsne', 
        group.by = metafeatures[1], 
        order = c('tex14wt',
                  'tex14ko'),
        cols = colswtko) + ggtitle('tex14 wt and ko germ cells and somatic niche') + NoAxes()

DimPlot(tex14.integrated,
        reduction = 'tsne', 
        group.by = metafeatures[1], 
        cols = colswtko) + ggtitle('tex14 wt and ko germ cells and somatic niche') + NoAxes()

# Identify Germ Cells and Subset PGCs ####
gcmarkers = c('Dazl','Stra8','Sycp1','Sycp3', 'Dppa3','Pou5f1')
FeaturePlot(tex14.integrated, gcmarkers, reduction = 'tsne', pt.size = 0.5)

tex14.integrated@meta.data$integrated_snn_res.1 #26 clusters
Idents(tex14.integrated) = 'integrated_snn_res.1'

DimPlot(tex14.integrated,
        reduction = 'tsne', 
        label = T,
        cols = cols30) + ggtitle('tex14 wt and ko germ cells and somatic niche') + NoAxes()

#sanity check graph
germcellids = WhichCells(tex14.integrated, idents = c(2,0,7,3,1,4,8,16,9,20,12,18))
DimPlot(tex14.integrated, 
        reduction = 'tsne', 
        label = T,
        cols = cols30,
        cells = germcellids) + ggtitle('tex14 wt and ko germ cells highlighted') + NoAxes()

pgc.integrated = subset(tex14.integrated,idents = c(2,0,7,3,1,4,8,16,9,20,12,18))

# PGC Only Clustering ####
pgc.integrated = RunPCA(object = pgc.integrated, npcs = 30)
pgc.integrated = JackStraw(object = pgc.integrated, num.replicate = 100)
pgc.integrated = ScoreJackStraw(pgc.integrated, dims = 1:15)
ElbowPlot(pgc.integrated) + ggtitle('PGCs from tex14 wt and ko')
saveRDS(pgc.integrated, 'pgc.integrated.rds')

#both the elbow plot and the jackstaw make me think 1:11 should be used
dims.used = 1:11
pgc.integrated = FindNeighbors(object = pgc.integrated, dims = dims.used)
pgc.integrated = FindClusters(pgc.integrated, resolution = seq(0.3,1,0.1))

metafeatures = colnames(pgc.integrated@meta.data)
metafeatures
# [1] "orig.ident"            
# [2] "nCount_RNA"            
# [3] "nFeature_RNA"          
# [4] "percent.blood"         
# [5] "percent.mito"          
# [6] "percent.ribo"          
# [7] "integrated_snn_res.0.4"
# [8] "integrated_snn_res.0.6"
# [9] "integrated_snn_res.0.8"
# [10] "integrated_snn_res.1"  
# [11] "seurat_clusters"       
# [12] "integrated_snn_res.0.3"
# [13] "integrated_snn_res.0.5"
# [14] "integrated_snn_res.0.7"
# [15] "integrated_snn_res.0.9"

pgc.integrated = RunTSNE(object = pgc.integrated, reduction = 'pca', dims = dims.used)
pgc.integrated = RunUMAP(object = pgc.integrated, reduction = 'pca', dims = dims.used)
saveRDS(pgc.integrated, 'pgc.integrated.rds')

#After careful comparison by eye, integrated_snn_res.0.5 was selected as the best clustering for PGC only
Idents(pgc.integrated) = 'integrated_snn_res.0.5'

DimPlot(pgc.integrated, 
        reduction = 'tsne', 
        cols = cols30,
        label = T) + ggtitle('tex14 wt and ko germ cells') + NoAxes()

# Split Dot Plot of Genes of Interest ####
Idents(pgc.integrated) <- 'integrated_snn_res.0.5'
Idents(pgc.integrated)
pgc.integrated <- RenameIdents(pgc.integrated,
                               `0` = 'Stra8+ / Sycp3+ Cells',
                               `1` = 'Sycp3+ Cells',
                               `2` = 'Stra8+ / Sycp3+ Cells',
                               `3` = 'Stra8+ Cells',
                               `4` = 'Pluripotent',
                               `5` = 'Cell Cycling / Pluripotent',
                               `6` = 'Pluripotent',
                               `7` = 'Outlier Germ Cell',
                               `8` = 'Cell Cycling / Pluripotent',
                               `9` = 'Outlier Germ Cell')

DimPlot(pgc.integrated, 
        reduction = 'tsne', 
        cols = cols30,
        label = T) + ggtitle('tex14 wt and ko germ cells') + NoAxes()

genelist = c('Rhox5','Actb','Dppa3','Phlda2','Rec8','Stra8','Sycp3')
DotPlot(pgc.integrated, 
        features = rev(genelist), 
        cols = c('blue', 'black'), 
        dot.scale = 8, 
        split.by = 'orig.ident') + RotatedAxis()

genelist = c('Dppa3','Sox2','Pou5f1','Nanog')
DotPlot(pgc.integrated, 
        features = rev(genelist), 
        cols = c('blue', 'black'), 
        dot.scale = 8, 
        split.by = 'orig.ident') + RotatedAxis()

# Investigation of RA Receptors in PGCs ####
FeaturePlot(tex14.integrated,
            c('Dppa3','Sycp3','Rara', 'Rarb', 'Rarg', 'Rxra', 'Rxrb', 'Rxrg'),
            reduction = 'tsne', 
            split.by = 'orig.ident')