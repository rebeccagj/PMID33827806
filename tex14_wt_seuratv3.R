library(Seurat)
library(dplyr)
library(rlang)
library(ggplot2)
library(cowplot)

tex14_wt_data = Read10X(data.dir = "tex14wt/outs/filtered_gene_bc_matrices/mm10/")
tex14_wt = CreateSeuratObject(counts = tex14_wt_data, min.cells = 3, project = "tex14wt", assay = "RNA")
rm(list=ls(pattern="_data"))

# Inital QC ####
#create metadata columns
tex14_wt[["percent.blood"]] = PercentageFeatureSet(object = tex14_wt, pattern = "^Hb")
tex14_wt[["percent.mito"]] = PercentageFeatureSet(object = tex14_wt, pattern = "^mt-")
tex14_wt[["percent.ribo"]] = PercentageFeatureSet(object = tex14_wt, pattern = "^Rps")

metafeatures = colnames(tex14_wt@meta.data)
VlnPlot(tex14_wt, features = metafeatures[2:6], group.by = "orig.ident", ncol = 5)

#pre filtered scatter and violin plots
preplot1a = FeatureScatter(object = tex14_wt, feature1 = metafeatures[2], feature2 = metafeatures[3])
preplot1b = FeatureScatter(object = tex14_wt, feature1 = metafeatures[2], feature2 = metafeatures[4])
preplot1c = FeatureScatter(object = tex14_wt, feature1 = metafeatures[2], feature2 = metafeatures[5])
CombinePlots(list(preplot1a, preplot1b, preplot1c), ncol = 3)

# Filter Cells ####
tex14_wt = subset(x = tex14_wt, subset = nCount_RNA < 25000 & nFeature_RNA > 200 & percent.mito < 7.0 & percent.blood < 0.5)

#post plots
postplot1a = FeatureScatter(object = tex14_wt, feature1 = metafeatures[2], feature2 = metafeatures[3])
postplot1b = FeatureScatter(object = tex14_wt, feature1 = metafeatures[2], feature2 = metafeatures[4])
postplot1c = FeatureScatter(object = tex14_wt, feature1 = metafeatures[2], feature2 = metafeatures[5])

VlnPlot(tex14_wt, features = metafeatures[2:6], group.by = "orig.ident", ncol = 5)

CombinePlots(list(postplot1a, postplot1b, postplot1c), ncol = 3)

# Normalize Data and Find Variable Genes ####
tex14_wt = NormalizeData(object = tex14_wt, normalization.method = "LogNormalize", scale.factor = 10000)
#Normalized values are stored in  pbmc[["RNA"]]@data.

tex14_wt = FindVariableFeatures(object = tex14_wt, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10_wt = head(VariableFeatures(object = tex14_wt), 10)

# plot variable features with and without labels
plot1d = VariableFeaturePlot(object = tex14_wt)
plot1e = LabelPoints(plot = plot1d, points = top10_wt) + ggtitle("WT")
CombinePlots(list(plot1d, plot1e))

# Scale Data ####
all_genes = rownames(tex14_wt)
tex14_wt = ScaleData(tex14_wt, all_genes)

# PCA
tex14_wt = RunPCA(tex14_wt, features = VariableFeatures(object = tex14_wt))
VizDimLoadings(tex14_wt, dims = 1:9, reduction = "pca")
tex14_wt = JackStraw(object = tex14_wt, num.replicate = 100)
tex14_wt = ScoreJackStraw(tex14_wt, dims = 1:20)
JackStrawPlot(tex14_wt, dims = 1:20)
ElbowPlot(tex14_wt) + ggtitle('All Cells from tex14 wt')
saveRDS(tex14_wt, "tex14_wt.rds")

# Total Cells Clustering ####
dims.used = 1:15
tex14_wt = FindNeighbors(object = tex14_wt, dims = dims.used)
tex14_wt = FindClusters(tex14_wt, resolution = seq(0.4,1,0.2))
#tex14_wt@meta.data$integrated_snn_res.* 
#here is where different res's are stored now
metafeatures = colnames(tex14_wt@meta.data)
metafeatures
# [1] "orig.ident"     
# [2] "nCount_RNA"     
# [3] "nFeature_RNA"   
# [4] "percent.blood"  
# [5] "percent.mito"   
# [6] "percent.ribo"   
# [7] "RNA_snn_res.0.4"
# [8] "RNA_snn_res.0.6"
# [9] "RNA_snn_res.0.8"
# [10] "RNA_snn_res.1"  
# [11] "seurat_clusters"

tex14_wt = RunTSNE(object = tex14_wt, reduction = "pca", dims = dims.used)
tex14_wt = RunUMAP(object = tex14_wt, reduction = "pca", dims = dims.used)
saveRDS(tex14_wt, "tex14_wt.rds")

# Generate TSNE/UMAP for clustering of total cells ####
cols30 = c('#339966','#800000','#808080','#000080','#99CCFF',
           '#9bc39b','#008000','#9999FF','#3366FF','#FF99CC',
           '#00FF00','#FF9900','#993366','#FF6600','#CC99FF',
           '#0000FF','#808000','#660066','#008080','#FFCC99',
           '#FF00FF','#666699','#FF8080','#c3c374','#3366FF',
           '#00FFFF','#008080','#0066CC','#000000','#FF0000')
colswtko = c('black','#007FFF')

DimPlot(tex14_wt, 
        reduction = "tsne", 
        group.by = metafeatures[10], 
        label = F,
        cols = cols30) + ggtitle('tex14 wt germ cells and somatic niche') + NoAxes()

DimPlot(tex14_wt, 
        reduction = "umap", 
        group.by = metafeatures[10], 
        label = F,
        cols = cols30) + ggtitle('tex14 wt germ cells and somatic niche') + NoAxes()

# Identify Germ Cells and Subset PGCs ####
gcmarkers = c('Dazl','Stra8','Sycp1','Sycp3', 'Dppa3','Pou5f1')
FeaturePlot(tex14_wt, gcmarkers, reduction = 'tsne')

tex14_wt@meta.data$RNA_snn_res.1 #21
Idents(tex14_wt) = "RNA_snn_res.1"
FeaturePlot(tex14_wt, "Dazl", reduction = "tsne")
#sanity check graph
germcellids = WhichCells(tex14_wt, idents = c(12,18,11,8,9,7,1,2,0,3))
DimPlot(tex14_wt, 
        reduction = "tsne", 
        label = T,
        cols = cols30,
        cells = germcellids) + ggtitle('tex14 wt germ cells highlighted') + NoAxes()

pgc_wt = subset(tex14_wt,idents = c(12,18,11,8,9,7,1,2,0,3))

# PGC Only Clustering ####
pgc_wt = RunPCA(object = pgc_wt, npcs = 30)
pgc_wt = JackStraw(object = pgc_wt, num.replicate = 100)
pgc_wt = ScoreJackStraw(pgc_wt, dims = 1:15)
JackStrawPlot(pgc_wt, dims = 1:15)
ElbowPlot(pgc_wt) + ggtitle('PGCs from tex14 wt')
saveRDS(pgc_wt, "pgc_wt.rds")

#the elbow plot makes me think 1:10 should be used
dims.used = 1:10
pgc_wt = FindNeighbors(object = pgc_wt, dims = dims.used)
pgc_wt = FindClusters(pgc_wt, resolution = seq(0.3,1,0.1))

metafeatures = colnames(pgc_wt@meta.data)
metafeatures
# [1] "orig.ident"     
# [2] "nCount_RNA"     
# [3] "nFeature_RNA"   
# [4] "percent.blood"  
# [5] "percent.mito"   
# [6] "percent.ribo"   
# [7] "RNA_snn_res.0.4"
# [8] "RNA_snn_res.0.6"
# [9] "RNA_snn_res.0.8"
# [10] "RNA_snn_res.1"  
# [11] "seurat_clusters"
# [12] "RNA_snn_res.0.3"
# [13] "RNA_snn_res.0.5"
# [14] "RNA_snn_res.0.7"
# [15] "RNA_snn_res.0.9"

pgc_wt = RunTSNE(object = pgc_wt, reduction = "pca", dims = dims.used)
pgc_wt = RunUMAP(object = pgc_wt, reduction = "pca", dims = dims.used)
pgc_wt = FindVariableFeatures(pgc_wt, selection.method = "vst", nfeatures = 5000)
saveRDS(pgc_wt, "pgc_wt.rds")

#After careful comparison by eye, RNA_snn_res.0.6 was selected as the best clustering for PGC only
Idents(pgc_wt) = "RNA_snn_res.0.6"
DimPlot(pgc_wt, 
        reduction = "tsne", 
        cols = cols30,
        label = T) + ggtitle('tex14 wt germ cells') + NoAxes()

pgc_wt = RenameIdents(pgc_wt,
                       `0` = "Oct4+ Pluripotent",
                       `1` = "Stra8+ Cells",
                       `2` = "Stra8+ Cells",
                       `3` = "Stra8+ Cells",
                       `4` = "Sycp3+ Cells",
                       `5` = "Sycp3+ Cells",
                       `6` = "Outlier Germ Cell",
                       `7` = "Transition",
                       `8` = "Oct4+ Pluripotent",
                       `9` = "Outlier Germ Cell")
TSNEPlot(pgc_wt, label = T)

# Differentially Expressed Genes Between Pluriopotent and Meiotically Arrested Cells ####
# plu = 0,8
# mei = 4,5

# Identify differential expressed genes across conditions ####
Idents(pgc_wt) = "RNA_snn_res.0.6"
wtplumei = subset(pgc_wt, idents = c("0","8","4","5"))
wtplumei = RenameIdents(wtplumei,
                         `0` = "mei",
                         `4` = "plu",
                         `5` = "plu",
                         `8` = "mei")

Idents(wtplumei)
wtplumei = FindVariableFeatures(wtplumei)
avgwt = log1p(AverageExpression(wtplumei)$RNA)
avgwt$gene = rownames(avgwt)
avgwt$diff = avgwt$mei - avgwt$plu
avgwt = tibble::rownames_to_column(avgwt)
avgwt = avgwt %>% dplyr::arrange(avgwt$diff)
avgwt = tibble::column_to_rownames(avgwt)

n = 30
wtnames = c(head(avgwt$gene, n), tail(avgwt$gene, n))
wtgraph = ggplot(avgwt, aes(mei, plu)) + geom_point() + ggtitle("wt")
wtgraph = LabelPoints(wtgraph, points = wtnames, repel = T, xnudge = 0, ynudge = 0)
wtgraph

plumarkers_wt = FindMarkers(wtplumei, ident.1 = "plu", ident.2 = "mei", test.use = 'wilcox')
plumarkers_wt$diff = plumarkers_wt$pct.1 - plumarkers_wt$pct.2
plumarkers_wt = tibble::rownames_to_column(plumarkers_wt)
head(plumarkers_wt)
plumarkers_wt = plumarkers_wt %>% dplyr::arrange(desc(avg_logFC))
plumarkers_wt$rowname[1:10]
str(plumarkers_wt)
plumarkers_wt = tibble::column_to_rownames(plumarkers_wt)
write.csv(plumarkers_wt,quote = F,file = "tex14-wt_plu_vs_mei_genes.csv")