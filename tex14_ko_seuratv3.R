library(Seurat)
library(dplyr)
library(rlang)
library(ggplot2)
library(cowplot)

tex14_ko_data = Read10X(data.dir = "tex14ko/outs/filtered_gene_bc_matrices/mm10/")
tex14_ko = CreateSeuratObject(counts = tex14_ko_data, min.cells = 3, project = "tex14ko", assay = "RNA")
rm(list=ls(pattern="_data"))

#Inital QC ####
#create metadata columns
tex14_ko[["percent.blood"]] = PercentageFeatureSet(object = tex14_ko, pattern = "^Hb")
tex14_ko[["percent.mito"]] = PercentageFeatureSet(object = tex14_ko, pattern = "^mt-")
tex14_ko[["percent.ribo"]] = PercentageFeatureSet(object = tex14_ko, pattern = "^Rps")

metafeatures = colnames(tex14_ko@meta.data)
VlnPlot(tex14_ko, features = metafeatures[2:6], group.by = 'orig.ident', ncol = 5)

#pre filtered scatter and violin plots
preplot1a = FeatureScatter(object = tex14_ko, feature1 = metafeatures[2], feature2 = metafeatures[3])
preplot1b = FeatureScatter(object = tex14_ko, feature1 = metafeatures[2], feature2 = metafeatures[4])
preplot1c = FeatureScatter(object = tex14_ko, feature1 = metafeatures[2], feature2 = metafeatures[5])
CombinePlots(list(preplot1a, preplot1b, preplot1c), ncol = 3)

# Filter Cells ####
tex14_ko = subset(x = tex14_ko, subset = nCount_RNA < 25000 & nFeature_RNA > 200 & percent.mito < 7.0 & percent.blood < 0.5)

#post plots
postplot1a = FeatureScatter(object = tex14_ko, feature1 = metafeatures[2], feature2 = metafeatures[3])
postplot1b = FeatureScatter(object = tex14_ko, feature1 = metafeatures[2], feature2 = metafeatures[4])
postplot1c = FeatureScatter(object = tex14_ko, feature1 = metafeatures[2], feature2 = metafeatures[5])

VlnPlot(object = tex14_ko, features = metafeatures[2:6], group.by = "orig.ident", ncol = 5)

CombinePlots(list(postplot1a, postplot1b, postplot1c), ncol = 3)

# Normalize Data and Find Variable Genes ####
tex14_ko = NormalizeData(object = tex14_ko, normalization.method = "LogNormalize", scale.factor = 10000)
#Normalized values are stored in  pbmc[["RNA"]]@data.

tex14_ko = FindVariableFeatures(object = tex14_ko, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10_ko = head(VariableFeatures(object = tex14_ko), 10)

# plot variable features with and without labels
plot1d = VariableFeaturePlot(object = tex14_ko)
plot1e = LabelPoints(plot = plot1d, points = top10_ko) + ggtitle("ko")
CombinePlots(list(plot1d, plot1e))

# Scale Data ####
all_genes = rownames(tex14_ko)
tex14_ko = ScaleData(tex14_ko, all_genes)

# PCA
tex14_ko = RunPCA(tex14_ko, features = VariableFeatures(object = tex14_ko))
VizDimLoadings(tex14_ko, dims = 1:9, reduction = "pca")
tex14_ko = JackStraw(object = tex14_ko, num.replicate = 100)
tex14_ko = ScoreJackStraw(tex14_ko, dims = 1:20)
JackStrawPlot(tex14_ko, dims = 1:20)
ElbowPlot(tex14_ko) + ggtitle('All Cells from tex14 ko')
saveRDS(tex14_ko, "tex14_ko.rds")

# Total Cells Clustering ####
dims.used = 1:15
tex14_ko = FindNeighbors(object = tex14_ko, dims = dims.used)
tex14_ko = FindClusters(tex14_ko, resolution = seq(0.4,1,0.2))
#tex14_ko@meta.data$integrated_snn_res.* 
#here is where different res's are stored now
metafeatures = colnames(tex14_ko@meta.data)
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

tex14_ko = RunTSNE(object = tex14_ko, reduction = "pca", dims = dims.used)
tex14_ko = RunUMAP(object = tex14_ko, reduction = "pca", dims = dims.used)
saveRDS(tex14_ko, "tex14_ko.rds")

# Generate TSNE/UMAP for clustering of total cells ####
cols30 = c('#339966','#800000','#808080','#000080','#99CCFF',
           '#9bc39b','#008000','#9999FF','#3366FF','#FF99CC',
           '#00FF00','#FF9900','#993366','#FF6600','#CC99FF',
           '#0000FF','#808000','#660066','#008080','#FFCC99',
           '#FF00FF','#666699','#FF8080','#c3c374','#3366FF',
           '#00FFFF','#008080','#0066CC','#000000','#FF0000')
colswtko = c('black','#007FFF')

DimPlot(tex14_ko, 
        reduction = "tsne", 
        group.by = metafeatures[10], 
        label = F,
        cols = cols30) + ggtitle('tex14 ko germ cells and somatic niche') + NoAxes()

DimPlot(tex14_ko, 
        reduction = "umap", 
        group.by = metafeatures[10], 
        label = F,
        cols = cols30) + ggtitle('tex14 ko germ cells and somatic niche') + NoAxes()

# Identify Germ Cells and Subset PGCs ####
gcmarkers = c('Dazl','Stra8','Sycp1','Sycp3', 'Dppa3','Pou5f1')
FeaturePlot(tex14_ko, gcmarkers, reduction = 'tsne')

tex14_ko@meta.data$RNA_snn_res.1 #21
Idents(tex14_ko) = "RNA_snn_res.1"
FeaturePlot(tex14_ko, "Dazl", reduction = "tsne")
#sanity check graph
germcellids = WhichCells(tex14_ko, idents = c(10,15,7,8,3,13,6,1,4,2))
DimPlot(tex14_ko, 
        reduction = "tsne", 
        label = T,
        cols = cols30,
        cells = germcellids) + ggtitle('tex14 wt germ cells highlighted') + NoAxes()

pgc_ko = subset(tex14_ko,idents =  c(10,15,7,8,3,13,6,1,4,2))

# PGC Only Clustering ####
pgc_ko = RunPCA(object = pgc_ko, npcs = 30)
pgc_ko = JackStraw(object = pgc_ko, num.replicate = 100)
pgc_ko = ScoreJackStraw(pgc_ko, dims = 1:15)
JackStrawPlot(pgc_ko, dims = 1:15)
ElbowPlot(pgc_ko) + ggtitle('PGCs from tex14 wt')
saveRDS(pgc_ko, "pgc_ko.rds")

#the elbow plot makes me think 1:10 should be used
dims.used = 1:10
pgc_ko = FindNeighbors(object = pgc_ko, dims = dims.used)
pgc_ko = FindClusters(pgc_ko, resolution = seq(0.3,1,0.1))

metafeatures = colnames(pgc_ko@meta.data)
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
Idents(pgc_ko) = "RNA_snn_res.0.6"
DimPlot(pgc_ko, 
        reduction = "tsne",
        cols = cols30,
        label = T) + ggtitle('tex14 ko germ cells') + NoAxes()

pgc_ko = RenameIdents(pgc_ko,
                       `0` = "Stra8+ Cells",
                       `1` = "Sycp3+ / Stra8+ Cells",
                       `2` = "Stra8+ Cells",
                       `3` = "Sycp3+ Cells",
                       `4` = "Transition",
                       `5` = "Stra8+ Cells",
                       `6` = "Oct4+ Pluripotent",
                       `7` = "Outlier Germ Cell",
                       `8` = "Transition",
                       `9` = "Outlier Germ Cell")
TSNEPlot(pgc_ko, label = T)

# Differentially Expressed Genes Between Pluriopotent and Meiotically Arrested Cells ####
# plu = 6
# mei = 3

# Identify differential expressed genes across conditions
Idents(pgc_ko) = "RNA_snn_res.0.6"
koplumei = subset(pgc_ko, idents = c("3","6"))
koplumei = RenameIdents(koplumei,
                         `3` = "mei",
                         `6` = "plu")

Idents(koplumei)
koplumei = FindVariableFeatures(koplumei)
avgko = log1p(AverageExpression(koplumei)$RNA)
avgko$gene = rownames(avgko)
avgko$diff = avgko$mei - avgko$plu
avgko = tibble::rownames_to_column(avgko)
avgko = avgko %>% dplyr::arrange(avgko$diff)
avgko = tibble::column_to_rownames(avgko)

n = 30
konames = c(head(avgko$gene, n), tail(avgko$gene, n))
kograph = ggplot(avgko, aes(mei, plu)) + geom_point() + ggtitle("ko")
kograph = LabelPoints(kograph, points = konames, repel = T, xnudge = 0, ynudge = 0)
kograph

plumarkers_ko = FindMarkers(koplumei, ident.1 = "plu", ident.2 = "mei", test.use = 'wilcox')
plumarkers_ko$diff = plumarkers_ko$pct.1 - plumarkers_ko$pct.2
plumarkers_ko = tibble::rownames_to_column(plumarkers_ko)
head(plumarkers_ko)
plumarkers_ko = plumarkers_ko %>% dplyr::arrange(desc(avg_logFC))
plumarkers_ko$rowname[1:10]
str(plumarkers_ko)
plumarkers_ko = tibble::column_to_rownames(plumarkers_ko)
write.csv(plumarkers_ko,quote = F,file = "tex14-ko_plu_vs_mei_genes.csv")