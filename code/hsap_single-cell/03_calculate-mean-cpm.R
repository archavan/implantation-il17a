# Before running this script, 
# (1) 10x single cell raw data should be downloaded from:
# Suryawanshi et al Science Advances 2018, RNA-seq data: BioProject PRJNA492902
# decidua data: SRA id SRR7895954
# villi data:   SRA id SRR7895962
# (2) Map raw data using `cellranger` using code in `code/hsap_single-cell`

### packages ==================================================================
library(tidyverse)
library(Seurat)

### data ======================================================================
decidua.data <- Read10X(data.dir = "data/hsap_single-cell/02_cellranger-count/decidua/outs/raw_feature_bc_matrix/")
decidua <- CreateSeuratObject(counts = decidua.data, project = "implantation", min.cells = 0, min.features = 200)

villi.data <- Read10X(data.dir = "data/hsap_single-cell/02_cellranger-count/villi/outs/raw_feature_bc_matrix/")
villi <- CreateSeuratObject(counts = villi.data, project = 'implantation', min.cells = 0, min.features = 200)

### qc ========================================================================
# percent mt
decidua[["percent.mt"]] <- PercentageFeatureSet(decidua, pattern = "^MT-")
villi[["percent.mt"]] <- PercentageFeatureSet(villi, pattern = "^MT-")

# Visualize QC metrics 
VlnPlot(decidua, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(villi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(decidua, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(decidua, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

plot3 <- FeatureScatter(villi, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(villi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot3, plot4))

# filter
decidua <- subset(decidua, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
villi <- subset(villi, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

### normalization =============================================================
## normalization to cpm
decidua <- NormalizeData(decidua, normalization.method = "RC", scale.factor = 1e6) 
villi <- NormalizeData(villi, normalization.method = "RC", scale.factor = 1e6) 

### calculate means ===========================================================
# calculate mean cpm for decidua and villi cells to create a "bulk" representation of the data. 

# decidua
decidua.cpm <- as.matrix(decidua[["RNA"]]@data)
decidua.cpm.mean <- rowMeans(decidua.cpm, dims = 1) %>% 
  as.data.frame()
decidua.cpm.mean$gene.name <- row.names(decidua.cpm.mean)
colnames(decidua.cpm.mean) <- c("decidua", "hsap_gene_name")
decidua.cpm.mean <- decidua.cpm.mean[, c(2,1)]

write.csv(decidua.cpm.mean, "data/hsap_single-cell/decidua_cpm_mean.csv", row.names = FALSE)

# villi
villi.cpm <- as.matrix(villi[["RNA"]]@data)
villi.cpm.mean <- rowMeans(villi.cpm, dims = 1) %>% 
  as.data.frame()
villi.cpm.mean$gene.name <- row.names(villi.cpm.mean)
colnames(villi.cpm.mean) <- c("villi", "hsap_gene_name")
villi.cpm.mean <- villi.cpm.mean[, c(2,1)]

write.csv(villi.cpm.mean, "data/hsap_single-cell/villi_cpm_mean.csv", row.names = FALSE)

### end =======================================================================
