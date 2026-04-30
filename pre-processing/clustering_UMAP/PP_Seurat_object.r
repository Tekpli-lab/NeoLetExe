library(scDblFinder)
library(DoubletFinder)
library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (mixtools)
library (SoupX)
library (DropletUtils)
library (singleCellTK)
library (immunarch)
library (tidyr)

### Folder need to be organized according to PP_CellRanger
metadata <-data.frame(read.table("metadata.txt", header=T, sep="\t", stringsAsFactors=F))
## a file containing some simple metadata to be passed to the Seurat object

for (sample1 in metadata$samplename){
###### scRNA-seq expresssion data ######

a1 <- load10X(dataDir = paste (sample1, "/outs/count", sep=""))
a1 = autoEstCont(a1)
out = adjustCounts(a1)
a1 <- CreateSeuratObject(out,project = metadata$Name[metadata$samplename==sample1], min.cells = 1, min.features = 1)
a1[["percent.mt"]] <- PercentageFeatureSet(a1, pattern = "^MT-")
a1 <- subset(a1, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 100 & nCount_RNA > 500)

a2 <- Read10X(data.dir = paste (sample1, "/outs/count/raw_feature_bc_matrix", sep=""))
a <- emptyDrops (a2)
is.cell <- a$FDR <= 0.001
sum(is.cell, na.rm=TRUE)
cell.counts <- a2[,which(is.cell),drop=FALSE]
dim(cell.counts)
a2 <- CreateSeuratObject(counts=cell.counts, project = metadata$Name[metadata$samplename==sample1], min.cells = 1, min.features = 1)
a2[["percent.mt"]] <- PercentageFeatureSet(a2, pattern = "^MT-")
a2 <- subset(a2, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 15 & nCount_RNA > 500)

a <- subset(a1, cell = colnames(a2))
sce <- as.SingleCellExperiment(a)
sced <- scDblFinder(sce)
a@meta.data$doublets.score <- sced$scDblFinder.score
a@meta.data$doublets.class <- sced$scDblFinder.class

###### meta.data ######
a@meta.data$Condition <- metadata$Condition[metadata$samplename==sample1]
a@meta.data$Time_point <- metadata$Time_point[metadata$samplename==sample1]
a@meta.data$Trial <- metadata$Trial[metadata$samplename==sample1]
a@meta.data$patient_number <- metadata$patient_number[metadata$samplename==sample1]

###### summary data ######
assign(metadata$Name[metadata$samplename==sample1],a)
print (metadata$Name[metadata$samplename==sample1])
print (dim (get (metadata$Name[metadata$samplename==sample1])))
}

Merged <- merge(P01A, y = c(P01B, P01C,
P02A,  P02B,  P02C,
P03A,  P03B,  P03C,   
P04A,
P05A,  P05B,  P05C,  
P06A,  P06B,  P06C, 
P07A,  P07B,  P07C, 
P08A,  P08B,  P08C, 
P09A,  P09B,  P09C, 
P10A,  P10B,  P10C, 
P11A,  P11B,  P11C, 
P12A,  P12B,  P12C, 
P13A,  P13B,  P13C, 
P14A,  P14B,  P14C, 
P15A,  P15B,  P15C, 
P16A,  P16B,  P16C, 
P17A,  P17B,  P17C, 
P18A,  P18B,  P18C, 
P19A,  P19B,  P19C, 
P20A,  P20B,  P20C, 
P21A,  P21B,  P21C, 
P22A,  P22B,  P22C, 
P23A,  P23B,  P23C,
P24A,  P24B,  P24C,  
P25A,  P25B,  P25C),
add.cell.ids = c(
"P01A",  "P01B",  "P01C",
"P02A",  "P02B",  "P02C",
"P03A",  "P03B",  "P03C",   
"P04A",
"P05A",  "P05B",  "P05C",  
"P06A",  "P06B",  "P06C", 
"P07A",  "P07B",  "P07C", 
"P08A",  "P08B",  "P08C", 
"P09A",  "P09B",  "P09C", 
"P10A",  "P10B",  "P10C", 
"P11A",  "P11B",  "P11C", 
"P12A",  "P12B",  "P12C", 
"P13A",  "P13B",  "P13C", 
"P14A",  "P14B",  "P14C", 
"P15A",  "P15B",  "P15C", 
"P16A",  "P16B",  "P16C", 
"P17A",  "P17B",  "P17C", 
"P18A",  "P18B",  "P18C", 
"P19A",  "P19B",  "P19C", 
"P20A",  "P20B",  "P20C", 
"P21A",  "P21B",  "P21C", 
"P22A",  "P22B",  "P22C", 
"P23A",  "P23B",  "P23C",
"P24A",  "P24B",  "P24C",  
"P25A",  "P25B",  "P25C")
, project = "NeoLetExe")

print (data.frame(table(Merged@meta.data$orig.ident)))
print (data.frame(table(Merged@meta.data$orig.ident, Merged@meta.data$Time_point)))
print (data.frame(table(Merged@meta.data$Time_point,Merged@meta.data$orig.ident)))

Merged <- subset(Merged, subset = nFeature_RNA > 450 & nFeature_RNA < 5000 & percent.mt < 15 & nCount_RNA > 600 & doublets.score < 0.95)

print (data.frame(table(Merged@meta.data$orig.ident)))
print (data.frame(table(Merged@meta.data$Time_point)))
dim (Merged)

library (gdata)
keep(Merged, sure=TRUE)
save.image("Merged_NeoLetExe_final_git.RData")


##########################################
##### Scale and regress out variables ####
##########################################
module load  R/4.1.2-foss-2021b
taskset --cpu-list 50-80 R

library(Seurat)
library(plyr)
library(dplyr)

load("Merged_NeoLetExe_final_git.RData")

Merged <- NormalizeData(Merged, normalization.method = "LogNormalize", scale.factor = 10000)
print (dim (Merged))
lapply(Merged@meta.data, class)
all.genes <- rownames(Merged)

Scaled <- ScaleData(Merged, features = all.genes, vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA"), split.by="orig.ident") 


library (gdata)
keep(Scaled, sure=TRUE)
save.image("Scaled_NeoLetExe_final_git.RData")
q()


##########################################
##### Scale and regress out variables ####
#########################################

module load  R/4.1.2-foss-2021b
taskset --cpu-list 150-180 R

library(Seurat)
library(plyr)
library(dplyr)
library(future)

load ("Scaled_NeoLetExe_final_git.RData")

variable <- FindVariableFeatures(Scaled, selection.method = "vst", nfeatures = 8000)
feature <- VariableFeatures(variable)
print (length (feature))
print (dim (Scaled))

Scaled <- RunPCA(Scaled, features = feature, npcs = 300, ndims.print = 1:5, nfeatures.print = 30, weight.by.var=T)


####################################################################################
##### Choose the Number of principal component and the resolution for clusters #####
####################################################################################
Scaled <- FindNeighbors(Scaled, reduction = "pca", dims = 1:140, k.param = 30, compute.SNN = T, nn.method = "rann", nn.eps = 0)
Scaled <- FindClusters(Scaled, resolution = c(0.8,1,1.2)) 
Idents(Scaled)

resol = 'RNA_snn_res.1'
Scaled <- SetIdent(object = Scaled, value = resol)

Scaled <- RunUMAP(Scaled, dims = 1:300, n.neighbors = 100L)
pdf ('UMAP_100L.pdf')
DimPlot(Scaled, label = F, group.by='orig.ident')+NoLegend()
dev.off()

library (gdata)
keep(Scaled, sure=TRUE)
save.image("Scaled_NeoLetExe_final_git.RData")


##############################################################################
#### INTEGRATION PIPELINE ##
##############################################################################

ifnb.list <- SplitObject(Scaled, split.by = "orig.ident")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 8000)
})

features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 8000)

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE,  npcs = 65)
})

anchors <- FindIntegrationAnchors(object.list = ifnb.list, reference = c(1, 4, 7, 10, 13, 16, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59, 62, 65, 68, 71), reduction = "rpca", dims = 1:65, anchor.features = 4500) #125 NR11C has so few cells that we need NPC < 128
anchors <- FindIntegrationAnchors(object.list = ifnb.list, reduction = "rpca", dims = 1:65, anchor.features = 4500) 

Scaled <- IntegrateData(anchorset = anchors, dims = 1:65, k.weight=80) 

Scaled <- ScaleData(Scaled, verbose = FALSE)
Scaled <- RunPCA(Scaled, verbose = FALSE, npcs =200)
#Scaled <- RunUMAP(Scaled, dims = 1:125)

Scaled <- FindNeighbors(Scaled, reduction = "pca", dims = 1:125) 
DefaultAssay(Scaled ) <- "integrated"
Scaled <- FindClusters(Scaled, resolution = c(1,1.2, 1.4, 1.6)) 
Idents(Scaled)

resol = 'integrated_snn_res.1.2'
Scaled <- SetIdent(object = Scaled, value = resol)
Scaled <- RunUMAP(Scaled, dims = 1:125)

library (gdata)
keep(Scaled, sure=TRUE)
save.image("Scaled_NeoLetExe_final_git_integrated.RData")
