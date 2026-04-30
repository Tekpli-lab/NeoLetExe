######################################################################################################################################################
## Take cells from a cell type and cluster it to discover cell states  ###############################################################################
## creates the outputs to create phylogentic trees and clones ########################################################################################
######################################################################################################################################################

library(Seurat)
library(plyr)
library(dplyr)
library(future)

load("Scaled_NeoLetExe_final_git.RData")
#### Load Seurat object from PP_Seurat_object

Subset_Cells <-Scaled[, Scaled@meta.data$main.ids=="T_cells"] 
rm (Scaled)
NL <- SplitObject(Subset_Cells, split.by = "orig.ident")

Merged <- merge(NL$P01A, y = c(NL$P01B, NL$P01C,
NL$P02A,  NL$P02B,  NL$P02C,
NL$P03A,  NL$P03B,  NL$P03C,   
NL$P04A,
NL$P05A,  NL$P05B,  NL$P05C,  
NL$P06A,  NL$P06B,  NL$P06C, 
NL$P07A,  NL$P07B,  NL$P07C, 
NL$P08A,  NL$P08B,  NL$P08C, 
NL$P09A,  NL$P09B,  NL$P09C, 
NL$P10A,  NL$P10B,  NL$P10C, 
NL$P11A,  NL$P11B,  NL$P11C, 
NL$P12A,  NL$P12B,  NL$P12C, 
NL$P13A,  NL$P13B,  NL$P13C, 
NL$P14A,  NL$P14B,  NL$P14C, 
NL$P15A,  NL$P15B,  NL$P15C, 
NL$P16A,  NL$P16B,  NL$P16C, 
NL$P17A,  NL$P17B,  NL$P17C, 
NL$P18A,  NL$P18B,  NL$P18C, 
NL$P19A,  NL$P19B,  NL$P19C, 
NL$P20A,  NL$P20B,  NL$P20C, 
NL$P21A,  NL$P21B,  NL$P21C, 
NL$P22A,  NL$P22B,  NL$P22C, 
NL$P23A,  NL$P23B,  NL$P23C,
NL$P24A,  NL$P24B,  NL$P24C,  
NL$P25A,  NL$P25B,  NL$P25C))

print (table(Merged@meta.data$orig.ident))
dim (Merged)
Merged <- NormalizeData(Merged, normalization.method = "LogNormalize", scale.factor = 10000)

hypoxia.genes <- c("VEGFA", "SLC2A1", "PGAM1", "ENO1", "LDHA", "TPI1", "P4HA1", "MRPS17", "CDKN3", "ADM", "NDRG1", "TUBB6","ALDOA", "MIF", "ACOT7")
stress.genes <- c("FOS", "HSPA1A", "JUN", "FOSB", "JUNB", "EGR1", "HSPA1B", "UBC", "ZFP36", "HSPB1", "HSP90AA1", "MT2A", "DNAJB1", "BTG2", "NR4A1", "CEBPD", "HSPA8", "MT1", "IER2", "DNAJA1", "SOCS3", "ATF3", "JUND", "CEBPB", "ID3", "PPP1R15A", "HSPE1", "CXCL1", "DUSP1", "HSP90AB1", "NFKBIA", "HSPH1")
interferon.genes <- c("GBP1", "TYMP", "IFI44L", "IFI35", "CXCL11", "ISG15", "RABGAP1L", "GCH1", "MELK", "IFIT5", "OAS1", "CXCL10", "APOBEC3G", "UBE2L6", "GBP2", "IFI27", "STX11", "IFI6", "IDO1", "LMO2", "STAT1", "RBCK1", "SP110", "IFIT2", "CASP1", "EIF2AK2", "CYTH1", "IL15", "C1orf38", "TRIM14", "MX2", "IFI30", "BAZ1A", "IRF7", "SP100", "WARS",  "IFI44", "PLSCR1", "PSMB9", "OASL", "RGL1", "NMI", "SAMHD1", "RSAD2", "CCL8", "TNFSF10", "TDRD7", "ISG20", "IFIT3", "OAS2", "MYD88", "HLA-F", "KIF5C", "ZBTB20", "TRIM21", "TLR3", "MX1", "RARRES3", "GMPR", "IL15RA", "PDZD2", "TREX1", "TRIM22", "XAF1", "PSMB8", "TAP1", "BST2", "PML")

Merged <- AddModuleScore(object = Merged, features = list(hypoxia.genes), ctrl = length(hypoxia.genes), name = 'hypoxia.Score')
Merged <- AddModuleScore(object = Merged, features = list(stress.genes), ctrl = length(stress.genes), name = 'stress.Score')
Merged <- AddModuleScore(object = Merged, features = list(interferon.genes), ctrl = length(interferon.genes), name = 'interferon.Score')
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Merged <- CellCycleScoring(Merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


library(future)
plan("multicore", workers=10)
options(future.globals.maxSize = 500*1024^3)

lapply(Merged@meta.data, class)
all.genes <- rownames(Merged)
Scaled <- ScaleData(Merged, features = all.genes, vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA", "interferon.Score1"), split.by="orig.ident") # T / NK cells
#Scaled <- ScaleData(Merged, features = all.genes, vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA"), split.by="orig.ident")  # Fibroblast, B cells endothelial
#Scaled <- ScaleData(Merged, features = all.genes, vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA", "hypoxia.Score1", "stress.Score1", "S.Score", "G2M.Score"), split.by="orig.ident") # Myeloid


variable <- FindVariableFeatures(Scaled, selection.method = "vst", nfeatures = 8000)
feature <- VariableFeatures(variable)

Scaled <- RunPCA(Scaled, features = feature, npcs = 100, ndims.print = 1:5, nfeatures.print = 30, weight.by.var=T)

NPC=75 # T/ NK cells, Fibroblast, Myeloid
# NPC=45 #B cells
# NPC=65 #Endothlial
Scaled <- FindNeighbors(Scaled, reduction = "pca", dims = 1:NPC, k.param = 30, compute.SNN = TRUE, nn.method = "rann", nn.eps = 0)
Scaled <- FindClusters(Scaled, resolution = c(0.8,1,1.2, 1.4)) 
Idents(Scaled)

resol=0.8
pdf(file="Scaled_NPC32_loop_resol.pdf",  paper="a4", onefile=TRUE)
for (i in c(1:4)){
Seurat <- SetIdent(object = Scaled, value = paste('RNA_snn_res', resol, sep='.'))
Seurat <- RunUMAP(Seurat, dims = 1:NPC)
plot(UMAPPlot(Seurat, reduction = "umap", pt.size=0.25,label=T))
print (resol)
resol= resol + 0.2
}
dev.off()

resol = 'RNA_snn_res.1.2'
Scaled <- SetIdent(object = Scaled, value = resol)
Scaled <- RunUMAP(Scaled, dims = 1:NPC)

Scaled@meta.data$RNA_snn_res.1.4 <- factor(Scaled@meta.data$RNA_snn_res.1.4, levels = c(0:length (unique(Scaled@meta.data$RNA_snn_res.1.4))))
Idents(Scaled) <- "RNA_snn_res.1.4"

library (gdata)
keep(Scaled, sure=TRUE)
save.image("Scaled_NL_Tcell.RData")

################################################
###### Supp Figure 5A ########## T/NK UMAP #####
################################################

load("Scaled_NL_Tcell.RData")

Idents(Scaled) <- "final.ids"
u <- Scaled[["umap"]]@cell.embeddings
cl <- Idents(Scaled)

gg_color_hue <- function(n) {
  hues = seq(30, 500, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pal <- gg_color_hue(length(unique(cl)))

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
pal1 <- getPalette(length(unique(cl)))
palopal <- Palo(u,cl,pal1, rgb_weight=c(10,10,10))

palopal [1] <-"#977899"  # NK_cyto
palopal [2] <-"#e0dede" # UNA
palopal [3] <-"#5fb0e8" # CD8_EM
palopal [4] <-"#00C957" # CD4_EM
palopal [5] <-"#F57C7C" # CD4_low
palopal [6] <-"#afdea2" # CD4_Naive   
palopal [7] <-"#3385BB" # CD8_EMRA
palopal [8] <-"#e2cae3" # NK_rest
palopal [9] <-"#DE9E83" # Treg
palopal [10] <-"#7F9D55" # CD4_Ex
palopal [11] <-"#99caf7" # CD8_Naive
palopal [12] <- "#FE8D19" # GDT
palopal [13] <-"#3f5361" # CD8_Ex
palopal [14] <-"#d40f0f" # MAIT
palopal [15] <-"#8a5a46" # ILC3
palopal [16] <-"#FBB268" # CRTAM_GDT

pdf ('Tcells_palo.pdf')
DimPlot(Scaled, label = T, group.by='final.ids', cols=palopal, repel=F, label.box=F, label.size=4.5, raster.dpi = c(1024, 1024), raster=T, label.color= palopal) +NoLegend()
DimPlot(Scaled, label = F, group.by='final.ids', cols=palopal, repel=F, label.box=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=F, label.color= palopal) +NoLegend()
dev.off()


T_cell_annot_genes <- c(
'LEF1', 'TCF7', 'SELL', 'CCR7',  #4Naive
'LEF1', 'TCF7', 'SELL', 'NELL2', 'TOB1',#8Naive
'IL7R', 'CD40LG', 'FOS', 'JUN', #4EM
'GZMK', 'EOMES', 'ITM2C', 'CCL4',#8EM
'GZMH', 'GZMA', 'GZMB', #8EMRA
'CD200', 'NMB','CXCL13', 'TIGIT', 'BTLA', 'HAVCR2', 'TOX1',  #4EX
'PDCD1','CXCL13', 'LAG3', 'TIGIT', 'HAVCR2',#8EX
'FOXP3', 'SAT1', 'IL2RA', 'CTLA4', #REG
'AREG', 'XCL1', 'KLRC1', 'NCAM1', #NKREST 
'FCGR3A', 'FGFBP2', 'TYROBP', 'NKG7', 'KJRF1',  #NKCYTO
"SLC4A10" ,"IL23R", "IL4I1", "RORC", "CCR6", #MAIT
'TRDV2', 'TRGV9',  #gdt
'ANK3', 'NR3C2', 'PRKCA', 'THEMIS', #low
 "KIT", "TNFSF11", "PCDH9", "SPINK2", "SCN1B", "IL23R",#ILC3
"KLRC2", "KIR2DL3", "ZNF683", "KLRC3" #CRTAM_GDT
)

Idents(Scaled) <- "final.ids"

cluster.averages <- AverageExpression(Scaled, return.seurat=TRUE)
levels (cluster.averages) <-  c('CD4_Naive', 'CD8_Naive', 'CD4_EM',  'CD8_EM', 'CD8_EMRA', 'CD4_Ex', 'CD8_Ex', 'Treg' ,'NK_rest', 'NK_cyto', 'MAIT', 'GDT', 'CD4_low', 'ILC3', 'CRTAM_GDT', 'UNA')
palocol <- c("#afdea2", "#99caf7", "#00C957" , "#5fb0e8", "#3385BB", "#7F9D55", "#3f5361", "#DE9E83", "#e2cae3", "#977899", "#d40f0f", "#FE8D19", "#F57C7C", "#8a5a46", "#FBB268", "#e0dede")

pdf ('heatmap_marker_genes_tcell.pdf')
DoHeatmap(cluster.averages,features = T_cell_annot_genes, draw.lines = F)+ scale_fill_gradientn(colors = c("dodgerblue4","white","brown3")) + ggtitle ('main_cell_types_marker')
DoHeatmap(cluster.averages,features = T_cell_annot_genes, draw.lines = F, group.colors=palocol)+ scale_fill_gradientn(colors = c("dodgerblue4","white","brown3")) + ggtitle ('main_cell_types_marker')
dev.off()


################################################################
###### Supp Figure 5B ########## fibroblasts UMAP###############
################################################################

load("Scaled_NR_Fibro.RData")
library (Palo)
library (Seurat)
library (RColorBrewer)
Idents(Scaled) <- "final.ids"
u <- Scaled[["umap"]]@cell.embeddings
cl <- Idents(Scaled)

gg_color_hue <- function(n) {
  hues = seq(30, 500, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pal <- gg_color_hue(length(unique(cl)))

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
pal1 <- getPalette(length(unique(cl)))
palopal <- Palo(u,cl,pal1, rgb_weight=c(10,10,10))

palopal [1] <-"#cccc00"  # CAF_MIF
palopal [2] <-"#B294C7" # CAF_ECM
palopal [3] <-"#10c266" # CAF_COMP
palopal [4] <-"#233ed9" # iCAF_MSC
palopal [5] <-"#FE870D" # CAF_ALPL
palopal [6] <-"#ab3030" # PVL_diff   
palopal [7] <-"#e0dede" # UNA
palopal [8] <-"#fca2a2" # PVL_immature
palopal [9] <-"#00aeef" # CAF_transit

pdf ('Fibro_palo.pdf')
DimPlot(Scaled, label = T, group.by='final.ids', cols=palopal, repel=F, label.box=F, label.size=4.5, raster.dpi = c(1024, 1024), raster=T, label.color= palopal) +NoLegend()
DimPlot(Scaled, label = F, group.by='final.ids', cols=palopal, repel=F, label.box=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=F, label.color= palopal) +NoLegend()
dev.off()
DimPlot(Scaled, label = T, group.by='final.ids', cols=palopal, repel=F, label.box=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=T, label.color= palopal) +NoLegend()

annot_genes <- c(
'COX6C', 'CCND1', 'PDPF', 'RAB13', 'MIF', # CAF_MIF
'ROBO2', 'KCND2', 'MEG8', 'KAZN', 'FOXP2', #CAF ECM
'MMP11', 'COL10A1', 'COMP', 'SDC1', 'LRRC15', 'MATN3', # CAF_COMP
'APOD', 'CFD', 'GPX3', 'GPC3', 'CXCL14', 'ADH1B', #iCAF_MSC - APOD - prolif metasta -CXCL14 and 12 -- stem cell
'ALPL', 'C1QTNF3', 'SFRP4', 'CPZ', # CAF_ALPL
'PRG4', 'CLEC3B', 'PCOLCE2', 'PI16', 'ACKR3', #CAF_transit
'RERGL', 'CASQ2', 'MYH11', 'MUSTN1', 'ENTPD3', 'SUSD5', 'AKAP6', #PVL_diff
'COL25A1', 'LINC02237', 'PDE1C', 'CD36', 'ADGRF5', 'HIGD1B' #PVL_imm
 )

Idents(Scaled) <- "final.ids"
cluster.averages <- AverageExpression(Scaled, return.seurat=TRUE)
levels (cluster.averages) <-  c('CAF_MIF', 'CAF_ECM', 'CAF_COMP',  'iCAF_MSC', 'CAF_ALPL', 'CAF_transit', 'PV_diff', 'PV_immature' ,'UNA')
palocol <- c("#cccc00", "#B294C7", "#10c266" , "#233ed9", "#FE870D", "#00aeef", "#ab3030", "#fca2a2", "#e0dede")

pdf ('heatmap_marker_genes_fibro.pdf')
DoHeatmap(cluster.averages,features = annot_genes, draw.lines = F)+ scale_fill_gradientn(colors = c("dodgerblue4","white","brown3")) + ggtitle ('main_cell_types_marker')
DoHeatmap(cluster.averages,features = annot_genes, draw.lines = F, group.colors=palocol)+ scale_fill_gradientn(colors = c("dodgerblue4","white","brown3")) + ggtitle ('main_cell_types_marker')
dev.off()

#################################################
###### Supp Figure 5C ########## myeloid UMAP ###
#################################################

load("Scaled_NL_Macro.RData")

library (Palo)
library (Seurat)
library (RColorBrewer)

Scaled@meta.data$final.ids <- metadata$final.ids
Idents(Scaled) <- "final.ids"
u <- Scaled[["umap"]]@cell.embeddings
cl <- Idents(Scaled)

gg_color_hue <- function(n) {
  hues = seq(30, 500, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pal <- gg_color_hue(length(unique(cl)))

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
pal1 <- getPalette(length(unique(cl)))
palopal <- Palo(u,cl,pal1, rgb_weight=c(10,10,10))
                   
palopal [1] <-"#d4d124"  # MF_LYVE1
palopal [2] <-"#e0dede" # UNA
palopal [3] <-"#DE9E83" # MF_angio
palopal [4] <-"#3385BB" # MF_REG_M1
palopal [5] <-"#c046c7" # DC_CLEC9A
palopal [6] <-"#FE8D19" # Mono_CD16   
palopal [7] <-"#9f5da3" # DC_CD1c
palopal [8] <-"#FBB268" # Mono_CD14
palopal [9] <-"#E42622" # MF_LA_M2
palopal [10] <-"#7F9D55" # MF_BNC2
palopal [11] <-"#ab8990" # MF_IFN
palopal [12] <-"#c774cc" # DC_CCR7
palopal [13] <-"#6DBD57" # Granulocyte


pdf ('Myelo_palo.pdf')
DimPlot(Scaled, label = T, group.by='final.ids', cols=palopal, repel=F, label.box=F, label.size=4.5, raster.dpi = c(1024, 1024), raster=T, label.color= palopal) +NoLegend()
DimPlot(Scaled, label = F, group.by='final.ids', cols=palopal, repel=F, label.box=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=F, label.color= palopal) +NoLegend()
dev.off()
DimPlot(Scaled, label = T, group.by='final.ids', cols=palopal, repel=F, label.box=F, label.size=4.5, raster=T, label.color= palopal) +NoLegend()

annot_genes <- c(
'RNASE1', 'LYVE1', 'F13A1', 'SELENOP', 'MAF', #Mf_LYVE1 #### - LA TAM like - Interstitial -Kupffer
'PLDL4', 'APOE', 'APOC1', 'C1QA', 'C1QB', 'CTSD', #MF_LA_M2 #### Mf_CX3CR1_APOE - LA-TAMs - classical M2
'CXCL5', 'CCL7', 'CXCL1', 'ANGPTL4', 'IL1RN',  'VEGFA', 'MMP9',  #  MF_angio #### Mf_CX3CR1_CCL4 - iflam TAM
'CXCL10', 'CXCL9','IFIT3', 'CCL18', 'IFIT1', 'IFIT2', 'GBP1',    # MF_IFN #### Mf_MMP IFN TAM #Mf_SLC2A1 - angio TAM
'KCNQ3', 'DIRC3','BNC2', 'MERTK', 'MS4A4E',  'TPRG1', 'FTX', 'MEF2A','FRMD4A',  #MF_BNC2 #### 
'EGR2', 'PLAU','CCL4L2', 'CCL3L1', 'CCL4',  'HSPA1B',  #MF_REG #### 
'FABP5', 'COX6C', 'CD24', 'CRABP2', 'PPDPF',  #MF_rest
'FCER1A', 'CD1C', 'CLEC10A', 'CD1E', 'CST3', #DC-CD1C 
'IRF8', 'CLEC9A', 'IDO1', 'S100B', #DC_CLEC9A
'CCR7', 'CCL22', 'CCL19', 'LAMP3',  #DC_CCR7
'CD14', 'LYZ', 'S100A8', 'S100A9', 'VCAN', #Mono_CD14
'LST1', 'FCGR3A', 'CDKN1C', 'CD52', 'IFITM2' #Mono_CD16
)

cluster.averages <- AverageExpression(Scaled, return.seurat=TRUE)
levels (cluster.averages) <-  c('MF_LYVE1', 'MF_LA_M2', 'MF_angio',  'MF_IFN', 'MF_BNC2', 'MF_REG_M1', 'Granulocyte', 'DC_CD1c' ,'DC_CLEC9A' , 'DC_CCR7', 'Mono_CD14', 'Mono_CD16', 'UNA')
palocol <- c("#d4d124", "#E42622", "#ab8990" , "#DE9E83", "#7F9D55", "#3385BB", "#6DBD57", "#9f5da3", "#c046c7", "#c774cc", "#FBB268", "#FE8D19", '#e0dede')

pdf ('heatmap_marker_genes_myelo.pdf')
DoHeatmap(cluster.averages,features = annot_genes, draw.lines = F)+ scale_fill_gradientn(colors = c("dodgerblue4","white","brown3")) + ggtitle ('main_cell_types_marker')
DoHeatmap(cluster.averages,features = annot_genes, draw.lines = F, group.colors=palocol)+ scale_fill_gradientn(colors = c("dodgerblue4","white","brown3")) + ggtitle ('main_cell_types_marker')
dev.off()


#################################################
###### Supp Figure 5D ########## B cells UMAP ###
#################################################

load("Scaled_NL_Bcells_and_plasma.RData")
library (Palo)
library (Seurat)
library (RColorBrewer)
Scaled@meta.data$final.ids <- metadata$final.ids

Idents(Scaled) <- "final.ids"
u <- Scaled[["umap"]]@cell.embeddings
cl <- Idents(Scaled)

gg_color_hue <- function(n) {
  hues = seq(30, 500, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pal <- gg_color_hue(length(unique(cl)))

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
pal1 <- getPalette(length(unique(cl)))
palopal <- Palo(u,cl,pal1, rgb_weight=c(10,10,10))
                   
palopal [1] <-"#e0dede"  # UNA
palopal [2] <-"#FBB268" # Plasmablasts
palopal [3] <-"#3385BB" # B_cells_naive
palopal [4] <-"#6DBD57" # B_cells_memory

pdf ('Bcells_palo.pdf')
DimPlot(Scaled, label = T, group.by='final.ids', cols=palopal, repel=F, label.box=F, label.size=4.5, raster.dpi = c(1024, 1024), raster=T, label.color= palopal) +NoLegend()
DimPlot(Scaled, label = F, group.by='final.ids', cols=palopal, repel=F, label.box=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=F, label.color= palopal) +NoLegend()
dev.off()
DimPlot(Scaled, label = T, group.by='final.ids', cols=palopal, repel=F, label.box=F, label.size=4.5, raster=F, label.color= palopal) +NoLegend()

annot_genes <- c(
'RGS13', 'IGLC2', 'IGLC1', 'IGLC3',  #Plasmablast
'IGKC', 'MZB1','JCHAIN', 'MZB1',  #Plasmablast
'TCL1A', 'FCER2', 'IL4R', 'CD37',#B_cells_Naive
'NR4A2', 'MS4A1', 'CD69', 'CD83' #B_cells_memory
 )


cluster.averages <- AverageExpression(Scaled, return.seurat=TRUE)
levels (cluster.averages) <-  c('Plasmablasts', 'B_cells_naive', 'B_cells_memory', 'UNA')
palocol <- c( "#FBB268", "#3385BB", "#6DBD57", '#e0dede')

pdf ('heatmap_marker_genes_bcells.pdf')
DoHeatmap(cluster.averages,features = annot_genes, draw.lines = F)+ scale_fill_gradientn(colors = c("dodgerblue4","white","brown3")))
DoHeatmap(cluster.averages,features = annot_genes, draw.lines = F, group.colors=palocol)+ scale_fill_gradientn(colors = c("dodgerblue4","white","brown3"))
dev.off()

######################################################
###### Supp Figure 5E ########## Endothelial UMAP ####
######################################################

load("Scaled_NL_Endoth.RData"
library (Seurat)
library (Palo)
library(ggplot2)
library (RColorBrewer)
Idents(Scaled) <- "final.ids"
u <- Scaled[["umap"]]@cell.embeddings
cl <- Idents(Scaled)

gg_color_hue <- function(n) {
  hues = seq(30, 500, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pal <- gg_color_hue(length(unique(cl)))

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
pal1 <- getPalette(length(unique(cl)))
palopal <- Palo(u,cl,pal1, rgb_weight=c(10,10,10))

palopal [1] <-"#52AF43"  # Tip_cells
palopal [2] <-"#A6CEE3" # Endothelial_Venous
palopal [3] <-"#e0dede" # UNA
palopal [4] <-"#FF99CC" # Endothelial_Capillary
palopal [5] <-"#F06C45" # Endothelial_arterial

pdf ('Endothelial_palo.pdf')
DimPlot(Scaled, label = T, group.by='final.ids', cols=palopal, repel=F, label.box=F, label.size=4.5, raster.dpi = c(1024, 1024), raster=T, label.color= palopal) +NoLegend()
DimPlot(Scaled, label = F, group.by='final.ids', cols=palopal, repel=F, label.box=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=F, label.color= palopal) +NoLegend()
dev.off()

Idents(Scaled) <- "final.ids"
annot_genes <- c(
'ESM1', 'FSCN1', 'VWA1', 'NID2', #E_Tip_cells
'RASA4', 'TLL1', 'LRRC1', 'LYST', #E_Venous_B2M 
'ACKR1', 'CLU', 'CCL14', 'PLAT',#E_Venous_ACKR1
'CA4', 'TIMP4', 'FABP4', 'LPL',  #E_Capillary_CA4
'FBLN5', 'IGFBP3', 'CXCL12', 'JAG1'#E_Arterial
 )
 
cluster.averages <- AverageExpression(Scaled, return.seurat=TRUE)
levels (cluster.averages) <-  c('Tip_cells', 'Endothelial_Venous', 'Endothelial_Capillary',  'Endothelial_arterial', 'UNA')
palocol <- c("#52AF43", "#A6CEE3", "#FF99CC" , "#F06C45", "#e0dede")

pdf ('heatmap_marker_genes_endo.pdf')
DoHeatmap(cluster.averages,features = annot_genes, draw.lines = F)+ scale_fill_gradientn(colors = c("dodgerblue4","white","brown3")) + ggtitle ('main_cell_types_marker')
DoHeatmap(cluster.averages,features = annot_genes, draw.lines = F, group.colors=palocol)+ scale_fill_gradientn(colors = c("dodgerblue4","white","brown3")) + ggtitle ('main_cell_types_marker')
dev.off()
