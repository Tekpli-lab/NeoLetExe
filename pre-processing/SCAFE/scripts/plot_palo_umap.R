library(Seurat)
library(ggplot2)
library(Palo)
library(RColorBrewer)
library(optparse)


option_list = list(make_option(c("-o", "--output"), type = "character", default = NULL, help = "Path to output pdf", metavar = "character"),
		   make_option(c("-i", "--input"), type = "character", default = NULL, help = "Input seurat object .rds file path", metavar = "character"))

opt = parse_args(OptionParser(option_list = option_list))
input_path <- opt$input
output_path <- opt$output

seurat_object <- readRDS(input_path)


Idents(seurat_object) <- "main.ids"
u <- seurat_object[["umap"]]@cell.embeddings
cl <- Idents(seurat_object)

gg_color_hue <- function(n) {
  hues = seq(30, 500, length = n + 1)
  hues = seq(30, 500, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pal <- gg_color_hue(length(unique(cl)))

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
pal1 <- getPalette(length(unique(cl)))
palopal <- Palo(u,cl,pal1, rgb_weight=c(10,10,10))


palopal ["Epithelial_cells"] <-"#A6CEE3" # Epith
palopal ["CAF"] <-"#F06C45" # CAF
palopal ["T_cells"] <-"#B294C7" # Tcell
palopal ["B_cells"] <- "#AB9AFF" # Bcells
palopal ["Plasmablasts"] <- "#FF9933"# Plasma
palopal ["Endothelial_cells"] <- "#52AF43" # Endoth
palopal ["Myeloid_cells"] <- "#FF99CC" # Myeloid
palopal ["PV_Fibroblasts"] <- "#DAA70F" # PVFIBRO
palopal ["Mast_cells"] <- "#666666"# Mast
palopal ["Neutrophils"] <- "#006633" # Neutro
palopal ["Adipocytes"] <- "#996600" # Adipocy
palopal ["DC"] <- "#9A4D73" # DC

med <- aggregate(u,list(cl),median)
colnames(med)[1] <- 'main.ids'

d = cbind(u, seurat_object[["main.ids"]])

pdf(output_path)
ggplot() +
  geom_point(data=d,aes(x=UMAP_1,y=UMAP_2,col=main.ids),size=0.1,alpha=0.3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  scale_color_manual(values=palopal) +
  theme(legend.position = 'none')

ggplot() +
  geom_point(data=d,aes(x=UMAP_1,y=UMAP_2,col=main.ids),size=0.1,alpha=0.3) +
  theme_void() +
  scale_color_manual(values=palopal) +
  theme(legend.position = 'none')

ggplot() +
  geom_point(data=d,aes(x=UMAP_1,y=UMAP_2,col=main.ids),size=0.1,alpha=0.3) +
  geom_text(data=med,aes(x=UMAP_1,y=UMAP_2,label=main.ids),size=5) +
  theme_void() +
  scale_color_manual(values=palopal) +
  theme(legend.position = 'none')

ggplot() +
  geom_point(data=d,aes(x=UMAP_1,y=UMAP_2,col=main.ids),size=0.1,alpha=0.3) +
  geom_text(data=med,aes(x=UMAP_1,y=UMAP_2,label=main.ids),size=5) +
  theme_void() +
  scale_color_manual(values=palopal)
dev.off()



