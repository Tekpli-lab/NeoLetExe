##################################################
####### FIGURE 2 ##################################
###################################################

###################################################
###### Figure 2A ########## FISH PLOTS ############
###################################################

### example of fishplot with patient 08

### in bash
sed '/^norm_/d' <  P08.cell_groupings > P08i.cell_groupings 
##### .cell_groupings are obtained by running PP_InferCNV_per_sample.r

library (fishplot)
patient_number <- '08' 

load("metadata_git.RData")
#### load the metadata available at EGA
target <- data.frame(read.table(paste0("P", patient_number, ".cell_groupings"), header=TRUE, sep="\t", stringsAsFactors=F))
#### Need the output of the InferCNV per sample PP_InferCNV_per_sample.r

counts_tot <- as.data.frame.matrix(table(metadata$time_point, metadata$patient_number))
n_time <- counts_tot[, colnames(counts_tot) %in% as.numeric(patient_number)]

target$cell_group_name <- gsub("all_observations.all_observations.", "", target$cell_group_name)
target <- target[substr(target$cell, 1,1) == "P",]
target$pat_num <- gsub("_([^å]*)", "", target$cell)

tab <- as.matrix(as.data.frame.matrix(prop.table(table(target$cell_group_name, target$pat_num)) * 100))
counts <- as.matrix(as.data.frame.matrix(table(target$cell_group_name,target$pat_num)))

tab <- tab[rowSums(tab) > 5,]
timepoints = c(0, 200, 400)
n_time <- n_time[c(1, 2, 3)]
prop <- t(round((t(counts) / n_time) * 1000000))
prop <- (prop / sum(prop)) * 100
prop <- prop[rownames(prop) %in% rownames(tab),]

par1 <- colSums(prop[1:4,])
par2 <- prop[5,]
prop1 <- rbind(par1, par2, prop)
parents = c(0, 0, 1, 1, 1, 1, 2)
fish = createFishObject(prop1, parents, timepoints=timepoints)
fish = layoutClones(fish)
fish = setCol(fish, col=c("#9DABC5", "#81BED0", "#D89560", "#85593E", "#C55E7B", "#9F8272", "#365584"))


pdf(paste0('fishplot_NL', patient_number, '.pdf'))
fishPlot(fish, shape="spline", vlines=timepoints, vlab=c("A", "B", "C"), bg.type='solid', bg.col='white', col.vline="black",
  border=0.5, cex.title=0.8, cex.vlab=3, pad.left=0.1, col.border='black')
dev.off()


###################################################
###### Figure 2B ########## CIRCOS PLOTS ##########
###################################################

library(Seurat)
library(SeuratData)
library(Palo)
library(ggplot2)
library (RColorBrewer)
library (tabula)
library (stringi)
library (ggrepel)
library(ggsignif)
library(plyr)
library(dplyr)
library(future)
library (tidyverse)
library (stringr)
library(circlize)

test <-  data.frame(read.table("results_gene_status_stats.txt", header=T, sep="\t", stringsAsFactors=F))
#### Can be found in Supplementary Tables and results from PP_CNA_genes_level

test1 <- (test[test$n_gene.x > 50,])
amp_agress <- test1[(test1$fdr_agress_amp) < 0.05 & (test1$enrich_agress_amp) >2 & !is.na(test1$fdr_agress_amp) & !is.na(test1$enrich_agress_amp), ]
amp_sensit <- test1[(test1$fdr_sensit_amp) < 0.05 & (test1$enrich_sensit_amp) >2 & !is.na(test1$fdr_sensit_amp) & !is.na(test1$enrich_sensit_amp), ]
del_agress <- test1[(test1$fdr_agress_del) < 0.05 & (test1$enrich_agress_del) >2 & !is.na(test1$fdr_agress_del) & !is.na(test1$enrich_agress_del), ]
del_sensit <- test1[(test1$fdr_sensit_del) < 0.05 & (test1$enrich_sensit_del) >2 & !is.na(test1$fdr_sensit_del) & !is.na(test1$enrich_sensit_del), ]

circos.par("track.height" = 1)
chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                 "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                 "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"
                 )
chrom_lengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                   159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                   115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                   59128983, 63025520, 48129895, 51304566)

centromeres <- c(
  chr1=121535434,  chr2=92326171,  chr3=90504854,  chr4=49660117,
  chr5=46405641,   chr6=58830166,  chr7=58054331,  chr8=43838887,
  chr9=47367679,   chr10=39254935, chr11=51644205, chr12=34856694,
  chr13=16000000,  chr14=16000000, chr15=19000000, chr16=35335801,
  chr17=22263006,  chr18=15460898, chr19=24681782, chr20=26369569,
  chr21=11288129,  chr22=13000000
)

df_agress_amp = data.frame(sectors = amp_agress$chr, x = (as.numeric(amp_agress$start) + as.numeric(amp_agress$end))/2, y = 1, col= 'red')
df_agress_amp <- df_agress_amp[!df_agress_amp$sectors=="chrMT",]
					  
df_agress_del = data.frame(sectors = del_agress$chr,x = (as.numeric(del_agress$start) + as.numeric(del_agress$end))/2, y = 1,col= 'blue') 
df_agress_del <- df_agress_del[!df_agress_del$sectors=="chrMT",]

df_amp_sensit = data.frame(sectors = amp_sensit$chr,x = (as.numeric(amp_sensit$start) + as.numeric(amp_sensit$end))/2,  y = 1,col= 'red')
df_amp_sensit <- df_amp_sensit[!df_amp_sensit$sectors=="chrMT",]
					  
df_del_sensit = data.frame(sectors = del_sensit$chr, x = (as.numeric(del_sensit$start) + as.numeric(del_sensit$end))/2, y = 1,col= 'blue')  
df_del_sensit <- df_del_sensit[!df_del_sensit$sectors=="chrMT",]		
			  
pdf('Figure 2B Circos plot resitant clones.pdf')			  
circos.clear()
col_text <- "black"
circos.par("track.height"=0.1, gap.degree=2, cell.padding=c(0, 0, 0, 0), track.margin = c(0.04, 0.04))
circos.initialize(factors = chromosomes, xlim = matrix(c(rep(0, 22), chrom_lengths), ncol=2))

# genomes
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
chr=CELL_META$sector.index
xlim=CELL_META$xlim
ylim=CELL_META$ylim
circos.text(mean(xlim), mean(ylim), chr, cex=1, col=col_text, 
facing="bending.inside", niceFacing=TRUE, font=2)
}, bg.col=NA, bg.border=F, track.height=0.1)

# genomes x axis
brk <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)*10^8
circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
circos.axis(h="top", major.at=brk, labels=round(brk/10^8, 1), labels.cex=0.7, 
col=col_text, labels.col=col_text, lwd=0.7, labels.facing="clockwise",  labels.font=2)
}, bg.border=F)

circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
chr=CELL_META$sector.index
xlim=CELL_META$xlim
ylim=CELL_META$ylim  
cen <- centromeres[chr]
circos.lines(c(cen, cen), c(ylim[1], ylim[2]), sector.index = chr, col = "black", lty = 2, lwd = 0.5)
}, bg.border=T, track.height=0.27)
circos.trackLines(df_agress_amp$sectors, df_agress_amp$x, df_agress_amp$y, col = df_agress_amp$col, type='h')

circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
chr=CELL_META$sector.index
xlim=CELL_META$xlim
ylim=CELL_META$ylim
cen <- centromeres[chr]
circos.lines(c(cen, cen), c(ylim[1], ylim[2]), sector.index = chr, col = "black", lty = 2, lwd = 0.5)
}, bg.border=T, track.height=0.27)
circos.trackLines(df_agress_del$sectors, df_agress_del$x, df_agress_del$y, col = df_agress_del$col, type='h')
dev.off()


pdf('Figure 2B Circos plot sensitive clones.pdf')	
circos.clear()
col_text <- "black"
circos.par("track.height"=0.1, gap.degree=2, cell.padding=c(0, 0, 0, 0), track.margin = c(0.04, 0.04))
circos.initialize(factors = chromosomes, xlim = matrix(c(rep(0, 22), chrom_lengths), ncol=2))

# genomes
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
chr=CELL_META$sector.index
xlim=CELL_META$xlim
ylim=CELL_META$ylim
circos.text(mean(xlim), mean(ylim), chr, cex=1, col=col_text, 
facing="bending.inside", niceFacing=TRUE, font=2)
}, bg.col=NA, bg.border=F, track.height=0.1)

# genomes x axis
brk <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)*10^8
circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
circos.axis(h="top", major.at=brk, labels=round(brk/10^8, 1), labels.cex=0.7, 
col=col_text, labels.col=col_text, lwd=0.7, labels.facing="clockwise",  labels.font=2)
}, bg.border=F)

circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
chr=CELL_META$sector.index
xlim=CELL_META$xlim
ylim=CELL_META$ylim
cen <- centromeres[chr]
circos.lines(c(cen, cen), c(ylim[1], ylim[2]), sector.index = chr, col = "black", lty = 2, lwd = 0.5)

}, bg.border=T, track.height=0.27)
circos.trackLines(df_amp_sensit$sectors, df_amp_sensit$x, df_amp_sensit$y, col = df_amp_sensit$col, type='h')

circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
chr=CELL_META$sector.index
xlim=CELL_META$xlim
ylim=CELL_META$ylim
cen <- centromeres[chr]
circos.lines(c(cen, cen), c(ylim[1], ylim[2]), sector.index = chr, col = "black", lty = 2, lwd = 0.5)

}, bg.border=T, track.height=0.27)
circos.trackLines(df_del_sensit$sectors, df_del_sensit$x, df_del_sensit$y, col = df_del_sensit$col, type='h')

dev.off()



#################################################################################################
###### Figure 2C ########## Coorelation heatmap between genes sensitive and resistant ###########
#################################################################################################


library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (mixtools)
library (pheatmap)
library (amap)
library (gwasrapidd)
library (ggplot2)

mean_mat <-  data.frame(read.table("correlation_matrix.txt", header=T, sep="\t", stringsAsFactors=F))
#correlation_matrix output of PP_differential_expression_clones

gene_status <- data.frame(read.table("gene_sensitive_resistant_annotation.txt", header=TRUE, sep="\t", stringsAsFactors=F))
#gene_status found in supplementary Table 7

colnames(mean_mat) <- gsub("\\.(?![0-9])", "-", colnames(mean_mat), perl = TRUE)
colnames(mean_mat) <- gsub("IGKV3.11", "IGKV3-11", colnames(mean_mat))
colnames(mean_mat) <- gsub("IGKV3.20", "IGKV3-20", colnames(mean_mat))

breaksList = seq(-0.5, 0.5, by = 0.05)
col.pal <- colorRampPalette(c("dodgerblue4","white","brown3"))(length(breaksList))

dcols1 <- hcluster(t(mean_mat), method = 'correlation', link = 'ward')
drows1 <- hcluster(mean_mat, method = 'correlation', link = 'ward')
cutree1<-cutree(drows1, k=10)
cutree1<- cutree1[rownames(mean_mat)]


cutree1[cutree1==1]<- "cluster1"
cutree1[cutree1==2]<- "cluster2"
cutree1[cutree1==3]<- "cluster3"
cutree1[cutree1==4]<- "cluster4"
cutree1[cutree1==5]<- "cluster5"
cutree1[cutree1==6]<- "cluster6"
cutree1[cutree1==7]<- "cluster7"
cutree1[cutree1==8]<- "cluster8"
cutree1[cutree1==9]<- "cluster9"
cutree1[cutree1==10]<- "cluster10"
cutree1[cutree1==11]<- "cluster11"
cutree1[cutree1==12]<- "cluster12"

all (colnames (mean_mat)==rownames (gene_status))
annotation_row <- data.frame(cbind( cutree1, gene_status))
colnames (annotation_row) <- c('clusters', 'cutree', 'status')

annotation_col <- annotation_row
annotation_col$status <- gsub ('sensit', 'sensitive', annotation_col$status)
annotation_col$status <- gsub ('agress', 'resistant', annotation_col$status)


anno_colors = list(clusters = c("cluster1" = "darkgreen", 
								"cluster2" = "darkorange",
								"cluster3" = "blue",
								"cluster4" = "darkgoldenrod",
								"cluster5" = "pink",
								"cluster6" = "coral",
								"cluster7" = "brown",
								"cluster8" = "cyan",
								"cluster9" = "chartreuse",
								"cluster10" = "purple"),
status = c("sensitive" = "green", "resistant" = "red", "both"="yellow"))



pdf ('Figure 2C correalation heatmap res sens.pdf')
pheatmap(mean_mat, color = col.pal, show_rownames = F, show_colnames = F, main = "",  
cluster_rows = drows1, cluster_cols = dcols1, breaks = breaksList, annotation_col=annotation_col, annotation_colors = anno_colors, fontsize = 10, border_color = NA)
dev.off()


#################################################################################################
###### Figure 2D ########## Coorelation between scores ##########################################
#################################################################################################

library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (mixtools)
library (pheatmap)
library (amap)
library (gwasrapidd)
library (ggplot2)
library(reshape2)
library(dendextend)
library(ggdendro)

cor_sco <- as.matrix(data.frame(read.table("correlation_matrices_scores.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)))
#correlation matrix output for PP_gene_cluster_correlation

rownames(cor_sco) <- gsub('[_]{2}', '', rownames(cor_sco))
colnames(cor_sco) <- gsub('[_]{2}', '', colnames(cor_sco))
rownames(cor_sco) <- gsub('_score', '', rownames(cor_sco))
colnames(cor_sco) <- gsub('_score', '', colnames(cor_sco))

rownames(cor_sco) <- gsub("[._]", " ", rownames(cor_sco))
colnames(cor_sco) <- gsub("[._]", " ", colnames(cor_sco))
rownames(cor_sco) <- gsub("\\s+", " ", rownames(cor_sco))
colnames(cor_sco) <- gsub("\\s+", " ", colnames(cor_sco))
rownames(cor_sco) <- trimws(rownames(cor_sco))
colnames(cor_sco) <- trimws(colnames(cor_sco))

dcols <- hcluster(t(cor_sco), method = 'correlation', link = 'ward')
drows <- hcluster(cor_sco, method = 'correlation', link = 'ward')

get_lower_tri <- function(cormat) {
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
lower_tri <- get_lower_tri(cor_sco[drows$order, dcols$order])
melted_lower_tri <- melt(lower_tri, na.rm = TRUE)

p <- ggplot(melted_lower_tri, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "brown3", midpoint = 0) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7,
      color = c(rep("black", 4), "darkgreen", rep("black", 11), "brown", rep("black", 5), "coral", 
                "black", "cyan", rep("black", 5), "purple", rep("black", 3), "chartreuse", 
                "darkorange", "black", "black", "darkgoldenrod", rep("black", 14), "pink", 
                "black", "black", "blue", "black")),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = unit(c(1, 1, 1, 3), "lines"),
    panel.border = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_fixed()
print(p)
ggsave('Figure 2Da correlation triangle scores.pdf', plot = p, width = 10, height = 8)

library(dendextend)
library(ggdendro)
library(ggplot2)
custom_colors <- c("red", "blue", "green", "purple") 
dend_object <- hcluster(t(cor_sco), method = 'correlation', link = 'ward')
dend_data <- as.dendrogram(dend_object)
dend_colored <- dend_data %>% color_branches(k = 4, col = custom_colors)  
dend_colored <- dend_colored %>% set("branches_lwd", 1)  %>% set("labels_cex", 0.4)
pdf ('Figure 2Db dendrogram.pdf',height = 7 / 3, width = 7)
plot(dend_colored, lwd = 4, axes = T, cex.lab = 0.5, cex.axis = 0.5, cex = 0.2 )
dev.off()



#################################################################################################
###### Figure 2E ########## Scores in clones resistant and sensitive clones #####################
#################################################################################################
library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (mixtools)
library (pheatmap)
library (amap)
library (gwasrapidd)
library (ggplot2)
library(ggplot2)
library(reshape2)
library(Seurat)
library(ggplot2)
library(ggrepel)
library (msigdbr)
library (dplyr)
library (clusterProfiler)
library (dplyr)

scores <- as.data.frame(data.frame(read.table("scores_gene_cluster_clones.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)))
#scores generated  in PP_gene_cluster_correlation

load("metadata_git.RData")

 library (dplyr)
 df_averages <- data.frame (scores %>%
  group_by(clone_id) %>%
  summarise(across(starts_with("cluster"), mean, na.rm = TRUE))
)

rownames (df_averages) <- df_averages[,1]
df_averages <- t(df_averages[, -1])
df_averages <- t(scale (t(df_averages)))


clones_ord <-c( unique (metadata$clone_id[metadata$clone_cat2=="res" &!is.na(metadata$clone_cat2)]), unique (metadata$clone_id[metadata$clone_cat2=="sens" &!is.na(metadata$clone_cat2)]))
clones_cat <-c( rep ('res' , length (unique (metadata$clone_id[metadata$clone_cat2=="res" &!is.na(metadata$clone_cat2)]))), rep ('sens', length (unique (metadata$clone_id[metadata$clone_cat2=="sens" &!is.na(metadata$clone_cat2)]))))

colnames (df_averages) <- gsub ("-", "_", colnames (df_averages))
df_averages <- df_averages[, match(clones_ord, colnames(df_averages))]

mean_mat <- as.matrix (df_averages)

breaksList = seq(-3, 3, by = 0.05)
col.pal <- colorRampPalette(c("dodgerblue4","white","brown3"))(length(breaksList))

dcols <- hcluster(t(mean_mat), method = 'pearson', link = 'complete')
drows <- hcluster(mean_mat, method = 'euclidean', link = 'ward')



CL <- data.frame(cbind (clones_ord, clones_cat))
annotation_col <- data.frame(CL[,-1])
rownames (annotation_col) <- CL$clones_ord
colnames (annotation_col) <- 'clones_cat'
annotation_col$PAM50 <- ifelse (gsub("^.|_.*", "", rownames (annotation_col)) %in% c('02', '06', '09', '11', '20', '22',  '25'), "LumA",gsub("^.|_.*", "", rownames (annotation_col)))
annotation_col$PAM50 <- ifelse (annotation_col$PAM50%in% c('01', '03',  '05', '07', '08', '13', '15', '16', '17', '18', '19',  '21', '23', '24'), "LumB", annotation_col$PAM50)
annotation_col$PAM50 <- ifelse (annotation_col$PAM50%in% c('12', '04'), NA, annotation_col$PAM50)
annotation_col$PAM50 <- ifelse (annotation_col$PAM50%in% c('10', '14'), "Her2_enriched", annotation_col$PAM50)

annotation_col$patient_number <- gsub("P(.*?)_.*", "\\1", rownames (annotation_col))


anno_colors = list(ID = c("sens" = "green", "res" = "red"),
clones_cat = c("sens" = "green", "res" = "red"),   PAM50 = c("LumA" = "blue", "LumB" = "light blue", "Her2_enriched"="deep pink"),
patient_number= c(
"01" = "darkgreen", 
"02" = "darkorange",  
"03" = "blue", 
"04" = "darkgoldenrod",
"05" = "pink", 
"06" = "coral", 
"07" = "deeppink4", 
"08" = "brown", 
"09" = "cyan", 
"10" = "chartreuse", 
"11" = "purple", 
"12" = "black", 
"13" = "gray", 
"14" = "orange", 
"15" = "deeppink", 
"16" = "red", 
"17" = "chocolate",
"18" = "yellow",
"19" = "yellow3", 
"20" = "cadetblue", 
"21" = "darkblue", 
"22" = "khaki", 
"23" = "violet", 
"24" = "turquoise", 
"25" = "limegreen"))

pdf ('Figure 2E Modules average in clones.pdf')
pheatmap(mean_mat, color = col.pal, show_rownames = T, show_colnames = T, main = "",  cluster_rows = drows, cluster_cols = dcols, 
breaks = breaksList, annotation_col=annotation_col, annotation_colors = anno_colors, fontsize = 10, border_color = NA)
dev.off()

