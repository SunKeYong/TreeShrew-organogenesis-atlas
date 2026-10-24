# umap visualization of all Treeshrew cCREs in Fig2
# 
setwd("path/Peak_embedding")
#
TS_Peak_UMAP_meta <- read.csv("TreeShrew_peaks_adata_UMAP_metadata.csv")
head(TS_Peak_UMAP_meta)
#
# 创建一个伪Seurat文件，方便分析数据
#
x.num <- c(rep(1,25000000),rep(0,25000000),rep(2,25000000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(TS_Peak_UMAP_meta$cell)),
                      sample(x.num,length(TS_Peak_UMAP_meta$cell)),
                      sample(x.num,length(TS_Peak_UMAP_meta$cell)),
                      sample(x.num,length(TS_Peak_UMAP_meta$cell)),
                      sample(x.num,length(TS_Peak_UMAP_meta$cell)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- TS_Peak_UMAP_meta$cell
#
TS_Peak_UMAP <- CreateSeuratObject(counts = t(x.num.sample), project = "TS_Peak", min.cells = 0, min.features = 0)
TS_Peak_UMAP
TS_Peak_UMAP <- NormalizeData(TS_Peak_UMAP, normalization.method = "LogNormalize", scale.factor = 10000)
TS_Peak_UMAP <- FindVariableFeatures(TS_Peak_UMAP, selection.method = "vst", nfeatures = 4)
TS_Peak_UMAP <- ScaleData(TS_Peak_UMAP, features = rownames(TS_Peak_UMAP))
variablegene <- VariableFeatures(object = TS_Peak_UMAP)
TS_Peak_UMAP <- RunPCA(TS_Peak_UMAP, features = variablegene,npcs =2,reduction.name = "UMAP",reduction.key = "UMAP")
TS_Peak_UMAP
#
#
s.UMAP <- cbind(TS_Peak_UMAP_meta$x,
                TS_Peak_UMAP_meta$y)
#
rownames(s.UMAP ) <- TS_Peak_UMAP_meta$cell
colnames(s.UMAP ) <- c("UMAP_1","UMAP_2")
TS_Peak_UMAP@reductions$UMAP@cell.embeddings <- s.UMAP 
#
DimPlot(TS_Peak_UMAP,reduction = "UMAP",alpha = 0.03,raster = T)+NoLegend()
#
TS_Peak_UMAP$peak_clusters_0.8 <- TS_Peak_UMAP_meta$peak_clusters_0.8
TS_Peak_UMAP$peak_clusters_1 <- TS_Peak_UMAP_meta$peak_clusters_1
TS_Peak_UMAP$peak_clusters_2 <- TS_Peak_UMAP_meta$peak_clusters_2
#
DimPlot(TS_Peak_UMAP,reduction = "UMAP",alpha = 0.03,raster = T,label = T,group.by = "peak_clusters_0.8")+NoLegend()
DimPlot(TS_Peak_UMAP,reduction = "UMAP",alpha = 0.03,raster = T,label = T,group.by = "peak_clusters_1")+NoLegend()
DimPlot(TS_Peak_UMAP,reduction = "UMAP",alpha = 0.03,raster = T,label = T,group.by = "peak_clusters_2")+NoLegend()
#
#
TS_Peak_UMAP@active.ident <- as.factor(TS_Peak_UMAP$peak_clusters_1)
table(TS_Peak_UMAP@active.ident)
#
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==0),]$cell) <- "Universal_peaks"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==1),]$cell) <- "Epithelial"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==2),]$cell) <- "Endothelium"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==3),]$cell) <- "Mesoderm"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==4),]$cell) <- "Mesoderm"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==5),]$cell) <- "Universal_peaks"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==6),]$cell) <- "Kidney"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==7),]$cell) <- "Epithelial"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==8),]$cell) <- "Immune"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==9),]$cell) <- "Neurons"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==10),]$cell) <- "Neurons"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==11),]$cell) <- "PNS_glia"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==12),]$cell) <- "Gut"   
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==13),]$cell) <- "Immune"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==14),]$cell) <- "Erythroid"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==15),]$cell) <- "Retina"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==16),]$cell) <- "Hepatocytes"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==17),]$cell) <- "Neurons"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==18),]$cell) <- "PNS_neurons"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==19),]$cell) <- "Olfactory"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==20),]$cell) <- "Muscle"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==21),]$cell) <- "Mesoderm"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==22),]$cell) <- "Neurons"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==23),]$cell) <- "Retina"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==24),]$cell) <- "Autonomic_nervous_system"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==25),]$cell) <- "Adrenocortical"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==26),]$cell) <- "Ependymal"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==27),]$cell) <- "Epithelial"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==28),]$cell) <- "Neurons"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==29),]$cell) <- "Neurons"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==30),]$cell) <- "Cardiomyocytes"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==31),]$cell) <- "Lung"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==32),]$cell) <- "Mesoderm"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==33),]$cell) <- "Epithelial"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==34),]$cell) <- "Mesoderm"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==35),]$cell) <- "Lens"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==36),]$cell) <- "Universal_peaks"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_1==37),]$cell) <- "Universal_peaks"
#
DimPlot(TS_Peak_UMAP,reduction = "UMAP",alpha = 0.03,raster = T,label = T)+NoLegend()
#
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_2==53),]$cell) <- "Neurons"
Idents(TS_Peak_UMAP,TS_Peak_UMAP_meta[which(TS_Peak_UMAP_meta$peak_clusters_2==49),]$cell) <- "PNS_neurons"
#
#
DimPlot(TS_Peak_UMAP,reduction = "UMAP",alpha = 0.03,raster = T,label = T)+NoLegend()
#
TS_Peak_UMAP$Peak_ID <- as.character(TS_Peak_UMAP@active.ident)
#
TS_Peak_UMAP@active.ident <- as.factor(TS_Peak_UMAP$Peak_ID)
levels(TS_Peak_UMAP)
#
TS_peak_color <- c("#DDAD3B","#00BFFF","#9970AB","#FF7F00",
                   "#32CD32","#3288BD","#48D1CC","#B15928",
                   "#6495ED","#D6604D","#EE82EE","red",
                   "#00FF7F","#FDCDAC","#9E0142","#B3CDE3",
                   "#87CEFA","#CDE8B9","#0000FF","#0099FF",
                   "#DE77AE")
#
DimPlot(TS_Peak_UMAP,reduction = "UMAP",alpha = 0.03,raster = T,label = T,cols = TS_peak_color)+NoLegend()
DimPlot(TS_Peak_UMAP,reduction = "UMAP",alpha = 0.1,raster = F,label = T,cols = TS_peak_color)+NoLegend()
#
p <- DimPlot(TS_Peak_UMAP,cols = TS_peak_color,label = F,
             repel = T,pt.size = 0.00000000002,
             reduction = "UMAP",
             shuffle=T,raster=F,alpha = 0.05,
             raster.dpi = c(1312, 1312))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
#
ggsave("./TS_Peak_UMAP1.png",plot=p,device="png",dpi=300,units = "cm",width = 20,height = 20)
#
p <- DimPlot(TS_Peak_UMAP,cols = TS_peak_color,label = F,
             repel = T,pt.size = 0.00000000002,
             reduction = "UMAP",
             shuffle=T,raster=F,alpha = 0.03,
             raster.dpi = c(1312, 1312))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
#
ggsave("./TS_Peak_UMAP2.png",plot=p,device="png",dpi=300,units = "cm",width = 20,height = 20)
#
#
#
cistopic1 <- read.csv("../Topics_otsu/topic_62.bed",sep = "\t",header = F)
cistopic1_tmp <- paste(cistopic1$V1,cistopic1$V2+1,sep = ":");cistopic1_tmp <- paste(cistopic1_tmp,cistopic1$V3,sep="-");cistopic1_tmp <- intersect(cistopic1_tmp,Cells(TS_Peak_UMAP))
TS_Peak_UMAP$tmp <- "tmp";TS_Peak_UMAP@active.ident <- as.factor(TS_Peak_UMAP$tmp);Idents(TS_Peak_UMAP,cells=cistopic1_tmp) <- "cistopic1";levels(TS_Peak_UMAP) <- c("tmp","cistopic1")
#
p <- DimPlot(TS_Peak_UMAP,cols = c("gray","#3FAD32"),label = F,
             repel = T,pt.size = 0.00000000002,cells.highlight =cistopic1_tmp,
             reduction = "UMAP",  sizes.highlight = 0.3,  cols.highlight = "#407834",
             shuffle=T,raster=F,alpha = 0.03,
             raster.dpi = c(1312, 1312))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("./TS_Peak_UMAP.topic_62.celltypeC1.1.png",plot=p,device="png",dpi=300,units = "cm",width = 20,height = 20)
#
cistopic1 <- read.csv("../Topics_otsu/topic_106.bed",sep = "\t",header = F)
cistopic1_tmp <- paste(cistopic1$V1,cistopic1$V2+1,sep = ":");cistopic1_tmp <- paste(cistopic1_tmp,cistopic1$V3,sep="-");cistopic1_tmp <- intersect(cistopic1_tmp,Cells(TS_Peak_UMAP))
TS_Peak_UMAP$tmp <- "tmp";TS_Peak_UMAP@active.ident <- as.factor(TS_Peak_UMAP$tmp);Idents(TS_Peak_UMAP,cells=cistopic1_tmp) <- "cistopic1";levels(TS_Peak_UMAP) <- c("tmp","cistopic1")
#
p <- DimPlot(TS_Peak_UMAP,cols = c("gray","#3FAD32"),label = F,
             repel = T,pt.size = 0.00000000002,cells.highlight =cistopic1_tmp,
             reduction = "UMAP",  sizes.highlight = 0.3,  cols.highlight = "#407834",
             shuffle=T,raster=F,alpha = 0.03,
             raster.dpi = c(1312, 1312))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("./TS_Peak_UMAP.topic_106.celltypeC6.1.png",plot=p,device="png",dpi=300,units = "cm",width = 20,height = 20)
#
cistopic1 <- read.csv("../Topics_otsu/topic_134.bed",sep = "\t",header = F)
cistopic1_tmp <- paste(cistopic1$V1,cistopic1$V2+1,sep = ":");cistopic1_tmp <- paste(cistopic1_tmp,cistopic1$V3,sep="-");cistopic1_tmp <- intersect(cistopic1_tmp,Cells(TS_Peak_UMAP))
TS_Peak_UMAP$tmp <- "tmp";TS_Peak_UMAP@active.ident <- as.factor(TS_Peak_UMAP$tmp);Idents(TS_Peak_UMAP,cells=cistopic1_tmp) <- "cistopic1";levels(TS_Peak_UMAP) <- c("tmp","cistopic1")
#
p <- DimPlot(TS_Peak_UMAP,cols = c("gray","#3FAD32"),label = F,
             repel = T,pt.size = 0.00000000002,cells.highlight =cistopic1_tmp,
             reduction = "UMAP",  sizes.highlight = 0.3,  cols.highlight = "#407834",
             shuffle=T,raster=F,alpha = 0.03,
             raster.dpi = c(1312, 1312))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("./TS_Peak_UMAP.topic_134.celltypeC10.1.png",plot=p,device="png",dpi=300,units = "cm",width = 20,height = 20)
#
#
##################################################
# # umap visualization of gini score
# calculate gini index
Ts_peak_matrix <- read.csv("matrix_count.txt",header = T,sep = "\t")
dim(Ts_peak_matrix)
colnames(Ts_peak_matrix) <- unlist(strsplit(colnames(Ts_peak_matrix),".p"))[seq(1,2*197,2)]
Ts_peak_matrix[1:5,1:5]
#
peak_name <- read.csv("TreeShrew_keep_Final.matrix.peak_counts.bed",header = F,sep = "\t")
#
peak_name <- peak_name[,c(1,2,3,4)]
colnames(peak_name) <- c("chr","start","end","treeshrew_peak")
head(peak_name)
#
Ts_peak_matrix2 <- cbind(peak_name,Ts_peak_matrix)
Ts_peak_matrix2[1:5,1:15]
#
Ts_peak_matrix3 <- Ts_peak_matrix2[,-4]
head(Ts_peak_matrix3)
#
write.table(Ts_peak_matrix3,file = "Treeshrew.Raw.PeakMatrix.Per.CellType.bed", sep = "\t", col.names = F, row.names = F, quote = F)
#
count_mat <- as.matrix(Ts_peak_matrix2[, -c(1:4)])
libsize <- colSums(count_mat)
# CPM
cpm_mat <- t(t(count_mat) / libsize * 1e6)
cpm_mat[1:5, 1:5]
count_mat[1:5, 1:5]
#
######
# 2. 按文章公式计算 Gini
gini_paper <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  
  if (sum(x) == 0) {
    return(NA_real_)
  }
  
  x <- sort(x)              # 文章要求 ascending order
  n <- length(x)
  i <- seq_len(n)
  
  gini <- sum((2 * i - n - 1) * x) / (n * sum(x))
  
  return(gini)
}
#
#
# 3. 每个 cCRE 计算一个 Gini
gini_score <- apply(cpm_mat, 1, gini_paper)
head(gini_score)
#
gini_score2 <- as.data.frame(gini_score)
head(gini_score2)
hist(gini_score2$gini_score)
#
TS_Peak_UMAP$gini_score <- gini_score2$gini_score
FeaturePlot(TS_Peak_UMAP,features = "gini_score",reduction = "UMAP")+NoLegend()
#
library(RColorBrewer)
FeaturePlot(TS_Peak_UMAP,features = "gini_score",raster = T,alpha = 0.03,reduction = "UMAP",
                 pt.size = 0.00000000002,raster.dpi = c(1312, 1312))+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))[1:10])
##################
#
ggsave("TS_Gini_Score.pdf", plot = p, width = 8, height = 6,device="pdf")
##
p2 <- FeaturePlot(TS_Peak_UMAP,features = "gini_score",raster = F,alpha = 0.1,pt.size = 0.00000000002,reduction = "UMAP")+ 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))[1:10])+NoLegend()+labs(title = NULL)+theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+theme(axis.line.x=element_blank(),axis.line.y=element_blank())+theme(panel.background = element_rect( colour = "black", size = 0))+theme(panel.border = element_blank())
ggsave("TS_Gini_Score.png", plot = p2, device="png",dpi=400,units = "cm",width = 15,height = 15)
#
#
#