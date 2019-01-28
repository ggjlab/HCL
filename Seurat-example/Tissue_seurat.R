setwd("./FetalThymus2")

## Reading files(dges)
a<-read.table("./FetalThymus2_dge.txt.gz",row.names = 1,header = T)
FetalThymus2<-a

name<-"FetalThymus2"
###### bulid dge below 500 UMI and dge more than 500 UMI
FetalThymus2_500less <- FetalThymus2[,colSums(FetalThymus2)<500 & colSums(FetalThymus2)> 100]
FetalThymus2_500more<-FetalThymus2[,colSums(FetalThymus2)>=500]

FetalThymus2_Anno <- data.frame(Cell_barcode= colnames(FetalThymus2_500more),
                               Sample      = replicate("FetalThymus",n=ncol(FetalThymus2_500more)),
                               Batch       = replicate("FetalThymus2",n=ncol(FetalThymus2_500more)))
colnames(FetalThymus2_500more) <- paste("2",colnames(FetalThymus2_500more),sep = ".")                           
colnames(FetalThymus2_500more) <- paste("FetalThymus",colnames(FetalThymus2_500more),sep = "_")                           
FetalThymus2_Anno[,"Cell_id"]  <- colnames(FetalThymus2_500more)
FetalThymus2_Anno[,"Cluster_id"] = replicate("1",n=ncol(FetalThymus2_500more))
FetalThymus2_Anno$Ages<-"10W"
FetalThymus2_Anno$Development_stage<-"Fetus"
FetalThymus2_Anno$Method<-rep("Microwell-seq")
FetalThymus2_Anno$Gender<-"Male"
FetalThymus2_Anno$Source<-rep("HCL")
FetalThymus2_Anno$Biomaterial<-rep("FetalThymus")
FetalThymus2_Anno$Name<-rep("FetalThymus2_10W")



##  make background
name<-"FetalThymus2"
name_background <- paste(name,"background", sep="_")
name_500more  <-  paste(name,"500more", sep="_")
name_500less  <-  paste(name,"500less", sep="_")

## check the data condition
par(mfrow=c(2,1))
hist(colSums(FetalThymus2_500more),breaks = 200)
hist(colSums(FetalThymus2_500more>0),breaks = 200)
abline(v=300)
summary(colSums(FetalThymus2_500more>0))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#71.0   386.0   440.0   473.8   530.0  1561.0 

library(Seurat)
seqwell <- CreateSeuratObject(raw.data = Matrix(as.matrix(FetalThymus2_500more),sparse=T)
                              ,min.cells = 3,min.genes = 300,names.delim = "\\.")  ##no normarlize

dim(seqwell@data)
#19211  9801
mito.genes <- grep(pattern = "^MT-", x = rownames(x = seqwell@data), value = TRUE)
percent.mito <- colSums(seqwell @raw.data[mito.genes, ])/colSums(seqwell @raw.data)
seqwell <- AddMetaData(object = seqwell, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = seqwell , features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = seqwell , gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seqwell, gene1 = "nUMI", gene2 = "nGene")
seqwell<- FilterCells(object = seqwell , subset.names = c("nGene", "percent.mito"), 
                      low.thresholds = c(300, -Inf), high.thresholds = c(2500, 0.2))
seqwell <- NormalizeData(object = seqwell, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
seqwell<- ScaleData(seqwell,vars.to.regress=c("nUMI", "percent.mito"), do.par = TRUE, num.cores =8)
par(mfrow=c(1,1))
seqwell<- FindVariableGenes(object = seqwell, mean.function = ExpMean, dispersion.function = LogVMR 
                            ,x.low.cutoff = 0.01, 
                            x.high.cutoff = 6, y.cutoff = 0.5)
length(seqwell @var.genes)# 1893
#hv.genes <- head(rownames(seqwell@hvg.info), 2000)

pbmc<-seqwell
rm(seqwell)

var.gene<-pbmc@var.genes
var.gene<-var.gene[!grepl(pattern = "*RPS",x=var.gene)]
var.gene<-var.gene[!grepl(pattern = "*RPL",x=var.gene)]
var.gene<-var.gene[!grepl(pattern = "*MT",x=var.gene)]
length(var.gene)
#2076

# Perform linear dimensional reduction
pbmc <- RunPCA(object = pbmc, pc.genes = var.gene, pcs.compute = 50, do.print = TRUE, 
               pcs.print = 1:5, genes.print = 5)
# Determine statistically significant principal components
pbmc <- JackStraw(object = pbmc, num.replicate = 100, num.pc = 40, num.cores = 8,do.par = TRUE)
# The JackStrawPlot function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line). In this case it appears that PCs 1-10 are significant.
JackStrawPlot(object = pbmc, PCs = 1:40)#25
# A more ad hoc method for determining which PCs to use is to look at a plot of the standard deviations of the principle components and draw your cutoff where there is a clear elbow in the graph. This can be done with PCElbowPlot. In this example, it looks like the elbow would fall around PC 9.
PCElbowPlot(object = pbmc,num.pc = 50)#14
PCHeatmap(object = pbmc, pc.use = 1:15, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = pbmc, pc.use = 16:30, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = pbmc, pc.use = 31:50, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)


# Run Non-linear dimensional reduction (tSNE)
# Seurat continues to use tSNE as a powerful tool to visualize and explore these datasets. While we no longer advise clustering directly on tSNE components, cells within the graph-based clusters determined above should co-localize on the tSNE plot. This is because the tSNE aims to place cells with similar local neighborhoods in high-dimensional space together in low-dimensional space. As input to the tSNE, we suggest using the same PCs as input to the clustering analysis, although computing the tSNE based on scaled gene expression is also supported using the genes.use argument.
pbmc <- RunTSNE(object = pbmc, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = pbmc,do.label = T, pt.size = 1,label.size = 5)

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:15, save.SNN = TRUE, 
                     resolution =c(0.6,0.8,1,1.4,2,2.5,4),force.recalc = T,k.param=15)

pbmc <- RunTSNE(object = pbmc, reduction.use = "pca", dims.use = 1:12, tsne.method = "FIt-SNE", 
                nthreads = 8, reduction.name = "FItSNE", reduction.key = "FItSNE_", 
                fast_tsne_path = "/home/ggj/Documents/tools/FIt-SNE-master/bin/fast_tsne", 
                max_iter = 2000,perplexity=100)
DimPlot(object = pbmc, reduction.use = "FItSNE", no.legend = F, do.return = TRUE, 
        pt.size = 1,group.by = "res.0.6",do.label = T)+ggtitle("res.0.6") #17
DimPlot(object = pbmc, reduction.use = "FItSNE", no.legend = F, do.return = TRUE, 
        pt.size = 1,group.by = "res.0.8",do.label = T) +ggtitle("res.0.8")#17
DimPlot(object = pbmc, reduction.use = "FItSNE", no.legend = F, do.return = TRUE, 
        pt.size = 1,group.by = "res.1",do.label = T) +ggtitle("res.1")#19
DimPlot(object = pbmc, reduction.use = "FItSNE", no.legend = F, do.return = TRUE, 
        pt.size = 1,group.by = "res.1.4",do.label = T)+ggtitle("res.1.4") #19
DimPlot(object = pbmc, reduction.use = "FItSNE", no.legend = F, do.return = TRUE, 
        pt.size = 1,group.by = "res.2",do.label = T)+ggtitle("res.2") #23
DimPlot(object = pbmc, reduction.use = "FItSNE", no.legend = F, do.return = TRUE, 
        pt.size = 1,group.by = "res.2.5",do.label = T) +ggtitle("res.2.5")#27
DimPlot(object = pbmc, reduction.use = "FItSNE", no.legend = F, do.return = TRUE, 
        pt.size = 1,group.by = "res.4",do.label = T)+ggtitle("res.4") #32
DimPlot(object = pbmc, reduction.use = "FItSNE", no.legend = F, do.return = TRUE, 
        pt.size = 1,do.label = T) 




pbmc<-SetAllIdent(pbmc,id="res.0.6")
aa<-FindMarkers(pbmc,0,1)# merge 

pbmc<-SetAllIdent(pbmc,id="res.0.8")
aa<-FindMarkers(pbmc,2,1)# merge 


pbmc<-SetAllIdent(pbmc,id="res.1.4")

current.cluster.ids <- 0:13
new.cluster.ids <-c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)
new.cluster.ids <-c(1,1,1,1,1,1,1,1,1,2,3,4,5,6)
pbmc@ident <- plyr::mapvalues(pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
table(pbmc@ident)
DimPlot(object = pbmc, reduction.use = "FItSNE", no.legend = F, do.return = TRUE, 
        pt.size = 1,do.label = T) 


save.image("./FetalThymus2.RData")
pbmc.markers<-FindAllMarkers(pbmc, only.pos = TRUE,  thresh.use = 0.25,min.pct = 0.15) 
pbmc.markers <- pbmc.markers[with(pbmc.markers, order(cluster,-avg_logFC, p_val_adj)),]
table(pbmc.markers$cluster)
library(gdata)
WriteXLS::WriteXLS(pbmc.markers,"./markers.xlsx")
library(dplyr)
pbmc.markers %>% group_by(cluster) %>% top_n(20, avg_logFC) ->top20
DoHeatmap(pbmc, genes.use = top20$gene,  slim.col.label = TRUE, remove.key = TRUE,cex.row =3,group.cex=5)
save.image("./FetalThymus2.RData")
save(pbmc,pbmc.markers,FetalThymus2_Anno,file = "./FetalThymus2_pbmc.RData")

