
human<-readRDS("./HCL_all_v2_pse20.CPM.rds")

mouse<-readRDS("./MCA_V2.PSUDOCELL20.rds")

orth<-read.table("./Human_Mouse_one-one.orth",sep="\t")
#orth<-as.matrix(orth)
mouse<-as.data.frame(mouse)
mouse.orth<-mouse[as.character(orth[,4]),]
human<-as.data.frame(human)
human.orth<-human[as.character(orth[,2]),]


data<-cbind(mouse.orth,human.orth)
rownames(data)<-orth[,1]
data[is.na(data)]<-0

P<-read.table("./MCA_V2.psudocell20.phe",sep="\t",head=T)
P2<-read.table("./HCL_v2.pse20.phe",sep="\t",head=T)
P1<-rbind(P,P2)
colnames(P1)<-c("Sample_ID","Study_ID","Celltype")
data1<-data[,as.character(P1$Sample_ID)]


source("2017-08-28-runMN-US.R")
#library(gplots)
#library(RColorBrewer)

celltypes1 <-unique(as.character(P1$Celltype))


var.genes1=get_variable_genes(data1,P1)
length(var.genes1)
write.table(var.genes1,"var.genes_75.out",sep="\t",quote=F)#####--------
celltype.NV=run_MetaNeighbor_US(var.genes1,data1,celltypes1,P1)
write.table(celltype.NV,file="celltype.NV_SRS_75.out",sep="\t",quote=F)###---------

cols=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100))
breaks=seq(0,1,length=101)
pdf("celltype.NV_SRS_75.pdf")  #########--------------------------
heatmap.2(celltype.NV,trace="none",density.info="none",col=cols,breaks=breaks,cexRow=0.3,cexCol=0.3)
dev.off()
top_hits=get_top_hits(celltype.NV,P1,threshold=0.9,filename="top_hits_SRS_75.out") 
top_hits=get_top_hits(celltype.NV,P1,threshold=0.8,filename="top_hits_SRS_0.8_75.out")
top_hits=get_top_hits(celltype.NV,P1,threshold=0.7,filename="top_hits_SRS_0.7_75.out")
top_hits=get_top_hits(celltype.NV,P1,threshold=0.6,filename="top_hits_SRS_0.6_75.out")













