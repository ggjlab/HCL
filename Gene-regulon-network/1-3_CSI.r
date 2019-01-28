data<-read.table("Together.aucell.filter",sep="\t")
#data<-t(data)
Mcor<-cor(t(data))
n<-length(Mcor[,1])

CSI<-matrix(nrow=n,ncol=n)
for (i in 1:n){
	for(j in 1:n){
			r1<-length(which(Mcor[i,]>Mcor[i,j]))
			c1<-length(which(Mcor[,j]>Mcor[i,j]))
			CSI[i,j]<-1-(r1+c1)/((n-1)*2)
		}
	
	}
rownames(CSI)<-rownames(data)
 colnames(CSI)<-rownames(data)
 
#tf-tf network
r<-CSI	 
lower<-as.matrix(r[lower.tri(r)])
nn<-length(r[,1])
k<-0
lower.r<-matrix(nrow=length(lower))
for(i in 1:(nn-1)){
	for ( j in (i+1):nn){
		k=k+1
		tmp<-paste(rownames(r)[i],colnames(r)[j],lower[k],sep="\t")
		lower.r[k]=tmp
		
		}
	}
write.table(lower.r,file="TF-TFNetwork.CSI.celltype.lower.out",sep="\t",quote=F,row.names=F,col.names=F)

Lower<-read.table("TF-TFNetwork.CSI.celltype.lower.out",sep="\t")
out<-Lower[which(Lower[,3]>0.7),]
write.table(out,file="TFmodule.network.out",sep="\t",quote=F,row.names=F,col.names=F)


library(pheatmap)
library(RColorBrewer)

MS<-CSI
MS[MS<0.65]=0
col<-colorRampPalette(c("#FAF9DA","#28245F"))(100)
pdf("Total_TF-TF.065.pdf",width=11,height =10)
x=pheatmap(MS,
		color=col,
		clustering_method = "ward.D2",
		cex=0.6,
		#show_rownames=FALSE,
         show_colnames = FALSE
)
dev.off()

tree_order=x$tree_row$order
tree_order_name=rownames(data)[tree_order]
write.table(tree_order_name,"TF-TF.order.names",sep="\t",quote=F,col.names=F,row.names = F)
