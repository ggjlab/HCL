load("/home/jingjingw/Jingjingw/Project/2018-MH-new/Pseudocell/FetalStomach1_500more.RData")
name<-"FetalStomach1"
outfile1<-"Human_FetalStomach1_pseudocell20.Rds"
outfile2<-"Human_FetalStomach1_pseudocell20.pheno.out"



Inter<-get(paste(name,"pbmc",sep = "_"))
Inter[Inter<0]=0
idd<-get(paste(name,"Anno1",sep = "_"))
Inter.id<-cbind(rownames(idd),idd$Cluster_id)

rownames(Inter.id)<-rownames(idd)
colnames(Inter.id)<-c("CellID","Celltype")
Inter.id<-as.data.frame(Inter.id)
Inter1<-Inter[,Inter.id$CellID]
Inter<-as.matrix(Inter1)
pseudocell.size = 20 ## 10 test
new_ids_list = list()
for (i in 1:length(levels(Inter.id$Celltype))) {
	cluster_id = levels(Inter.id$Celltype)[i]
	cluster_cells <- rownames(Inter.id[Inter.id$Celltype == cluster_id,])
	cluster_size <- length(cluster_cells)		
	pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
	pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
	names(pseudo_ids) <- sample(cluster_cells)	
	new_ids_list[[i]] <- pseudo_ids		
	}
	
new_ids <- unlist(new_ids_list)
new_ids <- as.data.frame(new_ids)
new_ids_length <- table(new_ids)

new_colnames <- rownames(new_ids)  ###add
all.data<-Inter[,as.character(new_colnames)] ###add
all.data <- t(all.data)###add

new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
	list(name=new_ids[,1]),FUN=mean)
rownames(new.data)<-new.data$name
new.data<-new.data[,-1]

new_ids_length<-as.matrix(new_ids_length)##
short<-which(new_ids_length<10)##
new_good_ids<-as.matrix(new_ids_length[-short,])##
result<-t(new.data)[,rownames(new_good_ids)]
colnames(result)<-paste("Human",colnames(result),sep="")
rownames(result)<-rownames(Inter)
#saveRDS(result,file=outdir1[i]) ###
saveRDS(result,file=outfile1) ###
cellty<-gsub("[_]Cell[0-9]|[_]Cell[0-9][0-9]|[_]Cell[0-9][0-9][0-9]|[_]Cell[0-9][0-9][0-9][0-9]|[_]Cell[0-9][0-9][0-9][0-9][0-9]","",colnames(result))
new.phe<-paste(colnames(result),'HumanFetal',cellty,sep="\t")

#write.table(new.phe,file=outdir2[i],quote=F,row.names=F) ###

write.table(new.phe,file=outfile2,quote=F,row.names=F) ###




