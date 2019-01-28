data<-read.csv("aucell.csv")
rownames(data)<-data$Cell
data<-data[,-1]

#n<-gsub("[....][1-9]|[...]","",colnames(data))
#n<-as.matrix(n)
#data<-t(data)
#all.data <- data###add
#new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
#	list(name=n[,1]),FUN=mean)
#rownames(new.data)<-new.data$name
#new.data<-new.data[,-1]
#data1<-new.data#normalize

data1<-t(data)
coln<-colnames(data1)
coln1<-gsub("_Cell[0-9]|_Cell[0-9][0-9]|_Cell[0-9][0-9][0-9]","",coln)
coln1<-as.factor(coln1)
le<-levels(coln1)
nle<-length(le)
coln1<-as.matrix(coln1)
Result<-matrix(0,ncol=length(coln1),nrow=nle)
for (i in 1:nle){
	tmp<-which(coln1[,1] == as.character(le[i]))
	Result[i,tmp]<-1
}
colnames(Result)<-colnames(data1)
rownames(Result)<-le

Result1<-Result/(rowSums(Result))
#Result1[Result1==0]=0.000001

#write.table(Result1,file="JSD.mouse.celltype.input",sep="\t",quote=F)


#-------------------JSD
KLD=function(A,B){
	sum(A*log(A/B))
}
JSD=function(P,Q){
	M=(P+Q)/2
	jsd=0.5*KLD(P,M)+0.5*KLD(Q,M)
	return (jsd)
}

Input1<-data1
Input2<-Result1
Input1[Input1==0]<-0.0000001 #attention
Input1<-Input1/rowSums(Input1)
Input2[Input2==0]<-0.0000001#attention
Input2<-Input2/rowSums(Input2)
TFn<-length(Input1[,1])
Celln<-length(Input1[1,])
Celltypn<-length(Input2[,1])
JSD.result<-matrix(0,nrow=TFn,ncol=Celltypn)
for (i in 1:TFn){
	for(j in 1:Celltypn){
		jsd1<-JSD(Input1[i,],Input2[j,])
		JSD.result[i,j]=jsd1	
	}

}
rownames(JSD.result)<-rownames(Input1)
colnames(JSD.result)<-rownames(Input2)

JSDR<-1-sqrt(JSD.result)
write.table(JSDR,file="Human.RSS.total.out",sep="\t",quote=F)
