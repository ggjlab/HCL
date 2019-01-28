#load reference
load("/home/ggj/Rdata/201810/NewReference/ref-exp.RData")
#start scHCL
Test_dge<-read.table("./_dge_sample.csv",sep = ",",header = T,row.names = 1)#load the matirx of new cells,row is genename ,col is cellname
#log-normalized the uploaded dge
Test_dge<-as.matrix(t(t(Test_dge)/colSums(Test_dge))*100000) 
Test_dge<-log(Test_dge+1)
tst <- data.frame(matrix(nrow =length(ref[,1]),ncol = length(Test_dge[1,])))#the 104 is the test cell numbers 
rownames(tst)<-rownames(ref)
colnames(tst)<-colnames(Test_dge)
for (i in rownames(ref)) {ref
  if(i%in%rownames(Test_dge)) tst[i,]<- Test_dge[i,]
}
tst[is.na(tst)]<-0
ref<-log(ref_exp+1)
cors <- cor(ref,tst)
cors[is.na(cors)]<-0
cor1<-cors
cor1m<-apply(cor1,2,max) #
cor1S<-apply(cors,2,function(x) rownames(cors)[which.max(x)])
cor1r<-cbind(cor1S,cor1m)
#scHCL results is in scHCL dataframe
scHCL<-data.frame(cor1r)
colnames(scHCL)<-c("scHCL_result" ,  "cors_log")
