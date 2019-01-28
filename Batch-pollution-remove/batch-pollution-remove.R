#This scripy is used to remove polluted genes  in DGE.
################1.load Test DGE and selcet cell with 500 more UMI
setwd("/home/ggj/Rdata/JQ-Test/Data/Test/")
temp <- list.files(pattern="*txt.gz")
name <- character()
for(i in 1:length(temp)){
  message("loading DGE")
  name[i] <- unlist(strsplit(temp[i],"_dge"))[1]
  tmpvalue<-read.table(temp[i],header=T,row.names=1)
  assign(name[i],tmpvalue)
  message(paste(name[i],"is loaded"))
}
for(i in 1:length(temp)) {
  dge<-get(name[i]) 
  colnames( dge ) <- paste0(as.character(name[i]),".",colnames(get(name[i])))         
  assign(name[i], dge)
}

name_500less <- name
name_500less <- paste(name_500less,"500less", sep="_")

for(i in 1:length(name)){
  dge <- get(name[i]);
  temp <- dge[,colSums(dge)<500 & colSums(dge)> 50]
  assign(name_500less[i], temp)
}

### bulid dge more than 500 UMI
name_500more <- name
name_500more <- paste(name_500more,"500more", sep="_")

for(i in 1:length(name_500more)){
  dge <- get(name[i]);
  temp <- dge[,colSums(dge)>=500]
  assign(name_500more[i], temp)
  message(paste(name[i],"is done"))
}
################2.rmbatch in manual work
name<-"Test"
anno<-data.frame(matrix(unlist(strsplit(colnames(Test_500more),"\\.")),ncol = 2,byrow = T)[,2] )
colnames(anno)<-"Cell_barcode"
anno$Sample<-"Test"
anno$Batch<-"Test"
anno$Cell_id<-colnames(Test_500more)
anno$Cluster_id<-"1"
anno$Ages<-"age"
anno$Development_stage<-"Adult"
anno$Method<-"Microwell-seq"
anno$Gender<-"Female"
anno$Source<-"Test"
anno$Biomaterial<-"Test"
anno$Name<-"Test"
head(anno)
dim(anno)##k the annotations
# 2835   12

#raw<-get(name)
more500<-get(paste(name,"500more",sep = "_"))
less<-get(paste(name,"500less",sep = "_"))
raw<-merge(less,more500,by="row.names",all=T)
raw<-data.frame(raw[,-1],row.names = raw$Row.names)
raw[is.na(raw)]<-0
## determine the batch cells
par(mfrow=c(1,1))
hist(colSums(less),breaks = 2000,xlim = c(0,1000))
allumi<-data.frame(umi=colSums(less))
ssa<-allumi[with(allumi,order(umi,decreasing = F)),]## check the order
ssa[1:500]   ### check the UMI:50-500, use 500 cells
abline(v= 276 )
rm(ssa)
ss<-rownames(allumi)[with(allumi,order(umi,decreasing = F))][1:500]
less<-less[,ss]


## narrow down the gene
aa<-data.frame(gene=rowSums(less))
table(aa$gene>10)
usegene<-rownames(aa)[aa$gene>10]
more500<-more500[usegene,]
less<-less[usegene,]
raw<-raw[usegene,]
background <- data.frame(var=replicate(1,n = nrow(more500)),
                         cellnum_express =rowSums(more500>0),
                         rowMean_500more =rowMeans(more500),
                         row.names = rownames(more500)
                         ,rowMeans_all=rowMeans(raw)
)
temp <- merge(background,data.frame(rowMean_less =rowMeans(less)),all.x=F, by="row.names")
background <- data.frame(temp[,-1],row.names = temp[,1])
for (m in rownames(background)){
  background[m,"var"] <- var(as.numeric(more500[m,]))
  background[m,"sd"] <- sqrt(background[m,"var"])
}

background <- background[with(background,order(-rowMean_less,-rowMean_500more,-cellnum_express, -sd)),]
background$multi<-background$rowMean_less*background$sd
background<-background[background$multi>=1  ,]

#background<-background[!grepl(x=rownames(background),pattern = "*MT-"),]
#background<-background[!grepl(x=rownames(background),pattern = "*RPS"),]
#background<-background[!grepl(x=rownames(background),pattern = "*RPL"),]

summary(background$rowMean_500more/background$rowMean_less)
summary(background$rowMeans_all/background$rowMean_less)
plot(density(summary(colSums(less["HBB",]))))



## determine the coeffficient and select the med between 2-5  
med<-median(background$rowMeans_all/background$rowMean_less)
med
med<-median(background$rowMean_500more/background$rowMean_less)
med
# 2.761137



background[,"batchValue"] <- background[,"rowMean_less"]*med#  value to delete
background$batchValue <- round(background$batchValue) 
background <- background[background$batchValue>0,]


dge_m<-get(paste(name,"500more",sep = "_"))
m <- dge_m
for (i in rownames(background)) { m[i,] <- m[i,]-background[i,"batchValue"] }
sum(dge_m)# 7201740
sum(m)#4157436
m[m<0] <- 0  # 
sum(m)#5963150
(sum(dge_m)-sum(m))/sum(dge_m)
#0.1719848
rowSums(dge_m["HBB",]> 0)
rowSums(m["HBB",]>0)

par(mfrow=c(2,1))
plot(density(summary(colSums(dge_m["HBB",]))))
plot(density(summary(colSums(m["HBB",]))))

Test_500less<-Test_500less
Test_500more<-Test_500more
Test_Anno<-anno
Test_rm.batch <- m
Test_background <- background
Test_less<-less

save(Test_rm.batch,
     Test_background,
     Test_less,
     Test_500more,
     Test_Anno,
     file = "/home/ggj/Rdata/201901/Test_500more_rmbatch.RData")
