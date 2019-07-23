#############
library(Seurat)
setwd("/home/ggj/github/HCL/HCL/Scenic_R/example_data/")
#load your DGE
exprMat <- readRDS("./HCL_v2.pse20.SRS_17000.rds")
gene<-data.frame(colSums(exprMat>0))
#load the annotation for each cell
ident<-read.table("./HCL_v2.pse20.SRS_17000.phe",header =T,row.names = 1)
dir.create("SCENIC")
setwd("SCENIC/")
cellInfo <-merge(ident,gene,by="row.names",all=T)
cellInfo<-data.frame(cellInfo[,-1],row.names = cellInfo$Row.names)
cellInfo[is.na(cellInfo)]<-0
cellInfo<-cellInfo[,-1]
colnames(cellInfo)<- c('CellType','nGene')
cellInfo$cluster<-gsub("Human","",cellInfo$CellType)
dim(exprMat)
head(cellInfo)
setwd("/home/ggj/github/HCL/HCL/Scenic_R/example_data/SCENIC/")
dir.create("int")
saveRDS(cellInfo, file="./int/cellInfo.Rds")
table(cellInfo$cluster)
colVars <- list(CellType=c('Human1'='#ffff00',
                           'Human2'='#1ce6ff',
                           'Human3'='#ff34ff',
                           'Human4'='#ff4a46',
                           'Human5'='#008941',
                           'Human6'='#006fa6',
                           'Human7'='#a30059',
                           'Human8'='#ffdbe5',
                           'Human9'='#7a4900',
                           'Human10'='#0000a6',
                           'Human11'='#63ffac',
                           'Human12'='#b79762',
                           'Human13'='#004d43',
                           'Human14'='#8fb0ff',
                           'Human15'='#997d87',
                           'Human16'='#5a0007',
                           'Human17'='#809693',
                           'Human18'='#feffe6',
                           'Human19'='#1b4400',
                           'Human20'='#4fc601',
                           'Human21'='#3b5dff',
                           'Human22'='#4a3b53',
                           'Human23'='#ff2f80',
                           'Human24'='#61615a',
                           'Human25'='#ba0900',
                           'Human26'='#6b7900',
                           'Human27'='#00c2a0',
                           'Human28'='#ffaa92',
                           'Human29'='#ff90c9',
                           'Human30'='#b903aa',
                           'Human31'='#d16100',
                           'Human32'='#ddefff',
                           'Human33'='#000035',
                           'Human34'='#7b4f4b',
                           'Human35'='#a1c299',
                           'Human36'='#300018',
                           'Human37'='#0aa6d8',
                           'Human38'='#013349',
                           'Human39'='#00846f',
                           'Human40'='#372101',
                           'Human41'='#ffb500',
                           'Human42'='#c2ffed',
                           'Human43'='#a079bf',
                           'Human44'='#cc0744',
                           'Human45'='#c0b9b2',
                           'Human46'='#c2ff99',
                           'Human47'='#001e09',
                           'Human48'='#00489c',
                           'Human49'='#6f0062',
                           'Human50'='#0cbd66',
                           'Human51'='#eec3ff',
                           'Human52'='#456d75',
                           'Human53'='#b77b68',
                           'Human54'='#7a87a1',
                           'Human55'='#788d66',
                           'Human56'='#885578',
                           'Human57'='#fad09f',
                           'Human58'='#ff8a9a',
                           'Human59'='#d157a0',
                           'Human60'='#bec459',
                           'Human61'='#456648',
                           'Human62'='#0086ed',
                           'Human63'='#886f4c',
                           'Human64'='#34362d',
                           'Human65'='#b4a8bd',
                           'Human66'='#00a6aa',
                           'Human67'='#452c2c',
                           'Human68'='#636375',
                           'Human69'='#a3c8c9',
                           'Human70'='#ff913f',
                           'Human71'='#938a81',
                           'Human72'='#575329',
                           'Human73'='#00fecf',
                           'Human74'='#b05b6f',
                           'Human75'='#8cd0ff',
                           'Human76'='#3b9700',
                           'Human77'='#04f757',
                           'Human78'='#c8a1a1',
                           'Human79'='#1e6e00',
                           'Human80'='#7900d7',
                           'Human81'='#a77500',
                           'Human82'='#6367a9',
                           'Human83'='#a05837',
                           'Human84'='#6b002c',
                           'Human85'='#772600',
                           'Human86'='#d790ff',
                           'Human87'='#9b9700',
                           'Human88'='#549e79',
                           'Human89'='#fff69f',
                           'Human90'='#201625',
                           'Human91'='#72418f',
                           'Human92'='#bc23ff',
                           'Human93'='#99adc0',
                           'Human94'='#3a2465',
                           'Human95'='#922329',
                           'Human96'='#5b4534',
                           'Human97'='#fde8dc',
                           'Human98'='#404e55',
                           'Human99'='#0089a3',
                           'Human100'='#cb7e98',
                           'Human101'='#a4e804',
                           'Human102'='#324e72'))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="./int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))
library(SCENIC)
org="hgnc" # or hgnc, or dmel
dbDir="/home/ggj/Rdata/201906/Human/databases/" # RcisTarget databases location
myDatasetTitle="SCENIC example on Human1.1" # choose a name for your analysis
dbs <- c("hg19-500bp-upstream-10species.mc9nr.feather","hg19-tss-centered-5kb-10species.mc9nr.feather")
names(dbs)<-c("500bp","5kb")
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 

setwd("/home/ggj/github/HCL/HCL/Scenic_R/example_data/SCENIC/")
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
saveRDS(scenicOptions, file="./int/scenicOptions.Rds")



###
### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
##1 load AUCell matrix from PyScenic
setwd("/home/ggj/github/HCL/HCL/Scenic_R/example_data")
regulonAUC<-importAUCfromText("./aucell.csv")
regulonAUC
dim(regulonAUC)
saveRDS(regulonAUC, file="./int/3.4_regulonAUC.Rds")

##2 run SCENIC tSNE
setwd("/home/ggj/github/HCL/HCL/Scenic_R/example_data/SCENIC/")
nPcs <- c(50) 
scenicOptions@settings$seed <- 123 # same seed for all of them
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE,filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5)
# Using only "high-confidence" regulons (normally similar)
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 50
scenicOptions@settings$defaultTsne$perpl <- 50
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
scenicOptions<-readRDS("./int/scenicOptions.Rds")
# Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat) # default t-SNE
savedSelections <- shiny::runApp(aucellApp)
# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
auc<-read.csv("/home/ggj/github/HCL/HCL/Scenic_R/example_data/aucell.csv",row.names = 1)
usethrethold<-apply(auc,2,summary)
colnames(usethrethold)<-gsub("\\...","(+)",colnames(usethrethold))
dim(usethrethold)
A0.5<-0.5*usethrethold[6,]
names(A0.5)
auc[1:5,1:5]
setdiff(colnames(usethrethold),rownames(regulonAUC))
setdiff(rownames(regulonAUC),colnames(usethrethold))
#runSCENIC aucell binarize
setwd("/home/ggj/github/HCL/HCL/Scenic_R/example_data/SCENIC/")
#0.5*max
newThresholds<-A0.5
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
runSCENIC_4_aucell_binarize(scenicOptions)





