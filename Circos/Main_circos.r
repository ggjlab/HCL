library(circlize)
#library(migest)
library(dplyr)
library(gdata)
library(RColorBrewer)

color_species = structure(c("#2E8B57", "#FF4500"), names = c("HCL","MCA"))

DF<-read.table("top_hits_SRS_75.out.result",sep="\t",head=T)


#all_regions = unique(Phe$Cluster)

all_regions = unique(c(as.character(DF$Cluster1), as.character(DF$Cluster2)))
#color_regions = structure(rev(rainbow(length(all_regions))), names = as.character(all_regions))
# color_regions = structure(c("#E41A1C","#377EB8","#4DAF4A","#FCCDE5","#B3DE69","#A65628","#6A3D9A","#1B9E77","#CAB2D6","#66A61E","#D95F02","#A6761D",
#  "#E6AB02","#7570B3"),names = as.character(all_regions))

color.list<-read.table("color.list",head=T,sep="=")
rownames(color.list)<-color.list$name
color.list<-color.list[as.character(all_regions),]
#color_regions = structure(rev(rainbow(length(all_regions))), names = as.character(all_regions))
color_regions = structure(as.character(color.list$color), names = as.character(all_regions))




df2 = data.frame(from=paste(DF$Species1,DF$Cluster1,sep="|"),to=paste(DF$Species2,DF$Cluster2,sep="|"),value=DF$Mean_AUROC)
#df3<-factor(df2)
#df2<-data.frame(Phe$Cell)
#combined = unique(data.frame(regions = Phe$Cluster, species = Phe$Species, stringsAsFactors = FALSE))
combined = unique(data.frame(regions = c(as.character(DF$Cluster1), as.character(DF$Cluster2)), 
    species = c(as.character(DF$Species1), as.character(DF$Species2)), stringsAsFactors = FALSE))
combined = combined[order(combined$species, combined$regions), ]
order = paste(combined$species, combined$regions, sep = "|")
grid.col = structure(color_regions[combined$regions], names = order)
gap = rep(1, length(order))
gap[which(!duplicated(combined$species, fromLast = TRUE))] = 5

pdf("HM-circos-new_Cluster.pdf")
circos.par(gap.degree = gap,start.degree=270)
chordDiagram(df2, order = order, 
	annotationTrack = c("grid"),
    grid.col = grid.col, directional = FALSE,
    preAllocateTracks = list(
        track.height = 0.04,
        track.margin = c(0.05, 0)
    )
)
for(species in unique(combined$species)) {
    l = combined$species == species
    sn = paste(combined$species[l], combined$regions[l], sep = "|")
    highlight.sector(sn, track.index = 1, col = color_species[species], 
        #text = species, 
        niceFacing = TRUE)
}
circos.clear()

legend("bottomleft", pch = 15, col = color_regions, 
    legend = names(color_regions), cex = 0.3)
legend("bottomright", pch = 15, col = color_species, 
    legend = names(color_species), cex = 0.6)


dev.off()




