
## Author: Qi Yu
## Contact: dryuqi@gmail.com

library("lattice")
library("reshape2")
library(RColorBrewer)
library(gridExtra)
library(reshape2)

#Plot heat map  
defplot<-function(file,name=NA,chr=NA){
matrix1<-read.table(file,header=FALSE)
names(matrix1)<-c("x","y","z")
png(paste(chr,"_",name,"log.png",sep=""),width=900, height=900, res=200)
myplot<-levelplot(log(z) ~ x*y, data=matrix1, xlab="X" , col.regions = heat.colors(100)[length(heat.colors(100)):1] , main=paste(name,chr))
print(myplot)
dev.off()
return(matrix1)
}


pachy<-defplot("pachytene_oe_chr2_55_70_100kb",name="Pachytene",chr="chr2")

#secondway

name="lep"

data<- as.matrix(read.table("inter_30_Lep.hic_chr12_70_80_100kb_matrix.txt", header=FALSE, sep = "\t",
                            #row.names = 1,
                            as.is=TRUE))


png(paste(chr,name,"_matrix100kb.png",sep=""),width=1200, height=1200, res=200)
image(x=1:100, y=1:101, z=data,
      #zlim=c(-1,1),
      col=colorRampPalette(c("red","yellow","white"))(256), useRaster=TRUE)

dev.off()

#calculate pearson correaltion and plot

cm <- cor(data)

png(paste(chr,name,"_corr100kb.png",sep=""),width=1200, height=1200, res=200)
image(x=1:101, y=1:101, z=cm,
      zlim=c(-1,1),
      col=colorRampPalette(c("red", "white","blue"))(256), useRaster=TRUE)
dev.off()

#subtract first compartment

png(paste(chr,name,"_ppc1_100kb.png",sep=""),width=1200, height=800, res=200)
cm[is.na(cm)] <- 0
princp = princomp(cm)
plot(princp$loadings[,1],type='l')
dev.off()

#subtract second compartment

png(paste(chr,name,"_pc2_100kb.png",sep=""),width=1200, height=800, res=200)

plot(princp$loadings[,2],type='l')
dev.off()

#######################################################################################


