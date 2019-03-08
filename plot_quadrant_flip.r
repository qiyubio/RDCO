## Author: Qi Yu
## Contact: dryuqi@gmail.com

library("lattice")
library("reshape2")
library(RColorBrewer)
library(gridExtra)
library(reshape2)


col.l<-heat.colors(100)[length(heat.colors(100)):1]

defplot2flip<-function(file,chr=NA,name=NA){
  
  matrix1<-read.table(file,header=FALSE)
  names(matrix1)<-c("x","y","z")
  matrix1$y<-max(matrix1$y)-matrix1$y
  #return(normalized_count)
  #my.at <- seq(20, 55, 85,140)
  png(paste(chr,"_",name,".four_chr1_chr3.png",sep=""),width=1000, height=1000, res=200)
  #b <- c(["B6_chr1","CAST_chr1","B6_chr3","CAST_chr3",seq(133,139))
  myplot1<-levelplot(log(z) ~ x*y, data=matrix1, xlab=" ",ylab=" ", col.regions = col.l, colorkey=list(space="bottom"),scales=list(x=list(alternatin=2, at=c(20, 60, 95, 130),labels=c("B6_chr1","CAST_chr1","B6_chr3","CAST_chr3")),y=list(at=c(15, 50, 85,125),labels=rev(c("B6_chr1","CAST_chr1","B6_chr3","CAST_chr3")))), main=paste(name,chr))
  print(myplot1)
  
  dev.off()
  
}


defplot2flip("CRD018512_bowtie2_mask_DBA2J.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_5000000.matrix",name="E7.5",chr="chr1_chr3")


