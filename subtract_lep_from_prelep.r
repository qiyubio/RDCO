## Author: Qi Yu
## Contact: dryuqi@gmail.com

library("lattice")
library("reshape2")
library(RColorBrewer)
library(gridExtra)
library(reshape2)


col.l<-colorRampPalette(brewer.pal(100,'YlGnBu'))(100)
m1<-0.1
m2<-11
ss<-4
ll<-12

defplot2<-function(file,chr=NA,name=NA){
  
  matrix1<-read.table(file,header=FALSE)
  names(matrix1)<-c("x","y","z")
  #matrix1$x<-(matrix1$x-min(matrix1$x)+1)x5000
  center<-min(matrix1$x)+round((max(matrix1$x)-min(matrix1$x))/2)
 
  png(paste(chr,"_",name,".subset.png",sep=""),width=1600, height=1500, res=200)
  matrix2<-matrix1[matrix1[,"x"] <= (min(matrix1$x)+ll) & matrix1[,"x"] > (min(matrix1$x)+ss) & matrix1[,"y"] <= (min(matrix1$x)+ll) & matrix1[,"y"] > (min(matrix1$x)+ss) ,]
  myplot2<-levelplot(log(z) ~ x*y, data=matrix2, xlab="X" , at=seq(m1, m2, length.out=10),col.regions =col.l, main=paste(name,chr))
  
  matrix3<-matrix1[matrix1[,"x"] < (center+ll) & matrix1[,"x"] >= (center+ss) & matrix1[,"y"] < (center+ll) & matrix1[,"y"] >= (center+ss) ,]
  myplot3<-levelplot(log(z) ~ x*y, data=matrix3, xlab="X" , at=seq(m1, m2, length.out=10),col.regions =col.l, main=paste(name,chr))
  
  matrix4<-matrix1[matrix1[,"x"] < (min(matrix1$x)+ll) & matrix1[,"x"] >= (min(matrix1$x)+ss) & matrix1[,"y"] <= (center+ll) & matrix1[,"y"] > (center+ss) ,]
  myplot4<-levelplot(log(z) ~ x*y, data=matrix4, xlab="X" , at=seq(m1, m2, length.out=10),col.regions =col.l, main=paste(name,chr))
  
  matrix5<-matrix1[matrix1[,"x"] <= (center+ll) & matrix1[,"x"] > (center+ss) & matrix1[,"y"] < (min(matrix1$x)+ll) & matrix1[,"y"] >= (min(matrix1$x)+ss) ,]
  myplot5<-levelplot(log(z) ~ x*y, data=matrix5, xlab="X" , at=seq(m1, m2, length.out=10),col.regions =col.l, main=paste(name,chr))
  
  print(grid.arrange(myplot4, myplot3, myplot2,myplot5,ncol=2))
  
  
  #print(myplot1)
  
  dev.off()
  
}

defplot2("181016_9014_3_CTTGTA_GC_050218_Hi_C_B6XCAST_Preleptotene.mm10_genome_mask_CAST.bwt2pairs_VI_chr3_5000000.matrix",name="prelep_highres",chr="chr3")

defplot2("180530_prelep_B6XCAST_genome_mask_CAST.bwt2pairs_heterozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="prelep",chr="chr3")

##############compare two stages#######################################

col.l<-colorRampPalette(brewer.pal(100,'YlGnBu'))(100)
m1<-0.1
m2<-11
ss<-4
ll<-12

defplot3<-function(file,file2=NA,chr=NA,name=NA){
  
  matrix1<-read.table(file,header=FALSE)
  names(matrix1)<-c("x","y","z")
  #matrix1$x<-(matrix1$x-min(matrix1$x)+1)x5000
  center<-min(matrix1$x)+round((max(matrix1$x)-min(matrix1$x))/2)
  
  matrix12<-read.table(file2,header=FALSE)
  names(matrix12)<-c("x","y","z")
  #matrix1$x<-(matrix1$x-min(matrix1$x)+1)x5000
  center2<-min(matrix12$x)+round((max(matrix12$x)-min(matrix12$x))/2)

  matrix22<-matrix12[matrix12[,"x"] <= (min(matrix12$x)+ll) & matrix12[,"x"] > (min(matrix12$x)+ss) & matrix12[,"y"] <= (min(matrix12$x)+ll) & matrix12[,"y"] > (min(matrix12$x)+ss) ,]
   
  matrix32<-matrix12[matrix12[,"x"] < (center+ll) & matrix12[,"x"] >= (center+ss) & matrix12[,"y"] < (center+ll) & matrix12[,"y"] >= (center+ss) ,]
   
  matrix42<-matrix12[matrix12[,"x"] < (min(matrix12$x)+ll) & matrix12[,"x"] >= (min(matrix12$x)+ss) & matrix12[,"y"] <= (center+ll) & matrix12[,"y"] > (center+ss) ,]
   
  matrix52<-matrix12[matrix12[,"x"] <= (center+ll) & matrix12[,"x"] > (center+ss) & matrix12[,"y"] < (min(matrix12$x)+ll) & matrix12[,"y"] >= (min(matrix12$x)+ss) ,]
 
  png(paste(chr,"_",name,".subset_compare.png",sep=""),width=1600, height=1500, res=200)
  
  matrix21<-matrix1[matrix1[,"x"] <= (min(matrix1$x)+ll) & matrix1[,"x"] > (min(matrix1$x)+ss) & matrix1[,"y"] <= (min(matrix1$x)+ll) & matrix1[,"y"] > (min(matrix1$x)+ss) ,]
  matrix31<-matrix1[matrix1[,"x"] < (center+ll) & matrix1[,"x"] >= (center+ss) & matrix1[,"y"] < (center+ll) & matrix1[,"y"] >= (center+ss) ,]
  matrix41<-matrix1[matrix1[,"x"] < (min(matrix1$x)+ll) & matrix1[,"x"] >= (min(matrix1$x)+ss) & matrix1[,"y"] <= (center+ll) & matrix1[,"y"] > (center+ss) ,]
  matrix51<-matrix1[matrix1[,"x"] <= (center+ll) & matrix1[,"x"] > (center+ss) & matrix1[,"y"] < (min(matrix1$x)+ll) & matrix1[,"y"] >= (min(matrix1$x)+ss) ,]
  
  matrix2<- matrix21/matrix22
  matrix3<- matrix31/matrix32
  matrix4<- matrix41/matrix42
  matrix5<- matrix51/matrix52
  
  myplot2<-levelplot(log(z) ~ x*y, data=matrix2, xlab="X" , at=seq(m1, m2, length.out=10),col.regions =col.l, main=paste(name,chr))
  
  myplot3<-levelplot(log(z) ~ x*y, data=matrix3, xlab="X" , at=seq(m1, m2, length.out=10),col.regions =col.l, main=paste(name,chr))
  
  myplot4<-levelplot(log(z) ~ x*y, data=matrix4, xlab="X" , at=seq(m1, m2, length.out=10),col.regions =col.l, main=paste(name,chr))
  
  myplot5<-levelplot(log(z) ~ x*y, data=matrix5, xlab="X" , at=seq(m1, m2, length.out=10),col.regions =col.l, main=paste(name,chr))
  
  
  print(grid.arrange(myplot4, myplot3, myplot2,myplot5,ncol=2))
  
  
  #print(myplot1)
  
  dev.off()
  
}

defplot3("181016_9014_3_CTTGTA_GC_050218_Hi_C_B6XCAST_Preleptotene.mm10_genome_mask_CAST.bwt2pairs_VI_chr3_5000000.matrix",file2="180530_prelep_B6XCAST_genome_mask_CAST.bwt2pairs_heterozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="prelep_highres-prelep",chr="chr3")
