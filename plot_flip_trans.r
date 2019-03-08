## Author: Qi Yu
## Contact: dryuqi@gmail.com

library("lattice")
library("reshape2")
library(RColorBrewer)
library(gridExtra)
library(reshape2)


defplot2flip_trans<-function(file,chr=NA,name=NA){
  
  matrix1<-read.table(file,header=FALSE)
  names(matrix1)<-c("x","y","z")
  matrix1$y<-max(matrix1$y)-matrix1$y
  #return(normalized_count)
  #my.at <- seq(20, 55, 85,140)
  png(paste(chr,"_",name,".four_chr1_chr3_trans.png",sep=""),width=1000, height=1000, res=200)
  #b <- c(["B6_chr1","CAST_chr1","B6_chr3","CAST_chr3",seq(133,139))
  myplot1<-levelplot(z ~ x*y, data=matrix1, xlab=" ",ylab=" ", col.regions = col.l, colorkey=list(space="bottom"),scales=list(x=list(alternatin=2, at=c(20, 60, 95, 130),labels=c("B6_chr1","CAST_chr1","B6_chr3","CAST_chr3")),y=list(at=c(15, 50, 85,125),labels=rev(c("B6_chr1","CAST_chr1","B6_chr3","CAST_chr3")))), main=paste(name,chr))
  print(myplot1)
  dev.off()
  
}

defplot2flip_trans("181016_9014_3_CTTGTA_GC_050218_Hi_C_B6XCAST_Preleptotene.mm10_genome_mask_CAST.bwt2pairs_VI_chr1_chr3_5000000_trans.matrix",name="preleptotene_highres",chr="chr1_chr3")
defplot2flip_trans("180530_lep_B6XCAST_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_5000000_trans.matrix",name="leptotene",chr="chr1_chr3")
defplot2flip_trans("180901_9011_1_ACTGAT_GC_072318_Hi-C_B6XCAST_Zygotene.mm10_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_5000000_trans.matrix",name="zygotene",chr="chr1_chr3")
defplot2flip_trans("181016_9014_4_ACTTGA_GC_072318_Hi-C_B6XCAST_Pachytene.mm10_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_5000000_trans.matrix",name="pachytene_highres",chr="chr1_chr3")
defplot2flip_trans("180903_9013_1_GAGTGG_GC_072318_Hi-C_B6XCAST_Diplotene.mm10_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_5000000_trans.matrix",name="diplotene",chr="chr1_chr3")
defplot2flip_trans("180901_9011_1_TTAGGC_GC_072318_Hi-C_B6XCAST_2N.mm10_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_5000000_trans.matrix",name="2N",chr="chr1_chr3")
#########################################################
col.l<-heat.colors(100)[length(heat.colors(100)):1]
chr<-"chr1_chr3"
name<-"E7.5"
defplot3flip_trans<-function(file,chr=NA,name=NA){
  
  matrix1<-read.table(file,header=FALSE)
  names(matrix1)<-c("x","y","z")
  matrix1$y<-max(matrix1$y)-matrix1$y
  #return(normalized_count)
  #my.at <- seq(20, 55, 85,140)
  png(paste(chr,"_",name,".four_chr1_chr3_trans.png",sep=""),width=1000, height=1000, res=200)
  #b <- c(["B6_chr1","CAST_chr1","B6_chr3","CAST_chr3",seq(133,139))
  myplot1<-levelplot(z ~ x*y, data=matrix1, xlab=" ",ylab=" ", col.regions = col.l, colorkey=list(space="bottom"),scales=list(x=list(alternatin=2, at=c(20, 60, 95, 130),labels=c("B6_chr1","DBA2J_chr1","B6_chr3","DBA2J_chr3")),y=list(at=c(15, 50, 85,125),labels=rev(c("B6_chr1","DBA2J_chr1","B6_chr3","DBA2J_chr3")))), main=paste(name,chr))
  print(myplot1)
  dev.off()
  
}
defplot3flip_trans("CRD018512_bowtie2_mask_DBA2J.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_5000000_trans.matrix",name="E7.5",chr="chr1_chr3")

defplot3flip<-function(file,chr=NA,name=NA){
  
  matrix1<-read.table(file,header=FALSE)
  names(matrix1)<-c("x","y","z")
  matrix1$y<-max(matrix1$y)-matrix1$y
  #return(normalized_count)
  #my.at <- seq(20, 55, 85,140)
  png(paste(chr,"_",name,".four_chr1_chr3.png",sep=""),width=1000, height=1000, res=200)
  #b <- c(["B6_chr1","CAST_chr1","B6_chr3","CAST_chr3",seq(133,139))
  myplot1<-levelplot(log(z) ~ x*y, data=matrix1, xlab=" ",ylab=" ", col.regions = col.l, colorkey=list(space="bottom"),scales=list(x=list(alternatin=2, at=c(20, 60, 95, 130),labels=c("B6_chr1","DBA2J_chr1","B6_chr3","DBA2J_chr3")),y=list(at=c(15, 50, 85,125),labels=rev(c("B6_chr1","DBA2J_chr1","B6_chr3","DBA2J_chr3")))), main=paste(name,chr))
  print(myplot1)
  dev.off()
  
}

defplot3flip("CRD018512_bowtie2_mask_DBA2J.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_5000000.matrix",name="E7.5",chr="chr1_chr3")

chr<-"chr1_chr3"
name<-"ES cell line"
defplot4flip_trans<-function(file,chr=NA,name=NA){
  
  matrix1<-read.table(file,header=FALSE)
  names(matrix1)<-c("x","y","z")
  matrix1$y<-max(matrix1$y)-matrix1$y
  #return(normalized_count)
  #my.at <- seq(20, 55, 85,140)
  png(paste(chr,"_",name,".four_chr1_chr3_trans.png",sep=""),width=1000, height=1000, res=200)
  #b <- c(["B6_chr1","CAST_chr1","B6_chr3","CAST_chr3",seq(133,139))
  myplot1<-levelplot(z ~ x*y, data=matrix1, xlab=" ",ylab=" ", col.regions = col.l, colorkey=list(space="bottom"),scales=list(x=list(alternatin=2, at=c(20, 60, 95, 130),labels=c("CAST_chr1","129s_chr1","CAST_chr3","129s_chr3")),y=list(at=c(15, 50, 85,125),labels=rev(c("CAST_chr1","129s_chr1","CAST_chr3","129s_chr3")))), main=paste(name,chr))
  print(myplot1)
  dev.off()
  
}

defplot4flip_trans("SRR2240730_bowtie2_mask_CAST_129s.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_5000000_trans.matrix",name="ES",chr="chr1_chr3")
#defplot2flip_trans("CRD018512_bowtie2_mask_DBA2J.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_5000000.matrix",name="E7.5",chr="chr1_chr3")

defplot4flip<-function(file,chr=NA,name=NA){
  
  matrix1<-read.table(file,header=FALSE)
  names(matrix1)<-c("x","y","z")
  matrix1$y<-max(matrix1$y)-matrix1$y
  #return(normalized_count)
  #my.at <- seq(20, 55, 85,140)
  png(paste(chr,"_",name,".four_chr1_chr3.png",sep=""),width=1000, height=1000, res=200)
  #b <- c(["B6_chr1","CAST_chr1","B6_chr3","CAST_chr3",seq(133,139))
  myplot1<-levelplot(log(z) ~ x*y, data=matrix1, xlab=" ",ylab=" ", col.regions = col.l, colorkey=list(space="bottom"),scales=list(x=list(alternatin=2, at=c(20, 60, 95, 130),labels=c("CAST_chr1","129s_chr1","CAST_chr3","129s_chr3")),y=list(at=c(15, 50, 85,125),labels=rev(c("CAST_chr1","129s_chr1","CAST_chr3","129s_chr3")))), main=paste(name,chr))
  print(myplot1)
  dev.off()
  
}
defplot4flip("SRR2240730_bowtie2_mask_CAST_129s.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_5000000.matrix",name="ES",chr="chr1_chr3")


