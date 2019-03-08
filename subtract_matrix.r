
## Author: Qi Yu
## Contact: dryuqi@gmail.com

library("lattice")
library("reshape2")
library(RColorBrewer)
library(gridExtra)
library(reshape2)

#Plot heat map individually 
defplot<-function(file,name=NA,chr=NA){
matrix1<-read.table(file,header=FALSE)
names(matrix1)<-c("x","y","z")
png(paste(chr,"_",name,"log.png",sep=""),width=900, height=900, res=200)
myplot<-levelplot(log(z) ~ x*y, data=matrix1, xlab="X" , col.regions = heat.colors(100)[length(heat.colors(100)):1] , main=paste(name,chr))
print(myplot)
dev.off()
return(matrix1)
}


zygo<-defplot("180901_9011_1_ACTGAT_GC_072318_Hi-C_B6XCAST_Zygotene.mm10_genome_mask_CAST.bwt2pairs_heterozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="Zygotene",chr="chr3")
prelep<-defplot("180530_prelep_B6XCAST_genome_mask_CAST.bwt2pairs_heterozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="preleptotene",chr="chr3")
lep<-defplot("180530_lep_B6XCAST_genome_mask_CAST.bwt2pairs_heterozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="leptotene",chr="chr3")
N1<-defplot("180901_9011_1_TTAGGC_GC_072318_Hi-C_B6XCAST_2N.mm10_genome_mask_CAST.bwt2pairs_heterozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="2N",chr="chr3")
diplo<-defplot("180903_9013_1_GAGTGG_GC_072318_Hi-C_B6XCAST_Diplotene.mm10_genome_mask_CAST.bwt2pairs_heterozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="Diplotene",chr="chr3")
ES<-defplot("SRR2240730_bowtie2_mask_CAST_129s.bwt2pairs_homozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="EScell",chr="chr3")
ES75<-defplot("CRD018512_bowtie2_mask_DBA2J.bwt2pairs_interaction_homozygous_samechr.bam_ncc.txt_VI_chr3_5000000.matrix",name="ES7.5",chr="chr3")
prelep_highres<-defplot("181016_9014_3_CTTGTA_GC_050218_Hi_C_B6XCAST_Preleptotene.mm10_genome_mask_CAST.bwt2pairs_VI_chr3_5000000.matrix",name="prelep_highres",chr="chr3")

k<-list(pachy,zygo,prelep,lep,diplo,N1)

kname<-list("pachy","zygo","prelep","lep","diplo","2N")
newl<-list()
col.l<-rev(colorRampPalette(brewer.pal(9,'RdYlGn'))(100))

for (i in c(1:6)){
In<-k[[i]]-ES
normalized_pn<-melt(In)
names(normalized_pn)<-c("x","y","z")
normalized_pn$z<-normalized_pn$z-min(normalized_pn$z,na.rm=TRUE)
newl[[i]]<-normalized_pn
}



#####################################################################################################

png(paste(chr,"_",kname[[i]],"minus_six_log.png",sep=""),width=1500, height=2000, res=200)

#plot to the same figure with same scale

l<-c(newl[[1]]$z,newl[[2]]$z,newl[[3]]$z,newl[[4]]$z,newl[[5]]$z,newl[[6]]$z,na.rm=TRUE)
m1<-log(min(l[l>0],na.rm=TRUE))
m2<-log(max(newl[[1]]$z,newl[[2]]$z,newl[[3]]$z,newl[[4]]$z,newl[[5]]$z,newl[[6]]$z,na.rm=TRUE))

r1<-levelplot(log(z)~ x*y, data=newl[[1]], xlab="chr3(5Mb)" , ylab="chr3(5Mb)" , at=seq(m1, m2, length.out=10), col.regions =col.l , main=paste(kname[[1]],"-ES",sep=""))
r2<-levelplot(log(z)~ x*y, data=newl[[2]], xlab="chr3(5Mb)" , ylab="chr3(5Mb)" ,at=seq(m1, m2, length.out=10), col.regions =col.l , main=paste(kname[[2]],"-ES",sep=""))
r3<-levelplot(log(z)~ x*y, data=newl[[3]], xlab="chr3(5Mb)" , ylab="chr3(5Mb)",at=seq(m1, m2, length.out=10),col.regions =col.l , main=paste(kname[[3]],"-ES",sep=""))
r4<-levelplot(log(z)~ x*y, data=newl[[4]], xlab="chr3(5Mb)", ylab="chr3(5Mb)", at=seq(m1, m2, length.out=10),col.regions =col.l , main=paste(kname[[4]],"-ES",sep=""))
r5<-levelplot(log(z)~ x*y, data=newl[[5]], xlab="chr3(5Mb)",ylab="chr3(5Mb)", at=seq(m1, m2, length.out=10),col.regions =col.l , main=paste(kname[[5]],"-ES",sep=""))
r6<-levelplot(log(z)~ x*y, data=newl[[6]], xlab="chr3(5Mb)",ylab="chr3(5Mb)", at=seq(m1, m2, length.out=10),col.regions =col.l , main=paste(kname[[6]],"-ES",sep=""))

print(grid.arrange(r3, r4, r2,r1,r5,r6, ncol=2))

dev.off()
#title("", outer=TRUE)

