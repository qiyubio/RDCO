## Author: Qi Yu
## Contact: dryuqi@gmail.com

library("lattice")
library("reshape2")
library(RColorBrewer)
library(gridExtra)
library(reshape2)


defm<-function(file,chr=NA,name=NA){
  matrix1<-read.table(file,header=FALSE)
  names(matrix1)<-c("x","y","z")
  matrix1$y<-max(matrix1$y)-matrix1$y
  return(matrix1)
}

highres_p<-defm("181016_9014_4_ACTTGA_GC_072318_Hi-C_B6XCAST_Pachytene.mm10_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_5000000.matrix")
p<-defm("180903_9013_1_ACTTGA_GC_072318_Hi-C_B6XCAST_Pachytene.mm10_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_5000000.matrix")

mergetwo<-merge(highres_pre,pre,by=c("x","y"),all=T)
mergetwo[is.na(mergetwo)] <- 0
mergetwo$ratio<-mergetwo$z.x-mergetwo$z.y
mergetwo$ratio<-mergetwo$ratio-min(mergetwo$ratio)
chr<-"chr3"
name<-"Pachytene"


col.l<-colorRampPalette(brewer.pal(9,'YlGnBu'))(100)
png(paste(chr,"_",name,".subtract2.png",sep=""),width=1000, height=1000, res=200)
myplot1<-levelplot(log(ratio) ~ x*y, data=mergetwo, xlab=" ",ylab=" ", col.regions = col.l, colorkey=list(space="bottom"),scales=list(x=list(alternatin=2, at=c(20, 60, 95, 130),labels=c("B6_chr1","CAST_chr1","B6_chr3","CAST_chr3")),y=list(at=c(15, 50, 85,125),labels=rev(c("B6_chr1","CAST_chr1","B6_chr3","CAST_chr3")))), main="Pachy_highres-Pachy")
print(myplot1)
dev.off()

#########################################################
highres_pre<-defm("181016_9014_3_CTTGTA_GC_050218_Hi_C_B6XCAST_Preleptotene.mm10_genome_mask_CAST.bwt2pairs_VI_chr1_chr3_5000000.matrix")
pre<-defm("180530_prelep_B6XCAST_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_5000000.matrix")

mergetwo<-merge(highres_pre,pre,by=c("x","y"),all=T)
mergetwo[is.na(mergetwo)] <- 0
mergetwo$ratio<-mergetwo$z.x-mergetwo$z.y
mergetwo$ratio<-mergetwo$ratio-min(mergetwo$ratio)
chr<-"chr3"
name<-"Prelep"

col.l<-colorRampPalette(brewer.pal(9,'YlGnBu'))(100)
png(paste(chr,"_",name,".subtract2.png",sep=""),width=1000, height=1000, res=200)
myplot1<-levelplot(log(ratio) ~ x*y, data=mergetwo, xlab=" ",ylab=" ", col.regions = col.l, colorkey=list(space="bottom"),scales=list(x=list(alternatin=2, at=c(20, 60, 95, 130),labels=c("B6_chr1","CAST_chr1","B6_chr3","CAST_chr3")),y=list(at=c(15, 50, 85,125),labels=rev(c("B6_chr1","CAST_chr1","B6_chr3","CAST_chr3")))), main="Prel_highres-Prel")
print(myplot1)
dev.off()

############################################################
highres_pre<-defm("181016_9014_3_CTTGTA_GC_050218_Hi_C_B6XCAST_Preleptotene.mm10_genome_mask_CAST.bwt2pairs_VI_chr1_chr3_nocis_5000000.matrix")
pre<-defm("180530_lep_B6XCAST_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_nocis_5000000.matrix")

mergetwo<-merge(highres_pre,pre,by=c("x","y"),all=T)
mergetwo[is.na(mergetwo)] <- 0
mergetwo$ratio<-mergetwo$z.x-mergetwo$z.y
mergetwo$ratio<-mergetwo$ratio-min(mergetwo$ratio)
chr<-"chr3"
name<-"Prel_highres-Lep"

col.l<-colorRampPalette(brewer.pal(9,'YlGnBu'))(100)
png(paste(chr,"_",name,".subtract2.png",sep=""),width=1000, height=1000, res=200)
myplot1<-levelplot(log(ratio) ~ x*y, data=mergetwo, xlab=" ",ylab=" ", col.regions = col.l, colorkey=list(space="bottom"),scales=list(x=list(alternatin=2, at=c(20, 60, 95, 130),labels=c("B6_chr1","CAST_chr1","B6_chr3","CAST_chr3")),y=list(at=c(15, 50, 85,125),labels=rev(c("B6_chr1","CAST_chr1","B6_chr3","CAST_chr3")))), main=name)
print(myplot1)
dev.off()

highres_pre<-defm("181016_9014_3_CTTGTA_GC_050218_Hi_C_B6XCAST_Preleptotene.mm10_genome_mask_CAST.bwt2pairs_VI_chr1_chr3_nocis_5000000.matrix")
pre<-defm("180530_lep_B6XCAST_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_nocis_5000000.matrix")
name<-"Lep-Prel_highres"
#mergetwo<-mergetwo[which(mergetwo$x+mergetwo$y!=max(mergetwo$x)),]
mergetwo<-merge(highres_pre,pre,by=c("x","y"),all=T)
mergetwo[is.na(mergetwo)] <- 0
mergetwo$ratio<-mergetwo$z.y-mergetwo$z.x
#mergetwo$ratio<-mergetwo$ratio-min(mergetwo$ratio)
col.l<-colorRampPalette(brewer.pal(9,'YlGnBu'))(100)
png(paste(chr,"_",name,".subtract_no_dialog.png",sep=""),width=1000, height=1000, res=200)
myplot1<-levelplot(ratio ~ x*y, data=mergetwo, xlab=" ",ylab=" ", col.regions = col.l, colorkey=list(space="bottom"),scales=list(x=list(alternatin=2, at=c(20, 60, 95, 130),labels=c("B6_chr1","CAST_chr1","B6_chr3","CAST_chr3")),y=list(at=c(15, 50, 85,125),labels=rev(c("B6_chr1","CAST_chr1","B6_chr3","CAST_chr3")))), main=name)
print(myplot1)
dev.off()

def<-function(file,name){
highres_pre<-defm(file)
pre<-defm("180901_9011_1_TTAGGC_GC_072318_Hi-C_B6XCAST_2N.mm10_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_nocis_5000000.matrix")
name<-paste(name,"2N",sep='-')
#mergetwo<-mergetwo[which(mergetwo$x+mergetwo$y!=max(mergetwo$x)),]
mergetwo<-merge(highres_pre,pre,by=c("x","y"),all=T)
mergetwo[is.na(mergetwo)] <- 0
mergetwo$ratio<-mergetwo$z.x-mergetwo$z.y
mergetwo$ratio<-mergetwo$ratio-min(mergetwo$ratio)
col.l<-colorRampPalette(brewer.pal(9,'YlGnBu'))(100)
png(paste(chr,"_",name,".subtract_2N_no_dialog_log.png",sep=""),width=1000, height=1000, res=200)
myplot1<-levelplot(log(ratio) ~ x*y, data=mergetwo, xlab=" ",ylab=" ", col.regions = col.l, colorkey=list(space="bottom"),scales=list(x=list(alternatin=2, at=c(20, 60, 95, 130),labels=c("B6_chr1","CAST_chr1","B6_chr3","CAST_chr3")),y=list(at=c(15, 50, 85,125),labels=rev(c("B6_chr1","CAST_chr1","B6_chr3","CAST_chr3")))), main=name)
print(myplot1)
dev.off()
}

def("181016_9014_3_CTTGTA_GC_050218_Hi_C_B6XCAST_Preleptotene.mm10_genome_mask_CAST.bwt2pairs_VI_chr1_chr3_nocis_5000000.matrix","Prelep_highres")
def("180530_lep_B6XCAST_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_nocis_5000000.matrix","Leptotene")
def("180901_9011_1_ACTGAT_GC_072318_Hi-C_B6XCAST_Zygotene.mm10_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_nocis_5000000.matrix","Zygotene")
def("181016_9014_4_ACTTGA_GC_072318_Hi-C_B6XCAST_Pachytene.mm10_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_nocis_5000000.matrix",name="Pachytene_highres")
def("180903_9013_1_GAGTGG_GC_072318_Hi-C_B6XCAST_Diplotene.mm10_genome_mask_CAST.bwt2pairs_interaction_MD_CT_VI.txt_chr1_chr3_nocis_5000000.matrix",name="Diplotene")



########################################################################

