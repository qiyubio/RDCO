
## Author: Qi Yu
## Contact: dryuqi@gmail.com

library("lattice")
library("reshape2")
library(RColorBrewer)
library(gridExtra)
library(reshape2)


defplot2<-function(file,chr=NA,name=NA){
  
  matrix1<-read.table(file,header=FALSE)
  names(matrix1)<-c("x","y","z")
  matrix1$y<-max(matrix1$y)-matrix1$y
  #return(normalized_count)
  png(paste(chr,"_",name,".four.png",sep=""),width=1000, height=1000, res=200)

  myplot1<-levelplot(log(z) ~ x*y, data=matrix1, xlab="X" , col.regions = heat.colors(100)[length(heat.colors(100)):1], main=paste(name,chr))
  print(myplot1)
  
  dev.off()
  
}


defplot2("180903_9013_1_ACTTGA_GC_072318_Hi-C_B6XCAST_Pachytene.mm10_genome_mask_CAST.bwt2pairs_heterozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="Pachytene",chr="chr3")
defplot2("180901_9011_1_ACTGAT_GC_072318_Hi-C_B6XCAST_Zygotene.mm10_genome_mask_CAST.bwt2pairs_heterozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="Zygotene",chr="chr3")
defplot2("180530_prelep_B6XCAST_genome_mask_CAST.bwt2pairs_heterozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="preleptotene",chr="chr3")
defplot2("180530_lep_B6XCAST_genome_mask_CAST.bwt2pairs_heterozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="leptotene",chr="chr3")
defplot2("180901_9011_1_TTAGGC_GC_072318_Hi-C_B6XCAST_2N.mm10_genome_mask_CAST.bwt2pairs_heterozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="2N",chr="chr3")
defplot2("180903_9013_1_GAGTGG_GC_072318_Hi-C_B6XCAST_Diplotene.mm10_genome_mask_CAST.bwt2pairs_heterozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="Diplotene",chr="chr3")
defplot2("SRR2240730_bowtie2_mask_CAST_129s.bwt2pairs_homozygous.bam_ncc.txt_chr3.out_5000000.matrix",name="EScell",chr="chr3")
defplot2("CRD018512_bowtie2_mask_DBA2J.bwt2pairs_interaction_homozygous_samechr.bam_ncc.txt_VI_chr3_5000000.matrix",name="ES7.5",chr="chr3")
defplot2("181016_9014_3_CTTGTA_GC_050218_Hi_C_B6XCAST_Preleptotene.mm10_genome_mask_CAST.bwt2pairs_VI_chr3_5000000.matrix",name="prelep_highres",chr="chr3")

#############plot quadrant with nolog####################################

defplot2nolog<-function(file,chr=NA,name=NA){
  
  matrix1<-read.table(file,header=FALSE)
  names(matrix1)<-c("x","y","z")
  
  #return(normalized_count)
  png(paste(chr,"_",name,".four.png",sep=""),width=1000, height=1000, res=200)
  
  myplot1<-levelplot(z ~ x*y, data=matrix1, xlab="X" , col.regions = heat.colors(100)[length(heat.colors(100)):1], main=paste(name,chr))
  print(myplot1)
  
  dev.off()
  
}

defplot2nolog("181016_9014_4_ACTTGA_GC_072318_Hi-C_B6XCAST_Pachytene.mm10_genome_mask_CAST.bwt2pairs_interaction_diffchr.bam_ncc_MD_CT_RE_VI.txt_chr1_chr3_5000000.matrix",name="pachytene",chr="chr1_chr3")




