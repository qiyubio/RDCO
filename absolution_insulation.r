#library()
## Author: Qi Yu
## Contact: dryuqi@gmail.com


read_file<-function(filename,title=NA){
  homo<-read.table(filename,header=FALSE)
  names(homo)<-c("chr","start","end","insulation")
  homo<-homo[which(homo$insulation >-10),]
  new<-aggregate(by=list(homo[,1]), homo[,4], min)
  names(new)<-c("chr","minin")
  homonew<-join(homo,new,by = "chr")
  homonew$absolutein<-homonew$insulation-homonew$minin
  homonew$type<-title
  return(homonew)
}

fnew5<-read_file("./newdata/old/inter_30_preLep.hic_merge_100000.out_matrix_withhead.is500001.ids200001.insulation.bedGraph.bedGraph",title="PreLep")
fnew7<-read_file("./newdata/old/inter_30_Lep.hic_merge_100000.out_matrix_withhead.is500001.ids200001.insulation.bedGraph.bedGraph",title="Lep")
#fnew15<-read_file("4N_ZPD_mm10_40000_fragment1.out",title="ZPD")
fnewz<-read_file("./newdata/old/zygotene_30.hic_merge_insulation.bedGraph.bedGraph",title="Zygo")
fnewp<-read_file("./newdata/old/pachytene_30.hic_merge_insulation.bedGraph.bedGraph",title="Pachy")
fnewd<-read_file("./newdata/old/diplotene_30.hic_merge_insulation.bedGraph.bedGraph",title="Diplo")
fnews<-read_file("./newdata/old/spermatid_30.hic_merge_insulation.bedGraph.bedGraph",title="Spermatid")
fnewsd<-read_file("./boundaries/sperm_inter_30.hic_merge_100000.out_matrix_withhead.is500001.ids200001.insulation.bedGraph.bedGraph",title="Sperm")

f<-rbind(fnew5,fnew7,fnewz,fnewp,fnewd,fnews,fnewsd)

f$type<-factor(f$type,levels=c("PreLep","Lep","Zygo","Pachy","Diplo","Spermatid","Sperm"))

test <- vector(mode="numeric", length=0)
winn<-sapply(2:60,function(i){
  test[i+1]<-round(40000*1.12^i)
})
tmp<-data.frame(winn)

#tmp$winny<-200*(tmp$winn^(-0.5))
tmp$winny<-400000*(tmp$winn^(-1))

f$mean_f <- ave(f$absolutein, as.factor(f$type), FUN=median)

  
  #geom_boxplot(width=0.4, notch=TRUE, size=1.5, aes(fill=mean_f))+scale_fill_gradient(low="white",high = "blue")+
#fnewa<-aggregate(fnew5,by=list(fnew5$xp),FUN=mean,na.rm=TRUE)

d3=ggplot(f,aes(x=type,y=absolutein,fill=mean_f))+
  theme_bw()+
  #geom_violin(scale = "width",size=1.5,trim = TRUE,position=position_dodge(0.85))+
  geom_boxplot(width=0.4, notch=TRUE, position=position_dodge(0.85),size=1.5)+
  scale_fill_gradient(low="white",high = "blue")+
  #scale_x_discrete(expand=c(0.2,0))+
  #annotation_logticks(sides="l")+
  #geom_dotplot(binaxis = "y", stackdir = "center",alpha=0.5)+
  
  ylab("Absolute Insulation score")+
  xlab("")+
  #scale_y_log10()+
  #ylim(0,3)+
  theme(text = element_text(size=20,family="sans"),axis.text.x=element_text(size=20,color="black"),legend.position="top",)+
  #scale_colour_discrete(name = "")+
  guides(colour=guide_legend( keywidth=2, keyheight=2, title=""))+
  #scale_color_manual(values=c("#C35817","grey"),name="")+
  #scale_color_manual(values=c("#cc0000","#009900","#ff9900"),breaks=c("Dog1","Dog2","Dog3"))+
  ggtitle("")
  #annotation_logticks(sides="l")
ggsave(plot=d3,filename=paste("prel_lep_Z_p_D_absolute_insulation_histone.tiff", sep=''),height=5,width=9)

aggregate(by=list(f$type), f$absolutein, median)
#fnew<-fnew[which(fnew$chr!="chrX" &fnew$chr!="chrY"),]
wilcox.test(fnew5$absolutein, fnew7$absolutein, alternative = "greater")$p.value
wilcox.test(fnew7$absolutein, fnewz$absolutein, alternative = "greater")$p.value
wilcox.test(fnewz$absolutein, fnewp$absolutein, alternative = "less")$p.value
wilcox.test(fnewp$absolutein, fnewd$absolutein, )$p.value
wilcox.test(fnewd$absolutein, fnews$absolutein, alternative = "less")$p.value
wilcox.test(fnews$absolutein, fnewsd$absolutein, alternative = "less")$p.value

#d3=ggplot(fnew)+
#  geom_line(aes(x=xp,y=yp,color=chr),size=1)+
#  geom_point(aes(x=xp,y=yp,color=chr),size=1)+
#  xlab("Genomic Distance")+
#  ylab("Interaction Properties")+
#  theme_bw()+
#  theme(text = element_text(size=25))+
#   ggtitle(title)+
#   scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
#   scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
#   annotation_logticks(sides="lb")
# ggsave(plot=d3,filename=paste(title,"_",genome,"_",winsize,"fragment_pos_per_chr_noXY.png", sep=''),height=10,width=16)
# 
# fnew1<-aggregate(fnew,by=list(fnew$xp),FUN=mean,na.rm=TRUE)
# 
# d2=ggplot(fnew1)+
#   geom_line(aes(x=xp,y=yp),size=1)+
#   geom_point(aes(x=xp,y=yp),size=1)+
#   xlab("Genomic Distance")+
#   ylab("Interaction Properties")+
#   theme_bw()+
#   theme(text = element_text(size=25))+
#   ggtitle(title)+
#   scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
#   scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
#   annotation_logticks(sides="lb")
# ggsave(plot=d2,filename=paste(title,"_",genome,"_",winsize,"fragment_pos.png", sep=''),height=10,width=16)
