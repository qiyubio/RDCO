## Author: Qi Yu
## Contact: dryuqi@gmail.com

library(ggplot2)
library(wesanderson)
args<-(commandArgs(TRUE))

#filename = args[1]
#title = args[2]

readfile<-function(filename,title=NA){
  homo<-read.table(filename)
  homo<-homo[c(2,5,7)]
  names(homo)<-c("start","count","end")
  homo$distance<-abs(homo$end-homo$start)/1000
  homo$type<-title
  return(homo)
}#end of function readfile

table<-readfile("/home/qiyu/data/hiC/190110_9018_1_CGTACG_GC_121718_HiChip_B6_4C_H3K4me3.intra.loop_counts.bedpe_at_B6_wt_DMC1_peaks_no_blacklist.bed_both_noself","4C_noself_HSs")
table2<-readfile("/home/qiyu/data/hiC/190110_9018_1_CGATGT_GC_121718_HiChip_B6_Lep_H3K4me3.intra.loop_counts.bedpe_at_B6_wt_DMC1_peaks_no_blacklist.bed_both_noself","Lep_noself_HSs")

new2<-rbind(table,table2)

#new2$type<-factor(new2$type,levels=c("PreL_Boris","P_Boris","P_CTCF","Z_Boris","Z_CTCF","P_Boris_DSB","P_CTCF_DSB"))
options(scipen=10000)

d2<-ggplot(new2,aes(x=distance,fill=type))+
  #geom_density(linesize=2)+
  geom_histogram(color="black")+
  facet_wrap( ~ type, ncol=1)+
  #facet_grid(type ~ .)+
  #geom_freqpoly()+
  #geom_density(size=1.5)+
  #scale_fill_manual(values=c(rev(wes_palette("Rushmore",5)),"blue","purple"),name="")+
  theme_bw()+
  xlab("loop Length (kb)")+
  scale_x_log10()+
  annotation_logticks(sides="b")+
  theme(text = element_text(size=25),legend.title=element_blank(),legend.justification=c(1,01))
#scale_x_log10()+scale_y_log10()+

ggsave(plot=d2,filename="noself_HSs_4C_Lep.png",height=7,width=10,dpi=200)

d2<-ggplot(new2,aes(x=count,color=type))+
  #geom_density(linesize=2)+
  stat_ecdf(size=2)+
  facet_wrap( ~ type, ncol=1)+
  #facet_grid(type ~ .)+
  #geom_freqpoly()+
  #geom_density(size=1.5)+
  #scale_fill_manual(values=c(rev(wes_palette("Rushmore",5)),"blue","purple"),name="")+
  theme_bw()+
  xlab("PET")+
  scale_x_log10()+
  scale_y_log10()+
  annotation_logticks(sides="lb")+
  theme(text = element_text(size=25),legend.title=element_blank(),legend.justification=c(1,01))
#scale_x_log10()+scale_y_log10()+

ggsave(plot=d2,filename="PET_count_noself_HSs_4C_Lep.png",height=7,width=10,dpi=200)

new3<-new2[which(new2$count>4),]
d2<-ggplot(new3,aes(x=distance,fill=type))+
  #geom_density(linesize=2)+
  geom_histogram(color="black")+
  facet_wrap( ~ type, ncol=1)+
  ggtitle(">=5")+
  #facet_grid(type ~ .)+
  #geom_freqpoly()+
  #geom_density(size=1.5)+
  #scale_fill_manual(values=c(rev(wes_palette("Rushmore",5)),"blue","purple"),name="")+
  theme_bw()+
  xlab("loop Length (kb)")+
  scale_x_log10()+
  annotation_logticks(sides="b")+
  theme(text = element_text(size=25),legend.title=element_blank(),legend.justification=c(1,01))
#scale_x_log10()+scale_y_log10()+
#scale_x_log10()+scale_y_log10()+

ggsave(plot=d2,filename="noself_HSs_4C_Lep_5.png",height=7,width=10,dpi=200)


