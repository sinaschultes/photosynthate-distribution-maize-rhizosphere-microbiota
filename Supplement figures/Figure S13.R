rm(list=ls())

########## load libraries 
  library(cowplot)
  library(ggh4x)
  library(writexl)
  library(readxl)
  library(dplyr)
  library(phyloseq)
  library(phyloseqCompanion)
  library(microbiome)
  library(rcompanion)
  library(plyr)
  library(FSA)
  library(foreach)
  library(doParallel)
  library(data.table)
  library(ggplot2)
  library(mirlyn)
  library(microViz)
  library(gridExtra)
  library(grid)
  library(car)
  library(tibble)
  library(ggpubr)

########## load HR-SIP output data file (supplied)
DATA <- readxl::read_excel("./13C label incorporators for plotting.xlsx")
DATA <- as.data.frame(DATA, na.rm=TRUE)

DATA$Origin <- as.factor(DATA$Origin)
DATA$Group <- as.factor(DATA$Group)
DATA$X13C_label_cat  <- as.factor(DATA$X13C_label_cat )
DATA$Method <- as.factor(DATA$Method)
DATA$weightedlog2fc <- as.numeric(DATA$weightedlog2fc)
DATA$plotgroup <-  paste0(DATA$Group, DATA$X13C_label_cat)
DATA$plotgroup <- as.factor(DATA$plotgroup)

sub1 <- subset(DATA, Method=="X13C_label_cat")
summary(sub1)

sub2 <- subset(DATA, Method=="Origin")
summary(sub2)


######################################## Figure S13

  ########## Fungi
  strip <- strip_themed(background_x = elem_list_rect(fill = c('grey90')))
  sub4 <- subset(DATA, Method=="Origin" & Group =="Fungi")

  p <- ggplot(sub4, aes(y=Class_Genus_ASV_no, x=weightedlog2fc, fill=Origin, shape=Origin)) 
  plot2 <- p + geom_point(aes(size = abundance),alpha=0.9,stat="identity") + 
    guides(fill = guide_legend(override.aes = list(size = 10))) + theme_bw()+
    scale_y_discrete(limits=rev,name="Unique incorporators") + 
    scale_x_continuous(name="Log2 fold change")+
    scale_shape_manual(values=c(21,22,23,24))+
    scale_size_continuous(name="Abs. abundance",breaks=c(100,1000,10000,100000),range=c(5,30),limits=c(0,max(ceiling(sub4$abundance)))) + 
    scale_fill_manual(values=c("Seedborne base"="rosybrown4","Seedborne tip"="salmon3","Shootborne base"="aquamarine4","Shootborne tip"="paleturquoise3"))+
    facet_wrap2("Group",nrow=1,strip = strip )+
    theme_linedraw()+    
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          plot.title=element_blank(),
          plot.margin=unit(c(0.5, 1, 0, 1), "cm"),
          strip.text=element_text(size=16, color="black"),
          axis.text.y = element_text(size=16, color="black",margin= margin(r=3)),
          axis.text.x = element_text(size=16, color="black",margin= margin(t=3)),
          legend.text = element_text(size=16,color="black"),
          legend.title = element_text(size=18,face="bold",color="black"),
          legend.position = ("none"))
  plot2



  ########## Bacteria
  strip <- strip_themed(background_x = elem_list_rect(fill = c('grey90')))
  sub3 <- subset(DATA, Method=="Origin" & Group =="Bacteria")

  p <- ggplot(sub3, aes(y=Class_Genus_ASV_no, x=weightedlog2fc, fill=Origin, shape=Origin)) 
  plot1 <- p + geom_point(aes(size = abundance),alpha=0.9,stat="identity") + 
    guides(fill = guide_legend(override.aes = list(size = 10))) + theme_bw()+
    scale_y_discrete(limits=rev,name="Unique incorporator ASVs") + 
    scale_x_continuous(name="Log2 fold change")+
    scale_shape_manual(name="Origin",values=c(21,22,23,24))+
    scale_size_continuous(name="Abs. abundance",breaks=c(100,1000,10000,100000),range=c(5,30),limits=c(0,max(ceiling(sub4$abundance)))) + 
    scale_fill_manual(name="Origin",values=c("Seedborne base"="rosybrown4","Seedborne tip"="salmon3","Shootborne base"="aquamarine4","Shootborne tip"="paleturquoise3"))+ 
    facet_wrap2("Group",nrow=1,strip = strip )+
    theme_linedraw()+    
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          plot.title=element_blank(),
          plot.margin=unit(c(0.5, 1, 0, 1), "cm"),
          strip.text=element_text(size=16, color="black"),
          axis.text.y = element_text(size=16, color="black",margin= margin(r=3)),
          axis.text.x = element_text(size=16, color="black",margin= margin(t=3)),
          legend.text = element_text(size=16,color="black"),
          legend.title = element_text(size=18,face="bold",color="black"),
          legend.position = ("right"))
  plot1


  ########## Cercozoa
  strip <- strip_themed(background_x = elem_list_rect(fill = c('grey90')))
  sub5 <- subset(DATA, Method=="Origin" & Group =="Protists")

  p <- ggplot(sub5, aes(y=Class_Genus_ASV_no, x=weightedlog2fc, fill=Origin, shape=Origin)) 
  plot3 <- p + geom_point(aes(size = abundance),alpha=0.9,stat="identity") + 
    guides(fill = guide_legend(override.aes = list(size = 10))) + theme_bw()+
    scale_y_discrete(limits=rev,name="Unique incorporator ASVs") + 
    scale_x_continuous(name="Log2 fold change")+
    scale_shape_manual(values=c(21,22,23,24))+
    scale_size_continuous(name="Abs. abundance",breaks=c(100,1000,10000,100000),range=c(5,30),limits=c(0,max(ceiling(sub4$abundance)))) + 
    scale_fill_manual(values=c("Seedborne base"="rosybrown4","Seedborne tip"="salmon3","Shootborne base"="aquamarine4","Shootborne tip"="paleturquoise3"))+
    facet_wrap2("Group",nrow=1,strip = strip )+
    theme_linedraw()+    
    theme(axis.title.y=element_blank(),
          axis.title.x=element_text(size=18,face="bold",margin= margin(r=3)),
          plot.title=element_blank(),
          plot.margin=unit(c(0.5, 1, 0, 1), "cm"),
          strip.text=element_text(size=16, color="black"),
          axis.text.y = element_text(size=16, color="black",margin= margin(r=3)),
          axis.text.x = element_text(size=16, color="black",margin= margin(t=3)),
          legend.position = ("none"))
  plot3


########## combine plot
allplot2 <- plot_grid(
  plot1, plot2,plot3,
  labels = "AUTO",nrow=3, align="v", axis="lr", label_size = 25,rel_heights = c(1.1, 0.5,1)
)
allplot2

png("Figure S13.png", width = 500, height = 550, units = 'mm', res = 1000)
grid.arrange(allplot2,ncol=1) 
dev.off()

