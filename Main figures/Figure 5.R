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

{DATA$Origin <- as.factor(DATA$Origin)
DATA$Group <- as.factor(DATA$Group)
DATA$X13C_label_cat  <- as.factor(DATA$X13C_label_cat )
DATA$Method <- as.factor(DATA$Method)
DATA$weightedlog2fc <- as.numeric(DATA$weightedlog2fc)
DATA$plotgroup <-  paste0(DATA$Group, DATA$X13C_label_cat)
DATA$plotgroup <- as.factor(DATA$plotgroup)}

sub1 <- subset(DATA, Method=="X13C_label_cat")
summary(sub1)

sub2 <- subset(DATA, Method=="Origin")
summary(sub2)


######################################## Figure 5ab

  ########## Fig. 5a
  strip <- strip_themed(background_x = elem_list_rect(fill = c("royalblue3" ,"#9C179EFF" ,"indianred1", "goldenrod1"))) 
  
  a = ggplot(sub1, 
             aes(x = Group,
                 y = weightedlog2fc,
                 fill = X13C_label_cat,
                 shape = X13C_label_cat))+ 
    geom_boxplot(outlier.size=-1)+
    stat_summary(fun=mean, geom="point", shape=4, size=3, color="black", fill="black") +
    scale_fill_manual(values=c("royalblue3" ,"#9C179EFF" ,"indianred1" ,"goldenrod1"))+ 
    scale_y_continuous(name="Log2 fold change")+
    geom_jitter(position=position_jitter(0.2), shape = 1,size=3)+     
    stat_boxplot(geom = "errorbar", width = 0.2)+           
    theme_linedraw()+    
    facet_wrap2("X13C_label_cat",nrow=1,strip = strip )+
    theme_linedraw()+    
    theme(axis.title.y=element_text(size=18,face="bold",margin= margin(r=3)),
          axis.title.x=element_blank(),
          plot.title=element_blank(),
          plot.margin=unit(c(1, 1, 1, 1), "cm"),
          strip.text=element_text(size=16, color="black"),
          axis.text.y = element_text(size=16, color="black",margin= margin(r=3)),
          axis.text.x = element_text(size=16, color="black",angle=40, hjust=1,vjust=1,margin= margin(t=3)),
          legend.position = ("none"))
  print(a)


  ########## Fig. 5b
  strip <- strip_themed(background_x = elem_list_rect(fill = c('rosybrown4','salmon3','aquamarine4', 'paleturquoise3')))
  
  b = ggplot(sub2, 
             aes(x = Group,
                 y = weightedlog2fc,
                 fill = Origin,
                 shape = Origin))+ 
    geom_boxplot(outlier.size=-1)+
    stat_summary(fun=mean, geom="point", shape=4, size=3, color="black", fill="black") +
    scale_fill_manual(values=c('rosybrown4','salmon3','aquamarine4', 'paleturquoise3'))+
    scale_y_continuous(name="Log2 fold change")+
    geom_jitter(position=position_jitter(0.2), shape = 1,size=3)+         
    stat_boxplot(geom = "errorbar", width = 0.2)+           
    theme_linedraw()+    
    facet_wrap2("Origin",nrow=1,strip = strip )+
    theme_linedraw()+    
    theme(axis.title.y=element_text(size=18,face="bold",margin= margin(r=3)),
          axis.title.x=element_blank(),
          plot.title=element_blank(),
          plot.margin=unit(c(1, 1, 1, 1), "cm"),
          strip.text=element_text(size=16, color="black"),
          axis.text.y = element_text(size=16, color="black",margin= margin(r=3)),
          axis.text.x = element_text(size=16, color="black",angle=40, hjust=1,vjust=1,margin= margin(t=3)),
          legend.position = ("none"))
  print(b)



  ########## Fig. 5a 5b combined
  plot5ab <- plot_grid(
    a, b,
    labels = "AUTO", ncol = 2, label_size = 40
  )
  plot5ab
  


######################################## Figure 5c

  ########## Fungi
  strip <- strip_themed(background_x = elem_list_rect(fill = c("grey90" ))) 
  sub4 <- subset(DATA, Method=="X13C_label_cat" & Group =="Fungi")

  p <- ggplot(sub4, aes(y=Class_Genus_ASV_no, x=weightedlog2fc, fill=X13C_label_cat, shape=X13C_label_cat)) 
  plot2 <- p + geom_point(aes(size = abundance),alpha=0.9,stat="identity") + 
    guides(fill = guide_legend(override.aes = list(size = 10))) + theme_bw()+
    scale_y_discrete(limits=rev,name="Unique incorporators") + 
    scale_x_continuous(name="Log2 fold change")+
    scale_shape_manual(name="13C label category",values=c(21,22,23,24))+
    scale_size_continuous(name="Abs. abundance",breaks=c(100,1000,10000,100000),range=c(5,30),limits=c(0,max(ceiling(sub4$abundance)))) + 
    scale_fill_manual(name="13C label category",values=c("Category 1"="royalblue3","Category 2"="#9C179EFF","Category 3"="indianred1","Category 4"="goldenrod1"))+
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
  plot2
  
  
  ########## Bacteria
  strip <- strip_themed(background_x = elem_list_rect(fill = c('grey90')))
  sub3 <- subset(DATA, Method=="X13C_label_cat" & Group =="Bacteria")

  p <- ggplot(sub3, aes(y=Class_Genus_ASV_no, x=weightedlog2fc, fill=X13C_label_cat, shape=X13C_label_cat)) 
  plot1 <- p + geom_point(aes(size = abundance),alpha=0.9,stat="identity") + 
    guides(fill = guide_legend(override.aes = list(size = 10))) + theme_bw()+
    scale_y_discrete(limits=rev,name="Unique incorporator ASVs") + 
    scale_x_continuous(name="Log2 fold change")+
    scale_shape_manual(values=c(21,22,23,24))+
    scale_size_continuous(name="Abs. abundance",breaks=c(100,1000,10000,100000),range=c(5,30),limits=c(0,max(ceiling(sub4$abundance)))) + 
    scale_fill_manual(values=c("Category 1"="royalblue3","Category 2"="#9C179EFF","Category 3"="indianred1","Category 4"="goldenrod1"))+ 
    facet_wrap2("Group",nrow=1,strip = strip )+
    theme_linedraw()+    
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          plot.title=element_blank(),
          plot.margin=unit(c(0.5, 1, 0, 1), "cm"),
          strip.text=element_text(size=16, color="black"),
          axis.text.y = element_text(size=16, color="black",margin= margin(r=3)),
          axis.text.x = element_text(size=16, color="black",margin= margin(t=3)),
          legend.position = ("none"))
  plot1
  
  

  ########## Cercozoa
  strip <- strip_themed(background_x = elem_list_rect(fill = c('grey90')))
  sub5 <- subset(DATA, Method=="X13C_label_cat" & Group =="Protists")
  
  p <- ggplot(sub5, aes(y=Class_Genus_ASV_no, x=weightedlog2fc, fill=X13C_label_cat, shape=X13C_label_cat)) 
  plot3 <- p + geom_point(aes(size = abundance),alpha=0.9,stat="identity") + 
    guides(fill = guide_legend(override.aes = list(size = 10))) + theme_bw()+
    scale_y_discrete(limits=rev,name="Unique incorporator ASVs") + 
    scale_x_continuous(name="Log2 fold change")+
    scale_shape_manual(values=c(21,22,23,24))+
    scale_size_continuous(name="Abs. abundance",breaks=c(100,1000,10000,100000),range=c(5,30),limits=c(0,max(ceiling(sub4$abundance)))) + 
    scale_fill_manual(values=c("Category 1"="royalblue3","Category 2"="#9C179EFF","Category 3"="indianred1","Category 4"="goldenrod1"))+
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


  ########## combine to form Fig. 5c
  plot5c <- plot_grid(
    plot1, plot2,plot3,
    labels = c("C"),nrow=3, align="v", axis="lr", label_size = 40, rel_heights = c(1, 0.75,1)
  )
  plot5c


########## combine plots a, b and c to form Fig. 5
allplot3 <- plot_grid(
  plot5ab, plot5c,rel_heights = c(0.6, 2),
  nrow = 2)
allplot3

png("Figure 5.png", width = 490, height = 700, units = 'mm', res = 1000) 
grid.arrange(allplot3,ncol=1) 
dev.off()


######################################## statistics

  ########## statistics Fig. 5a
  sub1 <- subset(DATA, Method=="X13C_label_cat" & X13C_label_cat =="Category 1")
  sub2 <- subset(DATA, Method=="X13C_label_cat" & X13C_label_cat =="Category 2")
  sub3 <- subset(DATA, Method=="X13C_label_cat" & X13C_label_cat =="Category 3")
  sub4 <- subset(DATA, Method=="X13C_label_cat" & X13C_label_cat =="Category 4")
  
  shapiro.test(sub1$weightedlog2fc)
  shapiro.test(sub2$weightedlog2fc)
  shapiro.test(sub3$weightedlog2fc)
  shapiro.test(sub4$weightedlog2fc)
  
  ggqqplot(sub1, x="weightedlog2fc")
  ggqqplot(sub2, x="weightedlog2fc")
  ggqqplot(sub3, x="weightedlog2fc")
  ggqqplot(sub4, x="weightedlog2fc")
  
  leveneTest(weightedlog2fc~Group, data=sub1) 
  leveneTest(weightedlog2fc~Group, data=sub2)
  leveneTest(weightedlog2fc~Group, data=sub3)
  leveneTest(weightedlog2fc~Group, data=sub4)
  
  
  ########## Kruskal Wallis rank sum test and Dunn’s post hoc tests with Benjamini-Hochberg correction
  sub1$Group <- as.factor(sub1$Group)
  kruskal.test(weightedlog2fc ~ Group, data = sub1) 
  
  sub2$Group <- as.factor(sub2$Group)
  kruskal.test(weightedlog2fc ~ Group, data = sub2) 
  
  sub3$Group <- as.factor(sub3$Group)
  kruskal.test(weightedlog2fc ~ Group, data = sub3) 
  
  sub4$Group <- as.factor(sub4$Group)
  kruskal.test(weightedlog2fc ~ Group, data = sub4) 

  #posthoc on category 1 subset
  Phocdunn <- dunnTest(weightedlog2fc ~ Group,
                       data=sub1,
                       method="bh") 
  Phocdunns <- Phocdunn$res
  cld <- cldList(comparison = Phocdunns$Comparison,
                 p.value    = Phocdunns$P.adj,
                 threshold  = 0.05)[1:2]
  cld
  
  #posthoc on category 4 subset
  Phocdunn <- dunnTest(weightedlog2fc ~ Group,
                       data=sub4,
                       method="bh") 
  Phocdunns <- Phocdunn$res
  cld <- cldList(comparison = Phocdunns$Comparison,
                 p.value    = Phocdunns$P.adj,
                 threshold  = 0.05)[1:2]
  cld
  


  ########## statistics Fig. 5b
  sub1 <- subset(DATA, Method=="Origin" & Origin =="Seedborne base")
  sub2 <- subset(DATA, Method=="Origin" & Origin =="Seedborne tip")
  sub3 <- subset(DATA, Method=="Origin" & Origin =="Shootborne base")
  sub4 <- subset(DATA, Method=="Origin" & Origin =="Shootborne tip")
  
  shapiro.test(sub1$weightedlog2fc)
  shapiro.test(sub2$weightedlog2fc)
  shapiro.test(sub3$weightedlog2fc)
  shapiro.test(sub4$weightedlog2fc)
  
  ggqqplot(sub1, x="weightedlog2fc")
  ggqqplot(sub2, x="weightedlog2fc")
  ggqqplot(sub3, x="weightedlog2fc")
  ggqqplot(sub4, x="weightedlog2fc")
  
  leveneTest(weightedlog2fc~Group, data=sub1) 
  leveneTest(weightedlog2fc~Group, data=sub2) 
  leveneTest(weightedlog2fc~Group, data=sub3) 
  leveneTest(weightedlog2fc~Group, data=sub4) 
  
  
  ########## Kruskal Wallis rank sum test and Dunn’s post hoc tests with Benjamini-Hochberg correction
  sub1$Group <- as.factor(sub1$Group)
  kruskal.test(weightedlog2fc ~ Group, data = sub1) 
  
  sub2$Group <- as.factor(sub2$Group)
  kruskal.test(weightedlog2fc ~ Group, data = sub2) 
  
  sub3$Group <- as.factor(sub3$Group)
  kruskal.test(weightedlog2fc ~ Group, data = sub3) 
  
  sub4$Group <- as.factor(sub4$Group)
  kruskal.test(weightedlog2fc ~ Group, data = sub4) 
  
  #posthoc on seedborne tip subset
  Phocdunn <- dunnTest(weightedlog2fc ~ Group,
                       data=sub2,
                       method="bh")    
  Phocdunns <- Phocdunn$res
  cld <- cldList(comparison = Phocdunns$Comparison,
                 p.value    = Phocdunns$P.adj,
                 threshold  = 0.05)[1:2]
  cld
  
