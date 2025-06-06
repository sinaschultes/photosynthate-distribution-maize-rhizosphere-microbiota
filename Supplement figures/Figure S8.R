rm(list = ls())

########## load libraries
library(writexl)
library(readxl)
library(dplyr)
library(plyr)
library(data.table)
library(ggplot2)
library(microViz)
library(gridExtra)
library(grid)
library(tibble)
library(ggpubr)
library(cowplot)


########## load and prep data (supplied)
  data <- read_excel("iCAMP results.xlsx")
  data <- as.data.frame(data)
  
  data$Mean <- as.numeric(data$Mean)
  data$Stdev <- as.numeric(data$Stdev)
  data$Meanpercent <- data$Mean*100
  data$Stdevpercent <- data$Stdev*100


########## plot1
  sub1 <- subset(data,Organism=="Bacteria" & Study=="I")
  sub1
  p1 <-
    ggplot(sub1, aes(x = Process, y = Meanpercent, fill = Group)) +
    geom_col(position = position_dodge()) +
    geom_errorbar(aes(x = Process, ymin = Meanpercent - Stdevpercent, ymax = Meanpercent + Stdevpercent), position = position_dodge(), colour = "black", alpha = 0.9, size = 0.4) +
    scale_y_continuous(name = "Relative importance [%]", limits = c(0, 102)) +
    scale_x_discrete(name = "Community assembly process", labels = c("Heterogeneous.Selection"="Het. select.","Homogeneous.Selection"="Hom. select.","Dispersal.Limitation"="Disp. limit.","Drift.and.Others"="Drift","Homogenizing.Dispersal"="Hom. disp.")) +
    scale_fill_manual(name = "Root section", values = c(
      'darkslateblue','chocolate2'
    )) + 
    geom_text(aes(label = CLD, y = Meanpercent + Stdevpercent + 3),position = position_dodge(0.9),vjust = 0,size = 6) +
    ggtitle("Bacteria")+
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_blank(),
      plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
      plot.title = element_text(size=24, hjust=0.5, color="black", margin = margin(b=15)),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = ("none"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16)
    )
  
  newSTorder = c("Heterogeneous.Selection","Homogeneous.Selection","Dispersal.Limitation","Drift.and.Others","Homogenizing.Dispersal")
  p1$data$Process <- as.character(p1$data$Process)
  p1$data$Process <- factor(p1$data$Process, levels=newSTorder)
  
  p1



########## plot2
  sub1 <- subset(data,Organism=="Bacteria" & Study=="II" & Factor=="Carbon_cat_4")
  sub1
  
  p2 <-
    ggplot(sub1, aes(x = Process, y = Meanpercent, fill = Group)) +
    geom_col(position = position_dodge()) +
    geom_errorbar(aes(x = Process, ymin = Meanpercent - Stdevpercent, ymax = Meanpercent + Stdevpercent), position = position_dodge(), colour = "black", alpha = 0.9, size = 0.4) +
    scale_y_continuous(name = "Relative importance [%]", limits = c(0, 102)) +
    scale_x_discrete(name = "Community assembly process", labels = c("Heterogeneous.Selection"="Het. select.","Homogeneous.Selection"="Hom. select.","Dispersal.Limitation"="Disp. limit.","Drift.and.Others"="Drift","Homogenizing.Dispersal"="Hom. disp.")) +
    scale_fill_manual(name = "Root section", values=c("royalblue3" ,"#9C179EFF" ,"indianred1" ,"goldenrod1"
    )) + 
    geom_text(aes(label = CLD, y = Meanpercent + Stdevpercent + 3),position = position_dodge(0.9),vjust = 0,size = 6) +
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 10), size = 18, face = "bold"),
      plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
      plot.title = element_text(size=24, hjust=0.5, color="black", margin = margin(b=15)),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = ("none"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16)
    )
  
  newSTorder = c("Heterogeneous.Selection","Homogeneous.Selection","Dispersal.Limitation","Drift.and.Others","Homogenizing.Dispersal")
  p2$data$Process <- as.character(p2$data$Process)
  p2$data$Process <- factor(p2$data$Process, levels=newSTorder)
  
  p2



########## plot3
  sub1 <- subset(data,Organism=="Bacteria" & Study=="II" & Factor=="Root.section")
  sub1
  
  p3 <-
    ggplot(sub1, aes(x = Process, y = Meanpercent, fill = Group)) +
    geom_col(position = position_dodge()) +
    geom_errorbar(aes(x = Process, ymin = Meanpercent - Stdevpercent, ymax = Meanpercent + Stdevpercent), position = position_dodge(), colour = "black", alpha = 0.9, size = 0.4) +
    scale_y_continuous(name = "Relative importance [%]", limits = c(0, 102)) +
    scale_x_discrete(name = "Community assembly process", labels = c("Heterogeneous.Selection"="Het. select.","Homogeneous.Selection"="Hom. select.","Dispersal.Limitation"="Disp. limit.","Drift.and.Others"="Drift","Homogenizing.Dispersal"="Hom. disp.")) +
    scale_fill_manual(name = "Root section", values = c(
      'darkslateblue','chocolate2'
    )) + 
    geom_text(aes(label = CLD, y = Meanpercent + Stdevpercent + 3),position = position_dodge(0.9),vjust = 0,size = 6) +
    ggtitle("Bacteria")+
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_blank(),
      plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
      plot.title = element_text(size=24, hjust=0.5, color="black", margin = margin(b=15)),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = ("none"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16)
    )
  
  newSTorder = c("Heterogeneous.Selection","Homogeneous.Selection","Dispersal.Limitation","Drift.and.Others","Homogenizing.Dispersal")
  p3$data$Process <- as.character(p3$data$Process)
  p3$data$Process <- factor(p3$data$Process, levels=newSTorder)
  
  p3



########## plot4
  sub1 <- subset(data,Organism=="Fungi" & Study=="I")
  sub1
  p4 <-
    ggplot(sub1, aes(x = Process, y = Meanpercent, fill = Group)) +
    geom_col(position = position_dodge()) +
    geom_errorbar(aes(x = Process, ymin = Meanpercent - Stdevpercent, ymax = Meanpercent + Stdevpercent), position = position_dodge(), colour = "black", alpha = 0.9, size = 0.4) +
    scale_y_continuous(name = "Relative importance [%]", limits = c(0, 102)) +
    scale_x_discrete(name = "Community assembly process", labels = c("Heterogeneous.Selection"="Het. select.","Homogeneous.Selection"="Hom. select.","Dispersal.Limitation"="Disp. limit.","Drift.and.Others"="Drift","Homogenizing.Dispersal"="Hom. disp.")) +
    scale_fill_manual(name = "Root section", values = c(
      'darkslateblue','chocolate2'
    )) + 
    geom_text(aes(label = CLD, y = Meanpercent + Stdevpercent + 3),position = position_dodge(0.9),vjust = 0,size = 6) +
    ggtitle("Fungi")+
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_blank(),
      plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
      plot.title = element_text(size=24, hjust=0.5, color="black", margin = margin(b=15)),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = ("none"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16)
    )
  
  newSTorder = c("Heterogeneous.Selection","Homogeneous.Selection","Dispersal.Limitation","Drift.and.Others","Homogenizing.Dispersal")
  p4$data$Process <- as.character(p4$data$Process)
  p4$data$Process <- factor(p4$data$Process, levels=newSTorder)
  
  p4



########## plot5
  sub1 <- subset(data,Organism=="Fungi" & Study=="II" & Factor=="Carbon_cat_4")
  sub1
  
  p5 <-
    ggplot(sub1, aes(x = Process, y = Meanpercent, fill = Group)) +
    geom_col(position = position_dodge()) +
    geom_errorbar(aes(x = Process, ymin = Meanpercent - Stdevpercent, ymax = Meanpercent + Stdevpercent), position = position_dodge(), colour = "black", alpha = 0.9, size = 0.4) +
    scale_y_continuous(name = "Relative importance [%]", limits = c(0, 102)) +
    scale_x_discrete(name = "Community assembly process", labels = c("Heterogeneous.Selection"="Het. select.","Homogeneous.Selection"="Hom. select.","Dispersal.Limitation"="Disp. limit.","Drift.and.Others"="Drift","Homogenizing.Dispersal"="Hom. disp.")) +
    scale_fill_manual(name = "13C label\ncategory", values=c("royalblue3" ,"#9C179EFF" ,"indianred1" ,"goldenrod1"
    )) + 
    geom_text(aes(label = CLD, y = Meanpercent + Stdevpercent + 3),position = position_dodge(0.9),vjust = 0,size = 6) +
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 10), size = 18, face = "bold"),
      plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
      plot.title = element_blank(),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = ("none"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16)
    )
  
  newSTorder = c("Heterogeneous.Selection","Homogeneous.Selection","Dispersal.Limitation","Drift.and.Others","Homogenizing.Dispersal")
  p5$data$Process <- as.character(p5$data$Process)
  p5$data$Process <- factor(p5$data$Process, levels=newSTorder)
  
  p5

  

########## plot6
  sub1 <- subset(data,Organism=="Fungi" & Study=="II" & Factor=="Root.section")
  sub1
  
  p6 <-
    ggplot(sub1, aes(x = Process, y = Meanpercent, fill = Group)) +
    geom_col(position = position_dodge()) +
    geom_errorbar(aes(x = Process, ymin = Meanpercent - Stdevpercent, ymax = Meanpercent + Stdevpercent), position = position_dodge(), colour = "black", alpha = 0.9, size = 0.4) +
    scale_y_continuous(name = "Relative importance [%]", limits = c(0, 102)) +
    scale_x_discrete(name = "Community assembly process", labels = c("Heterogeneous.Selection"="Het. select.","Homogeneous.Selection"="Hom. select.","Dispersal.Limitation"="Disp. limit.","Drift.and.Others"="Drift","Homogenizing.Dispersal"="Hom. disp.")) +
    scale_fill_manual(name = "Root section", values = c(
      'darkslateblue','chocolate2'
    )) + 
    geom_text(aes(label = CLD, y = Meanpercent + Stdevpercent + 3),position = position_dodge(0.9),vjust = 0,size = 6) +
    ggtitle("Fungi")+
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_blank(),
      plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
      plot.title = element_text(size=24, hjust=0.5, color="black", margin = margin(b=15)),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = ("none"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16)
    )
  
  newSTorder = c("Heterogeneous.Selection","Homogeneous.Selection","Dispersal.Limitation","Drift.and.Others","Homogenizing.Dispersal")
  p6$data$Process <- as.character(p6$data$Process)
  p6$data$Process <- factor(p6$data$Process, levels=newSTorder)
  
  p6


########## combine to create Figure S8
allplot <- plot_grid(
  p1,p4,p3,p6,p2,p5, 
  labels = c("A","","B", "","",""), ncol = 2, label_size = 30, scale=1
)
allplot


png("Figure S6.png", width = 380, height = 400, units = 'mm', res = 1000)
grid.arrange(allplot,ncol=1) 
  dev.off()

