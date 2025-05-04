rm(list=ls())

########## load libraries
library(vegan)
library(readxl)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(scales)
library(grid)
library(gridExtra)


########## load and prep data (supplied)
data <- read.csv2("fraction selection example.csv") 
data

data$copy.number <- as.numeric(data$copy.number)
data$copy.number.norm <- as.numeric(data$copy.number.norm)
data$density <- as.numeric(data$density)
data$sample <- as.character(data$sample)
data$DNA.conc <- as.numeric(data$DNA.conc)


  ########## plot A total DNA concentration
  df <-subset(data, organism=="16S")
  df
  
  sup1 <- expression(''^"12"*'C-control sample            ')
  sup2 <- expression(''^"13"*'C-labeled sample')
  
  total <- ggplot(df, aes(x = density, y = DNA.conc, label = fraction)) + 
    geom_line(aes(color = label), linewidth=1) +
    geom_point(aes(color = label),size=3)+
    scale_color_manual(values = c("springgreen3","tomato3","dodgerblue1"),labels=c("12C"=sup1, "13C"=sup2)) + 
    scale_y_continuous(name="Total DNA [ng/µl]")+
    scale_x_continuous((name="Buoyant density [g/ml]"), limits=c(1.685,1.75),breaks = seq(1.69,1.75,0.01))+
    theme_linedraw() +                                         
    theme(axis.title.y=element_text(margin= margin(r=3)),
          plot.title=element_blank(),
          strip.text=element_blank(),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          axis.title = element_text(size=18, face="bold"),
          axis.text.y = element_text(size=16, color="black"),
          axis.text.x = element_text(size=16, angle=90, vjust=0.4,color="black"),
          legend.title = element_blank(),
          legend.position = "none",
          legend.text = element_text(size=16))
  total
  
  
  ########## plot B copy nr bacteria
  df <-subset(data, organism=="16S")
  df
  
  sup1 <- expression(''^"12"*'C-control sample            ')
  sup2 <- expression(''^"13"*'C-labeled sample')
  
  p16S <- ggplot(df, aes(x = density, y = copy.number.norm, label = fraction)) + 
    geom_line(aes(color = label), linewidth=1) +
    geom_point(aes(color=label),size=3)+
    scale_color_manual(values = c("springgreen3","tomato3","dodgerblue1"),labels=c("12C"=sup1, "13C"=sup2)) + 
    scale_y_continuous(name="Copy number 16S [norm.]")+
    scale_x_continuous((name="Buoyant density [g/ml]"), limits=c(1.685,1.75),breaks = seq(1.69,1.75,0.01))+
    theme_linedraw() +                                          
    theme(axis.title.y=element_text(margin= margin(r=3)),
          plot.title=element_blank(),
          strip.text=element_blank(),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          axis.title = element_text(size=18, face="bold"),
          axis.text.y = element_text(size=16, color="black"),
          axis.text.x = element_text(size=16, angle=90, vjust=0.4,color="black"),
          legend.title = element_blank(),
          legend.position = "none",
          legend.text = element_text(size=16))
  p16S
  
  
  
  ########## plot C copy nr fungi
  df <-subset(data, organism=="ITS")
  df
  
  sup1 <- expression(''^"12"*'C-control sample            ')
  sup2 <- expression(''^"13"*'C-labeled sample')
  
  pITS <- ggplot(df, aes(x = density, y = copy.number.norm, label = fraction)) + 
    geom_line(aes(color = label), linewidth=1) +
    geom_point(aes(color=label),size=3)+
    scale_color_manual(values = c("springgreen3","tomato3","dodgerblue1"),labels=c("12C"=sup1, "13C"=sup2)) + 
    scale_y_continuous(name="Copy number ITS [norm.]")+
    scale_x_continuous((name="Buoyant density [g/ml]"), limits=c(1.685,1.75),breaks = seq(1.69,1.75,0.01))+
    theme_linedraw() +                                         
    theme(axis.title.y=element_text(margin= margin(r=3)),
          plot.title=element_blank(),
          strip.text=element_blank(),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          axis.title = element_text(size=18, face="bold"),
          axis.text.y = element_text(size=16, color="black"),
          axis.text.x = element_text(size=16, angle=90, vjust=0.4,color="black"),
          legend.title = element_blank(),
          legend.position = "none",
          legend.text = element_text(size=16))
  pITS
  
  
  
  ########## plot D copy nr cercozoa
  df <-subset(data, organism=="18S")
  df
  
  sup1 <- expression(''^"12"*'C-control sample            ')
  sup2 <- expression(''^"13"*'C-labeled sample')
  
  p18S <- ggplot(df, aes(x = density, y = copy.number.norm, label = fraction)) + 
    geom_line(aes(color = label), linewidth=1) +
    geom_point(aes(color=label),size=3)+
    scale_color_manual(values = c("springgreen3","tomato3","dodgerblue1"),labels=c("12C"=sup1, "13C"=sup2)) + 
    scale_y_continuous(name="Copy number 18S [norm.]")+
    scale_x_continuous((name="Buoyant density [g/ml]"), limits=c(1.685,1.75),breaks = seq(1.69,1.75,0.01))+
    theme_linedraw() +                                        
    theme(axis.title.y=element_text(margin= margin(r=3)),
          plot.title=element_blank(),
          strip.text=element_blank(),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          axis.title = element_text(size=18, face="bold"),
          axis.text.y = element_text(size=16, color="black"),
          axis.text.x = element_text(size=16, angle=90, vjust=0.4,color="black"),
          legend.title = element_blank(),
          legend.position = "none",
          legend.text = element_text(size=16))
  p18S
    


########## export Figure S18
library(cowplot)
combiplot <- plot_grid(
  total, p16S, pITS, p18S,
  labels = "AUTO", ncol = 2, label_size = 25
)
combiplot

png("Figure S18.png", width = 450, height = 300, units = "mm", res = 1000)
grid.arrange(combiplot, ncol = 1)
dev.off()

