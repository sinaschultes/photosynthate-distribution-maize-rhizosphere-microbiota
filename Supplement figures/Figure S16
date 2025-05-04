rm(list=ls())

########## load libraries
library(vegan)
library(readxl)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(lubridate)
library(scales)


  ########## load and prep data (supplied)
  data <- read.csv2("gasflow 13CO2.csv") 
  data
  
  df <- data %>%
    select(Time, ch1.12CO2, ch2.12CO2, ch2.13CO2) %>%
    gather(key = "variable", value = "value", -Time)
  df$value <- as.numeric(df$value)
  
  df$Time <- hms(df$Time)

  ########## plot
  sup1 <- expression(''^"12"*'CO'['2']*'  chamber I    ')
  sup2 <- expression(''^"12"*'CO'['2']*' chamber II    ')
  sup3 <- expression(''^"13"*'CO'['2']*' chamber II')
  
  g <- ggplot(df, aes(x = Time, y = value, group=variable)) + 
    geom_line(aes(color = variable), linewidth=.8) + 
    scale_color_manual(values = c("dodgerblue1", "springgreen3","tomato3"),labels=c("ch1.12CO2"=sup1, "ch2.12CO2"=sup2,"ch2.13CO2"=sup3,"H1"="rH [%] ch. 1", "H2"= "rH [%] ch. 2", "T1"="T [°C] ch. 1", "T2"="T [°C] ch. 2")) + 
    scale_y_continuous(name="Gas concentration [ppm]",limits=c(0,500), expand=c(0,0), minor_breaks = seq(0 , 500, 25))+
    scale_x_time(name="A labelling day [time]",breaks = breaks_width("120 min"),minor_breaks = breaks_width("60 min"))+
    theme_linedraw() +                                          
    theme(axis.title.y=element_text(margin= margin(r=3)),
          plot.title=element_blank(),
          strip.text=element_blank(),
          axis.title = element_text(size=18, face="bold"),
          axis.text.y = element_text(size=16, color="black"),
          axis.text.x = element_text(size=16, angle=90, vjust=0.4,color="black"),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size=16))
  g
  
########## export Figure S16
png("Figure S16.png", width = 200, height = 170, units = 'mm', res = 1000)
grid.arrange(g,ncol=1) 
dev.off()  
