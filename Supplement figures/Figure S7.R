rm(list=ls())

########## load libraries
  library(dplyr) 
  library(ggplot2)
  library(data.table)
  library(scales)
  library(writexl)
  library(readxl)
  library(dplyr)
  library(plyr)
  library(foreach)
  library(data.table)
  library(ggplot2)
  library(mirlyn)
  library(tidyr)
  library(microViz)
  library(gridExtra)
  library(grid)
  library(tibble)
  library(ggpubr)

         
########## load and prep data root type
  DATA <- read_xlsx("stamp output root type.xlsx")
  
  DATA <- as.data.frame(DATA)
  DATA$genus <- as.character(DATA$genus)
  DATA[, 8:19] <- lapply(DATA[, 8:19], as.numeric)
  
  heatmap_data <- DATA %>%
    select(Organism,genus, matches("mean_rel_freq_\\(perc\\)"))
  heatmap_data
  
  heatmap_long <- heatmap_data %>%
    pivot_longer(
      cols = -c(Organism,genus),
      names_to = "Root_type",
      values_to = "MeanRelFreq"
    )
  heatmap_long
  
  heatmap_long <- heatmap_long %>%
    mutate(Root_type = gsub("_mean_rel_freq_\\(perc\\)", "", Root_type),
           Root_type = gsub("_", " ", Root_type)) 
  heatmap_long

  ########## plot root type
  mine.heatmap.r <- ggplot(data = heatmap_long, mapping = aes(x = Root_type,
                                                         y = reorder(genus,-MeanRelFreq),
                                                         fill = MeanRelFreq)) +
    geom_tile() +
    theme(panel.background = element_blank())+
    xlab(label = "Sample_ID")+
    facet_wrap(~Organism, ncol=1, scales="free_y")+
    scale_fill_gradient2(name = "rel. abundance (%)", trans="log10", labels = scales::label_number(accuracy = 0.01, big.mark = ","),na.value=muted("dodgerblue3"),
                         low = muted("dodgerblue3"),
                         mid = "lightgoldenrod1",
                         high = "red",midpoint=0.01,limits = c(0.008, 65)) +
    scale_x_discrete(name = "Root type")+
    scale_y_discrete(name="Responsive genera >0.25% total rel. abn.", position="right")+
    theme(axis.text.y = element_text(face="italic"))
  newSTorder = c("primary","seminal", "crown 1","crown 2","crown 3","crown 4")
  mine.heatmap.r$data$Root_type <- as.character(mine.heatmap.r$data$Root_type)
  mine.heatmap.r$data$Root_type <- factor(mine.heatmap.r$data$Root_type, levels=newSTorder)
  
  mine.heatmap.r


  ########## load and prep data 11C level
  DATA <- read_xlsx("stamp output 11Clevel.xlsx")
  
  DATA <- as.data.frame(DATA)
  DATA$genus <- as.character(DATA$genus)
  head(DATA)
  DATA[, 8:13] <- lapply(DATA[, 8:13], as.numeric)
  
  heatmap_data <- DATA %>%
    select(Organism,genus, matches("mean_rel_freq_\\(perc\\)"))
  heatmap_data
  
  heatmap_long <- heatmap_data %>%
    pivot_longer(
      cols = -c(Organism,genus),
      names_to = "X11C_level",
      values_to = "MeanRelFreq"
    )
  heatmap_long
  
  heatmap_long <- heatmap_long %>%
    mutate(X11C_level = gsub("_mean_rel_freq_\\(perc\\)", "", X11C_level),
           X11C_level = gsub("_", " ", X11C_level)) 
  heatmap_long

  ########## plot root type
  mine.heatmap.c <- ggplot(data = heatmap_long, mapping = aes(x = X11C_level,
                                                              y = reorder(genus,-MeanRelFreq),
                                                              fill = MeanRelFreq)) +
    geom_tile() +
    theme(panel.background = element_blank())+
    xlab(label = "Sample_ID")+
    facet_wrap(~Organism, ncol=1, scales="free_y")+
    scale_fill_gradient2(name = "rel. abundance (%)", trans="log10", labels = scales::label_number(accuracy = 0.01, big.mark = ","),na.value=muted("dodgerblue3"),
                         low = muted("dodgerblue3"),
                         mid = "lightgoldenrod1",
                         high = "red",midpoint=0.01,limits = c(0.008, 65)) +
    scale_x_discrete(name = "Root type")+
    scale_y_discrete(name="Responsive genera >0.25% total rel. abn.", position="left")+
    theme(axis.text.y = element_text(face="italic"))
  newSTorder = c("low","medium", "high")
  mine.heatmap.c$data$X11C_level <- as.character(mine.heatmap.c$data$X11C_level)
  mine.heatmap.c$data$X11C_level <- factor(mine.heatmap.c$data$X11C_level, levels=newSTorder)
  
  mine.heatmap.c

########## combine and export
png("Figure S7.png", width = 360, height = 140, units = 'mm', res = 600)
grid.arrange(mine.heatmap.c, mine.heatmap.r, ncol=2, nrow=1) 
dev.off()
