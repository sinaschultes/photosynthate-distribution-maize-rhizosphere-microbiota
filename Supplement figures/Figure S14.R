rm(list = ls())

########## load libraries
library(cowplot)
library(ggh4x)
library(writexl)
library(readxl)
library(dplyr)
library(microbiome)
library(rcompanion)
library(plyr)
library(data.table)
library(ggplot2)
library(mirlyn)
library(microViz)
library(gridExtra)
library(grid)
library(car)
library(tibble)
library(ggpubr)

########## load data (supplied)
DATA <- read_excel("./origin incorporators for stackedbar plot.xlsx")
DATA <- as.data.frame(DATA, na.rm = TRUE)
DATA1 <- subset(DATA, Group=="Bacteria")
DATA2 <- subset(DATA, Group=="Fungi")
DATA3 <- subset(DATA, Group=="Cercozoa")


######################################## Figure S14

  ########## plot bacteria
  myownpalette <- c("coral2","palevioletred","darkolivegreen","lightgoldenrod2",
                    "slateblue","plum3","darkgreen","powderblue","darksalmon",
                    "darkgoldenrod2","palegreen1","steelblue4",
                    "mediumorchid4","navajowhite","cornflowerblue","thistle3","aquamarine4","coral4",
                    "firebrick3", "grey80","tan3","lightpink3","cyan3","peachpuff3",
                    "salmon1", "midnightblue","darkolivegreen3","skyblue1")
  
  bac <- ggplot(DATA1, aes(x = Location, y = abundance, fill = Genus.x)) +
    geom_bar(width = 0.9, stat = "identity", colour = "black", position = "fill") +
    theme_classic() +
    guides(fill = guide_legend(nrow = 28)) +
    scale_fill_manual(name = "Bacterial incorporator genera", values = myownpalette) +
    scale_y_continuous(name = "Rel. abundance in 13C heavy fractions") +
    scale_x_discrete(name = "Location", labels = c(
      "prim_tip" = "Primary tip", "prim_base" = "Primary base", "sem_tip" = "Seminal tip", "sem_base" = "Seminal base",
      "cr1_tip" = "Crown 1 tip", "cr1_base" = "Crown 1 base", "cr2_tip" = "Crown 2 tip", "cr2_base" = "Crown 2 base", "cr3_tip" = "Crown 3 tip"
    )) +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_blank(),
      plot.margin = unit(c(0.5, 0, 0, 1.5), "cm"),
      plot.title = element_blank(),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", angle = 40, hjust = 1, vjust = 1, margin = margin(t = 3)),
      legend.position = ("right"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.location = "plot",
      legend.text = element_text(size = 16)
    )
  newSTorder <- c("prim_tip", "prim_base", "sem_tip", "sem_base", "cr1_tip", "cr1_base", "cr2_tip", "cr2_base", "cr3_tip")
  bac$data$Location <- as.character(bac$data$Location)
  bac$data$Location <- factor(bac$data$Location, levels = newSTorder)
  bac



  ########## plot fungi
  myownpalette <- c("coral2","palevioletred","cornflowerblue","lightgoldenrod2","powderblue",
                    "slateblue","salmon1","plum3","navajowhite",
                    "skyblue1","steelblue4",
                    "aquamarine4","coral4","darkolivegreen","thistle3","grey80",
                    "lightpink3", "cyan3", "peachpuff3","darkgoldenrod2","darkolivegreen3","mediumorchid4",
                    "tan3", "palegreen1","midnightblue","firebrick3","darkgreen","darksalmon")
  
  fung <- ggplot(DATA2, aes(x = Location, y = abundance, fill = Genus.x)) +
    geom_bar(width = 0.9, stat = "identity", colour = "black", position = "fill") +
    theme_classic() +
    guides(fill = guide_legend(nrow = 28)) +
    scale_fill_manual(name = "Fungal incorporator genera", values = myownpalette) +
    scale_y_continuous(name = "Rel. abundance in 13C heavy fractions") +
    scale_x_discrete(name = "Location", labels = c(
      "prim_tip" = "Primary tip", "prim_base" = "Primary base", "sem_tip" = "Seminal tip", "sem_base" = "Seminal base",
      "cr1_tip" = "Crown 1 tip", "cr1_base" = "Crown 1 base", "cr2_tip" = "Crown 2 tip", "cr2_base" = "Crown 2 base", "cr3_tip" = "Crown 3 tip"
    )) +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_blank(),
      plot.margin = unit(c(0.5, 0, 0, 1.5), "cm"),
      plot.title = element_blank(),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", angle = 40, hjust = 1, vjust = 1, margin = margin(t = 3)),
      legend.position = ("right"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.location = "plot",
      legend.text = element_text(size = 16)
    )
  newSTorder <- c("prim_tip", "prim_base", "sem_tip", "sem_base", "cr1_tip", "cr1_base", "cr2_tip", "cr2_base", "cr3_tip")
  fung$data$Location <- as.character(fung$data$Location)
  fung$data$Location <- factor(fung$data$Location, levels = newSTorder)
  fung


  ########## plot cercozoa
  myownpalette <- c("cornflowerblue","peachpuff3","powderblue","grey80","lightgoldenrod2",
                    "darkgoldenrod2","steelblue4","slateblue","salmon1","skyblue1","navajowhite","darkolivegreen",
                    "mediumorchid4",
                    "lightpink3","coral4","coral2","thistle3","plum3",
                    "aquamarine4", "cyan3", "palevioletred","salmon1","darkolivegreen3",
                    "tan3", "palegreen1","midnightblue","firebrick3","darkgreen","darksalmon")
  
  prot <- ggplot(DATA3, aes(x = Location, y = abundance, fill = Genus.x)) +
    geom_bar(width = 0.9, stat = "identity", colour = "black", position = "fill") +
    theme_classic() + 
    guides(fill = guide_legend(nrow = 28)) + 
    scale_fill_manual(name = "Protistan incorporator genera", values = myownpalette) + 
    scale_y_continuous(name = "Rel. abundance in 13C heavy fractions") +
    scale_x_discrete(name = "Location", labels = c(
      "prim_tip" = "Primary tip", "prim_base" = "Primary base", "sem_tip" = "Seminal tip", "sem_base" = "Seminal base",
      "cr1_tip" = "Crown 1 tip", "cr1_base" = "Crown 1 base", "cr2_tip" = "Crown 2 tip", "cr2_base" = "Crown 2 base", "cr3_tip" = "Crown 3 tip"
    )) +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
      plot.margin = unit(c(0.5, 0, 0, 1.5), "cm"),
      plot.title = element_blank(),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", angle = 40, hjust = 1, vjust = 1, margin = margin(t = 3)),
      legend.position = ("right"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.location = "plot",
      legend.text = element_text(size = 16)
    )
  newSTorder <- c("prim_tip", "prim_base", "sem_tip", "sem_base", "cr1_tip", "cr1_base", "cr2_tip", "cr2_base", "cr3_tip")
  prot$data$Location <- as.character(prot$data$Location)
  prot$data$Location <- factor(prot$data$Location, levels = newSTorder)
  prot


########## save Figure S14
allplot <- plot_grid(
  bac, fung, prot,
  labels = "AUTO", align = "v", axis = "lr", ncol = 1, label_size = 40
)
allplot

png("Figure S14.png", width = 450, height = 600, units = "mm", res = 1000)
grid.arrange(allplot, ncol = 1)
dev.off()
