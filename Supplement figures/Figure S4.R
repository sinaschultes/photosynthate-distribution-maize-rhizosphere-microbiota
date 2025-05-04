rm(list = ls())

# load libraries
library(vegan)
library(dplyr)
library(agricolae)
library(rcompanion)
library(car)
library(readxl)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggExtra)
library(FSA)

########## load data

data <- read.csv2("IRMS results.csv")

data$Location <- as.factor(data$Location)
data$Root_section <- as.factor(data$Root_section)
data$Root_type <- as.factor(data$Root_type)

data$massfrac_13C_rhizo_percent <- as.numeric(data$massfrac_13C_rhizo_percent)
data$massfrac_13C_root_percent <- as.numeric(data$massfrac_13C_root_percent)

head(data)


######################################## Figure S4
data <- subset(data, Label == "12C")

  ########## Fig. S4a
  a <- ggplot(
    data,
    aes(
      x = Location,
      y = massfrac_13C_rhizo_percent,
      fill = Location,
      shape = NA,
      group = Location
    )
  ) +
    geom_boxplot(outlier.size = -1) +
    scale_fill_manual(values = c("gray30", "brown", "brown", "chocolate2", "chocolate2", "darkgoldenrod2", "darkgoldenrod2", "springgreen4", "springgreen4", "darkslateblue", "darkslateblue", "darkorchid3", "darkorchid3")) +
    ggtitle("Root type") +
    scale_y_continuous(name = "13C in rhizosphere [%]") +
    scale_x_discrete(name = "Location") + 
    geom_jitter(position = position_jitter(0.2), shape = 1) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 4, color = "black", fill = "black") +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3)),
      axis.title.x = element_text(margin = margin(t = 3)),
      plot.title = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      strip.text = element_blank(),
      axis.title = element_text(size = 18, face = "bold"), 
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", angle = 40, hjust = 1, vjust = 1, margin = margin(t = 3)), 
      legend.position = ("none")
    )
  newSTorder <- c("Bulk soil", "Primary tip", "Primary base", "Seminal tip", "Seminal base", "Crown 1 tip", "Crown 1 base", "Crown 2 tip", "Crown 2 base", "Crown 3 tip")
  a$data$Location <- as.character(a$data$Location)
  a$data$Location <- factor(a$data$Location, levels = newSTorder)
  print(a)

  
  ########## Fig. S4b
  data2 <- subset(data, Location != "Bulk soil")
  
  b <- ggplot(
    data2, 
    aes(
      x = Location,
      y = massfrac_13C_root_percent,
      fill = Location,
      shape = NA, 
      group = Location
    )
  ) +
    geom_boxplot(outlier.size = -1) +
    scale_fill_manual(values = c("brown", "brown", "chocolate2", "chocolate2", "darkgoldenrod2", "darkgoldenrod2", "springgreen4", "springgreen4", "darkslateblue", "darkslateblue", "darkorchid3", "darkorchid3")) +
    ggtitle("Root type") +
    scale_y_continuous(name = "13C in root tissue [%]") +
    scale_x_discrete(name = "Location") + 
    geom_jitter(position = position_jitter(0.2), shape = 1) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 4, color = "black", fill = "black") +
    stat_boxplot(geom = "errorbar", width = 0.2) + 
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3)),
      axis.title.x = element_text(margin = margin(t = 3)),
      plot.title = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      strip.text = element_blank(),
      axis.title = element_text(size = 18, face = "bold"), 
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", angle = 40, hjust = 1, vjust = 1, margin = margin(t = 3)),
      legend.position = ("none")
    )
  newSTorder <- c("Primary tip", "Primary base", "Seminal tip", "Seminal base", "Crown 1 tip", "Crown 1 base", "Crown 2 tip", "Crown 2 base", "Crown 3 tip")
  b$data$Location <- as.character(b$data$Location)
  b$data$Location <- factor(b$data$Location, levels = newSTorder)
  print(b)

  
  
########## combine plots to generate Figure S4
library(cowplot)
allplot <- plot_grid(
  a, b,
  labels = "AUTO", ncol = 2, label_size = 25
)
allplot

png("Figure S4.png", width = 350, height = 150, units = "mm", res = 1000)
grid.arrange(allplot, ncol = 1)
dev.off()
