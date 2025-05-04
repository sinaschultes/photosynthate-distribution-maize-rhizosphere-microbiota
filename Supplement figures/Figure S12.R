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

data <- subset(data, Label == "13C" & Location != "Bulk soil")

head(data)


######################################## Figure S11

  #make a plot to assign groups based on 13C mass fraction in rhizosphere samples
  head(data)
  data$Location <- as.factor(data$Location)
  data$massfrac_13C_rhizo_percent <- as.numeric(data$massfrac_13C_rhizo_percent)
  data$ID <- 1:36
  data

  p <- ggplot(data, aes(ID, massfrac_13C_rhizo_percent, color = factor(Location), shape = factor(Root_section))) +
    geom_hline(yintercept = 6.9, linetype = "solid", color = "cornflowerblue", linewidth = 1.5) + 
    geom_label(label = "  Category 1 (1.5% - 6.9%)  ", x = 47, y = 4, label.padding = unit(1, "lines"), label.size = 2, color = "white", fill = "cornflowerblue") +
    geom_hline(yintercept = 12.3, linetype = "solid", color = "cornflowerblue", linewidth = 1.5) + 
    geom_label(label = " Category 2 (6.9% - 12.3%) ", x = 47, y = 9, label.padding = unit(1, "lines"), label.size = 2, color = "white", fill = "cornflowerblue") +
    geom_hline(yintercept = 17.7, linetype = "solid", color = "cornflowerblue", linewidth = 1.5) + 
    geom_label(label = "Category 3 (12.3% - 17.7%)", x = 47, y = 14, label.padding = unit(1, "lines"), label.size = 2, color = "white", fill = "cornflowerblue") +
    geom_label(label = "Category 4 (17.7% - 23.1%)", x = 47, y = 19, label.padding = unit(1, "lines"), label.size = 2, color = "white", fill = "cornflowerblue") +
    geom_point(aes(colour = factor(Location)), size = 3) +
    scale_color_manual(name = "Location", values = c("brown", "brown", "chocolate2", "chocolate2", "darkgoldenrod2", "darkgoldenrod2", "springgreen4", "springgreen4", "darkslateblue")) +
    scale_shape_manual(name = "Root section", values = c(19, 17)) +
    scale_x_continuous(name = "Sample number", limits = c(0, 50)) +
    scale_y_continuous(name = "13C in rhizosphere [%]") +
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3)),
      axis.title.x = element_text(margin = margin(t = 3)),
      plot.title = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      strip.text = element_blank(),
      axis.title = element_text(size = 18, face = "bold"), 
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = ("right")
    )
  
  newSTorder <- c("Primary tip", "Primary base", "Seminal tip", "Seminal base", "Crown 1 tip", "Crown 1 base", "Crown 2 tip", "Crown 2 base", "Crown 3 tip")
  p$data$Location <- as.character(p$data$Location)
  p$data$Location <- factor(p$data$Location, levels = newSTorder)
  p
  
  png("Figure S11.png", width = 325, height = 200, units = "mm", res = 1000)
  grid.arrange(p, ncol = 1)
  dev.off()

