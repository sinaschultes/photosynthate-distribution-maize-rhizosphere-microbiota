rm(list = ls())

########## load libraries
library(vegan)
library(dplyr)
library(tibble)
library(multcompView)
library(tidyr)
library(readxl)
library(ggh4x)
library(agricolae)
library(broom)
library(ggpubr)
library(car)
library(gridExtra)


  ########## load and prep data
  data <- read_excel("alpha diversity.xlsx")
  data <- as.data.frame(data)
  
  data$Organism <- as.factor(data$Organism)
  data$X11C_level <- as.factor(data$X11C_level)
  data$Copynr <- as.numeric(data$Copynr)
  
  data$logCopynr <- log10(data$Copynr)
  data <- subset(data, Root_section != "bulk")


  ########## exploratory tests
  sub1 <- subset(data, Organism == "Bacteria")
  sub2 <- subset(data, Organism == "Fungi")
  
  shap.test <- shapiro.test(sub1$logCopynr)
  shap.test
  shap.test <- shapiro.test(sub1$Shannon)
  shap.test
  shap.test <- shapiro.test(sub2$logCopynr)
  shap.test
  shap.test <- shapiro.test(sub2$Shannon)
  shap.test
  
  ggqqplot(sub1, x = "logCopynr")
  ggqqplot(sub2, x = "logCopynr")
  
  leveneTest(Shannon ~ X11C_level, data = sub1)
  leveneTest(Shannon ~ X11C_level, data = sub2)



  ########## plot bacteria
  sub1$X11C_level <- factor(sub1$X11C_level, levels = c("low", "medium", "high"))
  
  p <- ggscatter(sub1,
    x = "logCopynr", y = "Shannon",
    add = "reg.line", conf.int = TRUE,
    cor.coef = FALSE, cor.method = "kendall",
    color = "X11C_level", shape = "X11C_level",
    palette = c("darkslateblue", "chocolate2", "brown"),
    add.params = list(color = "gray40", fill = "lightgrey")
  ) +
    labs(
      x = "Copy number/g dry soil (log10)",
      y = "Shannon index",
      color = "11C level",
      shape = "11C level"
    ) +
    theme_linedraw() +
    annotate("text", x = 10.9, y = 5.0, label = "Pearson correlation\nr = -0.23, p = 2.9e-02", size = 5, hjust = 0) +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
      plot.title = element_blank(),
      panel.spacing.y = unit(1.5, "lines"),
      panel.spacing.x = unit(1.5, "lines"),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      strip.text = element_text(size = 18, color = "black"),
      strip.placement = "outside",
      strip.background = element_rect(size = 16, fill = "black", color = "grey90"),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16),
      legend.position = ("right")
    )
  
  p

  ########## stats bacteria
  corr <- cor.test(sub1$logCopynr, sub1$Shannon, method = "pearson")
  corr

  
  
  ########## plot fungi
  sub2$X11C_level <- factor(sub2$X11C_level, levels = c("low", "medium", "high"))
  
  p2 <- ggscatter(sub2,
    x = "logCopynr", y = "Shannon",
    add = "reg.line", conf.int = TRUE,
    cor.coef = FALSE, cor.method = "kendall",
    color = "X11C_level", shape = "X11C_level",
    palette = c("darkslateblue", "chocolate2", "brown"),
    add.params = list(color = "gray40", fill = "lightgrey")
  ) +
    labs(
      x = "Copy number/g dry soil (log10)",
      y = "Shannon index",
      color = "11C level",
      shape = "11C level"
    ) +
    theme_linedraw() +
    annotate("text", x = 7.2, y = 3.3, label = "Kendall correlation\ntau = -0.24, p = 9.9e-04", size = 5, hjust = 0) +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
      plot.title = element_blank(),
      panel.spacing.y = unit(1.5, "lines"),
      panel.spacing.x = unit(1.5, "lines"),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      strip.text = element_text(size = 18, color = "black"),
      strip.placement = "outside",
      strip.background = element_rect(size = 16, fill = "black", color = "grey90"),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16),
      legend.position = ("right")
    )
  
  p2
  
  ########## stats fungi
  corr <- cor.test(sub2$logCopynr, sub1$Shannon, method = "kendall")
  corr

  
########## combine plots to generate Figure S11
library(cowplot)
combiplot <- plot_grid(
  p, p2,
  labels = "AUTO", ncol = 1, label_size = 25
)
combiplot

png("Figure S11.png", width = 300, height = 320, units = "mm", res = 1000)
grid.arrange(combiplot, ncol = 1)
dev.off()
