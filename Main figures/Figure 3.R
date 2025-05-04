rm(list = ls())

########## load libraries
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

data <- subset(data, Label == "13C")
head(data)


######################################## Figure 3

  ########## Fig. 3a
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
    scale_y_continuous(name = "13C in rhizosphere [%]", limits = c(0, 25)) +
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
  
  
  ########## Fig. 3b
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
    scale_y_continuous(name = "13C in root tissue [%]", limits = c(0, 100)) +
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

  

  ########## Stats for Fig. 3a and 3b
  shap.test <- shapiro.test(data$massfrac_13C_rhizo_percent)
  shap.test
  shap.test <- shapiro.test(data2$massfrac_13C_root_percent)
  shap.test
  
  ggqqplot(data, x = "massfrac_13C_rhizo_percent")
  ggqqplot(data2, x = "massfrac_13C_root_percent")
  
  leveneTest(massfrac_13C_rhizo_percent ~ Location, data = data) 
  leveneTest(massfrac_13C_root_percent ~ Location, data = data2) 
  
  # ANOVA + posthoc 3a
  type.aov <- aov(massfrac_13C_rhizo_percent ~ Location, data = data)
  summary(type.aov) 
  TukeyHSD(type.aov)$Location
  Tukey <- HSD.test(type.aov, "Location", group = TRUE, console = TRUE)#CLD
  
  # ANOVA + posthoc 3b
  type.aov <- aov(massfrac_13C_root_percent ~ Location, data = data2)
  summary(type.aov) 
  TukeyHSD(type.aov)
  Tukey <- HSD.test(type.aov, "Location", group = TRUE, console = TRUE)#CLD
  

  
  ########## Fig. 3c
  model <- lm(massfrac_13C_rhizo_percent ~ massfrac_13C_root_percent, data = data2)
  par(mfrow = c(2, 2))
  plot(model)
  
  coef(lm(data2$massfrac_13C_rhizo_percent ~ data2$massfrac_13C_root_percent))[2]#get slope
  
  c <- ggplot(data2, aes(x = massfrac_13C_root_percent, y = massfrac_13C_rhizo_percent, color = Root_type)) +
    geom_point(aes(size = 3, shape = Root_section, color = Root_type)) +
    geom_rug() +
    geom_smooth(
      method = lm,
      se = TRUE, fullrange = FALSE, color = "black"
    ) +
    geom_label(label = "all samples", x = 31, y = 12.5, label.padding = unit(0.3, "lines"), size = 5, label.size = 0.5, color = "white", fill = "black") +
    scale_y_continuous(name = "13C in rhizosphere [%]") +
    scale_x_continuous(name = "13C in root tissue [%]") +
    scale_color_manual(name = "Root type", values = c("brown", "chocolate2", "darkgoldenrod2", "springgreen4", "darkslateblue", "darkorchid3")) +
    scale_fill_manual(values = c("grey70")) +
    scale_shape_manual(values = c(19, 17)) +
    ggpubr::stat_cor(size = 6, color = "black") +
    guides(size = "none") +
    guides(color = "none") +
    guides(shape = "none") +
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
      legend.position = ("right"),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16)
    )
  newSTorder <- c("Primary", "Seminal", "Crown 1", "Crown 2", "Crown 3")
  c$data$Root_type <- as.character(c$data$Root_type)
  c$data$Root_type <- factor(c$data$Root_type, levels = newSTorder)
  c


  ########## Fig. 3d
  datatips <- subset(data2, Root_section == "tip")
  database <- subset(data2, Root_section == "base")
  coef(lm(datatips$massfrac_13C_rhizo_percent ~ datatips$massfrac_13C_root_percent))[2]#get slope
  coef(lm(database$massfrac_13C_rhizo_percent ~ database$massfrac_13C_root_percent))[2]#get slope
  
  d <- ggplot(data2, aes(x = massfrac_13C_root_percent, y = massfrac_13C_rhizo_percent, color = Root_section)) +
    geom_point(aes(color = Root_section, shape = Root_section, size = 3)) +
    geom_rug(aes(color = Root_section)) +
    geom_smooth(aes(fill = Root_section),
      method = lm,
      se = TRUE, fullrange = FALSE
    ) +
    geom_label(label = "root tips", x = 31, y = 17, label.padding = unit(0.3, "lines"), size = 5, label.size = 0.5, color = "white", fill = "maroon") +
    geom_label(label = "root bases", x = 12, y = 7, label.padding = unit(0.3, "lines"), size = 5, label.size = 0.5, color = "white", fill = "royalblue4") +
    scale_y_continuous(name = "13C in rhizosphere [%]") +
    scale_x_continuous(name = "13C in root tissue [%]") +
    scale_color_manual(name = "Root section", values = c("royalblue4", "maroon")) +
    scale_fill_manual(name = "Root section", values = c("grey70", "grey70")) +
    scale_shape_manual(name = "Root section", values = c(19, 17)) +
    ggpubr::stat_cor(aes(color = Root_section), label.x = 3, size = 6) +
    guides(size = "none") +
    guides(fill = "none") +
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
      legend.position = ("none"),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16)
    )
  d



########## combine plots to generate Figure 3
library(cowplot)
combiplot <- plot_grid(
  b, a, c, d,
  labels = "AUTO", ncol = 2, label_size = 25
)
combiplot

png("Figure 3.png", width = 350, height = 300, units = "mm", res = 1000)
grid.arrange(combiplot, ncol = 1)
dev.off()
