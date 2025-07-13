rm(list = ls())

########## load libraries
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(data.table)
  library(scales)

  ########## load and prep data
  DATA <- read.delim("MRI root length.txt")
  DATA <- as.data.table(DATA)
  DATA2 <- DATA %>%
    group_by(Sample_ID, Sample_name)

  newSTorder <- c("prim", "sem", "cr1", "cr2", "cr3", "cr4")
  DATA2$Sample_name <- factor(DATA2$Sample_name, levels = newSTorder)
  DATA2$Sample_ID <- factor(DATA2$Sample_ID, levels = c("0", "6", "13", "20"))


  ########## plot
  rootlength_plot <- ggplot(DATA2, aes(
    x = Sample_ID,
    y = Mean,
    group = Sample_name,
    colour = Sample_name
  )) +
    geom_errorbar(
      aes(ymin = Mean - SD, ymax = Mean + SD),
      width = 0.3,
      size = 1,
     ) +
    geom_line(size = 1, linetype = "longdash") +
    geom_point(size = 2.5, position = position_dodge()) +
    scale_color_manual(
      name = "Root type",
      values = c(
        "prim" = "brown",
        "sem" = "chocolate2",
        "cr1" = "darkgoldenrod2",
        "cr2" = "springgreen4",
        "cr3" = "darkslateblue",
        "cr4" = "darkorchid3"
      ),
      labels = c(
        "prim" = "Primary",
        "sem" = "Seminal",
        "cr1" = "Crown 1",
        "cr2" = "Crown 2",
        "cr3" = "Crown 3",
        "cr4" = "Crown 4"
      )
    ) +
    scale_x_discrete(labels = c("0" = "day 0", "6" = "day 6", "13" = "day 13", "20" = "day 20")) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab("Day after sowing") +
    ylab("Mean root length in mm") +
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
      plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
      plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15)),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = "right",
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
    )
  
  rootlength_plot

########## create Figure S17
png("Figure S17.png", width = 190, height = 115, units = "mm", res = 1000)
grid.arrange(rootlength_plot, ncol = 1)
dev.off()
