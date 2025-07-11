rm(list = ls())

########## load libraries
  library(phyloseq)
  library(phyloseqCompanion)
  library(microbiome)
  library(cowplot)
  library(writexl)
  library(readxl)
  library(dplyr)
  library(qiime2R)
  library(plyr)
  library(foreach)
  library(doParallel)
  library(data.table)
  library(ggplot2)
  library(mirlyn)
  library(microViz)
  library(gridExtra)
  library(ggConvexHull)
  library(grid)
  library(tibble)
  library(ggpubr)


  ########## load and prep data
  ps_16S <- readRDS("study1_ps_16S.rds")
  ps_16S <- subset_samples(ps_16S, Root_section != "bulk")
  
  ps_ITS <- readRDS("study1_ps_ITS.rds")
  ps_ITS <- subset_samples(ps_ITS, Root_section != "bulk")
  
  ps_18S <- readRDS("study1_ps_18S.rds")
  ps_18S <- subset_samples(ps_18S, Root_section != "bulk")

  ########## rarefy to even depth
  set.seed(20)
  ps_r_16S <- ps_16S %>%
    rarefy_even_depth(sample.size = min(sample_sums(ps_16S)), rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
  
  set.seed(20)
  ps_r_ITS <- ps_ITS %>%
    rarefy_even_depth(sample.size = min(sample_sums(ps_ITS)), rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
  
  set.seed(20)
  ps_r_18S <- ps_18S %>%
    rarefy_even_depth(sample.size = min(sample_sums(ps_18S)), rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)



  ######################################## Figure S5
  
    ########## Bacteria root type
    ps1_16S <- transform_sample_counts(ps_r_16S, function(x) 1E6 * x / sum(x))
    ps1.ord.16S <- ordinate(ps1_16S, "NMDS", "bray",k=3)
    
    plot1 <- plot_ordination(ps1_16S, ps1.ord.16S, type = "samples", shape = "Root_section", color = "Root_type")
    
    p1 <- plot1 + geom_point(size = 4.0, stroke = 2.0) +
      scale_shape_manual(name = "Root section", values = c(1, 17)) +
      scale_color_manual(name = "Root type", values = c("brown", "chocolate2", "darkgoldenrod1", "springgreen4", "darkslateblue", "darkorchid3"), labels = c("prim" = "Primary", "sem" = "Seminal", "cr1" = "Crown 1", "cr2" = "Crown 2", "cr3" = "Crown 3", "cr4" = "Crown 4")) +
      geom_convexhull(alpha = 0, size = 0.3, show.legend = FALSE, inherit.aes = FALSE, aes(x = NMDS1, y = NMDS2, color = Root_type, fill = Root_type)) +
      theme_linedraw() +
      theme(
        axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15)),
        axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
        axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
        legend.position = ("none"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text.align = 0,
        legend.text = element_text(size = 16)
      )
    
    newSTorder <- c("prim", "sem", "cr1", "cr2", "cr3", "cr4")
    p1$data$Root_type <- as.character(p1$data$Root_type)
    p1$data$Root_type <- factor(p1$data$Root_type, levels = newSTorder)
    
    plot(p1)
    p1$layers <- p1$layers[-1]
    plot(p1)
  
  
  
    ########## Bacteria 11C
    plot2 <- plot_ordination(ps1_16S, ps1.ord.16S, type = "samples", shape = "Root_section", color = "X11C_level")
    
    p2 <- plot2 + geom_point(size = 4.0, stroke = 2.0) +
      scale_shape_manual(name = "Root section", values = c(1, 17, 17)) +
      scale_color_manual(name = "11C level", values = c("darkslateblue", "chocolate2", "brown")) +
      geom_convexhull(alpha = 0, size = 0.3, show.legend = FALSE, inherit.aes = FALSE, aes(x = NMDS1, y = NMDS2, color = X11C_level, fill = X11C_level)) +
      theme_linedraw() +
      theme(
        axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15)),
        axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
        axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
        legend.position = ("none"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text.align = 0,
        legend.text = element_text(size = 16)
      )
    newSTorder <- c("low", "medium", "high")
    p2$data$X11C_level <- as.character(p2$data$X11C_level)
    p2$data$X11C_level <- factor(p2$data$X11C_level, levels = newSTorder)
    
    plot(p2)
    p2$layers <- p2$layers[-1]
    plot(p2)
  
  
  
    ########## Fungi root type
    ps1_ITS <- transform_sample_counts(ps_r_ITS, function(x) 1E6 * x / sum(x))
    ps1.ord.ITS <- ordinate(ps1_ITS, "NMDS", "bray", k=3)
    
    plot3 <- plot_ordination(ps1_ITS, ps1.ord.ITS, type = "samples", shape = "Root_section", color = "Root_type")
    
    p3 <- plot3 + geom_point(size = 4.0, stroke = 2.0) +
      scale_shape_manual(name = "Root section", values = c(1, 17)) +
      scale_color_manual(name = "Root type", values = c("brown", "chocolate2", "darkgoldenrod1", "springgreen4", "darkslateblue", "darkorchid3"), labels = c("prim" = "Primary", "sem" = "Seminal", "cr1" = "Crown 1", "cr2" = "Crown 2", "cr3" = "Crown 3", "cr4" = "Crown 4")) +
      geom_convexhull(alpha = 0, size = 0.3, show.legend = FALSE, inherit.aes = FALSE, aes(x = NMDS1, y = NMDS2, color = Root_type, fill = Root_type)) +
      theme_linedraw() +
      theme(
        axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15)),
        axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
        axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
        legend.position = ("none"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text.align = 0,
        legend.text = element_text(size = 16)
      )
    
    newSTorder <- c("prim", "sem", "cr1", "cr2", "cr3", "cr4")
    p3$data$Root_type <- as.character(p3$data$Root_type)
    p3$data$Root_type <- factor(p3$data$Root_type, levels = newSTorder)
    
    plot(p3)
    p3$layers <- p3$layers[-1]
    plot(p3)
  
  
    ########## Fungi 11C
    plot4 <- plot_ordination(ps1_ITS, ps1.ord.ITS, type = "samples", shape = "Root_section", color = "X11C_level")
    
    p4 <- plot4 + geom_point(size = 4.0, stroke = 2.0) +
      scale_shape_manual(name = "Root section", values = c(1, 17, 17)) +
      scale_color_manual(name = "11C level", values = c("darkslateblue", "chocolate2", "brown")) +
      geom_convexhull(alpha = 0, size = 0.3, show.legend = FALSE, inherit.aes = FALSE, aes(x = NMDS1, y = NMDS2, color = X11C_level, fill = X11C_level)) +
      theme_linedraw() +
      theme(
        axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15)),
        axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
        axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
        legend.position = ("none"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text.align = 0,
        legend.text = element_text(size = 16)
      )
    
    newSTorder <- c("low", "medium", "high")
    p4$data$X11C_level <- as.character(p4$data$X11C_level)
    p4$data$X11C_level <- factor(p4$data$X11C_level, levels = newSTorder)
    
    plot(p4)
    p4$layers <- p4$layers[-1]
    plot(p4)
    
  
  
    ########## Cercozoa root type
    ps1_18S <- transform_sample_counts(ps_r_18S, function(x) 1E6 * x / sum(x))
    ps1.ord.18S <- ordinate(ps1_18S, "NMDS", "bray",k=3)
    
    plot5 <- plot_ordination(ps1_18S, ps1.ord.18S, type = "samples", shape = "Root_section", color = "Root_type")
    
    p5 <- plot5 + geom_point(size = 4.0, stroke = 2.0) +
      scale_shape_manual(name = "Root section", values = c(1, 17)) +
      scale_color_manual(name = "Root type", values = c("brown", "chocolate2", "darkgoldenrod1", "springgreen4", "darkslateblue", "darkorchid3"), labels = c("prim" = "Primary", "sem" = "Seminal", "cr1" = "Crown 1", "cr2" = "Crown 2", "cr3" = "Crown 3", "cr4" = "Crown 4")) +
      geom_convexhull(alpha = 0, size = 0.3, show.legend = FALSE, inherit.aes = FALSE, aes(x = NMDS1, y = NMDS2, color = Root_type, fill = Root_type)) +
      theme_linedraw() +
      theme(
        axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15)),
        axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
        axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
        legend.position = ("none"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text.align = 0,
        legend.text = element_text(size = 16)
      )
    
    newSTorder <- c("prim", "sem", "cr1", "cr2", "cr3", "cr4")
    p5$data$Root_type <- as.character(p5$data$Root_type)
    p5$data$Root_type <- factor(p5$data$Root_type, levels = newSTorder)
    
    plot(p5)
    p5$layers <- p5$layers[-1]
    plot(p5)
    
  
  
    ########## Cercozoa 11C
    plot6 <- plot_ordination(ps1_18S, ps1.ord.18S, type = "samples", shape = "Root_section", color = "X11C_level")
    
    p6 <- plot6 + geom_point(size = 4.0, stroke = 2.0) +
      scale_shape_manual(name = "Root section", values = c(1, 17, 17)) +
      scale_color_manual(name = "11C level", values = c("darkslateblue", "chocolate2", "brown")) +
      geom_convexhull(alpha = 0, size = 0.3, show.legend = FALSE, inherit.aes = FALSE, aes(x = NMDS1, y = NMDS2, color = X11C_level, fill = X11C_level)) +
      theme_linedraw() +
      theme(
        axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
        plot.title = element_text(size = 24, hjust = 0.5, color = "black", margin = margin(b = 15)),
        axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
        axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
        legend.position = ("none"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text.align = 0,
        legend.text = element_text(size = 16)
      )
  
      newSTorder <- c("low", "medium", "high")
      p6$data$X11C_level <- as.character(p6$data$X11C_level)
      p6$data$X11C_level <- factor(p6$data$X11C_level, levels = newSTorder)
      
      plot(p6)
      p6$layers <- p6$layers[-1]
      plot(p6)


  ########## combine plots to generate Figure S5
  combiplot <- plot_grid(
    p2, p1, p4, p3, p6, p5,
    ncol = 2, label_size = 25, scale = 1
  )
  combiplot
  
  png("Figure S5.png", width = 250, height = 320, units = "mm", res = 1000)
  grid.arrange(combiplot, ncol = 1)
  dev.off()
