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


########## data loading and prep
  #load data bacteria
  ps_16S <- readRDS("study2_ps_16S.rds")

  #load data fungi
  ps_ITS <- readRDS("study2_ps_ITS.rds")

  #load data cercozoa
  ps_18S <- readRDS("study2_ps_18S.rds")

  #rarefy
  set.seed(20)
  ps_r_16S <- ps_16S %>%
    rarefy_even_depth(sample.size = min(sample_sums(ps_16S)), rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
  
  set.seed(20)
  ps_r_ITS <- ps_ITS %>%
    rarefy_even_depth(sample.size = min(sample_sums(ps_ITS)), rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
  
  set.seed(20)
  ps_r_18S <- ps_18S %>%
    rarefy_even_depth(sample.size = min(sample_sums(ps_18S)), rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

  #add column SampleID to metadata
  sample_data(ps_r_16S)$SampleID <- rownames(sample_data(ps_16S)) 
  sample_data(ps_r_ITS)$SampleID <- rownames(sample_data(ps_ITS)) 
  sample_data(ps_r_18S)$SampleID <- rownames(sample_data(ps_18S)) 
  
  #make 13C-heavy subsets
  ps_heavy_16S <- subset_samples(ps_r_16S, Group == "labeled_heavy")
  ps_heavy_ITS <- subset_samples(ps_r_ITS, Group == "labeled_heavy")
  ps_heavy_18S <- subset_samples(ps_r_18S, Group == "labeled_heavy")



######################################## Figure 4
  
  ########## A Bacteria
  sample_data(ps_heavy_16S)$Mass_fraction_13C_rhizosphere <- as.numeric(sample_data(ps_heavy_16S)$Mass_fraction_13C_rhizosphere)
  ps1_16S <- transform_sample_counts(ps_heavy_16S, function(x) 1E6 * x / sum(x)) 
  get_taxa_unique(ps1_16S, "Phylum")
  ps1.ord.16S <- ordinate(ps1_16S, "NMDS", "bray") 
  
  plot1 <- plot_ordination(ps1_16S, ps1.ord.16S, type = "samples", shape = "X11C_level", color = "Root_type")
  
  p1 <- plot1 + geom_point(size = 4.0, stroke = 2.0) +
    scale_shape_manual(name="PET-based\n11C level",values = c(1, 19, 17)) + 
    scale_color_manual(name="Root type",values = c("brown", "chocolate2", "darkgoldenrod1", "springgreen4", "darkslateblue"),labels=c("prim"="Primary", "sem"="Seminal","cr1"="Crown 1","cr2"="Crown 2","cr3"="Crown 3")) +
    geom_convexhull(alpha = 0, size = 0.3, show.legend = FALSE, inherit.aes = FALSE, aes(x = NMDS1, y = NMDS2, color = Root_type, fill = Root_type)) + 
    #annotate("text", x = 0.75, y = -0.95, label = paste("Stress =", round(ps1.ord.16S$stress, 3))) +
    ggtitle("Bacteria")+
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
      plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
      plot.title = element_text(size=24, hjust=0.5, color="black", margin = margin(b=15)),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = ("none"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16)
    )
  
  newSTorder <- c("prim", "sem", "cr1", "cr2", "cr3")
  p1$data$Root_type <- as.character(p1$data$Root_type)
  p1$data$Root_type <- factor(p1$data$Root_type, levels = newSTorder)
  
  newSTorder2 <- c("low", "medium", "high")
  p1$data$X11C_level <- as.character(p1$data$X11C_level)
  p1$data$X11C_level <- factor(p1$data$X11C_level, levels = newSTorder2)
  
  plot(p1)
  p1$layers <- p1$layers[-1]
  plot(p1)



  ########## B Bacteria
  plot2 <- plot_ordination(ps1_16S, ps1.ord.16S, type = "samples", shape = "Root_type", color = "Mass_fraction_13C_rhizosphere") 
  
  p2 <- plot2 + geom_point(size = 4.0, stroke = 2.0) +
    scale_shape_manual(name="Root type",values = c(1, 3, 10, 19, 15),labels=c("prim"="Primary", "sem"="Seminal","cr1"="Crown 1","cr2"="Crown 2","cr3"="Crown 3")) +  
    scale_color_gradient(name="Mass fraction of 13C\nin rhizosphere",low = "blue", high = "red") + 
    #annotate("text", x = 0.75, y = -0.95, label = paste("Stress =", round(ps1.ord.16S$stress, 3))) +
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
      plot.margin = unit(c(0.8,0.5,0.5,0.5), "cm"),
      plot.title = element_blank(),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = ("none"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16)
    )
  newSTorder <- c("prim", "sem", "cr1", "cr2", "cr3")
  p2$data$Root_type <- as.character(p2$data$Root_type)
  p2$data$Root_type <- factor(p2$data$Root_type, levels = newSTorder)
  
  plot(p2)
  p2$layers <- p2$layers[-1]
  plot(p2)

  
  
  ########## A Fungi
  sample_data(ps_heavy_ITS)$Mass_fraction_13C_rhizosphere <- as.numeric(sample_data(ps_heavy_ITS)$Mass_fraction_13C_rhizosphere)
  ps1_ITS <- transform_sample_counts(ps_heavy_ITS, function(x) 1E6 * x / sum(x)) 
  get_taxa_unique(ps1_ITS, "Phylum") 
  ps1.ord.ITS <- ordinate(ps1_ITS, "NMDS", "bray") 
  
  plot3 <- plot_ordination(ps1_ITS, ps1.ord.ITS, type = "samples", shape = "X11C_level", color = "Root_type") 
  
  p3 <- plot3 + geom_point(size = 4.0, stroke = 2.0) +
    scale_shape_manual(name="PET-based\n11C level",values = c(1, 19, 17)) + 
    scale_color_manual(name="Root type",values = c("brown", "chocolate2", "darkgoldenrod1", "springgreen4", "darkslateblue"),labels=c("prim"="Primary", "sem"="Seminal","cr1"="Crown 1","cr2"="Crown 2","cr3"="Crown 3")) +
    geom_convexhull(alpha = 0, size = 0.3, show.legend = FALSE, inherit.aes = FALSE, aes(x = NMDS1, y = NMDS2, color = Root_type, fill = Root_type)) + 
    #annotate("text", x = 0.75, y = -0.95, label = paste("Stress =", round(ps1.ord.ITS$stress, 3))) +
    ggtitle("Fungi")+
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
      plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
      plot.title = element_text(size=24,hjust=0.5, color="black", margin = margin(b=15)),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = ("none"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16)
    )
  
  newSTorder <- c("prim", "sem", "cr1", "cr2", "cr3")
  p3$data$Root_type <- as.character(p3$data$Root_type)
  p3$data$Root_type <- factor(p3$data$Root_type, levels = newSTorder)
  
  newSTorder2 <- c("low", "medium", "high")
  p3$data$X11C_level <- as.character(p3$data$X11C_level)
  p3$data$X11C_level <- factor(p3$data$X11C_level, levels = newSTorder2)
  
  plot(p3)
  p3$layers <- p3$layers[-1]
  plot(p3)



  ########## B Fungi
  plot4 <- plot_ordination(ps1_ITS, ps1.ord.ITS, type = "samples", shape = "Root_type", color = "Mass_fraction_13C_rhizosphere") 
  
  p4 <- plot4 + geom_point(size = 4.0, stroke = 2.0) +
    scale_shape_manual(name="Root type",values = c(1, 3, 10, 19, 15),labels=c("prim"="Primary", "sem"="Seminal","cr1"="Crown 1","cr2"="Crown 2","cr3"="Crown 3")) +  
    scale_color_gradient(name="Mass fraction\n[%] of 13C\n in rhizosphere",low = "blue", high = "red") + 
    #annotate("text", x = 0.75, y = -0.95, label = paste("Stress =", round(ps1.ord$stress, 3))) +
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
      plot.margin = unit(c(0.8,0.5,0.5,0.5), "cm"),
      plot.title = element_blank(),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = ("none"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16)
    )
  
  newSTorder <- c("prim", "sem", "cr1", "cr2", "cr3")
  p4$data$Root_type <- as.character(p4$data$Root_type)
  p4$data$Root_type <- factor(p4$data$Root_type, levels = newSTorder)
  
  plot(p4)
  p4$layers <- p4$layers[-1]
  plot(p4)



  ########## A Cercozoa
  sample_data(ps_heavy_18S)$Mass_fraction_13C_rhizosphere <- as.numeric(sample_data(ps_heavy_18S)$Mass_fraction_13C_rhizosphere)
  ps1_18S <- transform_sample_counts(ps_heavy_18S, function(x) 1E6 * x / sum(x)) 
  get_taxa_unique(ps1_18S, "Phylum") 
  ps1.ord.18S <- ordinate(ps1_18S, "NMDS", "bray") 
  
  plot5 <- plot_ordination(ps1_18S, ps1.ord.18S, type = "samples", shape = "X11C_level", color = "Root_type") 
  
  p5 <- plot5 + geom_point(size = 4.0, stroke = 2.0) +
    scale_shape_manual(name="PET-based\n11C level",values = c(1, 19, 17)) +
    scale_color_manual(name="Root type",values = c("brown", "chocolate2", "darkgoldenrod1", "springgreen4", "darkslateblue"),labels=c("prim"="Primary", "sem"="Seminal","cr1"="Crown 1","cr2"="Crown 2","cr3"="Crown 3")) +
    geom_convexhull(alpha = 0, size = 0.3, show.legend = FALSE, inherit.aes = FALSE, aes(x = NMDS1, y = NMDS2, color = Root_type, fill = Root_type)) + 
    #annotate("text", x = 0.75, y = -0.95, label = paste("Stress =", round(ps1.ord.18S$stress, 3))) +
    ggtitle("Cercozoa")+
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
      plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
      plot.title = element_text(size=24,hjust=0.5, color="black", margin = margin(b=15)),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = ("none"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16)
    )
  
  newSTorder <- c("prim", "sem", "cr1", "cr2", "cr3")
  p5$data$Root_type <- as.character(p5$data$Root_type)
  p5$data$Root_type <- factor(p5$data$Root_type, levels = newSTorder)
  
  newSTorder2 <- c("low", "medium", "high")
  p5$data$X11C_level <- as.character(p5$data$X11C_level)
  p5$data$X11C_level <- factor(p5$data$X11C_level, levels = newSTorder2)
  
  plot(p5)
  p5$layers <- p5$layers[-1]
  plot(p5)
  
  
  
  ########## B Cercozoa
  plot6 <- plot_ordination(ps1_18S, ps1.ord.18S, type = "samples", shape = "Root_type", color = "Mass_fraction_13C_rhizosphere") 
  
  p6 <- plot6 + geom_point(size = 4.0, stroke = 2.0) +
    scale_shape_manual(name="Root type",values = c(1, 3, 10, 19, 15),labels=c("prim"="Primary", "sem"="Seminal","cr1"="Crown 1","cr2"="Crown 2","cr3"="Crown 3")) +  
    scale_color_gradient(name="Mass fraction of 13C\nin rhizosphere",low = "blue", high = "red") + 
    #annotate("text", x = 0.75, y = -0.95, label = paste("Stress =", round(ps1.ord$stress, 3))) +
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
      plot.margin = unit(c(0.8,0.5,0.5,0.5), "cm"),
      plot.title = element_blank(),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = ("none"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16)
    )
  
  newSTorder <- c("prim", "sem", "cr1", "cr2", "cr3")
  p6$data$Root_type <- as.character(p6$data$Root_type)
  p6$data$Root_type <- factor(p6$data$Root_type, levels = newSTorder)
  
  plot(p6)
  p6$layers <- p6$layers[-1]
  plot(p6)


########## combine plots to generate Figure 4
combiplot <- plot_grid(
  p1,p3,p5,p2,p4,p6, 
  labels = c("A","","", "B","",""), ncol = 3, label_size = 25, scale=1
)
combiplot

png("Figure 4.png", width = 400, height = 250, units = 'mm', res = 1000)
grid.arrange(combiplot,ncol=1) 
dev.off()

