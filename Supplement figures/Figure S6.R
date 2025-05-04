rm(list = ls())

########## load libraries
library(dplyr)
library(phyloseq)
library(phyloseqCompanion)
library(microbiome)
library(plyr)
library(foreach)
library(data.table)
library(ggplot2)
library(mirlyn)
library(microViz)
library(gridExtra)
library(grid)
library(tibble)
library(phyloseqCompanion)
library(ggpubr)
library(vegan)
library(phyloseq)
library(vegan)
library(pairwiseAdonis)
library(tidyr)
library(dplyr)
library(reshape2)


########## data loading and data prep
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

  sample_data(ps_heavy_16S)$Mass_fraction_13C_rhizosphere <- as.numeric(sample_data(ps_heavy_16S)$Mass_fraction_13C_rhizosphere)
  sample_data(ps_heavy_ITS)$Mass_fraction_13C_rhizosphere <- as.numeric(sample_data(ps_heavy_ITS)$Mass_fraction_13C_rhizosphere)
  sample_data(ps_heavy_18S)$Mass_fraction_13C_rhizosphere <- as.numeric(sample_data(ps_heavy_18S)$Mass_fraction_13C_rhizosphere)
  
  
  ########## perform pairwise adonis
  
  #define the function
  veganotu <- function(ps_obj) {
    OTU <- otu_table(ps_obj)
    if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    return(as(OTU, "matrix"))
  }
  
  ps_objects <- list("ps_heavy_16S","ps_heavy_ITS","ps_heavy_18S")
  
  #loop through each ps_ object
  for (ps_name in ps_objects) {
    cat("Processing:", ps_name, "\n")
    
    ps_obj <- get(ps_name)
    
    #prepare OTU and Bray-Curtis distance
    vegtab <- veganotu(ps_obj)
    otu.dist <- vegdist(vegtab, method='bray') 
    sample_data_df <- data.frame(sample_data(ps_obj))
    
    #run pairwise adonis
    pairwise_result <- pairwise.adonis(otu.dist, 
                                       factors = sample_data_df$Root_type, 
                                       p.adjust.m = "bonferroni") 
    
    pairwise_result <- as.data.frame(pairwise_result)
    
    #format results
    pairwise_result <- pairwise_result %>%
      separate(pairs, into = c("Group1", "Group2"), sep = " vs ") %>%
      mutate(sig_star = case_when(
        p.adjusted <= 0.001 ~ "***",
        p.adjusted <= 0.01  ~ "**",
        p.adjusted <= 0.05  ~ "*",
        TRUE ~ ""
      ))
    
    #create R² and significance matrices
    groups <- unique(c(pairwise_result$Group1, pairwise_result$Group2))
    mat <- matrix(NA, nrow = length(groups), ncol = length(groups),
                  dimnames = list(groups, groups))
    sig_mat <- mat
    
    for (i in 1:nrow(pairwise_result)) {
      g1 <- pairwise_result$Group1[i]
      g2 <- pairwise_result$Group2[i]
      r2 <- round(pairwise_result$R2[i], 2)
      star <- pairwise_result$sig_star[i]
      
      mat[g1, g2] <- r2
      mat[g2, g1] <- r2
      sig_mat[g1, g2] <- star
      sig_mat[g2, g1] <- star
    }
    
    #melt matrices and merge
    mat_melt <- melt(mat, na.rm = TRUE)
    sig_melt <- melt(sig_mat, na.rm = TRUE)
    colnames(sig_melt)[3] <- "sig"
    plot_df <- left_join(mat_melt, sig_melt, by = c("Var1", "Var2"))
    plot_df$label <- paste0(plot_df$value, plot_df$sig)
    plot_df$vars <- paste(plot_df$Var1, plot_df$Var2)
    
    #subset to specific pairs
    newdf <- subset(plot_df, vars %in% c("prim sem", "prim cr1", "prim cr2", "prim cr3", 
                                         "sem cr1", "sem cr2", "sem cr3", 
                                         "cr1 cr2", "cr1 cr3", "cr2 cr3"))
    
    #rename groups for readability
    newdf$Var1 <- recode(newdf$Var1, 
                         "prim" = "Primary", "sem" = "Seminal", 
                         "cr1" = "Crown 1", "cr2" = "Crown 2", "cr3" = "Crown 3")
    newdf$Var2 <- recode(newdf$Var2, 
                         "prim" = "Primary", "sem" = "Seminal", 
                         "cr1" = "Crown 1", "cr2" = "Crown 2", "cr3" = "Crown 3")
    
    #save each dataframe
    out_name <- paste0("pairwise_df_", ps_name)
    assign(out_name, newdf, envir = .GlobalEnv)
  }
  
  ########## pairwise adonis results
  pairwise_df_ps_heavy_16S
  pairwise_df_ps_heavy_ITS
  pairwise_df_ps_heavy_18S
  
  
  ########## plotting
  p1 <- ggplot(pairwise_df_ps_heavy_16S, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), size = 4) +
    scale_fill_gradient2(low = "white", high = "steelblue4", mid = "lightblue",
                         midpoint = 0.15, limit = c(0, 0.33), 
                         name = "R²") +
    theme_minimal() +
    labs(x = "", y = "", title = "Bacteria")+
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
      plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
      plot.title = element_text(size=24, hjust=0.5, color="black", margin = margin(b=15)),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, angle=45, hjust=1, color = "black", margin = margin(t = 3)),
      legend.position = ("none"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16),
       panel.grid = element_blank()) 
    
  newSTorder = c("Primary", "Seminal", "Crown 1","Crown 2","Crown 3")
  p1$data$Var1 <- as.character(p1$data$Var1)
  p1$data$Var1 <- factor(p1$data$Var1, levels=newSTorder)
  
  newSTorder2 = c("Primary", "Seminal", "Crown 1","Crown 2","Crown 3")
  p1$data$Var2 <- as.character(p1$data$Var2)
  p1$data$Var2 <- factor(p1$data$Var2, levels=newSTorder2)
  
  p1
  
  
  
  p2 <- ggplot(pairwise_df_ps_heavy_ITS, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), size = 4) +
    scale_fill_gradient2(low = "white", high = "steelblue4", mid = "lightblue",
                         midpoint = 0.15, limit = c(0, 0.33), 
                         name = "R²") +
    theme_minimal() +
    labs(x = "", y = "", title = "Fungi")+
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
      plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
      plot.title = element_text(size=24, hjust=0.5, color="black", margin = margin(b=15)),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, angle=45, hjust=1, color = "black", margin = margin(t = 3)),
      legend.position = ("none"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16),
      panel.grid = element_blank()) 
  
  newSTorder = c("Primary", "Seminal", "Crown 1","Crown 2","Crown 3")
  p2$data$Var1 <- as.character(p2$data$Var1)
  p2$data$Var1 <- factor(p2$data$Var1, levels=newSTorder)
  
  newSTorder2 = c("Primary", "Seminal", "Crown 1","Crown 2","Crown 3")
  p2$data$Var2 <- as.character(p2$data$Var2)
  p2$data$Var2 <- factor(p2$data$Var2, levels=newSTorder2)
  
  p2
  
  
  
  p3 <- ggplot(pairwise_df_ps_heavy_18S, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), size = 4) +
    scale_fill_gradient2(low = "white", high = "steelblue4", mid = "lightblue",
                         midpoint = 0.15, limit =  c(0, 0.33), 
                         name = "R²") +
    theme_minimal() +
    labs(x = "", y = "", title = "Cercozoa")+
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 3), size = 18, face = "bold"),
      plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
      plot.title = element_text(size=24, hjust=0.5, color="black", margin = margin(b=15)),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, angle=45, hjust=1, color = "black", margin = margin(t = 3)),
      legend.position = ("right"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text.align = 0,
      legend.text = element_text(size = 16),
      panel.grid = element_blank()) 
  
  newSTorder = c("Primary", "Seminal", "Crown 1","Crown 2","Crown 3")
  p3$data$Var1 <- as.character(p3$data$Var1)
  p3$data$Var1 <- factor(p3$data$Var1, levels=newSTorder)
  
  newSTorder2 = c("Primary", "Seminal", "Crown 1","Crown 2","Crown 3")
  p3$data$Var2 <- as.character(p3$data$Var2)
  p3$data$Var2 <- factor(p3$data$Var2, levels=newSTorder2)
  
  p3
  
  
  
########## combine plots to generate Figure S6
library(cowplot)
  
combiplot <- plot_grid(
  p1,p2,p3, 
  labels = c("","",""), ncol = 2, label_size = 25, scale=1
)
combiplot


png("Figure S6.png", width = 250, height = 250, units = 'mm', res = 1000)
grid.arrange(combiplot,ncol=1) 
dev.off()

