rm(list = ls())

########## load libraries
library(dplyr)
library(HTSSIP)
library(DESeq2)
library(phyloseq)
library(phyloseqCompanion)
library(plyr)
library(writexl)
library(mirlyn)
library(foreach)
library(doParallel)
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(tibble)
library(ggpubr)

########## load studyII data of respective organism. run complete script once for each organism in a separate directory (16S, ITS, 18S).
ps <- readRDS("study2_ps_16S.RDS") #run bacteria
sample_data(ps)

ps <- readRDS("study2_ps_ITS.RDS") #run fungi
sample_data(ps)

ps <- readRDS("study2_ps_18S.RDS") #run cercozoa
sample_data(ps)

group_column <- "Group"
padj_cutoff <- 0.05
ncores <- 6
doParallel::registerDoParallel(ncores)

m <- sample_data(ps)
m


######################################## identify label incorporators in samples of diff. 13C_label_category (cat1 - cat4)
######################################## 

  ########## define the function
  run_HRSIP_analysis <- function(X13C_group, ps, m, padj_cutoff = 0.1, output_file = "output.xlsx") {
    # prune samples based on Label and 13C category
    physeq_group <- prune_samples(
      (m$Label == "unlabeled" & m$X13C_label_cat == X13C_group) |
        (m$Label == "labeled" & m$X13C_label_cat == X13C_group),
      ps
    )
    
    # pet Label levels
    sample_data(physeq_group)$Label <- factor(
      sample_data(physeq_group)$Label,
      levels = c("unlabeled", "labeled")
    )
    
    # check sample no per group
    print(table(sample_data(physeq_group)$Group))
    
    # define sparsity threshold
    thresh <- seq(0.05, 0.95, 0.05)
    
    # run HRSIP
    df_l2fc <- HRSIP(
      physeq_group,
      design = ~Label,
      density_windows = data.frame(density_min = 1.72, density_max = 1.75),
      l2fc_threshold = 0.138,
      padj_cutoff = padj_cutoff,
      padj_method = "BH",
      sparsity_threshold = thresh,
      sparsity_apply = "heavy"
    )
    
    # how many incorporators?
    n_incorp <- df_l2fc %>%
      filter(padj < padj_cutoff) %>%
      summarize(n_incorp_OTUs = n_distinct(OTU))
    
    print(n_incorp)
    
    # filter padj for sign
    df_filt <- df_l2fc %>%
      filter(padj < padj_cutoff)
    
    # export filtered results
    write_xlsx(df_filt, path = output_file)
    
    return(df_filt)
  }
  
  
  ########## run function for each 13C-level group
  df_cat1 <- run_HRSIP_analysis(
    X13C_group = "1",
    ps = ps,
    m = m,
    padj_cutoff = 0.05,
    output_file = "./cat1.xlsx"
  )
  df_cat1$X13C_label_cat <- "1"
  
  df_cat2 <- run_HRSIP_analysis(
    X13C_group = "2",
    ps = ps,
    m = m,
    padj_cutoff = 0.05,
    output_file = "./cat2.xlsx"
  )
  df_cat2$X13C_label_cat <- "2"
  
  df_cat3 <- run_HRSIP_analysis(
    X13C_group = "3",
    ps = ps,
    m = m,
    padj_cutoff = 0.05,
    output_file = "./cat3.xlsx"
  )
  df_cat3$X13C_label_cat <- "3"
  
  df_cat4 <- run_HRSIP_analysis(
    X13C_group = "4",
    ps = ps,
    m = m,
    padj_cutoff = 0.05,
    output_file = "./cat4.xlsx"
  )
  df_cat4$X13C_label_cat <- "4"
  
  
  ########## combine output
  df_combined <- bind_rows(df_cat1, df_cat2, df_cat3, df_cat4)
  df_combined <- as.data.frame(df_combined)
  df_combined
  str(df_combined)

  
  ########## combine l2fc with abundance data from phyloseq object
  phy <- phyloseq_to_df(ps)
  phy$X13C_label_cat <- factor(trimws(as.character(phy$X13C_label_cat)))#X13C_label_cat should be clean and a factor (phyloseq_to_df can mess this up sometimes)
  
  subphy <- subset(phy, Group == "labeled_heavy" & X13C_label_cat == "1")
  subphy$OTU <- subphy$Row.names
  subdata <- subset(df_combined, X13C_label_cat == "1")
  incorps_1 <- subphy %>% right_join(subdata, by = "OTU")
  
  subphy <- subset(phy, Group == "labeled_heavy" & X13C_label_cat == "2")
  subphy$OTU <- subphy$Row.names
  subdata <- subset(df_combined, X13C_label_cat == "2")
  incorps_2 <- subphy %>% right_join(subdata, by = "OTU")
  
  subphy <- subset(phy, Group == "labeled_heavy" & X13C_label_cat == "3")
  subphy$OTU <- subphy$Row.names
  subdata <- subset(df_combined, X13C_label_cat == "3")
  incorps_3 <- subphy %>% right_join(subdata, by = "OTU")
  
  subphy <- subset(phy, Group == "labeled_heavy" & X13C_label_cat == "4")
  subphy$OTU <- subphy$Row.names
  subdata <- subset(df_combined, X13C_label_cat == "4")
  incorps_4 <- subphy %>% right_join(subdata, by = "OTU")
  
  
  incorps <- rbind(incorps_1, incorps_2, incorps_3, incorps_4)
  incorps$Genus_ASV <- paste(incorps$Genus.x, incorps$Row.names, sep = " _ ")
  incorps
  

  ########## get number of incorporator ASVs
  unique(incorps$OTU)
  unique(incorps$Genus.x)
  unique(incorps$Phylum.x)
  
  sub <- subset(incorps, X13C_label_cat.x == "1")
  unique(sub$OTU)
  sub <- subset(incorps, X13C_label_cat.x == "2")
  unique(sub$OTU)
  sub <- subset(incorps, X13C_label_cat.x == "3")
  unique(sub$OTU)
  sub <- subset(incorps, X13C_label_cat.x == "4")
  unique(sub$OTU)
  
  
  ########## prepare data for enrichment plots (aggregate abundance and log2fc by genusxASV)
  incorps$Genus_ASV <- paste(incorps$Genus.x, incorps$Row.names, sep = " _ ")
  incorps$newcol <- incorps$abundance * incorps$log2FoldChange
  
  newdata <- aggregate(incorps$newcol, by = list(incorps$Genus_ASV, incorps$X13C_label_cat.x), FUN = sum)
  colnames(newdata) <- c("Genus_ASV", "X13C_label_cat", "newcol")
  newdata
  
  newdata2 <- aggregate(incorps$abundance, by = list(incorps$Genus_ASV, incorps$X13C_label_cat.x), FUN = sum)
  colnames(newdata2) <- c("Genus_ASV", "X13C_label_cat", "abundance")
  newdata2
  
  DATA <- cbind(newdata, newdata2)
  idx <- which(duplicated(names(DATA)))
  DATA <- DATA[, -idx]
  
  DATA$weightedlog2fc <- DATA$newcol / DATA$abundance

  DATA$Method <- "X13C_label_cat"
  DATA$Group <- "Bacteria"#adapt to respective organism
  DATA

  write_xlsx(DATA, "./X13C incorporators for enrichment plots.xlsx")








######################################## identify label incorporators in samples of diff. origin (seed-borne tips and bases,
######################################## shoot-borne tips and bases) in the root system

  ########## define the function
  run_HRSIP_analysis <- function(origin_group, ps, m, padj_cutoff = 0.1, output_file = "output.xlsx") {
    # prune samples based on Label and Origin
    physeq_group <- prune_samples(
      (m$Label == "unlabeled" & m$Origin == origin_group) |
        (m$Label == "labeled" & m$Origin == origin_group),
      ps
    )
  
    # pet Label levels
    sample_data(physeq_group)$Label <- factor(
      sample_data(physeq_group)$Label,
      levels = c("unlabeled", "labeled")
    )
  
    # check sample no per group
    print(table(sample_data(physeq_group)$Group))
  
    # define sparsity threshold
    thresh <- seq(0.05, 0.95, 0.05)
  
    # run HRSIP
    df_l2fc <- HRSIP(
      physeq_group,
      design = ~Label,
      density_windows = data.frame(density_min = 1.72, density_max = 1.75),
      l2fc_threshold = 0.138,
      padj_cutoff = padj_cutoff,
      padj_method = "BH",
      sparsity_threshold = thresh,
      sparsity_apply = "heavy"
    )
  
    # how many incorporators?
    n_incorp <- df_l2fc %>%
      filter(padj < padj_cutoff) %>%
      summarize(n_incorp_OTUs = n_distinct(OTU))
  
    print(n_incorp)
  
    # filter padj for sign
    df_filt <- df_l2fc %>%
      filter(padj < padj_cutoff)
  
    # export filtered results
    write_xlsx(df_filt, path = output_file)
  
    return(df_filt)
  }
  
  
  ########## run function for each origin group
  df_seed_base <- run_HRSIP_analysis(
    origin_group = "seedborne_base",
    ps = ps,
    m = m,
    padj_cutoff = 0.05,
    output_file = "./seedbase.xlsx"
  )
  df_seed_base$Origin <- "seedborne_base"
  
  df_seed_tip <- run_HRSIP_analysis(
    origin_group = "seedborne_tip",
    ps = ps,
    m = m,
    padj_cutoff = 0.05,
    output_file = "./seedtip.xlsx"
  )
  df_seed_tip$Origin <- "seedborne_tip"
  
  df_shoot_tip <- run_HRSIP_analysis(
    origin_group = "shootborne_tip",
    ps = ps,
    m = m,
    padj_cutoff = 0.05,
    output_file = "./shoottip.xlsx"
  )
  df_shoot_tip$Origin <- "shootborne_tip"
  
  df_shoot_base <- run_HRSIP_analysis(
    origin_group = "shootborne_base",
    ps = ps,
    m = m,
    padj_cutoff = 0.05,
    output_file = "./shootbase.xlsx"
  )
  df_shoot_base$Origin <- "shootborne_base"
  
  
  ########## combine output 
  df_combined <- bind_rows(df_seed_base, df_seed_tip, df_shoot_tip, df_shoot_base)
  df_combined <- as.data.frame(df_combined)
  df_combined
  str(df_combined)


  ########## combine l2fc with abundance data from phyloseq object
  phy <- phyloseq_to_df(ps)
  
  subphy <- subset(phy, Group == "labeled_heavy" & Origin == "seedborne_tip")
  subphy$OTU <- subphy$Row.names
  subdata <- subset(df_combined, Origin == "seedborne_tip")
  incorps_seedtip <- subphy %>% right_join(subdata, by = "OTU")
  
  subphy <- subset(phy, Group == "labeled_heavy" & Origin == "seedborne_base")
  subphy$OTU <- subphy$Row.names
  subdata <- subset(df_combined, Origin == "seedborne_base")
  incorps_seedbase <- subphy %>% right_join(subdata, by = "OTU")
  
  subphy <- subset(phy, Group == "labeled_heavy" & Origin == "shootborne_tip")
  subphy$OTU <- subphy$Row.names
  subdata <- subset(df_combined, Origin == "shootborne_tip")
  incorps_shoottip <- subphy %>% right_join(subdata, by = "OTU")
  
  subphy <- subset(phy, Group == "labeled_heavy" & Origin == "shootborne_base")
  subphy$OTU <- subphy$Row.names
  subdata <- subset(df_combined, Origin == "shootborne_base")
  incorps_shootbase <- subphy %>% right_join(subdata, by = "OTU")
  
  
  incorps <- rbind(incorps_seedtip, incorps_seedbase, incorps_shoottip, incorps_shootbase)
  incorps$Genus_ASV <- paste(incorps$Genus.x, incorps$Row.names, sep = " _ ")
  incorps
  
  
  ########## get number of incorporator ASVs
  unique(incorps$OTU)
  unique(incorps$Genus.x)
  unique(incorps$Phylum.x)
  
  sub <- subset(incorps, Origin.x == "seedborne_tip")
  unique(sub$OTU)
  sub <- subset(incorps, Origin.x == "seedborne_base")
  unique(sub$OTU)
  sub <- subset(incorps, Origin.x == "shootborne_tip")
  unique(sub$OTU)
  sub <- subset(incorps, Origin.x == "shootborne_base")
  unique(sub$OTU)


  ########## prepare data for enrichment plots (aggregate abundance and log2fc by genusxASV)
  incorps$Genus_ASV <- paste(incorps$Genus.x, incorps$Row.names, sep = " _ ")
  incorps$newcol <- incorps$abundance * incorps$log2FoldChange
  
  newdata <- aggregate(incorps$newcol, by = list(incorps$Genus_ASV, incorps$Origin.x), FUN = sum)
  colnames(newdata) <- c("Genus_ASV", "Origin.x", "newcol")
  newdata
  
  newdata2 <- aggregate(incorps$abundance, by = list(incorps$Genus_ASV, incorps$Origin.x), FUN = sum)
  colnames(newdata2) <- c("Genus_ASV", "Origin.x", "abundance")
  newdata2
  
  DATA <- cbind(newdata, newdata2)
  idx <- which(duplicated(names(DATA)))
  DATA <- DATA[, -idx]
  
  DATA$weightedlog2fc <- DATA$newcol / DATA$abundance
  
  DATA$Method <- "Origin"
  DATA$Group <- "Fungi"#adapt to respective organism
  DATA
  write_xlsx(DATA, "./origin incorporators for enrichment plots.xlsx")
  

  ########## prepare data for stacked barplot
  incorps$Genus_ASV <- paste(incorps$Genus.x, incorps$Row.names, sep = " _ ")
  incorps
  
  incorps$newcol <-incorps$abundance*incorps$log2FoldChange
  
  newdata <- aggregate(incorps$newcol, by=list(incorps$Genus_ASV,incorps$Genus.x,incorps$Origin.x, incorps$Location), FUN=sum)
  colnames(newdata) <- c("Genus_ASV", "Genus.x","Origin.x","Location","newcol")
  newdata
  
  newdata2 <- aggregate(incorps$abundance, by=list(incorps$Genus_ASV,incorps$Genus.x,incorps$Origin.x,incorps$Location), FUN=sum)
  colnames(newdata2) <- c("Genus_ASV", "Genus.x","Origin.x","Location","abundance")
  newdata2
  
  DATA <- cbind(newdata,newdata2)
  idx = which(duplicated(names(DATA)))
  DATA = DATA[,-idx]
  
  DATA$weightedlog2fc <- DATA$newcol/DATA$abundance 
  DATA
  write_xlsx(DATA,"./origin incorporators for stackedbar plots.xlsx")

