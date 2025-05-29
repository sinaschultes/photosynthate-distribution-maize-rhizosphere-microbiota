rm(list = ls())

########## load libraries
library(vegan)
library(dplyr)
library(tidyr)
library(readxl)
library(phyloseq)
library(microbiome)
library(phyloseqCompanion)



######################################## create the alpha diversity dataset

  ########## load and prep data
  ps_16S <- readRDS("study1_ps_16S.rds")
  ps_ITS <- readRDS("study1_ps_ITS.rds")
  ps_18S <- readRDS("study1_ps_18S.rds")


  ########## rarefy
  set.seed(20)
  ps_r_16S <- ps_16S %>%
    rarefy_even_depth(sample.size = min(sample_sums(ps_16S)), rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
  
  set.seed(20)
  ps_r_ITS <- ps_ITS %>%
    rarefy_even_depth(sample.size = min(sample_sums(ps_ITS)), rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
  
  set.seed(20)
  ps_r_18S <- ps_18S %>%
    rarefy_even_depth(sample.size = min(sample_sums(ps_18S)), rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


  ########## add column SampleID to metadata
  sample_data(ps_r_16S)$SampleID <- rownames(sample_data(ps_16S)) 
  sample_data(ps_r_ITS)$SampleID <- rownames(sample_data(ps_ITS)) 
  sample_data(ps_r_18S)$SampleID <- rownames(sample_data(ps_18S)) 

  
  ########## write the function
  process_phyloseq_alpha <- function(ps_obj, organism_label) {
    # extract metadata 
    meta <- sample.data.frame(ps_obj) 
   
    # estimate alpha diversity indices
    alpha <- estimate_richness(ps_obj)
    alpha$SampleID <- rownames(alpha)
    
    # pielou's evenness
    even <- evenness(abundances(ps_obj)[], 'pielou')
    even$SampleID <- rownames(even)
    
    # combine
    merged <- right_join(meta, alpha, by = "SampleID")
    merged2 <- right_join(merged, even, by="SampleID")
    
    # add organism info
    merged2$Organism <- organism_label
    return(merged2)
  }

  
  ########## run for each ps object
  df_16S <- process_phyloseq_alpha(ps_r_16S, "16S")
  df_ITS <- process_phyloseq_alpha(ps_r_ITS, "ITS")
  df_18S <- process_phyloseq_alpha(ps_r_18S, "18S")

  
  ########## remove useless cols from output
  df_16S <- df_16S[ , !(names(df_16S) %in% c("Sample_number","Plant", "Lab_code","Appeared")) ]
  df_ITS <- df_ITS[ , !(names(df_ITS) %in% c("Sample_number", "Plant","Lab_code","Appeared")) ]
  df_18S <- df_18S[ , !(names(df_18S) %in% c("Sample_number","Plant","Sample_ID")) ]

  
  ########## merge into one data frame
  alpha_combined <- bind_rows(df_16S, df_ITS, df_18S)
  
  
########## export results
writexl::write_xlsx(alpha_combined, "alpha diversity.xlsx")

