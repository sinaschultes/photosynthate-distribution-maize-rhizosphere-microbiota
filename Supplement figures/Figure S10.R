rm(list = ls())

########## load libraries
library(vegan)
library(purrr)
library(dplyr)
library(tidyr)
library(readxl)
library(ggh4x)
library(dplyr)
library(agricolae)
library(broom)
library(ggpubr)
library(car)
library(tidyverse)
library(gridExtra)
library(writexl)


########## define function to extract p values from posthoc tests
get.pvalues <- function(pw_result) {
  pvals <- pw_result$p.value
  out <- as.vector(pvals)
  names(out) <- apply(expand.grid(rownames(pvals), colnames(pvals)), 1, function(x) paste(sort(x), collapse = "-"))
  out <- out[!is.na(out)]
  return(out)
}


  ########## load and prep data
  data <- read_excel("alpha diversity.xlsx")
  data <- as.data.frame(data)
  
  data$Organism <- as.factor(data$Organism)
  data$Root_section <- as.factor(data$Root_section)
  data$X11C_level <- as.factor(data$X11C_level)
  data$Root_type <- as.factor(data$Root_type)
  data$Copynr <- as.numeric(data$Copynr)
  data$logCopynr <- log10(data$Copynr)
  data
  
  
  ########## subset and pivot to long format
  sub1 <- subset(data, Organism != "Cercozoa")
  
  data_long <- sub1 %>%
    pivot_longer(
      cols = c(X11C_level, Root_type, Root_section),
      names_to = "GroupVar",
      values_to = "Level"
    )
  data_long$GroupVar <- factor(data_long$GroupVar, levels = c("X11C_level", "Root_type", "Root_section"))
  
  head(data_long)
  str(data_long)
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(agricolae)
  library(tibble)
  
  ########## define function to perform statistical tests
  run_anova_tukey <- function(sub_df) {
    response <- "logCopynr"
    group <- "Level"
  
    formula <- as.formula(paste(response, "~", group))
    aov_model <- aov(formula, data = sub_df)
  
    aov_summary <- summary(aov_model)[[1]]
    anova_p_value <- aov_summary["Pr(>F)"][1, 1]
    anova_df <- aov_summary["Df"][1, 1]
    anova_resid <- aov_summary["Df"][2, 1]
    anova_f <- aov_summary["F value"][1, 1]
    anova_result <- sprintf("F(%d,%d) = %.2f, p = %.2e", anova_df, anova_resid, anova_f, anova_p_value)
  
    tukey <- agricolae::HSD.test(aov_model, group, group = TRUE)
    cld <- tukey$groups %>%
      rownames_to_column(var = group) %>%
      mutate(
        Organism = unique(sub_df$Organism),
        GroupVar = unique(sub_df$GroupVar),
        Anova_result = anova_result,
        Anova_P_value = anova_p_value
      )
    return(cld)
  }
  
  ########## extract results
  tukey_results <- data_long %>%
    group_by(Organism, GroupVar) %>%
    group_modify(~ run_anova_tukey(.x)) %>%
    ungroup()
  
  ########## attach results to plotting dataset
  plot_annotations <- tukey_results %>%
    rename(Level = Level, CL_label = groups) %>%
    select(Organism, GroupVar, Level, CL_label, Anova_result)
  
  data_plot <- left_join(data_long, plot_annotations, by = c("Organism", "GroupVar", "Level"))
  
  ########## label data for CLD letters
  label_data_clean <- data_plot %>%
    group_by(Organism, GroupVar, Level, CL_label) %>%
    summarise(y_pos = max(logCopynr, na.rm = TRUE) * 1.01, .groups = "drop")
  
  ########## label data for ANOVA p values
  label_data_pval <- data_plot %>%
    group_by(Organism, GroupVar) %>%
    summarise(
      p_value_label = first(Anova_result),
      y_pos = max(logCopynr, na.rm = TRUE) * 1.09,
      .groups = "drop"
    )
  
  ########## annotation positions
  label_data_pval <- as.data.frame(label_data_pval)
  label_data_pval$x_pos <- c(0.5, 0.6, 0.5, 0.5, 0.6, 0.5)
  label_data_pval
  label_data_clean$GroupVar <- factor(label_data_clean$GroupVar, levels = c("X11C_level", "Root_type", "Root_section"))
  
  
  ########## plot
  y <- ggplot(
    data_long,
    aes(
      x = Level,
      y = logCopynr,
      fill = Level,
      shape = NA
    )
  ) +
    geom_boxplot(outlier.size = -1) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 4, color = "black", fill = "black") +
    scale_fill_manual(values = c(
      "bulk" = "grey40", "low" = "darkslateblue", "medium" = "chocolate2", "high" = "brown", "base" = "darkslateblue", "tip" = "chocolate2",
      "prim" = "brown", "sem" = "chocolate2", "cr1" = "darkgoldenrod2", "cr2" = "springgreen4", "cr3" = "darkslateblue", "cr4" = "darkorchid3"
    )) +
    scale_y_continuous(name = "Copy number/g dry soil (log10)", expand = expansion(mult = c(0.05, 0.1))) +
    geom_jitter(position = position_jitter(0.2), shape = 1, size = 3) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    theme_linedraw() +
    geom_text(data = label_data_clean, aes(x = Level, y = y_pos, label = CL_label), vjust = -0.8, size = 6, color = "black", na.rm = TRUE) +
    geom_text(data = label_data_pval, aes(x = x_pos, y = y_pos, label = p_value_label), inherit.aes = FALSE, hjust = 0, size = 6) +
    facet_grid2(GroupVar ~ Organism, scales = "free", independent = "all", switch = "y", labeller = labeller(GroupVar = c("X11C_level" = "11C level", "Root_type" = "Root type", "Root_section" = "Root section"))) +
    theme_linedraw() +
    theme(
      axis.title.y = element_text(margin = margin(r = 3), size = 18, face = "bold"),
      axis.title.x = element_blank(),
      plot.title = element_blank(),
      panel.spacing.y = unit(1.5, "lines"),
      panel.spacing.x = unit(1.5, "lines"),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      strip.text = element_text(size = 18, color = "black"),
      strip.placement = "outside",
      strip.background = element_rect(size = 16, fill = "black", color = "grey90"),
      axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 3)),
      axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 3)),
      legend.position = ("none")
    )
  newSTorder <- c("bulk", "base", "tip", "prim", "sem", "cr1", "cr2", "cr3", "cr4", "low", "medium", "high")
  y$data$Level <- as.character(y$data$Level)
  y$data$Level <- factor(y$data$Level, levels = newSTorder)
  
  print(y)

########## generate Figure S10
png("Figure S10.png", width = 300, height = 300, units = "mm", res = 1000)
grid.arrange(y, ncol = 1)
dev.off()
