rm(list = ls())

########## load libraries
library(vegan)
library(purrr)
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

get.pvalues <- function(pw_result) {
  pvals <- pw_result$p.value
  out <- as.vector(pvals)
  names(out) <- apply(expand.grid(rownames(pvals), colnames(pvals)), 1, function(x) paste(sort(x), collapse = "-"))
  out <- out[!is.na(out)]
  return(out)
}


######################################## stats and plot

  ########## load and prep data
  data <- read_excel("alpha diversity.xlsx")
  data <- as.data.frame(data)
  
  data$Organism <- as.factor(data$Organism)
  data$Root_section <- as.factor(data$Root_section)
  data$X11C_level <- as.factor(data$X11C_level)
  data$Root_type <- as.factor(data$Root_type)
  
  ########## subset and pivot to long format
  sub1 <- subset(data, Organism == "Cercozoa")
  
  data_long <- sub1 %>%
    pivot_longer(cols = c(Observed, Pielou), names_to = "Metric", values_to = "Value") %>%
    pivot_longer(
      cols = c(X11C_level, Root_type, Root_section),
      names_to = "Factor", values_to = "Level"
    )
  
  
  ########## exploratory tests
  shap.test <- shapiro.test(sub1$Observed)
  shap.test
  shap.test <- shapiro.test(sub1$Pielou)
  shap.test
  
  ggqqplot(sub1, x = "Observed")
  ggqqplot(sub1, x = "Pielou")
  
  leveneTest(Observed ~ Root_type, data = sub1)
  leveneTest(Observed ~ Root_section, data = sub1)
  leveneTest(Observed ~ X11C_level, data = sub1)
  leveneTest(Pielou ~ Root_type, data = sub1)
  leveneTest(Pielou ~ Root_section, data = sub1)
  leveneTest(Pielou ~ X11C_level, data = sub1)
  
  
  ########## make variable lists
  response_vars <- c("Observed", "Pielou")
  group_vars <- c("X11C_level", "Root_type", "Root_section")
  df <- sub1
  
  
  ########## write function to extract p-values and CLDs
  run_kruskal_wilcox <- function(response, group, df) {
    formula <- as.formula(paste(response, "~", group))
    kruskal <- kruskal.test(formula, data = df)
  
    #get p-values
    kruskal_stat <- kruskal$statistic
    kruskal_df <- kruskal$parameter
    kruskal_p_value <- kruskal$p.value
    kruskal_result <- sprintf(
      "Chisq(%d) = %.2f, p = %.2e",
      kruskal_df, kruskal_stat, kruskal_p_value
    )
  
    #get CLDs
    pwc <- pairwise.wilcox.test(df[[response]], df[[group]], p.adjust.method = "BH")
    out.p <- get.pvalues(pwc)
    pwc_letters <- multcompLetters(out.p, threshold = 0.05, reversed = FALSE)$Letters

    cld <- tibble::tibble(
      !!group := names(pwc_letters),
      groups = unname(pwc_letters),
      Response = response,
      Grouping = group,
      Kruskal_stat = as.numeric(kruskal_stat),
      Kruskal_df = as.numeric(kruskal_df),
      Kruskal_p_value = kruskal_p_value,
      Kruskal_result = kruskal_result
    )
  
    return(cld)
  }


  ########## results
  wilcox_results <- cross_df(list(response = response_vars, group = group_vars)) %>%
    mutate(results = map2(response, group, ~ run_kruskal_wilcox(.x, .y, df))) %>%
    unnest(results)
  wilcox_results <- as.data.frame(wilcox_results)

  ########## summarize results for plotting
  plot_annotations <- wilcox_results %>%
    dplyr::rename(CL_label = groups) %>%
    mutate(
      Metric = response,
      Factor = Grouping,
      Level = coalesce(X11C_level, Root_type, Root_section),
    ) %>%
    select(Metric, Factor, Level, CL_label, Kruskal_result) %>%
    distinct()
 
   data_plot <- left_join(data_long, plot_annotations,
    by = c("Metric", "Factor", "Level")
  )
  data_plot <- as.data.frame(data_plot)
  
  label_data_clean <- data_plot %>%
    select(Level, Metric, Factor, CL_label, Value) %>%
    distinct(Level, Metric, Factor, CL_label) %>%
    left_join(
      data_plot %>%
        group_by(Level, Metric, Factor) %>%
        summarise(y_pos = max(Value, na.rm = TRUE) * 1.05, .groups = "drop"),
      by = c("Level", "Metric", "Factor")
    )
  
  label_data_pval <- data_plot %>%
    select(Level, Metric, Factor, Kruskal_result, Value) %>%
    group_by(Metric, Factor) %>%
    summarise(
      p_value_label = shQuote(first(Kruskal_result)),
      y_pos = max(Value, na.rm = TRUE) * 1.2,
      .groups = "drop"
    )
  label_data_pval <- as.data.frame(label_data_pval)
  label_data_pval$x_pos <- c(0.5, 0.6, 0.5, 0.5, 0.6, 0.5)

  ########## remove CLDs from non-sign facets
  label_data_clean <- label_data_clean %>%
    mutate(CL_label = ifelse(Factor == "X11C_level" & Metric == "Observed", NA, CL_label))
  label_data_clean <- label_data_clean %>%
    mutate(CL_label = ifelse(Factor == "Root_section" & Metric == "Observed", NA, CL_label))
  label_data_clean <- label_data_clean %>%
    mutate(CL_label = ifelse(Factor == "Root_section" & Metric == "Pielou", NA, CL_label))
  
  ########## reassign letters so letter of first level is a
  label_data_clean <- label_data_clean %>%
    mutate(CL_label = case_when(
      Factor == "X11C_level" & Metric == "Pielou" & Level == "medium" ~ "a",
      Factor == "X11C_level" & Metric == "Pielou" & Level == "high" ~ "b",
      TRUE ~ CL_label
    ))
  

  ########## plot
  strip <- strip_themed(background_x = elem_list_rect(fill = c("grey90")))
  
  newSTorder2 <- c("X11C_level", "Root_type", "Root_section")
  data_plot$Factor <- factor(data_plot$Factor, levels = newSTorder2)
  label_data_clean$Factor <- factor(label_data_clean$Factor, levels = newSTorder2)
  label_data_pval$Factor <- factor(label_data_pval$Factor, levels = newSTorder2)


  y <- ggplot(
    data_plot,
    aes(
      x = Level,
      y = Value,
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
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    geom_jitter(position = position_jitter(0.2), shape = 1, size = 3) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    theme_linedraw() +
    geom_text(data = label_data_clean, aes(x = Level, y = y_pos, label = CL_label), vjust = -0.8, size = 6, color = "black", na.rm = TRUE) +
    geom_text(data = label_data_pval, aes(x = x_pos, y = y_pos, label = p_value_label), inherit.aes = FALSE, parse = TRUE, hjust = 0, size = 6) +
    facet_grid2(Factor ~ Metric, scales = "free", independent = "all", switch = "y", labeller = labeller(Factor = c("X11C_level" = "11C level", "Root_type" = "Root type", "Root_section" = "Root section"),Metric = c("Observed" = "ASV richness", "Pielou" = "Pielou evenness")))+
    theme_linedraw() +
    theme(
      axis.title.y = element_blank(),
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
  
png("Figure S9d.png", width = 300, height = 300, units = "mm", res = 1000)
grid.arrange(y, ncol = 1)
dev.off()
