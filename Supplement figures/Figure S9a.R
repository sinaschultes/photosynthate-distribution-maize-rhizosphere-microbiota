rm(list = ls())

########## load libraries
  library(multcompView)
  library(tibble)
  library(vegan)
  library(purrr)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(ggh4x)
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


  ########## exploratory tests
  sub1 <- subset(data, Organism == "Bacteria")
  sub2 <- subset(data, Organism == "Fungi")
  sub3 <- subset(data, Organism == "Cercozoa")
  
  shap.test <- shapiro.test(sub1$Shannon)
  shap.test
  shap.test <- shapiro.test(sub2$Shannon)
  shap.test
  shap.test <- shapiro.test(sub2$Shannon)
  shap.test
  
  ggqqplot(sub1, x = "Shannon")
  ggqqplot(sub2, x = "Shannon")
  ggqqplot(sub3, x = "Shannon")
  
  leveneTest(Shannon ~ Root_type, data = sub1)
  leveneTest(Shannon ~ Root_section, data = sub1)
  leveneTest(Shannon ~ X11C_level, data = sub1)
  leveneTest(Shannon ~ Root_type, data = sub2)
  leveneTest(Shannon ~ Root_section, data = sub2)
  leveneTest(Shannon ~ X11C_level, data = sub2)
  leveneTest(Shannon ~ Root_type, data = sub3)
  leveneTest(Shannon ~ Root_section, data = sub3)
  leveneTest(Shannon ~ X11C_level, data = sub3)

  ########## pivot to long format
  data_long <- data %>%
    pivot_longer(
      cols = c(X11C_level, Root_type, Root_section),
      names_to = "GroupVar",
      values_to = "Level"
    )
  data_long$GroupVar <- factor(data_long$GroupVar, levels = c("X11C_level", "Root_type", "Root_section"))
  data_long$Organism <- factor(data_long$Organism, levels = c("Bacteria", "Fungi", "Cercozoa"))
  head(data_long)
  data_long <- as.data.frame(data_long)


  ########## define function to preform statistical tests
  run_stat_test <- function(df, keys) {
    org <- keys$Organism
    groupvar <- keys$GroupVar
  
    if (length(org) != 1) stop("Multiple Organisms found")
  
    formula <- Shannon ~ Level
  
    if (org == "Bacteria") {
      aov_model <- aov(formula, data = df)
      aov_summary <- summary(aov_model)[[1]]
      p_val <- aov_summary$`Pr(>F)`[1]
  
      if (p_val > 0.05) {
        cld <- tibble::tibble(
          Level = unique(df$Level),
          CL_label = NA,
          Organism = org,
          GroupVar = groupvar,
          Anova_result = sprintf(
            "F(%d,%d) = %.2f, p = %.2e",
            aov_summary$Df[1],
            aov_summary$Df[2],
            aov_summary$`F value`[1],
            p_val
          ),
          Kruskal_result = NA_character_
        )
      } else {
        tukey <- agricolae::HSD.test(aov_model, "Level", group = TRUE)
        cld <- tukey$groups %>%
          tibble::rownames_to_column(var = "Level") %>%
          dplyr::rename(CL_label = groups) %>%
          mutate(
            Organism = org,
            GroupVar = groupvar,
            Anova_result = sprintf(
              "F(%d,%d) = %.2f, p = %.2e",
              aov_summary$Df[1],
              aov_summary$Df[2],
              aov_summary$`F value`[1],
              p_val
            ),
            Kruskal_result = NA_character_
          )
      }
    } else {
      kruskal <- kruskal.test(formula, data = df)
      p_val <- kruskal$p.value
  
      if (p_val > 0.05) {
        cld <- tibble::tibble(
          Level = unique(df$Level),
          CL_label = NA,
          Organism = org,
          GroupVar = groupvar,
          Anova_result = NA_character_,
          Kruskal_result = sprintf(
            "Chisq(%d) = %.2f, p = %.2e",
            kruskal$parameter,
            kruskal$statistic,
            p_val
          )
        )
      } else {
        pwc <- pairwise.wilcox.test(df$Shannon, df$Level, p.adjust.method = "BH")
        out.p <- get.pvalues(pwc)
        pwc_letters <- multcompLetters(out.p, threshold = 0.05, reversed = FALSE)$Letters
  
        cld <- tibble::tibble(
          Level = names(pwc_letters),
          CL_label = unname(pwc_letters),
          Organism = org,
          GroupVar = groupvar,
          Anova_result = NA_character_,
          Kruskal_result = sprintf(
            "Chisq(%d) = %.2f, p = %.2e",
            kruskal$parameter,
            kruskal$statistic,
            p_val
          )
        )
      }
    }
  
    return(cld)
  }


  ########## extract results
  results <- data_long %>%
    group_by(Organism, GroupVar) %>%
    group_map(run_stat_test) %>%
    bind_rows()


  ########## attach results to plotting dataset
  plot_annotations <- results %>%
    select(Organism, GroupVar, Level, CL_label, Anova_result, Kruskal_result) %>%
    mutate(
      stat_result = coalesce(Anova_result, Kruskal_result)
    )

  data_plot <- left_join(data_long, plot_annotations, by = c("Organism", "GroupVar", "Level"))


  ########## label data for CLD letters
  label_data_clean <- data_plot %>%
    group_by(Organism, GroupVar, Level, CL_label) %>%
    summarise(y_pos = max(Shannon, na.rm = TRUE) * 1.01, .groups = "drop")

  ########## label data for ANOVA p values
  label_data_pval <- data_plot %>%
    group_by(Organism, GroupVar) %>%
    summarise(
      p_value_label = first(stat_result),
      y_pos = max(Shannon, na.rm = TRUE) * 1.09,
      .groups = "drop"
    )

  ########## annotation positions
  label_data_pval <- as.data.frame(label_data_pval)
  label_data_pval$x_pos <- c(0.5, 0.6, 0.5, 0.5, 0.6, 0.5, 0.5, 0.5, 0.5)
  label_data_pval$y_pos <- 6.6
  label_data_clean$Organism <- factor(label_data_clean$Organism, levels = c("Bacteria", "Fungi", "Cercozoa"))
  label_data_clean$GroupVar <- factor(label_data_clean$GroupVar, levels = c("X11C_level", "Root_type", "Root_section"))

  min(data_long$Shannon)
  max(data_long$Shannon)
  
  ########## plot

  y <- ggplot(
    data_long,
    aes(
      x = Level,
      y = Shannon,
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
    scale_y_continuous(name = "Shannon index", expand = expansion(mult = c(0.05, 0.1))) +
    geom_jitter(position = position_jitter(0.2), shape = 1, size = 3) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    theme_linedraw() +
    geom_text(data = label_data_clean, aes(x = Level, y = y_pos, label = CL_label), vjust = -0.8, size = 6, color = "black", na.rm = TRUE) +
    geom_text(data = label_data_pval, aes(x = x_pos, y = y_pos, label = p_value_label), inherit.aes = FALSE, hjust = 0, size = 6) +
    facet_grid2(GroupVar ~ Organism, scales = "free", independent = "all", switch = "y", labeller = labeller(GroupVar = c("X11C_level" = "11C level", "Root_type" = "Root type", "Root_section" = "Root section"))) +
    coord_cartesian(ylim = c(0.8, 6.7)) +
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

########## generate Figure S9a
png("Figure S9a.png", width = 400, height = 300, units = "mm", res = 1000)
grid.arrange(y, ncol = 1)
dev.off()
