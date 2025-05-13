###################################
# REQUIRED LIBRARIES ###############

library(ggplot2)
library(ggpubr)
library(dplyr)

plot_variable_barplot <- function(data, variable, group_col, group_colors = NULL) {
  # Ensure group column is a factor
  data[[group_col]] <- as.factor(data[[group_col]])
  
  # Drop NA values
  df_plot <- data %>%
    filter(!is.na(.data[[variable]]), !is.na(.data[[group_col]]))
  
  # Perform Kruskal-Wallis test
  formula <- as.formula(paste(variable, "~", group_col))
  kruskal_result <- kruskal.test(formula, data = df_plot)
  kruskal_p <- kruskal_result$p.value
  
  # Generate pairwise comparisons
  group_levels <- levels(df_plot[[group_col]])
  comparisons <- combn(group_levels, 2, simplify = FALSE)
  
  # Get y-axis range
  max_y <- max(df_plot[[variable]], na.rm = TRUE)
  min_y <- min(df_plot[[variable]], na.rm = TRUE)
  padding <- (max_y - min_y) * 0.1
  label_y_positions <- seq(from = max_y + padding, by = padding, length.out = length(comparisons))
  
  # Create the base plot
  p <- ggplot(df_plot, aes_string(x = group_col, y = variable, fill = group_col)) +
    stat_summary(fun = mean, geom = "bar", color = "black", width = 0.6) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    stat_compare_means(
      comparisons = comparisons,
      label = "p.signif",
      method = "wilcox.test",
      label.y = label_y_positions
    ) +
    labs(
      title = paste0("Mean ", variable, " by ", group_col, "\nKruskal-Wallis p = ", format(kruskal_p, digits = 3)),
      x = group_col,
      y = variable
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    ) +
    expand_limits(y = c(min_y - padding, max_y + padding))
  
  # Apply user-defined colors if provided
  if (!is.null(group_colors)) {
    if (!all(group_levels %in% names(group_colors))) {
      stop("Some group levels are missing in the provided group_colors vector.")
    }
    p <- p + scale_fill_manual(values = group_colors)
  }
  
  print(p)
}

