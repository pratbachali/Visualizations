# ------------------ Libraries ------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(rstatix)
  library(ggpubr)
  library(gridExtra)
  library(grid)
  library(plyr)
  library(scales)
})

# ------------------ USER INPUTS ------------------

# Ensure `df` is loaded before running
# Example:
# df <- read.csv("your_data.csv")

grouping_vars <- c("cluster", "SUBSTANC")     # Group by cluster and treatment
timepoint_var <- "timepoint"                  # Timepoint variable
gene_vars <- c("Cell_Cycle", "IG_Chains", "Plasma_Cells")  # GSVA/gene columns

# Auto-generate timepoint color map
time_levels <- unique(df[[timepoint_var]])
timepoint_colors <- setNames(scales::hue_pal()(length(time_levels)), time_levels)

# Output folders
dir.create("plots_individual", showWarnings = FALSE)
dir.create("plots_grid", showWarnings = FALSE)

# ------------------ Checks ------------------

required_cols <- c(grouping_vars, timepoint_var, gene_vars)
missing <- setdiff(required_cols, colnames(df))
if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse = ", ")))

# Subset to required columns
df_vp <- df[, required_cols]

# ------------------ Functions ------------------

# Helper function to save plots in both tiff and pdf formats
save_plots <- function(plots, plot_name, format = c("tiff", "pdf")) {
  for (format_type in format) {
    if (format_type == "tiff") {
      filename <- file.path("plots_individual", paste0(plot_name, ".tiff"))
      tiff(filename, units = "in", width = 5, height = 5, res = 300, type = "cairo")
      print(plots)
      dev.off()
    } else if (format_type == "pdf") {
      filename <- file.path("plots_individual", paste0(plot_name, ".pdf"))
      pdf(filename, width = 5, height = 5)
      print(plots)
      dev.off()
    }
  }
}

# Helper function to save grid plots
save_grid_plots <- function(plots, grid_name) {
  layout <- matrix(1:length(plots), nrow = 1)
  file_name <- file.path("plots_grid", paste0(grid_name, "_grid.tiff"))
  
  tiff(file_name, units = "in", width = 15, height = 5, res = 300, type = "cairo")
  grid.arrange(grobs = plots, layout_matrix = layout,
               top = textGrob(grid_name, gp = gpar(fontsize = 18)))
  dev.off()
  
  # Save as PDF
  file_name_pdf <- file.path("plots_grid", paste0(grid_name, "_grid.pdf"))
  pdf(file_name_pdf, width = 15, height = 5)
  grid.arrange(grobs = plots, layout_matrix = layout,
               top = textGrob(grid_name, gp = gpar(fontsize = 18)))
  dev.off()
}

# Subset dataframe by grouping variables
subset_df_per_cluster_and_treatment <- function(df, cols2group) {
  df %>% group_nest(across(all_of(cols2group)))
}

# Calculate p-values for each gene
calculate_p_vals <- function(df, genes, timepoint_var) {
  signi_list <- list()
  for (gene in genes) {
    res <- df %>%
      rstatix::pairwise_t_test(as.formula(paste0(gene, " ~ ", timepoint_var)),
                               p.adjust.method = "bonferroni", paired = FALSE) %>%
      mutate(
        p.adj.signif = case_when(
          p.adj < 0.001 ~ "***",
          p.adj < 0.01  ~ "**",
          p.adj < 0.05  ~ "*",
          TRUE          ~ "ns"
        )
      )
    signi_list[[gene]] <- res
  }
  return(signi_list)
}

# Create violin plots
plot_violin_plot <- function(df, genes, timepoint_var, pvals_list, color_map) {
  names(pvals_list) <- genes
  df <- df[, c(genes, timepoint_var)]
  melted_df <- melt(df, id.vars = timepoint_var)
  melted_df <- plyr::rename(melted_df, c("value" = "GSVA_score"))
  
  plot_list <- list()
  for (gene in genes) {
    data_gene <- melted_df %>% filter(variable == gene)
    stats <- as_tibble(pvals_list[[gene]])
    
    if (nrow(stats) > 0) {
      upper_bound <- 2  # y-axis upper limit from coord_cartesian
      step <- 0.2       # space between significance bars
      num_comparisons <- nrow(stats)
      
      # Ensure the annotations fit within the y-axis limit
      max_position <- upper_bound - 0.1
      min_position <- upper_bound - step * num_comparisons
      
      # Assign descending positions to prevent overlap
      stats <- stats %>%
        arrange(p.adj) %>%
        mutate(y.position = seq(from = min_position, to = max_position, length.out = nrow(stats)))
    }
    
    p <- ggplot(data_gene, aes(x = .data[[timepoint_var]], y = GSVA_score, color = .data[[timepoint_var]])) +
      scale_color_manual(values = color_map) +
      facet_wrap(~variable, scales = "free") +
      theme_minimal() +
      coord_cartesian(ylim = c(-2, 2)) +
      theme(legend.position = "none") +
      geom_violin(position = position_dodge(0.2), trim = FALSE, size = 0.85) +
      geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, binwidth = 0.1)
    
    if (nrow(stats) > 0) {
      p <- p + stat_pvalue_manual(
        stats,
        label = "p.adj.signif",
        y.position = "y.position",
        tip.length = 0.01,
        label.size = 4.0,
        step.increase = 0.07,
        hide.ns = FALSE  # Show 'ns' when not significant
      )
    }
    
    plot_list[[gene]] <- p
  }
  
  return(plot_list)
}

# ------------------ Main Workflow ------------------
subset_dfs <- subset_df_per_cluster_and_treatment(df_vp, grouping_vars)

pvals_list <- list()
for (i in seq_along(subset_dfs$data)) {
  pvals_list[[i]] <- calculate_p_vals(subset_dfs$data[[i]], gene_vars, timepoint_var)
}
names(pvals_list) <- paste(subset_dfs[[grouping_vars[1]]], subset_dfs[[grouping_vars[2]]], sep = "_")

plot_list_all <- list()

# Process each subset of the data
for (i in seq_along(subset_dfs$data)) {
  key_name <- paste(subset_dfs[[grouping_vars[1]]][i], subset_dfs[[grouping_vars[2]]][i], sep = "_")
  message("Processing: ", key_name)
  
  gene_plots <- plot_violin_plot(
    df = subset_dfs$data[[i]],
    genes = gene_vars,
    timepoint_var = timepoint_var,
    pvals_list = pvals_list[[i]],
    color_map = timepoint_colors
  )
  
  plot_list_all[[key_name]] <- gene_plots
  
  # Save individual plots in both PDF and TIFF formats (ensuring no duplicates)
  for (gene in names(gene_plots)) {
    save_plots(gene_plots[[gene]], paste0(key_name, "_", gene))
  }
}

# ------------------ Save Grouped Grid Plots ------------------
for (group_name in names(plot_list_all)) {
  gene_set <- plot_list_all[[group_name]]
  save_grid_plots(gene_set, group_name)
}
