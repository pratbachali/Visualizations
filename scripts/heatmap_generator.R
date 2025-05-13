# ============================
# Required Libraries
# ============================
library(ComplexHeatmap)
library(openxlsx)
library(circlize)
library(grid)

# ============================
# Function: Generate Annotations
# ============================
generate_annotations <- function(data, variables, type = c("top", "bottom"), custom_colors = list()) {
  type <- match.arg(type)
  anno_list <- list()
  col_list <- list()
  
  for (var in variables) {
    val <- data[[var]]
    if (is.null(val)) next
    
    if (is.numeric(val)) {
      if (var == "SLEDAI") {
        anno_list[[var]] <- anno_points(val)
        col_list[[var]] <- colorRamp2(c(6, 10, 20), c("thistle2", "plum2", "hotpink"))
      } else if (type == "top" && var %in% c("ComplementC3", "ComplementC4")) {
        anno_list[[var]] <- anno_barplot(val)
      } else if (type == "top" && var == "ANTIdsDNA") {
        val <- log(as.numeric(val) + 1)
        anno_list[[var]] <- anno_lines(val)
      } else {
        anno_list[[var]] <- anno_points(val)
      }
    } else {
      val <- as.factor(val)
      if (var %in% names(custom_colors)) {
        pal <- custom_colors[[var]]
      } else {
        pal <- circlize::rand_color(length(levels(val)), luminosity = "bright")
        names(pal) <- levels(val)
      }
      anno_list[[var]] <- val
      col_list[[var]] <- pal
    }
  }
  
  do.call(HeatmapAnnotation, c(
    list(annotation_name_side = "right",
         gap = unit(1, "mm"),
         show_legend = (type == "bottom"),
         col = col_list),
    anno_list
  ))
}

# ============================
# Function: Generate Heatmap
# ============================
generate_clinical_heatmap <- function(input_file,
                                      bottom_anno_vars,
                                      top_anno_vars,
                                      output_png,
                                      custom_colors = list()) {
  
  message("Reading input file: ", input_file)
  df_hm <- read.xlsx(input_file)
  rownames(df_hm) <- df_hm$friendlyName
  
  # Identify cluster column
  cluster_col <- grep("cluster|vae", colnames(df_hm), value = TRUE)[1]
  df_hm[[cluster_col]] <- factor(df_hm[[cluster_col]])
  column_split <- df_hm[[cluster_col]]
  
  # Create Annotations
  message("Generating annotations...")
  ha1 <- generate_annotations(df_hm, bottom_anno_vars, type = "bottom", custom_colors = custom_colors)
  ha2 <- generate_annotations(df_hm, top_anno_vars, type = "top", custom_colors = custom_colors)
  
  # Prepare heatmap data
  excluded_cols <- unique(c("friendlyName", cluster_col, bottom_anno_vars, top_anno_vars))
  heatmap_vars <- setdiff(colnames(df_hm), excluded_cols)
  df_clinical <- df_hm[, heatmap_vars]
  df_clinical_t <- as.data.frame(t(df_clinical))
  
  col_fun_heatmap = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  # Generate heatmap
  message("Drawing heatmap...")
  set.seed(133)
  ht <- Heatmap(df_clinical_t,
                name = "clinicaldata",
                row_names_side = "left",
                cluster_rows = TRUE,
                row_km = 3,
                cluster_columns = FALSE,
                top_annotation = ha2,
                bottom_annotation = ha1,
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 10),
                row_gap = unit(1, "mm"),
                column_split = column_split,
                border = TRUE,
                show_column_names = FALSE,
                row_title = NULL,
                show_row_dend = FALSE)
  
  # Save to PNG
  message("Saving heatmap to: ", output_png)
  png(output_png, width = 18, height = 10, units = "in", res = 300)
  draw(ht)
  dev.off()
  
  message("✅ Heatmap saved as: ", output_png)
}
