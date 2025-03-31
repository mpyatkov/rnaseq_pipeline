library(tidyverse) #2.0.0
library(RColorBrewer) #1.1.3
library(pheatmap) #1.0.13
library(tools) #4.4.3
library(gridExtra) #2.3
library(openxlsx) #4.2.8
library(grid)

read_data <- function(input) {
  if (!file.exists(input)) {
    stop("Error: The file does not exist. Please check the file path.")
  }
  
  if (grepl("\\.txt$", input, ignore.case = TRUE)) {
    data <- read.table(input, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else if (grepl("\\.csv$", input, ignore.case = TRUE)) {
    data <- read.csv(input, header = TRUE, stringsAsFactors = FALSE)
  } else if (grepl("\\.(xls|xlsx)$", input, ignore.case = TRUE)) {
    data <- read_excel(input)
  } else {
    stop("Error: Unsupported file format. Please provide a .txt, .csv, or .xls/.xlsx file.")
  }
  
  return(data)
}

# Match column order for noncluster
match_column <- function(reference, col_need_match){
  target_order <- reference$Condition_Name
  matched_columns <- sapply(target_order, function(cond) {
    matched <- grep(cond, col_need_match, value = TRUE)
    if (length(matched) == 0) return(NA)
    matched[1]
  })
  col_need_match <- as.vector(na.omit(matched_columns))
}

# Select color
color_selected <- function(color_length) {
  color_total <- c(
    "#e6194b","#ffe119","#46f0f0","#f58231","#bcf60c", 
    "#ff00ff","#9a6324","#fffac8","#e6beff","#00bfff", 
    "#ffd8b1","#00ff7f","#f5a9bc","#1e90ff","#ffa500",
    "#98fb98","#911eb4","#afeeee","#fa8072","#9acd32",
    "#3cb44b","#000075","#808000","#cd5c5c","#dda0dd",
    "#40e0d0","#ff69b4","#8a2be2","#c71585","#5f9ea0",
    "#dc143c","#87cefa","#ff6347","#9932cc","#00ced1",
    "#ff4500","#6a5acd","#b0e0e6","#d2691e","#a9a9f5",
    "#adff2f","#8b0000","#7fffd4","#00fa9a","#ba55d3",
    "#2e8b57","#ffdab9","#b22222","#ffe4e1","#7b68ee"
  )
  
  if (color_length <= length(color_total)) {
    return(color_total[1:color_length])
  }
  
  warning("groups is larger than 50, color will randomly select")
  palette_fn <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
  extra_colors <- palette_fn(color_length - length(color_total))
  colors <- c(color_total, extra_colors)
  return(colors)
}

# Get tpm in any group, return gene_to_keep
get_genes_with_high_tpm_in_any_group <- function(data, annotation_col, tpm_cutoff) {
  # Get all column names starting with "tpm_"
  tpm_cols <- grep("^tpm_", colnames(data), value = TRUE)
  
  # Extract condition names by removing "tpm_" prefix
  condition_names <- sub("^tpm_", "", tpm_cols)
  
  # Create a mapping between tpm column names and their corresponding group
  tpm_info <- data.frame(
    tpm_col = tpm_cols,
    Condition_Name = condition_names,
    stringsAsFactors = FALSE
  )
  tpm_info$Group <- annotation_col[tpm_info$Condition_Name, "Group"]
  
  # Create a TPM matrix and rename columns with group names
  tpm_matrix <- data[, tpm_info$tpm_col]
  colnames(tpm_matrix) <- tpm_info$Group
  tpm_matrix$id <- data$id
  
  # Convert wide matrix to long format for group-wise processing
  tpm_long <- tpm_matrix %>%
    pivot_longer(-id, names_to = "Group", values_to = "TPM")
  
  # For each gene and group, check if all samples in the group have TPM > cutoff
  group_check <- tpm_long %>%
    group_by(id, Group) %>%
    summarise(all_above = all(TPM >= tpm_cutoff), .groups = "drop")
  
  # Keep genes where any group satisfies the condition
  gene_to_keep <- group_check %>%
    group_by(id) %>%
    summarise(any_group_all_above = any(all_above), .groups = "drop") %>%
    filter(any_group_all_above) %>%
    pull(id)
  
  return(gene_to_keep)
}

pair_heatmap <- function(input,
                         condition_1,
                         condition_2,
                         output_folder,
                         save_name,
                         custom_color,
                         using_package_f,
                         scale_function,
                         genic_region,
                         log2FC,
                         FDR,
                         tpm_cutoff,
                         cluster_number_row,
                         cluster_number_col,
                         order_for_noclustering,
                         annotation_colors,
                         folder_name){
  if (is.character(input) && length(input) == 1 && input != "") {
    if (!using_package_f %in% c("EdgeR", "DESeq2")) {
      stop("Error: 'using_package_f' must be either 'EdgeR' or 'DESeq2'.")
    }
    
    #READ DATA
    data <- read_data(input)
    
    group_data <- data %>%
      select(id,starts_with("tpm_mean"))
    
    data <- data %>%
      select(all_of(colnames(data)[
        which(colnames(data) == "id") : (which(grepl("^rpkm_", colnames(data)))[1] - 1)
      ]),
      starts_with("tpm"), -starts_with("tpm_mean"),
      starts_with(using_package_f),
      )
    annotation_col <- order_for_noclustering
    rownames(annotation_col) <- annotation_col$Condition_Name
    annotation_col$Condition_Name <- NULL
    
    if(tpm_cutoff!= 0){
    gene_to_keep <- get_genes_with_high_tpm_in_any_group(data = data, annotation_col = annotation_col, tpm_cutoff = tpm_cutoff)
    data_filtered <- data %>%
      filter(id %in% gene_to_keep) %>%
      filter(
        (abs(pull(select(., ends_with("log2FoldChange")))) > log2FC) & 
          (pull(select(., ends_with("FDR"))) < FDR)
      )}else if (tpm_cutoff == 0){
        data_filtered <- data %>%
          filter(
            (abs(pull(select(., ends_with("log2FoldChange")))) > log2FC) & 
              (pull(select(., ends_with("FDR"))) < FDR)  
          )}
  
    if (nrow(data_filtered) == 0) {
      stop("Error: No data left after filtering. Consider adjusting log2FC, FDR or tpm_cutoff thresholds.")
    }
    
    if (scale_function == "maxscale") {
      norm_factor <- apply(select(data_filtered, starts_with("tpm")), 1, function(x) max(abs(x), na.rm = TRUE))
      data_filtered <- data_filtered %>%
        mutate(across(starts_with("tpm"), 
                      ~ ifelse(norm_factor == 0, 0, .x / norm_factor), 
                      .names = "scale_{.col}")) %>%
        as.data.frame()  
    } 
    
    else if (scale_function == "zscale") {
      mean_x <- apply(select(data_filtered, starts_with("tpm")), 1, mean, na.rm = TRUE)
      sd_x <- apply(select(data_filtered, starts_with("tpm")), 1, sd, na.rm = TRUE)
      data_filtered <- data_filtered %>%
        mutate(across(starts_with("tpm"), 
                      ~ ifelse(sd_x == 0, 0, (.x - mean_x) / sd_x), 
                      .names = "scale_{.col}")) %>%
        as.data.frame() 
      
      min_scale <- apply(select(data_filtered, starts_with("scale_")), 1, min, na.rm = TRUE)
      max_scale <- apply(select(data_filtered, starts_with("scale_")), 1, max, na.rm = TRUE)
      max_norm <- pmax(abs(min_scale), max_scale)
      data_filtered <- data_filtered %>%
        mutate(across(starts_with("scale_"), 
                      ~ case_when(
                        .x != 0  ~ (.x / max_norm),
                        .x == 0 ~ 0
                      )))
    }
    
    upregulated_count <- sum(data_filtered %>% select(starts_with(using_package_f)) %>% apply(1, function(x) any(x > 1)))
    downregulated_count <- sum(data_filtered %>% select( starts_with(using_package_f)) %>% apply(1, function(x) any(x < -1)))
    
    
    df <- data_filtered %>%
      select("id", starts_with("scale_"))
    rownames(df) <- df$id
    df$id <- NULL
    df <- as.matrix(df)
    df[is.na(df)] <- 0 
    colnames(df) <- gsub("^scale_tpm_", "", colnames(df))
    
    valid_conditions <- rownames(annotation_col) %in% colnames(df)
    annotation_col <- annotation_col[valid_conditions, , drop = FALSE]
    annotation_colors$Group <- annotation_colors$Group[unique(annotation_col$Group)]
    
    HM_origin <- pheatmap(df, 
                          cellwidth = 20, 
                          cellheight = 2 , 
                          fontsize_col = 10,
                          fontsize_row = 2,
                          border_color = NA, 
                          cutree_cols = cluster_number_col, 
                          cutree_rows = cluster_number_row, 
                          cluster_rows = TRUE, 
                          cluster_cols = TRUE,
                          annotation_col = annotation_col,
                          annotation_colors = annotation_colors,
                          treeheight_row = 200,
                          treeheight_col = 100,
                          angle_col = 90, 
                          color = custom_color
    )
    
    row_order_origin <- HM_origin$tree_row$order
    col_order_origin <- HM_origin$tree_col$order
    
    #save col order
    col_order_table <- data.frame(
      COL_NUMBER = seq_along(col_order_origin),  
      COL_NAME = colnames(df)[col_order_origin],
      stringsAsFactors = FALSE
    )
    col_order_table$Group <- annotation_col[col_order_table$COL_NAME, "Group"]
    col_order_table <- col_order_table %>%
      select(COL_NUMBER,Group,COL_NAME)
    
    heatmap_width_ratio <- 0.75  
    heatmap_height_ratio <- 0.9 
    
    cellwidth <- (8.27 * heatmap_width_ratio * 72 - 275- 5 * cluster_number_col) / ncol(df)    
    cellheight <- if(ncol(df)>30){(11.69 * heatmap_height_ratio * 72-95- 5 * cluster_number_row) / nrow(df)}
    else{(11.69 * heatmap_height_ratio * 72-105- 5 * cluster_number_row) / nrow(df)}
    
    df <- df[, col_order_origin] 
    colnames(df) <- seq_len(ncol(df))
    annotation_col <- annotation_col[col_order_origin, , drop = FALSE]
    rownames(annotation_col) <- as.character(seq_len(ncol(df)))
    
    HM_scale <- pheatmap(df, 
                         cellwidth = cellwidth, 
                         cellheight = cellheight,
                         fontsize_col = 100/ncol(df),
                         fontsize_row = 2,
                         border_color = NA, 
                         cutree_cols = cluster_number_col, 
                         cutree_rows = cluster_number_row, 
                         annotation_col = annotation_col,
                         annotation_colors = annotation_colors,
                         cluster_rows = TRUE, 
                         cluster_cols = TRUE,
                         treeheight_row = 100,
                         treeheight_col = 100,
                         angle_col = 90, 
                         show_colnames = if(ncol(df)>50){FALSE}else{TRUE},
                         show_rownames = FALSE,
                         color = custom_color)
    
    origin_pdf_width <- (ncol(df) * 20) / 72 + 6    
    origin_pdf_height <- (nrow(df) * 2) / 72 + 7 + if(!is.na(cluster_number_row)){cluster_number_row}else{0}
    
    #add row order to data
    data_filtered <- data_filtered %>%
      mutate(heatmap_order = match(1:nrow(data_filtered), row_order_origin)) %>%
      arrange(heatmap_order)
    first_columns <- c("heatmap_order", "id")
    tpm_columns <- grep("^tpm_", colnames(data_filtered), value = TRUE)  
    scale_columns <- grep("^scale_", colnames(data_filtered), value = TRUE)
    
    method_columns <- grep(paste0("^", using_package_f), colnames(data_filtered), value = TRUE)
    all_other_columns <- colnames(data_filtered)
    counts_columns <- setdiff(all_other_columns, c(first_columns, scale_columns, tpm_columns, method_columns))
    data_filtered <- data_filtered %>%
      rename_with(~ paste0("counts_", .), all_of(counts_columns))
    counts_columns <- grep("^counts_", colnames(data_filtered), value = TRUE)
    
    if (scale_function == "zscale") {
      data_filtered <- data_filtered %>%
        rename_with(~ paste0("z", .), all_of(scale_columns))
      scale_columns <- paste0("z", scale_columns)
    } else if (scale_function == "maxscale") {
      data_filtered <- data_filtered %>%
        rename_with(~ paste0("max", .), all_of(scale_columns))
      scale_columns <- paste0("max", scale_columns)
    }
    
    
    data_filtered <- data_filtered %>%
      select(any_of(first_columns), any_of(scale_columns), any_of(counts_columns), any_of(tpm_columns), any_of(method_columns))
    
    
    # Remove fc, log2fc, FDR, pvalue, ...
    data <- data %>%
      select(-starts_with(using_package_f))
    
    tmp_path <- file.path(output_folder, "tmp/")
    long_filepath <- paste0(tmp_path, file_path_sans_ext(basename(input)), "_long.pdf")
    short_filepath <- paste0(tmp_path, file_path_sans_ext(basename(input)), "_short.pdf")
    
    
    # Save a
    pdf(long_filepath, width = origin_pdf_width, height = origin_pdf_height)
    grid.newpage()
    title <- paste0(folder_name,"\n",file_path_sans_ext(basename(input)), 
                    "\nTREATMENT: ", condition_2, 
                    ",\tCONTROL: ", condition_1,
                    ",\nUPREGULATED: ", upregulated_count,
                    ",\tDOWNREGULATED: ", downregulated_count,
                    "\nSCALE FUNCTION: ", scale_function, ",\tTPM CUTOFF: ", tpm_cutoff, 
                    "\nLOG2FC: ", log2FC, ",\tFDR: ", FDR)
    grid.text(title, x = 0.5, y = 0.999, gp = gpar(fontsize = 10, fontface = "bold"), just = c("center","top"))
    pushViewport(viewport(x = 0.5, y = 0.49, width = 0.9, height = 0.7)) 
    grid.draw(HM_origin$gtable) 
    popViewport()
    
    dev.off()
    
    
    # Save b
    pdf(short_filepath, width = 8.27, height = 11.69)
    
    title <- paste0(folder_name,"\n",file_path_sans_ext(basename(input)))
    grid.newpage()
    grid.text(title, x = 0.5, y = 0.99, gp = gpar(fontsize = 10, fontface = "bold"), just = c("center","top"))
    
    pushViewport(viewport(x = 0.35, y = 0.5, width = heatmap_width_ratio, height = heatmap_height_ratio))
    print(HM_scale, newpage = FALSE)
    
    table_theme <- ttheme_default(
      core = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                  padding = unit(c(1.5, 1.5), "mm")), 
      colhead = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                     padding = unit(c(1.5, 1.5), "mm")), 
      rowhead = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                     padding = unit(c(1.5, 1.5), "mm"))
    ) 
    info <- paste0("TREATMENT: ", condition_2, "\nCONTROL: ", condition_1,
                   "\n\nUPREGULATED: ", upregulated_count, "\nDOWNREGULATED: ",downregulated_count,
                   "\n\nSCALE FUNCTION: ", scale_function, "\nTPM CUTOFF: ", tpm_cutoff, "\nLOG2FC: ", log2FC,
                   "\nFDR: ", FDR)
    grid.text(info, x = 1, y = 0.9, gp = gpar(fontsize = 8), just = c("left","top"))
    pushViewport(viewport(x = 1, y = 0.45, width = 1-heatmap_width_ratio, height = 0.15, just = c("left","top")))
    grid.table(col_order_table, rows = NULL, cols = NULL, theme = table_theme)
    popViewport()
    
    dev.off()
    
    
    
    # Save c
    pattern <- "_(ExonCollapsed|FullGeneBody)_.*$"
    match <- str_extract(file_path_sans_ext(basename(input)), pattern)
    if (is.na(match)) {
      match <- "UnknownGeneType"
    } else {
      match <- substr(match, 2, nchar(match)) 
    }
    output_filename <- file.path(output_folder, paste0("HeatmapData_",folder_name, "_", match, "_log2FC_", log2FC, "_FDR_", FDR, ".xlsx"))
    write.xlsx(data_filtered, output_filename, rowNames = FALSE)
  }
  
  return(list(data = data,
              group_data = group_data,
              id_list = data_filtered$id,
              upregulated_count = upregulated_count,
              downregulated_count = downregulated_count,
              col_order_table = col_order_table,
              long_filepath = long_filepath,
              short_filepath = short_filepath))
  
}




heatmap_integrated <- function(union_id_list, 
                               all_data,
                               save_name,
                               log2FC,
                               FDR,
                               tpm_cutoff, 
                               scale_function, 
                               genic_region,
                               sample_type, 
                               output_folder, 
                               custom_color, 
                               cluster_number_row, 
                               cluster_number_col,
                               order_for_noclustering,
                               annotation_colors){
  filtered_data_list <- lapply(all_data, function(df) df[df$id %in% union_id_list, ])
  data_integrated <- filtered_data_list[[1]]
  
  for (i in 2:length(filtered_data_list)) {
    new_data <- filtered_data_list[[i]]
    new_cols <- setdiff(colnames(new_data), colnames(data_integrated))
    
    if (length(new_cols) > 0) { 
      data_integrated <- merge(data_integrated, new_data[, c("id", new_cols)], by = "id", all = TRUE)
    }
  }
  data_integrated <- data_integrated[, unique(colnames(data_integrated))]
  
  annotation_col <- order_for_noclustering
  rownames(annotation_col) <- annotation_col$Condition_Name
  annotation_col$Condition_Name <- NULL
 
  
  if (scale_function == "maxscale") {
    norm_factor <- apply(select(data_integrated, starts_with("tpm")), 1, function(x) max(abs(x), na.rm = TRUE))
    data_integrated <- data_integrated %>%
      mutate(across(starts_with("tpm"), 
                    ~ ifelse(norm_factor == 0, 0, .x / norm_factor), 
                    .names = "scale_{.col}")) %>%
      as.data.frame()  
  }else if (scale_function == "zscale") {
    mean_x <- apply(select(data_integrated, starts_with("tpm")), 1, mean, na.rm = TRUE)
    sd_x <- apply(select(data_integrated, starts_with("tpm")), 1, sd, na.rm = TRUE)
    data_integrated <- data_integrated %>%
      mutate(across(starts_with("tpm"), 
                    ~ ifelse(sd_x == 0, 0, (.x - mean_x) / sd_x), 
                    .names = "scale_{.col}")) %>%
      as.data.frame()  
    
    min_scale <- apply(select(data_integrated, starts_with("scale_")), 1, min, na.rm = TRUE)
    max_scale <- apply(select(data_integrated, starts_with("scale_")), 1, max, na.rm = TRUE)
    max_norm <- pmax(abs(min_scale), max_scale)
    data_integrated <- data_integrated %>%
      mutate(across(starts_with("scale_"), 
                    ~ case_when(
                      .x != 0  ~ (.x / max_norm),
                      .x == 0 ~ 0
                    )))
  }
  
  df <- data_integrated %>%
    select("id", starts_with("scale_"))
  rownames(df) <- df$id
  df$id <- NULL
  df <- as.matrix(df)
  df[is.na(df)] <- 0 
  colnames(df) <- gsub("^scale_tpm_", "", colnames(df))
  colnames(df) <- gsub("^mean_", "", colnames(df))
  #sort nonclustering col order
  if (!is.null(order_for_noclustering)) {
    target_order <- order_for_noclustering$Condition_Name
    matched_cols <- colnames(df)
    matched_index <- sapply(target_order, function(cond) {
      idx <- which(grepl(cond, matched_cols))
      if (length(idx) == 0) return(NA)
      idx[1]
    })
    matched_index <- matched_index[!is.na(matched_index)]
    df <- df[, matched_index, drop = FALSE]
  }
 
  valid_conditions <- rownames(annotation_col) %in% colnames(df)
  annotation_col <- annotation_col[valid_conditions, , drop = FALSE]
  annotation_colors$Group <- annotation_colors$Group[unique(annotation_col$Group)]
  
  
  #doing heatmap
  HM_origin <- pheatmap(df, 
                        cellwidth = 20, 
                        cellheight = 2 , 
                        fontsize_col = 10,
                        fontsize_row = 2,
                        border_color = NA, 
                        cutree_cols = if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]}, 
                        cutree_rows = if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}, 
                        cluster_rows = TRUE, 
                        cluster_cols = TRUE,
                        annotation_col = annotation_col,
                        annotation_colors = annotation_colors,
                        treeheight_row = 200,
                        treeheight_col = 80,
                        angle_col = 90, 
                        color = custom_color
  )
  
  row_order_origin <- HM_origin$tree_row$order
  col_order_origin <- HM_origin$tree_col$order
  #save col order
  col_order_table <- data.frame(
    COL_NUMBER = seq_along(col_order_origin),  
    COL_NAME = colnames(df)[col_order_origin],
    stringsAsFactors = FALSE
  )
  col_order_table$Group <- annotation_col[col_order_table$COL_NAME, "Group"]
  col_order_table <- col_order_table %>%
    select(COL_NUMBER,Group,COL_NAME)
  
  heatmap_width_ratio <- 0.75  
  heatmap_height_ratio <- 0.9 
  
  cellwidth <- (8.27 * heatmap_width_ratio * 72-275 - 5 *if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]}) / ncol(df)    
  cellheight <- if(ncol(df)>50){(11.69 * heatmap_height_ratio * 72-100 - 5*if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}) / nrow(df)}
  else{(11.69 * heatmap_height_ratio * 72- 110 - 5 * if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}) / nrow(df)}
  
  cellwidth_noclustering <- (8.27 * heatmap_width_ratio * 72-260) / ncol(df)    
  cellheight_noclustering <- if(ncol(df)>50){(11.69 * heatmap_height_ratio * 72 -15 - 5*if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}) / nrow(df)}
  else{(11.69 * heatmap_height_ratio * 72- 25 - 5 * if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}) / nrow(df)}
  
  
  HM_noclustering_origin <- pheatmap(df, 
                                     cellwidth = 20, 
                                     cellheight = 2 , 
                                     fontsize_col = 10,
                                     fontsize_row = 2,
                                     border_color = NA, 
                                     cutree_cols = if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]}, 
                                     cutree_rows = if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}, 
                                     cluster_rows = TRUE, 
                                     cluster_cols = FALSE,
                                     annotation_col = annotation_col,
                                     annotation_colors = annotation_colors,
                                     treeheight_row = 200,
                                     angle_col = 90, 
                                     color = custom_color
  )
  
  col_order_table_noclustering <- data.frame(
    COL_NUMBER = seq_along(colnames(df)),  
    COL_NAME = colnames(df),
    stringsAsFactors = FALSE
  )
  col_order_table_noclustering$Group <- annotation_col[col_order_table_noclustering$COL_NAME, "Group"]
  col_order_table_noclustering <- col_order_table_noclustering %>%
    select(COL_NUMBER,Group,COL_NAME)
  
  
  colnames(df) <- seq_len(ncol(df))
  rownames(annotation_col) <- as.character(seq_len(ncol(df)))
  
  HM_noclustering_scale <- pheatmap(df, 
                                    cellwidth = cellwidth_noclustering, 
                                    cellheight = cellheight_noclustering,
                                    fontsize_col = 100/ncol(df),
                                    fontsize_row = 2,
                                    border_color = NA, 
                                    cutree_cols = if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]}, 
                                    cutree_rows = if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}, 
                                    cluster_rows = TRUE, 
                                    cluster_cols = FALSE,
                                    annotation_col = annotation_col,
                                    annotation_colors = annotation_colors,
                                    treeheight_row = 130,
                                    angle_col = 90, 
                                    show_colnames = if(ncol(df)>50){FALSE}else{TRUE},
                                    show_rownames = FALSE,
                                    color = custom_color)
  
  df <- df[, col_order_origin] 
  annotation_col <- annotation_col[col_order_origin, , drop = FALSE]
  
  HM_scale <- pheatmap(df, 
                       cellwidth = cellwidth, 
                       cellheight = cellheight,
                       fontsize_col = 100/ncol(df),
                       fontsize_row = 2,
                       border_color = NA, 
                       cutree_cols = if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]}, 
                       cutree_rows = if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}, 
                       cluster_rows = TRUE, 
                       cluster_cols = TRUE,
                       annotation_col = annotation_col,
                       annotation_colors = annotation_colors,
                       treeheight_row = 130,
                       treeheight_col = 80,
                       angle_col = 90, 
                       show_colnames = if(ncol(df)>50){FALSE}else{TRUE},
                       show_rownames = FALSE,
                       color = custom_color)
  
  
  
  #PROCESSING FOR DATA_INTEGRATED
  if(sample_type == "individuals"){
    data_integrated <- data_integrated %>%
      mutate(heatmap_order = match(1:nrow(data_integrated), row_order_origin)) %>%
      arrange(heatmap_order)
    first_columns <- c("heatmap_order", "id")  
    scale_columns <- grep("^scale_", colnames(data_integrated), value = TRUE)  
    tpm_columns <- grep("^tpm_", colnames(data_integrated), value = TRUE)  
    all_other_columns <- colnames(data_integrated)
    counts_columns <- setdiff(all_other_columns, c(first_columns, scale_columns, tpm_columns))
    data_integrated <- data_integrated %>%
      rename_with(~ paste0("counts_", .), all_of(counts_columns))
    counts_columns <- grep("^counts_", colnames(data_integrated), value = TRUE)
    
    if (scale_function == "zscale") {
      data_integrated <- data_integrated %>%
        rename_with(~ paste0("z", .), all_of(scale_columns))
      scale_columns <- paste0("z", scale_columns)
    } else if (scale_function == "maxscale") {
      data_integrated <- data_integrated %>%
        rename_with(~ paste0("max", .), all_of(scale_columns))
      scale_columns <- paste0("max", scale_columns)
    }
    if (!is.null(order_for_noclustering)) {
    scale_columns <- match_column(order_for_noclustering, scale_columns)
    tpm_columns <- match_column(order_for_noclustering, tpm_columns)
    counts_columns <- match_column(order_for_noclustering, counts_columns)}
    data_integrated <- data_integrated %>%
      select(all_of(first_columns), all_of(scale_columns), all_of(counts_columns), all_of(tpm_columns))}
  
  else if(sample_type == "groups"){
    data_integrated <- data_integrated %>%
      mutate(heatmap_order = match(1:nrow(data_integrated), row_order_origin)) %>%
      arrange(heatmap_order)
    first_columns <- c("heatmap_order", "id")  
    scale_columns <- grep("^scale_", colnames(data_integrated), value = TRUE)  
    tpm_columns <- grep("^tpm_", colnames(data_integrated), value = TRUE) 
    
    if (scale_function == "zscale") {
      data_integrated <- data_integrated %>%
        rename_with(~ paste0("z", .), all_of(scale_columns))
      scale_columns <- paste0("z", scale_columns)
    } else if (scale_function == "maxscale") {
      data_integrated <- data_integrated %>%
        rename_with(~ paste0("max", .), all_of(scale_columns))
      scale_columns <- paste0("max", scale_columns)
    }
    if (!is.null(order_for_noclustering)) {
    scale_columns <- match_column(order_for_noclustering, scale_columns)
    tpm_columns <- match_column(order_for_noclustering, tpm_columns)}
    data_integrated <- data_integrated %>%
      select(all_of(first_columns), all_of(scale_columns), all_of(tpm_columns))

  }
  
  
  
  origin_pdf_width <- (ncol(df) * 20) / 72 + 6 + 0.3*if(!is.na(if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]})){if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}}else{0}
  origin_pdf_height <- (nrow(df) * 2) / 72 + 10 + 0.3*if(!is.na(if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]})){if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]}}else{0}
  origin_pdf_width_noclustering <- (ncol(df) * 20) / 72 + 6
  origin_pdf_height_noclustering <- (nrow(df) * 2) / 72 + 8
  
  if(sample_type == "individuals"){
    cluster_long_path <- paste0("Heatmap_AllIndividuals_SamplesClustered_Long_", save_name, ".pdf")
    nocluster_long_path <- paste0("Heatmap_AllIndividuals_NoSamplesClustering_Long_", save_name, ".pdf")
    cluster_short_path <- paste0("Heatmap_AllIndividuals_SamplesClustered_Short_", save_name, ".pdf")
    nocluster_short_path <- paste0("Heatmap_AllIndividuals_NoSamplesClustering_Short_", save_name, ".pdf")
    heatmapdata_path <- paste0("HeatmapData_AllIndividuals_",save_name,"_log2FC_",log2FC,"_FDR_",FDR,".xlsx")}
  else if(sample_type == "groups"){
    cluster_long_path <- paste0("Heatmap_AllGroups_SamplesClustered_Long_", save_name, ".pdf")
    nocluster_long_path <- paste0("Heatmap_AllGroups_NoSamplesClustering_Long_", save_name, ".pdf")
    cluster_short_path <- paste0("Heatmap_AllGroups_SamplesClustered_Short_", save_name, ".pdf")
    nocluster_short_path <- paste0("Heatmap_AllGroups_NoSamplesClustering_Short_", save_name, ".pdf")
    heatmapdata_path <- paste0("HeatmapData_AllGroups_",save_name,"_log2FC_",log2FC,"_FDR_",FDR,".xlsx")}
  
  # SAVE integrated long cluster
  pdf(file.path(output_folder, cluster_long_path), width = origin_pdf_width, height = origin_pdf_height)
  grid.newpage()
  title <- paste0(save_name, "_SampleClustered_",sample_type,
                  "\nSCALE FUNCTION: ", scale_function,",\tTPM CUTOFF: ", tpm_cutoff, 
                  "\nLOG2FC: ", log2FC,",\tFDR: ", FDR,
                  "\nGENE NUMBER: ", nrow(df))
  grid.text(title, x = 0.5, y = 0.999, gp = gpar(fontsize = 12, fontface = "bold"), just = c("center","top"))
  pushViewport(viewport(x = 0.5, y = 0.9, width = 0.9, height = 0.8, just = c("center","top"))) 
  grid.draw(HM_origin$gtable) 
  popViewport()
  dev.off()
  
  # SAVE integrated long noclustered
  pdf(file.path(output_folder, nocluster_long_path), width = origin_pdf_width_noclustering, height = origin_pdf_height_noclustering)
  grid.newpage()
  title <- paste0(save_name, "_NoSampleClustering_",sample_type,
                  "\nSCALE FUNCTION: ", scale_function,",\tTPM CUTOFF: ", tpm_cutoff, 
                  "\nLOG2FC: ", log2FC,",\tFDR: ", FDR,
                  "\nGENE NUMBER: ", nrow(df))
  grid.text(title, x = 0.5, y = 0.999, gp = gpar(fontsize = 12, fontface = "bold"), just = c("center","top"))
  pushViewport(viewport(x = 0.5, y = 0.9, width = 0.9, height = 0.8, just = c("center","top"))) 
  grid.draw(HM_noclustering_origin$gtable) 
  popViewport()
  dev.off()
  
  # SAVE integrated short noclustering
  pdf(file.path(output_folder, nocluster_short_path), width = 8.27, height = 11.69)
  title <- paste0(save_name,"_NoSampleClustering_",sample_type)
  info <- paste0("\nSCALE FUNCTION: ", scale_function,",\nTPM CUTOFF: ", tpm_cutoff, 
                 "\nLOG2FC: ", log2FC,",\nFDR: ", FDR,
                 "\nGENE NUMBER: ", nrow(df))
  grid.newpage()
  grid.text(title, x = 0.5, y = 0.98, gp = gpar(fontsize = 12, fontface = "bold"), just = c("center","top"))
  pushViewport(viewport(x = 0.35, y = 0.95, width = heatmap_width_ratio, height = heatmap_height_ratio, just = c("center","top")))
  print(HM_noclustering_scale, newpage = FALSE)
  table_theme <- ttheme_default(
    core = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                padding = unit(c(1.5, 1.5), "mm")), 
    colhead = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                   padding = unit(c(1.5, 1.5), "mm")), 
    rowhead = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                   padding = unit(c(1.5, 1.5), "mm"))
  )
  grid.text(info, x = 1, y = 1, gp = gpar(fontsize = 8), just = c("left","top"))
  pushViewport(viewport(x = 1.025, y = 0.5, width = 1-heatmap_width_ratio, height = 0.15, just = c("left","top")))
  grid.table(col_order_table_noclustering, rows = NULL, cols = NULL, theme = table_theme)
  popViewport()
  dev.off()
  
  
  # short clustered
  pdf(file.path(output_folder, cluster_short_path), width = 8.27, height = 11.69)
  title <- paste0(save_name,"_SampleClustered_",sample_type)
  info <- paste0("\nSCALE FUNCTION: ", scale_function,",\nTPM CUTOFF: ", tpm_cutoff, 
                 "\nLOG2FC: ", log2FC,",\nFDR: ", FDR,
                 "\nGENE NUMBER: ", nrow(df))
  grid.newpage()
  grid.text(title, x = 0.5, y = 0.98, gp = gpar(fontsize = 12, fontface = "bold"), just = c("center","top"))
  pushViewport(viewport(x = 0.35, y = 0.95, width = heatmap_width_ratio, height = heatmap_height_ratio, just = c("center","top")))
  print(HM_scale, newpage = FALSE)
  table_theme <- ttheme_default(
    core = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                padding = unit(c(1.5, 1.5), "mm")), 
    colhead = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                   padding = unit(c(1.5, 1.5), "mm")), 
    rowhead = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                   padding = unit(c(1.5, 1.5), "mm"))
  )
  grid.text(info, x = 1, y = 1, gp = gpar(fontsize = 8), just = c("left","top"))
  pushViewport(viewport(x = 1.025, y = 0.5, width = 1-heatmap_width_ratio, height = 0.15, just = c("left","top")))
  grid.table(col_order_table, rows = NULL, cols = NULL, theme = table_theme)
  popViewport()
  dev.off()
  
  write.xlsx(data_integrated, file.path(output_folder, heatmapdata_path), rowNames = FALSE)
  #RETURN
  return(list(col_order_table_noclustering = col_order_table_noclustering, col_order_table = col_order_table))
}




heatmap_analysis <- function(
    input_list,
    condition_1,
    condition_2,
    output_folder,
    folder_name_list,  # New parameter
    custom_color = NULL,
    using_package_f,
    scale_function,
    genic_region,
    log2FC,
    FDR,
    tpm_cutoff,
    cluster_number_row,
    cluster_number_col,
    groups_order_for_noclustering,
    individuals_order_for_noclustering,
    annotation_colors) {
  # CHECKING INPUT
  if (missing(input_list) || is.null(input_list)) {
    stop("Error: 'input_list' is mandatory and must be a file path or a list of file paths.")
  }
  if (!is.character(input_list)) {
    stop("Error: 'input' must be a list of character vectors.")
  }
  if (any(!file.exists(input_list))) {
    stop("Error: One or more input files do not exist. Please check the file paths.")
  }
  # CHECKING OUTPUT_FOLDER
  if (missing(output_folder) || is.null(output_folder)) {
    stop("Error: 'output_folder' is mandatory and must be a valid folder path.")
  }
  if (!dir.exists(output_folder)) {
    stop("Error: The specified 'output_folder' directory does not exist.")
  }
  # CHECKING log2FC
  if (!is.numeric(log2FC) || log2FC <= 0) {
    stop("Error: 'log2FC' must be a numeric value greater than 0.")
  }
  # CHECKING FDR
  if (!is.numeric(FDR) || FDR <= 0 || FDR > 1) {
    stop("Error: 'FDR' must be a numeric value between 0 and 1.")
  }
  # CHECKING USING_PACKAGE_F
  if (!using_package_f %in% c("EdgeR", "DESeq2")) {
    stop("Error: 'using_package_f' must be either 'EdgeR' or 'DESeq2'.")
  }
  # CHECK cluster_number_row
  if (is.null(cluster_number_row) || length(cluster_number_row) == 0) {
    cluster_number_row <- rep(1, length(input_list) + 2)
  } else if (length(cluster_number_row) != length(input_list) + 2) {
    stop("'cluster_number_row' must be a numeric vector of length equal to the number of input_list files +1.\n e.g: cluster_number_row = c(input1_ClusterNumber, input2_ClusterNumber, ..., integrated_ClusterNumber)")
  }
  # CHECK cluster_number_col
  if (is.null(cluster_number_col) || length(cluster_number_col) == 0) {
    cluster_number_col <- rep(1, length(input_list) + 2)
  } else if (length(cluster_number_col) != length(input_list) + 2) {
    stop("'cluster_number_col' must be a numeric vector of length equal to the number of input_list files +1.\n e.g: cluster_number_col = c(input1_ClusterNumber, input2_ClusterNumber, ..., integrated_ClusterNumber)")
  }
  # default value of custom_color
  if (is.null(custom_color) || length(custom_color) == 0 || length(custom_color) == 1) {
    custom_color <- c("red", "yellow", "#3fdc04")
  }
  custom_color <- colorRampPalette(custom_color)(200)
  # initialize storing value
  all_col_order_tables <- list()
  summary_data <- data.frame(File = character(), Upregulated = integer(), Downregulated = integer(), stringsAsFactors = FALSE)
  long_pdfs <- c()
  short_pdfs <- c()
  all_id_lists <- list()
  all_data <- list()
  all_group_data <- list()
  # Create folder
  tmp_path <- file.path(output_folder, "tmp/")
  dir.create(tmp_path, recursive = TRUE, showWarnings = FALSE)
  for (step in seq_along(input_list)) {
    
    print(folder_name_list[step])
    result <- pair_heatmap(
      input = input_list[step], 
      condition_1 = condition_1[step], 
      condition_2 = condition_2[step], 
      output_folder = output_folder, 
      custom_color = custom_color, 
      using_package_f = using_package_f, 
      scale_function = scale_function,
      genic_region = genic_region,
      log2FC = log2FC, 
      FDR = FDR, 
      tpm_cutoff = tpm_cutoff,
      cluster_number_row = cluster_number_row[step], 
      cluster_number_col = cluster_number_col[step],
      order_for_noclustering = individuals_order_for_noclustering,
      annotation_colors = annotation_colors,
      folder_name = folder_name_list[step] 
    )
    
    all_col_order_tables[[folder_name_list[step]]] <- result[[6]]
    long_pdfs <- c(long_pdfs, result[[7]])
    short_pdfs <- c(short_pdfs, result[[8]])
    all_id_lists[[step]] <- result[[3]]
    all_data[[step]] <- result[[1]]
    all_group_data[[step]] <- result[[2]]
    
    # get upregulated and downregulated number
    summary_data <- rbind(summary_data, data.frame(
      File = folder_name_list[step],
      `Condition_1(CONTROL)` = condition_1[step],
      `Condition_2(TREATMENT)` = condition_2[step],
      Downregulated = result[[4]],
      Upregulated = result[[5]]
    ))
  }
  
  
  union_id_list <- unique(unlist(all_id_lists))
  # GET G code
  all_g_codes <- unique(unlist(lapply(all_col_order_tables, function(tbl) {
    g_codes <- gsub(".*_(G[0-9]+)_.*", "\\1", tbl$COL_NAME)  
    g_codes[grepl("^G[0-9]+$", g_codes)]  
  })))
  save_name <- paste(sort(all_g_codes), collapse = "_")
  # Combined pdf long and short
  long_pdf_output <- file.path(output_folder, paste0("Heatmap_AllPairwiseComparisons_Long_", save_name, ".pdf"))
  short_pdf_output <- file.path(output_folder, paste0("Heatmap_AllPairwiseComparisons_Short_", save_name, ".pdf"))
  system(paste("pdftk", paste(long_pdfs, collapse = " "), "cat output", long_pdf_output))
  system(paste("pdftk", paste(short_pdfs, collapse = " "), "cat output", short_pdf_output))
  unlink(tmp_path, recursive = TRUE, force = TRUE)
  
  #HERE
  print("Do the individuals integrated analysis")
  individuals_result <- heatmap_integrated(union_id_list = union_id_list,
                                           all_data = all_data,
                                           save_name = save_name,
                                           tpm_cutoff = tpm_cutoff,
                                           log2FC = log2FC,
                                           FDR = FDR,
                                           scale_function = scale_function,
                                           genic_region = genic_region,
                                           sample_type = "individuals",
                                           output_folder = output_folder,
                                           custom_color = custom_color,
                                           cluster_number_row = cluster_number_row,
                                           cluster_number_col = cluster_number_col,
                                           annotation_colors = annotation_colors,
                                           order_for_noclustering = individuals_order_for_noclustering)
  print("Do the groups integrated analysis")
  groups_result <- heatmap_integrated(union_id_list = union_id_list, 
                                      all_data = all_group_data, 
                                      save_name = save_name,
                                      tpm_cutoff = tpm_cutoff,
                                      log2FC = log2FC,
                                      FDR = FDR,
                                      scale_function = scale_function, 
                                      genic_region = genic_region,
                                      sample_type = "groups", 
                                      output_folder = output_folder, 
                                      custom_color = custom_color, 
                                      cluster_number_row = cluster_number_row, 
                                      cluster_number_col = cluster_number_col,
                                      annotation_colors = annotation_colors,
                                      order_for_noclustering = groups_order_for_noclustering)
  
  summary_filepath <- file.path(output_folder, paste0("Summary_DEGs_SampleOrder_", save_name, "_log2FC_", log2FC, "_FDR_", FDR, ".xlsx"))
  excel_sheets <- list(
    summary_DEG = summary_data,  
    Individuals_NoCluster_ColOrder = individuals_result[["col_order_table_noclustering"]],  
    Individuals_Cluster_ColOrder = individuals_result[["col_order_table"]],
    Groups_NoCluster_ColOrder = groups_result[["col_order_table_noclustering"]],  
    Groups_Cluster_ColOrder = groups_result[["col_order_table"]]
  )
  
  shortened_sheet_name <- function(name) {
    name <- paste0(file_path_sans_ext(basename(name)),"_Order")
    if (nchar(name) > 31) {
      return(substr(name, nchar(name) - 30, nchar(name)))  
    } else {
      return(name)
    }
  }
  shortened_names <- lapply(names(all_col_order_tables), shortened_sheet_name)
  names(all_col_order_tables) <- shortened_names
  for (sheet_name in names(all_col_order_tables)) {
    excel_sheets[[sheet_name]] <- all_col_order_tables[[sheet_name]]
  }
  write.xlsx(excel_sheets, file = summary_filepath, rowNames = FALSE)
  return(summary_data)
}

do_heatmap <- function(
    input_folders,
    output_folder,
    custom_color = NULL,
    genic_region = "ExonCollapsed", # 'ExonCollapsed', 'FullGeneBody', or 'both'
    using_package_f = "EdgeR",  # 'EdgeR', or 'DESeq2'
    scale_function = "both", # "maxscale", "zscale", or "both"
    abs_log2FC = c(1,2),
    FDR = 0.05,
    tpm_cutoff = c(0,1),
    cluster_number_row = NULL,
    cluster_number_col = NULL) {
  
  if (is.character(input_folders) && length(input_folders) == 1) {
    if (!dir.exists(input_folders)) {
      stop("Error: The specified input folder does not exist.")
    }
    
    label_path <- file.path(input_folders, "00_Setup_Pipeline", "Sample_Labels.txt")
    sample_label <- read.table(label_path, header = TRUE, sep = ";", stringsAsFactors = FALSE)
    
    if (!all(c("Group", "Condition_Name", "Sample_ID") %in% colnames(sample_label))) {
      stop("Error: Sample_Labels.txt must contain columns: Group, Condition_Name, Sample_ID.")
    }
    # get color for each group
    group_list <- unique(sample_label$Group)
    group_colors <- color_selected(length(group_list))
    names(group_colors) <- group_list
    annotation_colors <- list(Group = group_colors)
    
    #get order for unclustering group 
    sample_label <- sample_label|>
      select(Group,Condition_Name,Sample_ID)
    groups_order_for_noclustering <- sample_label %>%
      distinct(Group, Condition_Name) %>%
      arrange(Group)
    
    #get order for unclustering individuals
    individuals_order_for_noclustering <- sample_label %>%
      mutate(Condition_Name = paste0(Condition_Name, "_", Sample_ID)) %>%
      select(Group, Condition_Name) %>%
      arrange(Group)
    
    #get all 09d folder
    all_subdirs <- list.dirs(input_folders, full.names = TRUE, recursive = FALSE)
    input_folders <- all_subdirs[grepl("^.*/09d", all_subdirs)]
    
    if (length(input_folders) == 0) {
      stop("Error: No subdirectories starting with '09d' found in the specified input folder.")
    }
    
    print(paste("Found", length(input_folders), "folders starting with '09d':"))
    print(input_folders)
    
  }
  
  else{groups_order_for_noclustering = NULL
  individuals_order_for_noclustering = NULL
  annotation_colors = NULL
  }
  
  if (!is.vector(input_folders) || !is.character(input_folders)) {
    stop("Error: 'input_folders' must be a character vector (list of folder paths).")
  }
  
  # Validate genic_region
  valid_genic_regions <- c("ExonCollapsed", "FullGeneBody", "both")
  if (!(genic_region %in% valid_genic_regions)) {
    stop("Error: 'genic_region' must be one of 'ExonCollapsed', 'FullGeneBody', or 'both'.")
  }
  
  # Validate scale_function
  valid_scale_functions <- c("maxscale", "zscale", "both")
  if (!(scale_function %in% valid_scale_functions)) {
    stop("Error: 'scale_function' must be one of 'maxscale', 'zscale', or 'both'.")
  }
  
  # Define the genic_region and scale_function combinations
  genic_regions <- if (genic_region == "both") c("ExonCollapsed", "FullGeneBody") else c(genic_region)
  scale_functions <- if (scale_function == "both") c("maxscale", "zscale") else c(scale_function)
  
  # Normalize output_folder path
  output_folder <- normalizePath(output_folder, winslash = "/", mustWork = FALSE)
  
  # Extract folder names from input_folders
  folder_name_list <- basename(input_folders)
  
  master_summary <- data.frame()
  
  # Iterate over each genic_region and scale_function
  for (gr in genic_regions) {
    for (sf in scale_functions) {
      for (lfc in abs_log2FC){
        for(tpm_co in tpm_cutoff){
        # Create output folder for this combination (only once)
        new_output_folder <- file.path(output_folder, paste0(gr, "_", sf, "_log2FC", lfc,"_FDR", FDR, "_TpmCutoff", tpm_co))
        if (!dir.exists(new_output_folder)) {
          dir.create(new_output_folder, recursive = TRUE)
        }
        
        # Initialize input and condition lists
        input_list <- c()
        condition_1_list <- c()
        condition_2_list <- c()
        folder_names <- c()  # Store corresponding folder names
        
        # Loop through each input_folder and collect matching files
        for (i in seq_along(input_folders)) {
          input_folder <- normalizePath(input_folders[i], winslash = "/", mustWork = FALSE)
          
          # Skip if input_folder doesn't exist
          if (!dir.exists(input_folder)) {
            warning(paste("Warning: The input folder", input_folder, "does not exist. Skipping."))
            next
          }
          
          # Find Condition_1.txt and Condition_2.txt
          condition_1_file <- file.path(input_folder, "Condition_1.txt")
          condition_2_file <- file.path(input_folder, "Condition_2.txt")
          
          # Ensure both condition files exist
          if (!file.exists(condition_1_file)) {
            warning(paste("Warning: Condition_1.txt not found in", input_folder, ". Skipping."))
            next
          }
          if (!file.exists(condition_2_file)) {
            warning(paste("Warning: Condition_2.txt not found in", input_folder, ". Skipping."))
            next
          }
          
          # Read condition files
          condition_1_data <- read.table(condition_1_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
          condition_2_data <- read.table(condition_2_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
          
          # Ensure 'Description' column exists
          if (!"Description" %in% colnames(condition_1_data)) {
            warning(paste("Warning: 'Condition_1.txt' in", input_folder, "is missing 'Description' column. Skipping."))
            next
          }
          if (!"Description" %in% colnames(condition_2_data)) {
            warning(paste("Warning: 'Condition_2.txt' in", input_folder, "is missing 'Description' column. Skipping."))
            next
          }
          
          # Extract conditions
          condition_1 <- condition_1_data$Description[1]
          condition_2 <- condition_2_data$Description[1]
          
          # Locate the genic_region-specific subdirectory
          matching_dirs <- list.dirs(input_folder, full.names = TRUE, recursive = FALSE)
          matching_dirs <- matching_dirs[grepl(paste0(gr, "$"), matching_dirs)]
          
          if (length(matching_dirs) == 0) {
            warning(paste("Warning: No subdirectories ending in", gr, "found in", input_folder))
            next
          }
          
          # Collect all DiffExp files from the matching directories
          for (dir_name in matching_dirs) {
            diffexp_files <- list.files(dir_name, pattern = "^DiffExp", full.names = TRUE, recursive = FALSE)
            
            # Store input files and corresponding conditions
            for (file in diffexp_files) {
              input_list <- c(input_list, file)
              condition_1_list <- c(condition_1_list, condition_1)
              condition_2_list <- c(condition_2_list, condition_2)
              folder_names <- c(folder_names, folder_name_list[i])  # Store folder name
            }
          }
        }
        
        # Skip processing if no valid DiffExp files were found
        if (length(input_list) == 0) {
          warning(paste("Warning: No 'DiffExp' files found for Gene Type:", gr, "and Scale Function:", sf))
          next
        }
        print(paste0("Processing: Gene Region = ", gr,", Scale Function = ", sf, ", log2FC = ", lfc, ", FDR = ", FDR, ", tpm_cutoff = ", tpm_co))
        # Call heatmap_analysis
        summary_data <- heatmap_analysis(
          input_list = input_list,
          condition_1 = condition_1_list,
          condition_2 = condition_2_list,
          output_folder = new_output_folder,
          folder_name_list = folder_names,  # Pass folder names
          custom_color = custom_color,
          using_package_f = using_package_f,
          genic_region = gr,
          scale_function = sf,
          log2FC = lfc,
          FDR = FDR,
          tpm_cutoff = tpm_co,
          cluster_number_row = cluster_number_row,
          cluster_number_col = cluster_number_col,
          groups_order_for_noclustering = groups_order_for_noclustering,
          individuals_order_for_noclustering = individuals_order_for_noclustering,
          annotation_colors = annotation_colors
        )
        summary_data$genic_region <- gr
        summary_data$abs_log2FC <- lfc
        summary_data$FDR <- FDR
        summary_data$tpm_cutoff <- tpm_co
        master_summary <- rbind(master_summary, summary_data)
        }
      }
    }
  }
    master_summary_path <- file.path(output_folder, "Master_Summary_DEG.xlsx")
    write.xlsx(master_summary, master_summary_path, rowNames = FALSE)
content <- "
Function: do_heatmap

Description:
This function processes differential expression data from input folders, applies 
log2 fold-change and FDR filtering, and optional TPM-based filtering, performs clustering, and normalizes the data 
using selected methods. The output includes heatmaps and summary reports.

Usage:
do_heatmap(
    input_folders,  
    output_folder,
    custom_color = NULL,
    genic_region = \"ExonCollapsed\",
    using_package_f = \"EdgeR\",
    scale_function = \"both\",
    abs_log2FC = c(1,2),
    FDR = 0.05,
    tpm_cutoff = c(0,1),
    cluster_number_row = NULL,
    cluster_number_col = NULL
)

Arguments:
- input_folders: Required. Path(s) to input folders. Can be a single path or a vector of multiple paths.
  If a single path is provided, all subfolders starting with '09d' will be processed and will automatically get Sample_Labels.txt file for Nonclustering's heatmap order.
  Example:
      input_folders <- \"/projectnb/wax-dk/max/G216_G221/Scripts/\"
      OR
      input_folders <- c(\"/projectnb/wax-dk/max/G216_G221/Scripts/09d_DE_29_RefSeqLncRNA76k\",
                         \"/projectnb/wax-dk/max/G216_G221/Scripts/09d_DE_30_RefSeqLncRNA76k\")

- output_folder: Required. Path where output files will be generated. Subdirectories will be 
  created under this folder formatted as output_folder/<genic_region>_<scale_function>_<log2FC>_<FDR>_<tpm_cutoff>.

- custom_color: Optional. Color gradient for heatmaps. Provide a vector specifying colors from low to high, 
  e.g., c(\"#123456\",\"#abcdef\",\"green\"). Default: red → yellow → green.

- genic_region: Optional. Region of genes to process. Options: \"ExonCollapsed\", \"FullGeneBody\", 
  or \"both\". Default: \"both\".

- using_package_f: Optional. Specifies which algorithm’s log2FC and FDR to use. Options: \"EdgeR\" or \"DESeq2\". 
  Default: \"EdgeR\".

- scale_function: Optional. Type of normalization. Options: \"maxscale\", \"zscale\", or \"both\". 
  Default: \"both\".
  \"maxscale\" = x/max(x)
  \"zscale\" = (x-mean(x))/std(x)/max(max((x-mean(x))/std(x)),-min((x-mean(x))/std(x)))

- abs_log2FC: Optional. Log2 fold change threshold. Must be greater than 0. Default: 1, 2.

- FDR: Optional. False Discovery Rate threshold. Must be between 0 and 1. Default: 0.05.

- tpm_cutoff: Optional. TPM cutoff. All samples with TPM larger than tpm_cutoff will be take into analysis. Default: 0, 1.

- cluster_number_row: Optional. Number of different parameters to set for row clustering. One parameter is needed for each pairwise comparison, one parameter for All Groups-based heatmaps, and one parameter for All Individuals-based heatmaps. The length of the vector that presents these parameters is equal to the number of pairwise comparison inputs plus 2. Format: c(input1_cluster_number,input2_cluster_number,input3_cluster_number,...,groups_cluster_number,individuals_cluster_number). Default: NULL.

- cluster_number_col: Optional. Number of column clusters. Must be similar to cluster_number_row. Default: NULL.

Details:
This function processes differential expression data from input folders, applies 
log2 fold-change and FDR filtering, performs clustering, and normalizes the data 
using selected methods. The output includes heatmaps and summary reports.

Output Files:
PDF Files:
- Heatmap_AllPairwiseComparisons_long_<GCODE>.pdf: Merged heatmaps for all pairs (original format).
- Heatmap_AllPairwiseComparisons_short_<GCODE>.pdf: Merged heatmaps for all pairs (formatted for printing).
- Heatmap_AllIndividuals_SamplesClustered_long_<GCODE>.pdf: Merged heatmap for all samples (original format).
- Heatmap_AllIndividuals_SamplesClustered_short_<GCODE>.pdf: Merged heatmap for all samples (formatted for printing).
- Heatmap_AllIndividuals_NoSampleClustering_long_<GCODE>.pdf: Merged heatmap for all samples without clustering (original format).
- Heatmap_AllIndividuals_NoSampleClustering_short_<GCODE>.pdf: Merged heatmap for all samples without clustering (formatted for printing).
- Heatmap_AllGroups_SamplesClustered_long_<GCODE>.pdf: Merged heatmap for all sample sets, using mean tpm (original format).
- Heatmap_AllGroups_SamplesClustered_short_<GCODE>.pdf: Merged heatmap for all sample sets, using mean tpm (formatted for printing).
- Heatmap_AllGroups_NoSampleClustering_long_<GCODE>.pdf: Merged heatmap for all sample sets, using mean tpm without clustering (original format).
- Heatmap_AllGroups_NoSampleClustering_short_<GCODE>.pdf: Merged heatmap for all samplesets, using mean tpm without clustering (formatted for printing).

Excel Files:
- HeatmapData_<input>_log2FC_<log2FC>_FDR_<FDR>.xlsx: Filtered information for each input pair, including 
  heatmap gene order, gene names, scaled TPM values, count values, TPM values, fold change, and FDR.
- HeatmapData_AllIndividuals_<GCODE>_log2FC_<log2FC>_FDR_<FDR>.xlsx: Merged filtered information for all input pairs, 
  including the union of all selected genes.
- HeatmapData_AllGroups_<GCODE>_log2FC_<log2FC>_FDR_<FDR>.xlsx: Merged filtered information for all input sample sets, 
  including the union of all selected genes, using mean tpm.
- summary_DEGs_SampleOrder_<GCODE>_log2FC_<log2FC>_FDR_<FDR>.xlsx: Summary of differential expression genes, 
  sample orders before and after clustering, and sample orders for each pair.

Examples:
do_heatmap(
    input_folders = c(\"/projectnb/wax-dk/max/G216_G221/Scripts\"),
    output_folder = \"/projectnb/wax-dk/max/G216_G221/Results\",
    custom_color = c(\"red\", \"yellow\", \"blue\"),
    gene_type = \"ExonCollapsed\",
    using_package_f = \"EdgeR\",
    scale_function = \"maxscale\",
    abs_log2FC = 1,
    FDR = 0.05
)
"
writeLines(content, file.path(output_folder,"README.txt"))
}

# Get variable from terminal
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
  input_folders <- args[1]
  output_folder <- args[2]
  
  do_heatmap(
    input_folders = input_folders,
    output_folder = output_folder
  )
} else {
  stop("Please provide input_folders and output_folder as command line arguments")
}

file.remove(file.path(output_folder, "Rplots.pdf"))
