# creating 3 types of Venn diagrams

# if VENN_INDIVIDUAL = 0 then
# create 1 venn diagramm
# draw_all: pairwise comparison of DiffExp_v2 files (ExonCollapsed only) 
# from different samples
# ex. G186_M2/G186_M1 vs G186_M3/G186_M1 
# ex. or 09a_1 vs 09a_2

# if VENN_INDIVIDUAL = 1 then
# create 2 venn diagrams for each individual folder (ex.09a_1)
# ex. venn diagram will be created for each feature (ExonCollapsed/IntronOnly...) set
# 1 type of venn - deseq/edger for each feature separately
# 2 type of venn - all feature together
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(RColorBrewer)
    library(VennDiagram)
    library(grid)
    library(gridExtra)
    library(stringr)
})

# convert file names to correct titles for venn diagrams
get_correct_names <- function(comparisons, samples, names) {
  library(tidyr)
  cmp <- read_delim(comparisons, delim = ";", col_names = T, trim_ws = T) %>%
    select(comparison = Comparison_Number, Condition_1, Condition_2)
  
  samples <- read_delim(samples, delim=";", col_names = T, trim_ws = T) %>%
    select(Group, Condition_Name) %>%
    distinct(Group, Condition_Name)
  
  cmp_samples <- left_join(cmp, samples, by=c("Condition_1"="Group")) %>%
    left_join(., samples, by=c("Condition_2"="Group")) %>% 
    select(comparison, Condition_1 = Condition_Name.x, Condition_2 = Condition_Name.y) %>% 
    mutate(old_name = paste0(Condition_1,"_",Condition_2),
           new_name = paste0(Condition_2," / \n", Condition_1)) %>% 
    select(comparison, old_name, new_name)
  
  left_join(tibble(old_name = names), cmp_samples) %>% pull(new_name)
}

# get the subset of up and down genes for deseq and edger
read_dataset <- function(dn, fc, padj) {
  df <- read_tsv(dn) %>% 
    select(id, DESeq2_padj_FDR, DESeq2_foldChange, EdgeR_padj_FDR, EdgeR_foldChange) %>% 
    mutate(edger.fc = ifelse(EdgeR_foldChange < 1, -1*(1/EdgeR_foldChange), EdgeR_foldChange),
           deseq.fc = ifelse(DESeq2_foldChange < 1, -1*(1/DESeq2_foldChange), DESeq2_foldChange))
  
  up_edger <- df %>% filter(edger.fc > fc & EdgeR_padj_FDR < padj) %>% pull(id)
  dn_edger <- df %>% filter(edger.fc < -fc & EdgeR_padj_FDR < padj) %>% pull(id)
  up_deseq <- df %>% filter(deseq.fc > fc & DESeq2_padj_FDR < padj) %>% pull(id)
  dn_deseq <- df %>% filter(deseq.fc < -fc & DESeq2_padj_FDR < padj) %>% pull(id)
  all_deseq <- df %>% filter(abs(deseq.fc) > fc & DESeq2_padj_FDR < padj) %>% pull(id)
  all_edger <- df %>% filter(abs(edger.fc) > fc & EdgeR_padj_FDR < padj) %>% pull(id)
  
  list(up_edger = up_edger,
       down_edger = dn_edger,
       up_deseq = up_deseq,
       down_deseq = dn_deseq,
       all_deseq = all_deseq,
       all_edger = all_edger)
}

# ex. we can use the same
# log2fc > 1 -- up genes EdgeR_padj_FDR < 0.05
# log2fc < -1 -- down genes EdgeR_padj_FDR < 0.05
# fc > 2 -- up genes EdgeR_padj_FDR < 0.05
# fc < -2 -- down genes EdgeR_padj_FDR < 0.05

# get list of colors from 2 to 4 by type all/up/down/updown/downup
get_colors <- function(n, type = "all") {
  if (type == "up") {
    c("orangered", "red", "salmon","pink2")[1:n]
  } else if (type == "down") {
    c("cornflowerblue", "cadetblue1", "deepskyblue1", "lightblue2")[1:n]
  } else if (type == "updown") {
    # only for two colors
    # c("orangered","cadetblue1")
    c("red","cadetblue1")
  } else if (type == "downup") {
    # only for two colors
    # c("cornflowerblue", "red")
    c("cadetblue1", "red")
  }
  else {
    colorRampPalette(brewer.pal(3, "Pastel2"))(n)
  }
}

# draw only one diagram
drawVenn <- function(data, names, col, pos, dist, main_title_height, main_title, main_col) {

  # do not create Rplots.pdf
  pdf(NULL)
  
  # switch off the logger
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  
  res <- venn.diagram(
    x = data,
    category.names = names,
    # filename = fname,
    filename = NULL,
    output=F,
    
    # Output features
    imagetype="png" ,
    height = 1024 ,
    width = 1024 ,
    resolution = 600,
    compression = "lzw",
    scaled = F,
    # inverted = T,
    # margin = 0.07,

    # main title
    main = main_title,
    main.fontface = "bold",
    main.cex = 1.5,
    main.fontfamily = "sans",
    main.pos = c(0.5, main_title_height),
    main.col = main_col,

    # Circles
    lwd = 1,
    # lty = 'blank',
    lty = 'solid',
    fill = col,
    
    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 1.1,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = pos,
    cat.dist = dist,
    cat.fontfamily = "sans"
    
    # diagram can be rotated, if it will required for 2 comparisons
    # rotation.degree = ifelse(length(data[[1]]) > length(data[[2]]),180,1)
  )
  
  # "close bracket" for pdf(null) 
  # dev.off()
  
  res
}

# draw two diagrams side-by-side by type (up/down/all)
# data <- list(x1, x2, x3)
# left diagram -> edger, right -> deseq
draw_pair <- function(data, names, type, parameters) {
  
  main_color<- if (type == "up") {
    "red"
  } else if (type == "down") {
    "blue"
  } else {
    "black"
  }
  
  num <- length(data)
  col <- get_colors(num, type)
  
  if (type == "updown") {
    edger_data <- list(data[[1]]$up_edger, data[[2]]$down_edger)
    deseq_data <- list(data[[1]]$up_deseq, data[[2]]$down_deseq)
    # new_names <- paste0(names, c("_up","_down"))
    new_names <- names
  } else if (type == "downup") {
    edger_data <- list(data[[1]]$down_edger, data[[2]]$up_edger)
    deseq_data <- list(data[[1]]$down_deseq, data[[2]]$up_deseq)
    # new_names <- paste0(names, c("_down","_up"))
    new_names <- names
  } else {
    edger_data <- lapply(data, function(x) {x[[paste0(type,"_edger")]]})
    deseq_data <- lapply(data, function(x) {x[[paste0(type,"_deseq")]]})
    # new_names <- paste0(names, "_", type)
    new_names <- names
  }
  
  print(length(unlist(edger_data)))
  
  label<-ifelse(type=="all", "all_DE", type)
  
  label <- if (type == "updown" || type == "downup") {
    paste0(label, " genes\n(red - up, blue -down)")
  } else {
    paste0(label, " genes")
  }
  
  edger <- drawVenn(edger_data, new_names, col, 
                    parameters[[num]]$pos, parameters[[num]]$dist, parameters[[num]]$main_title_height,
                    paste0("edgeR_", label), main_color)
  
  deseq <- drawVenn(deseq_data, new_names, col, 
                    parameters[[num]]$pos, parameters[[num]]$dist, parameters[[num]]$main_title_height,
                    paste0("DESeq_", label), main_color)
  
  grid.arrange(gTree(children=edger), gTree(children=deseq), ncol=2)
}

# draw full set diagrams (5 diagrams - for set of 2 comparisons, 3 diagrams - for 3 and 4 comparisons)
draw_all <- function(data, names, dataset_label, num_comparisons, parameters, output_name) {
  num <- length(data)
  
  all <- draw_pair(data, names, "all", parameters)
  up <- draw_pair(data, names, "up", parameters)
  down <- draw_pair(data, names, "down", parameters)
  
  if (num == 2) {
    updown <- draw_pair(data, names, "updown", parameters)
    downup <- draw_pair(data, names, "downup", parameters)
  }
  
  #--------------
  
  if (num == 2) {
    plot <- grid.arrange(gTree(children=gList(all)),
                         gTree(children=gList(up)),
                         gTree(children=gList(down)),
                         gTree(children=gList(updown)),
                         gTree(children=gList(downup)),
                         ncol=1, nrow=5)
  } else {
    plot <- grid.arrange(gTree(children=gList(all)),
                         gTree(children=gList(up)),
                         gTree(children=gList(down)),
                         ncol=1, nrow=3)
  }

  title <- paste0("Dataset: ", dataset_label,
                  " Comparisons: ", paste0(num_comparisons, collapse = " vs "))
  plot <- grid.arrange(gTree(children=gList(plot)), 
                       top=textGrob(title, 
                                    gp=gpar(fontsize=18, fontfamily="sans", fontface="bold")))
  
  # 6x6 one square, 6x12 one line
  height <- ifelse(num != 2,18,30) 
  width <- 12
  
  pdf(file=paste0(output_name,".pdf"), height = height, width = width)
  grid.draw(plot)
  dev.off()
}

# short version of the draw_pair, only for methods comparison
draw_pair_methods <- function(data, names, parameters, title) {
  
  num <- 2
  col_up <- get_colors(num, "up")
  col_down <- get_colors(num, "down")
  
  up_genes <- list( data[["up_edger"]], data[["up_deseq"]])
  down_genes <- list( data[["down_edger"]], data[["down_deseq"]])
  
  up_plot <- drawVenn(up_genes, c("up_edgeR", "up_DESeq"), col_up, parameters[[num]]$pos, 
                    parameters[[num]]$dist, parameters[[num]]$main_title_height,
                    title, "red")
  
  down_plot <- drawVenn(down_genes, c("down_edgeR","down_DESeq"), col_down, parameters[[num]]$pos, 
                      parameters[[num]]$dist, parameters[[num]]$main_title_height,
                      title, "blue")
  
  grid.arrange(gTree(children=up_plot), gTree(children=down_plot), ncol=2)
}

# combined plot of all features together for up/down and deseq/edger
draw_all_features_methods <- function(data, names, parameters) {
  num <- length(data)
  
  col_up <- get_colors(num, "up")
  col_down <- get_colors(num, "down")
  
  edger_data_up <- lapply(data, function(x) {x[[paste0("up","_edger")]]})
  deseq_data_up <- lapply(data, function(x) {x[[paste0("up","_deseq")]]})
  edger_data_down <- lapply(data, function(x) {x[[paste0("down","_edger")]]})
  deseq_data_down <- lapply(data, function(x) {x[[paste0("down","_deseq")]]})
  
  title <- "DESeq"
  up_deseq_plot <- drawVenn(deseq_data_up, paste0("up_", names), col_up, parameters[[num]]$pos, 
                            parameters[[num]]$dist, parameters[[num]]$main_title_height,
                            title, "red")
  
  title <- "edgeR"
  up_edger_plot <- drawVenn(edger_data_up, paste0("up_", names), col_up, parameters[[num]]$pos, 
                            parameters[[num]]$dist, parameters[[num]]$main_title_height,
                            title, "red")

  title <- "DESeq"
  down_deseq_plot <- drawVenn(deseq_data_down, paste0("down_", names), col_down, parameters[[num]]$pos, 
                        parameters[[num]]$dist, parameters[[num]]$main_title_height,
                        title, "blue")
  
  title <- "edgeR"
  down_edger_plot <- drawVenn(edger_data_down, paste0("down_", names), col_down, parameters[[num]]$pos, 
                              parameters[[num]]$dist, parameters[[num]]$main_title_height,
                              title, "blue")
  
  grid.arrange(gTree(children=up_edger_plot), gTree(children=down_edger_plot), 
               gTree(children=up_deseq_plot), gTree(children=down_deseq_plot),
               nrow = 2, ncol=2)
  
}

# wrapper for function "draw_all_features_methods"
# combined plot of all features together for up/down and deseq/edger
all_features_methods <- function(files, data, params, title, out_name) {
  features <- str_extract(files, "(?<=_)(ExonCollapsed|IntronicOnly|ExonOnly|IntronOnly|ExonicOnly|FullGeneBody)(?=_)")   
  daf <- draw_all_features_methods(data, features, params)
  
  plot <- grid.arrange(gTree(children=gList(daf)), 
                       top=textGrob(title,
                                    gp=gpar(fontsize=18, fontfamily="sans", fontface="bold")))
  
  
  pdf(file=paste0(out_name,".pdf"), height = 12, width = 12)
  grid.draw(plot)
  dev.off()
}

# wrapper for "draw_pair_methods"
# separate plots for each feature individually
# edger/deseq and up/down for each future
draw_pairwise_for_each_feature <- function(files, params, title, out_name) {

  all_plots <- lapply(files, function(fname) {
    tmp <- read_dataset(fname, 2, 0.05)
    feature <- str_extract(fname, "(?<=_)(ExonCollapsed|IntronicOnly|ExonOnly|IntronOnly|ExonicOnly|FullGeneBody)(?=_)")
    draw_pair_methods(tmp, c("edgeR","DESeq"), params, feature)
  })
  # grid.arrange(rectGrob(), rectGrob())
  plot <- marrangeGrob(all_plots, 
                       nrow=2, 
                       ncol=2, 
                       top=textGrob(title, gp=gpar(fontsize=18, fontfamily="sans", fontface="bold")))
  
  
  pdf(file=paste0(out_name,".pdf"), height = 10, width = 20)
  grid.draw(plot)
  dev.off()
}


#------- script start point 
args <- commandArgs(T)
DATASET_LABEL <- args[1]
OUTPUT_NAME <- args[2]
# VENN_INDIVIDUAL=0/1
VENN_INDIVIDUAL <- args[3]

# DATASET_LABEL <- "TEST_DATASET"
# OUTPUT_NAME <- "test"
# VENN_INDIVIDUAL <- "0"

# parameters for circles
params <- c()
# parameters for 2 circles
params[[2]] <- list(
  pos = c(-180, 0),
  dist = c(0.05, 0.05),
  main_title_height = 0.95
)
# parameters for 3 circles
params[[3]] <- list(
  pos = c(-24, 24, 180),
  dist = c(0.08, 0.08, 0.08),
  main_title_height = 1.05
)
# parameters for 4 circles
params[[4]] <- list(
  pos = c(-10, 5, -30, 30),
  dist = c(-0.38, -0.38, 0.15, 0.15),
  main_title_height = 1.0
)

# read files from the current directory
files <- list.files(pattern = "DiffExp[[_]|[:alnum:]]+.txt", full.names = T)

# read data
data <- lapply(files, function(x) {read_dataset(x, 2, 0.05)})

if (VENN_INDIVIDUAL == "0") {

   # extract group names (ex. *_ExonCollapsed_LONG..NAME_featureCounts)
   names <- str_extract(files, "(?<=(ExonCollapsed_))([[:alnum:]|[_]]+)(?=(_featureCounts|_htseq))")
   names <- get_correct_names("./Comparisons.txt", "./Sample_Labels.txt", names)

  # extract comparison numbers: 1_Diff..., 3_Diff... -> 1,3,..
  numbers <- str_extract(basename(files), "^(\\d)+")
  
  # draw all 
  draw_all(data,names, DATASET_LABEL, numbers, params, OUTPUT_NAME)
} else {
  
  # files <- list.files(pattern = "DiffExp[[_]|[:alnum:]]+.txt", full.names = T)
  number = str_extract(files[1], "(?<=./)(\\d)+")
  counter = str_extract(files[1], "(?<=_)(featureCounts|htseq)(?=\\.txt)")
  features <- str_extract(files, "(?<=_)(ExonCollapsed|IntronicOnly|ExonOnly|IntronOnly|ExonicOnly|FullGeneBody)(?=_)")
  comparison_name <- str_extract(files[1], paste0("(?<=(", features[1],"_))([[:alnum:]|[_]]+)(?=(_featureCounts|_htseq))"))

  title <- paste0(DATASET_LABEL, ", ", comparison_name)

  draw_pairwise_for_each_feature(files,params, title, paste0(OUTPUT_NAME,"_Individual"))

  # Draw all features together only if multiple features exist
  print(features)
  print(paste0("Number of features: ", length(unique(features))))

  if (length(unique(features)) > 1) {
     print("Draw all features together")
     all_features_methods(files, data, params, title, paste0(OUTPUT_NAME,"_All"))  
  }
}



