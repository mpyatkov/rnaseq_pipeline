library(readr)
library(dplyr)
library(RColorBrewer)
library(VennDiagram)
library(grid)
library(gridExtra)
library(stringr)

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
    scale = F,
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
    new_names <- paste0(names, c("_up","_down"))
  } else if (type == "downup") {
    edger_data <- list(data[[1]]$down_edger, data[[2]]$up_edger)
    deseq_data <- list(data[[1]]$down_deseq, data[[2]]$up_deseq)
    new_names <- paste0(names, c("_down","_up"))
  } else {
    edger_data <- lapply(data, function(x) {x[[paste0(type,"_edger")]]})
    deseq_data <- lapply(data, function(x) {x[[paste0(type,"_deseq")]]})
    new_names <- paste0(names, "_", type)
  }
  
  print(length(unlist(edger_data)))

  label<-ifelse(type=="all", "all_DE", type)

  edger <- drawVenn(edger_data, new_names, col, 
                    parameters[[num]]$pos, parameters[[num]]$dist, parameters[[num]]$main_title_height,
                    paste0("edgeR_", label, " genes"), main_color)
  
  deseq <- drawVenn(deseq_data, new_names, col, 
                    parameters[[num]]$pos, parameters[[num]]$dist, parameters[[num]]$main_title_height,
                    paste0("DESeq_", label, " genes"), main_color)

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
  
  pdf(file=output_name, height = height, width = width)
  grid.draw(plot)
  dev.off()
}

#------- script start point 
args <- commandArgs(T)
DATASET_LABEL <- args[1]
OUTPUT_NAME <- args[2]
# DATASET_LABEL <- "TEST_DATASET"
# OUTPUT_NAME <- "test.pdf"

# parameters for circles
params <- c()
# parameters for 2 circles
params[[2]] <- list(
  pos = c(-180, 0),
  dist = c(0.03, 0.03),
  main_title_height = 0.95
)
# parameters for 3 circles
params[[3]] <- list(
  pos = c(-24, 24, 135),
  dist = c(0.055, 0.055, 0.085),
  main_title_height = 1.05
)
# parameters for 4 circles
params[[4]] <- list(
  pos = c(0, 0, -30, 30),
  dist = c(-0.37, -0.33, 0.1, 0.1),
  main_title_height = 0.95
)

# read files from the current directory
files <- list.files(pattern = "DiffExp[[_]|[:alnum:]]+.txt", full.names = T)

# extract group names (ex. *_ExonCollapsed_LONG..NAME_featureCounts)
names <- str_extract(files, "(?<=(ExonCollapsed_))([[:alnum:]|[_]]+)(?=(_featureCounts|_htseq))")

# extract comparison numbers: 1_Diff..., 3_Diff... -> 1,3,..
numbers <- str_extract(basename(files), "^(\\d)+")

# read data
data <- lapply(files, function(x) {read_dataset(x, 2, 0.05)})

# draw all
draw_all(data,names, DATASET_LABEL, numbers, params, OUTPUT_NAME)
