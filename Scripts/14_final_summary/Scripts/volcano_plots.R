#!/usr/bin/env Rscript

DEBUG <- FALSE

library(argparser)
 
ParseArguments <- function() {
  p <- arg_parser('Volcano plots')
  p <- add_argument(p,'--segex_files_path', default=".", help="directory with segex output files")
  p <- add_argument(p, '--sample_labels', help='File with sample labels', default="Sample_Labels.txt")
  p <- add_argument(p, "--comparisons", help = "file with comparisons", default = "Comparisons.txt")
  p <- add_argument(p, "--output_prefix", help = "output prefix for each pdf file", default = "tmp")
  return(parse_args(p))
}

argv <- ParseArguments()

if (DEBUG){
  argv$segex_files_path = "/projectnb/wax-dk/max/20221110/HYPOX_76K/Scripts/14_final_summary/output/Segex_09d/Segex09d_ExonCollapsed/"
  argv$sample_labels <- "/projectnb/wax-dk/max/20221110/HYPOX_76K/Scripts/00_Setup_Pipeline/Sample_Labels.txt"
  argv$comparisons <- "/projectnb/wax-dk/max/20221110/HYPOX_76K/Scripts/00_Setup_Pipeline/Comparisons.txt"
}

print(argv)

## install lab libs
library(devtools)
# devtools::install_github("mpyatkov/NotationConverter", upgrade = "never", dependencies = FALSE)
# 
# devtools::install_github could produce the following error:
# Error: Failed to install 'NotationConverter' from GitHub:
# HTTP error 403.
# API rate limit exceeded for 192.12.187.131. (But here's the good news: Authenticated requests get a higher rate limit. Check out the documentation for more details.)
# 
#   Rate limit remaining: 0/60
#   Rate limit reset at: 2023-06-02 18:58:08 UTC
# 
#   To increase your GitHub API rate limit
#   - Use `usethis::create_github_token()` to create a Personal Access Token.
#   - Use `usethis::edit_r_environ()` and add the token as `GITHUB_PAT`.

devtools::install_local("/projectnb/wax-es/WAXMANLAB_SOFT/NotationConverter", upgrade = "never", dependencies = FALSE)
library(NotationConverter)

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(cowplot)
library(patchwork)

get_df <- function(fname, log2fc_thr, adjpval_thr, drop_not_signif = FALSE) {
  df <- read_tsv(fname, col_names = T, show_col_types = F) %>% 
    select(id = 1, ratio = 2, fc = 3, padj = 6) %>% 
    mutate(log2fc = log2(ratio))
  
  df1 <- df %>% 
    ## filter(!str_starts(id,"nc_")) %>% ## drop all lnc genes
    mutate(diffexpressed = case_when(log2fc > log2fc_thr & padj < adjpval_thr ~ "UP",
                                     log2fc < -log2fc_thr & padj < adjpval_thr ~ "DOWN",
                                     TRUE ~ "NO"),
           label = case_when(diffexpressed == "NO" ~ NA,
                             TRUE ~ id))
  
  if (drop_not_signif){
    df1 <- df1 %>% 
      filter(!is.na(label))
  }
  
  df1
}


## find min pvalue and max log2fc for all datasets to create correct limits for ggplot plots
get_minmax_for_all_df <- function(segex_path = "."){
  files <- list.files(path = segex_path, pattern = "txt", full.names = T)  
  
  tmp <- map_dfr(files, function(fname){
    inner_tmp <- get_df(fname, 1, 0.05)
    tibble(padj = min(inner_tmp$padj),
           log2fc = max(abs(inner_tmp$log2fc)))
  })
  
  min_pvalue <- min(tmp$padj)
  max_abs_x <- ceiling(max(tmp$log2fc))
  
  list(ylim = -log10(c(1, min_pvalue)),
       xlim = c(-max_abs_x, max_abs_x))
}

create_plot <- function(fname, log2fc_thr, adjpval_thr, title){
  
  # df_tmp <- get_df(fname, log2fc_thr = log2fc_thr, adjpval_thr = adjpval_thr)
  df_tmp <- get_df(fname, log2fc_thr = log2fc_thr, adjpval_thr = adjpval_thr) %>% 
    select(everything(),gname = label) %>% 
    notationConverter(., from = "segex", to = "mm10") %>% 
    select(everything(),label = gname) %>% 
    relocate(label, .after = everything())
  
  genes <- table(df_tmp$diffexpressed)
  up_genes <- genes['UP']
  down_genes <- genes['DOWN']
  
  up_genes <- ifelse(is.na(up_genes),0,up_genes)
  down_genes <- ifelse(is.na(down_genes),0,down_genes)
  
  final_subtitle <- str_glue("P-value < {adjpval_thr}. |Log2FC| > {log2fc_thr}. UP genes: {up_genes}, DOWN genes: {down_genes}")  
  
  p <- ggplot(data=df_tmp, aes(x=log2fc, y=-log10(padj), col=diffexpressed, label=label)) +
    geom_vline(xintercept=c(-log2fc_thr, log2fc_thr), col="gray") +
    geom_hline(yintercept=-log10(adjpval_thr), col="gray")+
    geom_point(size = 0.5) + 
    # geom_text_repel(show.legend  = FALSE, 
    #                 segment.color = 'darkgray') + ## fontface = 'bold', segment.linetype = 2, segment.color = 'black'
    geom_text_repel(show.legend  = FALSE,
                    box.padding = 0.3,
                    min.segment.length = 0,
                    segment.color = 'darkgray',
                    max.overlaps = 15,
                    size = 3) + ## fontface = 'bold', segment.linetype = 2, segment.color = 'black'
    #ylim(0,300)+
    ylab("-log10 (P-value)")+
    xlab("log2 (fold change)")+
    # scale_color_manual(values=c("blue", "gray", "red")) +
    scale_color_manual(values=c("blue", "red", "gray"), 
                       labels = c("DOWN","UP","Not significant"), 
                       breaks = c("DOWN","UP","Not significant"), 
                       limits = c("DOWN","UP","Not significant"))+
    scale_x_continuous(breaks = scales::pretty_breaks(10))+
    scale_y_continuous(breaks = scales::pretty_breaks(10))+
    guides(color=guide_legend(title = "Differential expression",
                              override.aes = list(size=5)))+
    theme_bw() +
    ggtitle(title, subtitle = final_subtitle) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = "bottom")
  
  p
}


## MAIN

cmp <- read_delim(argv$comparisons, delim = ";", col_names = T, trim_ws = T, comment = "#") %>%
  select(comparison = Comparison_Number, Condition_1, Condition_2)

samples <- read_delim(argv$sample_labels, delim=";", col_names = T, trim_ws = T, comment="#") %>%
  select(Group, Condition_Name) %>%
  distinct()

cmp_smp <- left_join(cmp, samples, by = c("Condition_1" = "Group")) %>% 
  select(everything(), Condition1_label = Condition_Name) %>% 
  left_join(., samples, by = c("Condition_2" = "Group")) %>% 
  select(everything(), Condition2_label = Condition_Name)  %>% 
  mutate(title = str_glue("{comparison}. {Condition2_label} vs {Condition1_label}")) %>% 
  select(comparison, title) %>% 
  deframe() ## convert to named vector

# ## Just for test
# tt <- list.files(path = argv$segex_files_path, pattern = "txt", full.names = T)
# tp <- create_plot(tt[[6]],
#             log2fc_thr = 1,
#             adjpval_thr = 0.05,
#             title = "final_title")
# 
# tp

plots <- map(list.files(path = argv$segex_files_path, pattern = "txt", full.names = T), function(full_fname){
  base_fname <- basename(full_fname)
  ix <- str_extract(base_fname,"^\\d+")
  final_title <- cmp_smp[[ix]]
  create_plot(full_fname,
              log2fc_thr = 1,
              adjpval_thr = 0.05,
              title = final_title)
})

## read all txt files and extract limits
limits <- get_minmax_for_all_df(segex_path = argv$segex_files_path)

## add limits
plots_same_scale <- map(plots, function(plot){
  plot+
    scale_x_continuous(breaks = scales::pretty_breaks(10), limits = c(limits$xlim[1],limits$xlim[2]))+
    scale_y_continuous(breaks = scales::pretty_breaks(10), limits = c(limits$ylim[1],limits$ylim[2]))
})

## save every thing
print("Saving individual scale plots...")
system.time(plots %>% 
  map(function(p){plot_grid(wrap_plots(p))}) %>% 
  marrangeGrob(nrow = 2, ncol = 1) %>%
  ggsave(filename = str_glue("{argv$output_prefix}_individual_scale_volcano.pdf"), width = 20, height = 15.50))

print("Saving common scale plots...")
system.time(plots_same_scale %>% 
  map(function(p){plot_grid(wrap_plots(p))}) %>% 
  marrangeGrob(nrow = 2, ncol = 1) %>%
  ggsave(filename = str_glue("{argv$output_prefix}_common_scale_volcano.pdf"), width = 20, height = 15.50))
