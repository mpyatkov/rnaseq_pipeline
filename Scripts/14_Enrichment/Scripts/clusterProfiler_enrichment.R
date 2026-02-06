#!/usr/bin/env Rscript

current_dir=getwd()
setwd("/projectnb/wax-es/routines/RENV_443_clusterProfiler")
source("renv/activate.R")
setwd(current_dir)

library(argparser)

p <- arg_parser('GO enrichment analysis')
p <- add_argument(p,'--input_path', default="./", help="Up/Down files path")
argv <- parse_args(p)

DEBUG <- F
if (DEBUG){
  #argv$input_file <- "/projectnb/wax-dk/max/SIA_NUS/KC14_RNASEQ/Scripts/09d_DE_1_RefSeqLncRNA76k/Output_DiffExp_1i_featureCounts_RefSeqLncRNA76k_ExonCollapsed/Up_Genes_EdgeR_RefSeqLncRNA76k_ExonCollapsed_MaleLiver_WT_KC14101WT_M01_KC14102WT_M02_KC14103WT_M03_FemaleLiver_WT_KC14119WT_M13_KC14120WT_M14_KC14121WT_M15_featureCounts.txt"
  #argv$input_path <- "/projectnb2/wax-dk/david/G224_Shashi_Pipe/Scripts/"
  argv$input_path <- "/projectnb/wax-dk/max/G251_GalunLab/Scripts//09d_DE_1_RefSeqLncRNA76k/Output_DiffExp_1i_featureCounts_RefSeqLncRNA76k_ExonCollapsed/Down_Genes_EdgeR_RefSeqLncRNA76k_ExonCollapsed_mir122_KO_Male_2wk_G251_M04M05M06_WT_Male_2wk_G251_M01M02M03_featureCounts.txt"
}
print(argv)

library(clusterProfiler)
library(org.Mm.eg.db)
library(writexl)
library(patchwork)

library(stringr)
library(ggplot2)
library(gridExtra)
library(grid)

library(purrr)
library(readr)
library(tidyr)
library(dplyr)

library(gridExtra)
library(cowplot)
library(enrichplot)

## FUNCS
shorten_term <- function(short_term) {
  if(str_length(short_term) > 58) {
    as.character(str_glue("{str_sub(short_term,1,55)}..."))
  } else {
    #str_pad(short_term,58,side = "right",pad="-")
    short_term
  }
}

construct_title <- function(filename, top, stats_title) {
  comparison_num <- str_extract(filename,"(?<=09d_DE_)([[:alnum:]]+)(?=_Ref)") %>%
    str_pad(., width = 2, side = "left", pad="0")
  
  filename <- basename(filename)
  
  updown <- str_extract(filename,"([[:alnum:]|[_]]+)(?<=_ExonCollapsed|_FullGeneBody)")
  body <- str_remove(filename,paste0(updown,"_")) %>% str_remove("_featureCounts.txt")
  str_glue("{comparison_num}: {top} (clusterProfiler {packageVersion('clusterProfiler')})\n",
           "{updown}, {stats_title}\n",
           "{body}")
}

extract_updown_and_body <- function(filename){
  fname <- basename(filename)
  str_extract(fname,"([[:alnum:]|[_]]+)(?<=_ExonCollapsed|_FullGeneBody)")
}


## treeplot
create_treeplot <- function(lst){
  go_obj <- lst$obj
  go_title <- lst$title
  
  #go_obj@result$Description <- as.character(str_glue("{go_obj@result$Description} ({go_obj@result$ID})"))
  go_obj@result$Description_ext <- as.character(str_glue("{go_obj@result$ID}: {go_obj@result$Description}"))
  
  p1 <- pairwise_termsim(go_obj) %>%
    treeplot(label_format = 10,
             offset.params = list(bar_tree = 100,tiplab = 5, hextend = 20),
             fontsize = 5,
             cluster.params = list(label_words_n = 4))+
    ggtitle(go_title)+
    theme(plot.title = element_text(size = 9))
  
  p1$layers[[7]]$aes_params$size <- 5
  
  ##
  new_labels <- tibble(lbl = p1$data$label) %>%
    left_join(., go_obj@result %>%
                select(lbl = Description, new_lbl = Description_ext), by = join_by(lbl)) %>%
    pull(new_lbl)
  
  #p1$data$label <- str_sub(p1$data$label,1,90)
  p1$data$label <- str_sub(new_labels,1,80)
  
  p1
}

## MAIN

mmdb <- org.Mm.eg.db

fname <- argv$input_path
print(paste0(">>> Processing:", fname))

df <- read_tsv(file = fname, col_names = T, show_col_types = F)
if (nrow(df) < 1){
  print(str_glue("ERROR: {fname} DOES NOT HAVE ANY SIGNIFICANT GENES"))
  stop("exit 1")
  q(save = "no", 0, FALSE)
}
  
## get_entrezids
gene_symbols <- df %>%
  pull(1) %>%
  discard(~str_detect(.x,"nc_"))

entrez_ids <- AnnotationDbi::select(mmdb, keys=gene_symbols, columns="ENTREZID", keytype="SYMBOL") %>%
  drop_na() %>% pull(ENTREZID)

# GO enrichment analysis
goterm_analysis <- enrichGO(gene = entrez_ids,
                            OrgDb         = mmdb,
                            ont           = "All",
                            pvalueCutoff  = 0.05,
                            readable      = TRUE)

goterm_analysis_mf <- enrichGO(gene = entrez_ids,
                               OrgDb         = mmdb,
                               ont           = "MF",
                               pvalueCutoff  = 0.05,
                               readable      = TRUE)

goterm_analysis_bp <- enrichGO(gene = entrez_ids,
                               OrgDb         = mmdb,
                               ont           = "BP",
                               pvalueCutoff  = 0.05,
                               readable      = TRUE)

goterm_analysis_cc <- enrichGO(gene = entrez_ids,
                               OrgDb         = mmdb,
                               ont           = "CC",
                               pvalueCutoff  = 0.05,
                               readable      = TRUE)

res <- bind_rows(goterm_analysis_mf@result %>% mutate(ONTOLOGY = "MF"),
                 goterm_analysis_bp@result %>% mutate(ONTOLOGY = "BP"),
                 goterm_analysis_cc@result %>% mutate(ONTOLOGY = "CC"))


genes_stat <- df %>%  pull(1) %>%
  list(total = length(.),
       pcg = discard(.,~str_detect(.x,"nc_")) %>% length(.),
       lncrna = keep(.,~str_detect(.x,"nc_"))%>% length(.))

genes_stat_title <- str_glue("total # of genes: {genes_stat$total}, PCG: {genes_stat$pcg}, lncRNA: {genes_stat$lncrna}")

title1 <- construct_title(fname, "GO enrichment top20 by ontologies", genes_stat_title)
title2 <- construct_title(fname, "GO enrichment top60 any ontology", genes_stat_title)
mf_title <- construct_title(fname, "MF Ontology (only PCG)", genes_stat_title)
bp_title <- construct_title(fname, "BP Ontology (only PCG)", genes_stat_title)
cc_title <- construct_title(fname, "CC Ontology (only PCG)", genes_stat_title)
all_ontologies_title <- construct_title(fname, "All Ontologies together (only PCG)", genes_stat_title)

## top20 by ONTOLOGY
res1 <- res %>%
  mutate(short_term = str_glue("{ID}: {Description}"),
         mnlog10padj = -log10(p.adjust)) %>%
  filter(mnlog10padj>3) %>%
  mutate(short_term = map_chr(short_term, shorten_term)) %>%
  arrange(desc(mnlog10padj)) %>%
  slice_head(n = 20, by = ONTOLOGY)

## top60 does not matter which ONTOLOGY
res2 <- res %>%
  mutate(short_term = str_glue("{ID}: {Description}"),
         mnlog10padj = -log10(p.adjust)) %>%
  filter(mnlog10padj>3) %>%
  mutate(short_term = map_chr(short_term, shorten_term)) %>%
  arrange(desc(mnlog10padj)) %>%
  slice_head(n = 60)

if(nrow(res1) < 1) {
  print(str_glue("ERROR: {fname} DOES NOT HAVE ANY SIGNIFICANT TERMS AFTER FILTRATION (-log10(p.adjust) > 3)."))
  stop("exit 2")
  q("no", 0, FALSE)
}

# if(nrow(res1) < 1) {
#   print(str_glue("{fname} DOES NOT HAVE ANY SIGNIFICANT TERMS AFTER FILTRATION (-log10(p.adjust) > 3)."))
#   stop(str_glue("{fname} DOES NOT HAVE ANY SIGNIFICANT TERMS AFTER FILTRATION (-log10(p.adjust) > 3)."))
# }

## plot
#cls<-c("#ccece6", "#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#006d2c", "#00441b")

p1 <- ggplot(res1, aes(x = mnlog10padj, y = reorder(short_term, mnlog10padj), fill = mnlog10padj))+
  #coord_flip() +
  scale_fill_gradient(low = "blue", high = "red")+
  #scale_fill_viridis_c(option = "magma")+
  #scale_fill_gradient()+
  #scale_fill_gradientn(colours = heat.colors(10))+
  #scale_fill_gradient(low = "grey", high = "brown")+
  #scale_fill_gradientn(colours = cls)+
  geom_bar(width=0.7, stat = "identity")+
  geom_text(aes(label = Count), vjust = 0.5, hjust = -0.3, color = "black", fontface = "bold")+
  facet_grid(rows = vars(ONTOLOGY), scales = "free_y", space = "free_y", switch = "y")+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 15, color = "black"), #hjust = 1, vjust = 0.5, family = "robotomono"
    panel.grid.minor = element_blank(),
    axis.line = element_line(),
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_text(size = 12, face = "bold"))+
  ylab("")+
  xlab("-log10(padj)") +
  #xlim(2.5, max(res1$mnlog10padj))+
  #scale_x_continuous(limits = c(2.5, max(res1$mnlog10padj)))+
  coord_cartesian(xlim = c(2.5, max(res1$mnlog10padj)))+
  labs(fill="-log10(p.adj)")

p1 <- p1+plot_annotation(title = title1)

p11 <- ggplot(res2, aes(x = mnlog10padj, y = reorder(short_term, mnlog10padj), fill = mnlog10padj))+
  #scale_fill_gradientn(colours = cls)+
  scale_fill_gradient(low = "blue", high = "red")+
  geom_bar(width=0.7, stat = "identity")+
  geom_text(aes(label = Count), vjust = 0.5, hjust = -0.3, color = "black", fontface = "bold")+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 15, color = "black"), #hjust = 1, vjust = 0.5, family = "robotomono"
    panel.grid.minor = element_blank(),
    axis.line = element_line(),
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_text(size = 12, face = "bold"))+
  ylab("")+
  xlab("-log10(padj)")+
  coord_cartesian(xlim = c(2.5, max(res2$mnlog10padj)))+
  labs(fill="-log10(p.adj)")

p11 <- p11+plot_annotation(title = title2)


p2 <- ggplot(res1, aes(x = mnlog10padj, y = reorder(short_term, mnlog10padj)))+
  geom_point(aes(size = Count, colour = mnlog10padj))+
  geom_point(aes(size = Count), shape = 21, colour = "black")+
  #scale_color_gradientn(colours = cls)+
  scale_color_gradient(low = "blue", high = "red")+
  facet_grid(rows = vars(ONTOLOGY), scales = "free_y", space = "free_y", switch = "y")+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 15, color = "black"), #hjust = 1, vjust = 0.5, family = "robotomono"
    panel.grid.minor = element_blank(),
    axis.line = element_line(),
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_text(size = 12, face = "bold"))+
  ylab("")+
  xlab("-log10(p.adj)")+
  labs(color="-log10(p.adj)", size="Genes\nnumber")+
  scale_size(range = c(1, 10)) +
  xlim(2.5, max(res1$mnlog10padj))

p2 <- p2+plot_annotation(title = title1)

p22 <- ggplot(res2, aes(x = mnlog10padj, y = reorder(short_term, mnlog10padj)))+
  geom_point(aes(size = Count, colour = mnlog10padj))+
  geom_point(aes(size = Count), shape = 21, colour = "black")+
  #geom_point(aes(size = Count),shape = 1, colour = "black")+
  #scale_color_gradientn(colours = cls)+
  scale_color_gradient(low = "blue", high = "red")+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 15, color = "black"), #hjust = 1, vjust = 0.5, family = "robotomono"
    panel.grid.minor = element_blank(),
    axis.line = element_line(),
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_text(size = 12, face = "bold"))+
  ylab("")+
  xlab("-log10(p.adj)")+
  labs(color="-log10(p.adj)", size="Genes\nnumber")+
  scale_size(range = c(1, 10)) +
  xlim(2.5, max(res2$mnlog10padj))

p22 <- p22+plot_annotation(title = title2)


plots <- list(plot_grid(p1),plot_grid(p11),plot_grid(p2),plot_grid(p22))

# pairwise_termsim(goterm_analysis_cc)
# sum(goterm_analysis_cc@result$p.adjust < 0.05) %>% filter(p.adjust < 0.05)

# output_fname <- basename(fname) %>% tools::file_path_sans_ext()  %>%  paste0(".pdf") %>% paste0(compar_num,"_",.)

compar_num <- str_extract(fname,"(?<=09d_DE_)([[:alnum:]]+)(?=_Ref)") %>%
  str_pad(., width = 2, side = "left", pad="0")

outputnames_enrichment <- paste0(str_pad(1:4, width = 2, side = "left", pad = "0"),
                      "_",
                      compar_num,
                      "_",
                      c("top20_bar", "top60_bar", "top20_dot", "top60_dot"),
                      extract_updown_and_body(fname),
                      ".pdf")

plots_with_name <- plots %>% set_names(outputnames_enrichment)
outputnames_ontologies <- paste0(str_pad(5:8, width = 2, side = "left", pad = "0"),
                      "_",
                      compar_num,
                      "_",
                      c("MF", "BP", "CC", "ALLONT"),
                      extract_updown_and_body(fname),
                      ".pdf")

treeplots <- list(
  list(obj = goterm_analysis_mf, title = mf_title),
  list(obj = goterm_analysis_bp, title = bp_title),
  list(obj = goterm_analysis_cc, title = cc_title),
  list(obj = goterm_analysis, title = all_ontologies_title)
) %>% set_names(outputnames_ontologies) %>% 
  keep(\(l){sum(l$obj@result$p.adjust<0.05) > 1}) %>%
  keep(\(l){zz <- pairwise_termsim(l$obj); nrow(zz@termsim) > 5}) %>%
  map(create_treeplot)
  
all_plots <- flatten(list(plots_with_name,treeplots))

## output into separate files
walk2(flatten(list(plots,treeplots)), names(all_plots), \(p,name) {
  ggsave(plot = p, filename = name, width = 13, height = 15.19)
})


output_fname_xlsx <- basename(fname) %>% tools::file_path_sans_ext()  %>%  paste0(".xlsx") %>% paste0(compar_num,"_",.)
write_xlsx(res, path = output_fname_xlsx, col_names = T)

# # list("aa" = c(1,2,3)) %>% enframe() %>% unnest_longer(col = "value")
