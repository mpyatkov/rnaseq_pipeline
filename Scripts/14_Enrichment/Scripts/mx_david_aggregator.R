suppressMessages(library(tidyverse))

remotes::install_cran("argparser", upgrade = "never", quiet = T)
library(argparser)

p <- arg_parser('DAVID GO enrichment combining')
p <- add_argument(p,'--input_path', default="./", help="DAVID files, usually txt files")
p <- add_argument(p,'--sample_labels', default="../00_Setup_Pipeline/Sample_Labels.txt", help="Sample_Labels.txt file")
p <- add_argument(p,'--comparisons', default="../00_Setup_Pipeline/Comparisons.txt", help="Comparisons.txt file")
argv <- parse_args(p)
print(argv)


## ONLY FOR DEBUGs
# argv$input_path <- "../Job_Summary/DAVID_results/ExonCollapsed/"
# argv$sample_labels <- "../../00_Setup_Pipeline/Sample_Labels.txt"
# argv$comparisons <- "../../00_Setup_Pipeline/Comparisons.txt"

## named vector with comparison names
simplified_comparisons <- function(argv_sample_labels, argv_comparisons){
  sample_labels <- read_delim(argv_sample_labels, delim = ";", col_names=T, show_col_types = F) %>% 
    select(Group,Condition_Name) %>% 
    distinct()
  
  ## named vector, v[['1']] -- extract first comparison name
  comparisons <- read_delim(argv_comparisons, delim = ";", col_names=T, show_col_types = F) %>% 
    pivot_longer(!Comparison_Number, names_to = "condition", values_to = "Group") %>% 
    left_join(., sample_labels, by = join_by(Group)) %>% 
    select(-Group) %>% 
    pivot_wider(names_from = condition, values_from = Condition_Name) %>% 
    mutate(Comparision_Name = str_glue("{Condition_2}_vs_{Condition_1}")) %>% 
    select(-Condition_1, -Condition_2) %>% 
    mutate(Comparison_Number = str_pad(Comparison_Number, width = 2, pad = "0")) %>% 
    deframe()
  
  comparisons
}

COMPAR_NAMES <- simplified_comparisons(argv$sample_labels, argv$comparisons)

extract_table_from_david <- function(line){
  header <- line %>% str_split("\\n", simplify = T) %>% .[[1]] 
  annotation_cluster <- str_extract(header, "([[:digit:]]+)(?=\t)") %>% as.numeric(.)
  annotation_es <- header %>% str_extract("([[:digit:]|\\.]+)$") %>% as.numeric(.) %>% round(.,digits = 5)

  body <- read_tsv(line, skip = 1, show_col_types = F) %>% 
    mutate(Annotation_Cluster = annotation_cluster,
           Enrichment_Score = annotation_es) %>% 
    relocate(Enrichment_Score,.before = Category) %>% 
    relocate(Annotation_Cluster,.before = Enrichment_Score)
  
  body
}

## aggregate unique genes
## c("a,b,c", "b,c,d") -> "a,b,c,d"
agg_uniq_genes_in_group <- function(ll){
  str_c(ll, collapse = ",") %>% 
    str_split(",") %>% 
    unlist(.) %>% 
    str_trim() %>% 
    unique() %>% 
    str_c(., collapse = ", ")
}

process_file <- function(f){
  
  fname <- basename(f)

  comparison_num_str <- str_extract(fname,"([[:alnum:]]+)(?=_)")
  #comparison_num_num <- str_extract(fname,"([[:alnum:]]+)(?=_)") %>% as.numeric()
  body <- ifelse(str_detect(fname, "FullGeneBody"), "FullGeneBody", "ExonCollapsed")
  updown <- ifelse(str_detect(fname, "Up_"), "Up","Down")
  
  comparison_name <- COMPAR_NAMES[[comparison_num_str]]
  
  tmp <- read_file(f)
  cluster_chunks <- str_split(tmp, "\\n\\n", simplify = T) %>% 
    discard(~str_length(.x) == 0)
  
  res <- map(cluster_chunks, \(x) extract_table_from_david(x)) %>% list_rbind()
  
  res %>% 
    mutate(Comparison_Num = comparison_num_str,
           Regulation = updown,
           Coverage = body,
           Comparison_Name = comparison_name) %>% 
    mutate(All_unique_genes_in_annotation_cluster = agg_uniq_genes_in_group(Genes), .by = c(Regulation, Comparison_Num, Comparison_Name, Annotation_Cluster)) %>% 
    mutate(Annotation_Position = row_number(), .by = c(Comparison_Num,Regulation,Coverage,Comparison_Name, Annotation_Cluster)) %>% 
    relocate(Annotation_Position,.after = Annotation_Cluster) %>% 
    relocate(Comparison_Name,.before =Annotation_Cluster) %>% 
    relocate(Coverage,.before =Comparison_Name) %>% 
    relocate(Regulation,.before =Coverage) %>% 
    relocate(Comparison_Num,.before = Regulation)
}

excollapsed <- list.files(path = argv$input_path, pattern = "ExonCollapsed", full.names = T) %>% 
  map(\(f) process_file(f)) %>% 
  list_rbind()

fullgenebody <- list.files(path = argv$input_path, pattern = "FullGeneBody", full.names = T) %>% 
  map(\(f) process_file(f)) %>% 
  list_rbind()

writexl::write_xlsx(list(ExonCollapsed = excollapsed,
                         FullGeneBody = fullgenebody,
                         Sample_labels = read_delim(argv$sample_labels, delim = ";", col_names=T, show_col_types = F),
                         Comparisons = read_delim(argv$comparisons, delim = ";", col_names=T, show_col_types = F)),
                    path = "Max_AggregatedDavidResults.xlsx")



