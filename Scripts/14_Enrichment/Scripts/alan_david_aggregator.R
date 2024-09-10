### DAVID output parsing script ###
# Created by Alan Downey-Wall
# July 22nd 2022

## Purpose
# Purpose of script is to parse DAVID output files into a more user friendly format,
# while also summarizing clusters for easy downstream analysis.

##########################################################################################################
########################################## Functions #####################################################
##########################################################################################################

# NOTE: Make sure to initialize all functions prior to running the main functions above.

##########################
#### Small functions #####
##########################

## Function for spliting term strings 
# pos = 1 - The term ID
# pos = 2 - The term name
term_split <- function(x,pos=1){
  if(pos==1){
    y <- ifelse(isTRUE(grepl('GO:',x)),unlist(strsplit(x,split = '~'))[1],unlist(strsplit(x,split = '~|:'))[1])
  }
  if(pos==2){
    y <- ifelse(isTRUE(grepl('GO:',x)),unlist(strsplit(x,split = '~'))[2],unlist(strsplit(x,split = '~|:'))[2])
  }
  return(y)
}

## Removes white spaces from strings within a vector
wsr <- function(x){
  gsub("\\s", "",x)
}

##########################
#### Main functions #####
##########################
multi_file_process <- function(x,delim='[_.]',pos=2){
  for(j in 1:length(x)){
    sfile <- file_process(x[j],delim = delim,pos = pos)
    if(!is_null(sfile)){
      if(j == 1){
        y <- sfile
      }else{y <- rbind(y,sfile)} 
    }
  }
  return(y)
}

#file <- files[2]
#temp <- file_process(file)

file_process <- function(file,delim='[_.]',pos=2){
    delim='[_.]'
    pos=1
    print(file)
    dat <- read_lines(file, progress = F)
    if(length(dat) > 1){
      fname <- basename(file)
      
      comparison_num_str <- str_extract(fname,"([[:alnum:]]+)(?=_)")
      #comparison_num_num <- str_extract(fname,"([[:alnum:]]+)(?=_)") %>% as.numeric()
      
      body <- ifelse(str_detect(fname, "FullGeneBody"), "FullGeneBody", "ExonCollapsed")
      updown <- ifelse(str_detect(fname, "Up_"), "Up","Down")
      
      comparison_name <- COMPAR_NAMES[[comparison_num_str]]
      
      ## Select cluster and enrichment information ##
      dat_annotation <- grep('Annotation',dat,value = T)
      cluster <- data.frame(str_split(dat_annotation,'\t| ',simplify = T))
      cluster_label <- data.frame(Annotation_Cluster=as.numeric(cluster$X3),Cluster_Score=as.numeric(cluster$X6))

      ## Processing category / term data
      c_sum <- clusters_process(dat,file)
      
      y <- cbind(cluster_label,c_sum) %>%
        mutate(All_Unique_Terms_Name_Copy=All_Unique_Terms_Name) %>%
        relocate(All_Unique_Terms_Name_Copy,.after = All_Unique_Terms_ID) %>%
        relocate(Top_Genes,.after = Cluster_Score) %>%
        relocate(Top_Count,.after = Cluster_Score) %>%
        relocate(Top_FDR,.after = Cluster_Score) %>%
        relocate(All_Unique_Terms_Name,.after = Cluster_Score) %>%
        relocate(Top_Term,.after = Cluster_Score)
      
      y <- y %>% mutate(Comparison_Num = comparison_num_str,
                  Regulation = updown,
                  Coverage = body,
                  Comparison_Name = comparison_name) %>% 
        select(Comparison_Num, Regulation, Coverage, Comparison_Name, everything())
      
    }else{y <- NULL}

    return(y)
}

clusters_process <- function(dat,file){

  block_start <- grep('Annotation',dat)
  block_end <- which(dat == "")-2
  
  header <- c("Category","Term","Count","Percent","PValue","Genes","List_Total",
              "Pop_Hits","Pop_Total","Fold_Enrichment","Bonferroni","Benjamini","FDR")
  top_header <- paste0("Top_",header)
  
  for(i in 1:length(block_start)){
    out <- read_tsv(file = file,skip = block_start[i],n_max = block_end[i] - block_start[i], show_col_types = F, progress = F)
    # Top process for each enrichment cluster   
    top_process <- out[1,] 
    # Relabel
    colnames(top_process) <- top_header
    top_process <- top_process %>%
      mutate(Top_Term_ID = term_split(Top_Term,pos=1),
             Top_Term_Name = str_to_sentence(term_split(Top_Term,pos=2))) %>%
      relocate(Top_Term_ID,.after=Top_Term) %>%
      relocate(Top_Term_Name,.after=Top_Term_ID)
    
    ## Add summary of terms within cluster ##
    top_process$All_Terms <- paste(out$Term,collapse=',')
    top_process$All_Terms_ID <- paste(sapply(out$Term,term_split,pos=1),collapse=',')
    top_process$All_Terms_Name <- paste(str_to_sentence(sapply(out$Term,term_split,pos=2)),collapse=',')
    top_process$All_Terms_Count <- length(out$Term)
    top_process$All_Unique_Terms_ID <- paste(unique(toupper(sapply(out$Term,term_split,pos=1))),collapse=',')
    top_process$All_Unique_Terms_Name <- paste(unique(str_to_sentence(sapply(out$Term,term_split,pos=2))),collapse=',')
    top_process$All_Unique_Terms_Name_Count <- length(unique(str_to_sentence(sapply(out$Term,term_split,pos=2))))
    
    ## Add summary of genes within cluster ##
    # wsr = white space remove (function above)
    
    # top_process$All_Genes <- paste(wsr(unlist(strsplit(out$Genes,','))),collapse = ',')
    # top_process$All_Genes_Count <- length(wsr(unlist(strsplit(out$Genes,','))))
    unique_genes <- unique(wsr(unlist(strsplit(out$Genes,','))))
    
    top_process$All_Unique_Genes_Count <- length(unique_genes)
    top_process$All_Unique_Genes <- paste(unique_genes,collapse = ',')
    
    if(i == 1){
     y <- top_process 
    }else{
      y <- rbind(y,top_process)
    }
  }
  return(y)
}

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
    deframe()
  
  comparisons
}

########################
##### User Input #######
########################

library(tidyverse)
remotes::install_cran("writexl", upgrade = "never")
library(writexl)
remotes::install_cran("argparser", upgrade = "never", quiet = T)
library(argparser)

p <- arg_parser('DAVID GO enrichment combining')
p <- add_argument(p,'--input_path', default="./", help="DAVID files, usually txt files")
p <- add_argument(p,'--sample_labels', default="../00_Setup_Pipeline/Sample_Labels.txt", help="Sample_Labels.txt file")
p <- add_argument(p,'--comparisons', default="../00_Setup_Pipeline/Comparisons.txt", help="Comparisons.txt file")
argv <- parse_args(p)
print(argv)

# argv$input_path <- "../Job_Summary/DAVID_results/ExonCollapsed/"
# argv$sample_labels <- "../../00_Setup_Pipeline/Sample_Labels.txt"
# argv$comparisons <- "../../00_Setup_Pipeline/Comparisons.txt"

COMPAR_NAMES <- simplified_comparisons(argv$sample_labels, argv$comparisons)

### Parse multiple files
# Provide vector of file paths for each DAVID output file
input_folder <- argv$input_path
files_exoncollapsed <- list.files(input_folder,pattern = "ExonCollapsed",full.names = T)
files_fullgenebody <- list.files(input_folder,pattern = "FullGeneBody", full.names = T)
df_exoncollapsed <- multi_file_process(files_exoncollapsed)
df_fullgenebody <- multi_file_process(files_fullgenebody)
write_xlsx(list(ExonCollapsed = df_exoncollapsed, FullGeneBody = df_fullgenebody), path = 'Alan_AggregatedDavidResults.xlsx')
