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
      GeneSetID <- str_split(basename(file),delim,simplify = T)[pos]
      
      ## Select cluster and enrichment information ##
      dat_annotation <- grep('Annotation',dat,value = T)
      cluster <- data.frame(str_split(dat_annotation,'\t| ',simplify = T))
      cluster_label <- data.frame(Annotation_Cluster=as.numeric(cluster$X3),Cluster_Score=as.numeric(cluster$X6))
      #cluster$X2 <- as.numeric(str_split(cluster$X2," ",simplify = T)[,3])
      
      ## Processing category / term data
      c_sum <- clusters_process(dat,file)
      
      y <- cbind(GeneSetID,cluster_label,c_sum)
      
      y <- y %>%
        mutate(All_Unique_Terms_Name_Copy=All_Unique_Terms_Name) %>%
        relocate(All_Unique_Terms_Name_Copy,.after = All_Unique_Terms_ID) %>%
        relocate(Top_Genes,.after = Cluster_Score) %>%
        relocate(Top_Count,.after = Cluster_Score) %>%
        relocate(Top_FDR,.after = Cluster_Score) %>%
        relocate(All_Unique_Terms_Name,.after = Cluster_Score) %>%
        relocate(Top_Term,.after = Cluster_Score)
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
argv <- parse_args(p)
print(argv)

### Parse multiple files
# Provide vector of file paths for each DAVID output file
input_folder <- argv$input_path
files_exoncollapsed <- list.files(input_folder,pattern = "ExonCollapsed",full.names = T)
files_fullgenebody <- list.files(input_folder,pattern = "FullGeneBody", full.names = T)
df_exoncollapsed <- multi_file_process(files_exoncollapsed)
df_fullgenebody <- multi_file_process(files_fullgenebody)
write_xlsx(list(ExonCollapsed = df_exoncollapsed, FullGeneBody = df_fullgenebody), path = 'Alan_AggregatedDavidResults.xlsx')
