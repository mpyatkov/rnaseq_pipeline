suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

args <- commandArgs(T)
INPUT_PREFIX <- args[1]
OUTPUT_PREFIX <- args[2]

prepare_for_separation <- function(l){
  str_replace(l,"^\\(","") %>% 
    str_replace(., "\\)$","") %>% 
    str_replace_all(.,"\\(","_") %>%
    str_replace_all(.,"\\)","_") %>% 
    str_replace(.,"__","_") 
}

read_tsv(paste0(INPUT_PREFIX,".tsv"), col_names = T) %>% 
  #select(transcript_id) %>% 
  mutate(transcript_id = prepare_for_separation(transcript_id)) %>% 
  separate(transcript_id, into = c("strand_tmp", "gene_id", "tr_id","expr","abs_frac"), sep = "_", convert = T) %>% 
  select(gene_id, transcript_id = tr_id, expr, abs_frac, everything(), -strand_tmp) %>% 
  write_tsv(paste0(OUTPUT_PREFIX,".tsv"))

# 
# orig <- read_tsv(paste0(INPUT_PREFIX,".tsv"), col_names = T)
# 
# new <- orig %>% 
#   select(transcript_id) %>% 
#   mutate(transcript_id = prepare_for_separation(transcript_id)) %>% 
#   separate(transcript_id, into = c("strand", "gene_id", "tr_id","expr","abs_frac"), sep = "_", convert = T) %>% 
#   write_tsv()

# mutate(strand = str_extract(transcript_id,"\\+|\\-"),
#        gene_id = str_extract(transcript_id,"(?<=\\))([[:alnum:]]+)(?=\\_)"),
#        transcript_id = str_extract(transcript_id,"(?<=\\_)([[:alnum:]]+)(?=\\()"),
#        expr = as.double(str_extract(transcript_id,"(?<=\\()([[:digit:]|\\.]+)(?=\\))")),
#        abs_frac = as.double(str_extract(transcript_id,"(?<=\\()([[:digit:]|\\.]+)(?=\\)$)"))) %>% 
