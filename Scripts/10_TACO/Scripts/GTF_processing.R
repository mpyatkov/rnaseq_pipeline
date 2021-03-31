suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))

args <- commandArgs(T)
INPUT_PREFIX <- args[1]
OUTPUT_PREFIX <- args[2]

print("Starting GTF postprocessing...")

gtf <- rtracklayer::import(paste0(INPUT_PREFIX,".gtf"))

transcripts <- as_tibble(gtf) %>% 
  filter(type == "transcript") %>% 
  rowwise() %>% 
  mutate(new_name = paste0("(",strand,")",
                          paste0(gene_id,"_",transcript_id),
                          "(",round(as.double(expr),1),")",
                          "(",round(as.double(abs_frac),3),")"),
         ix=paste0(gene_id,"|",transcript_id)) %>% 
  select(ix, new_name)

## GTF file
gtf_processed <- as_tibble(gtf) %>% 
  rowwise() %>% 
  mutate(ix=paste0(gene_id,"|",transcript_id)) %>% 
  left_join(., transcripts, by =c("ix")) %>% 
  mutate(transcript_id = new_name) %>% 
  select(-new_name, -ix)

rtracklayer::export(gtf_processed, paste0(OUTPUT_PREFIX,".gtf"))

print("Starting BED postprocessing...")

## BED file 
# read_tsv did not recognize column 11 correctly, this is the spec for the bed file
cols_spec <- cols(
  X1 = "c",
  X2 = "d",
  X3 = "d",
  X4 = "c",
  X5 = "d",
  X6 = "c",
  X7 = "d",
  X8 = "d",
  X9 = "d",
  X10 = "d",
  X11 = "c",
  X12 = "c"
)

bed_processed <- read_tsv(paste0(INPUT_PREFIX, ".bed"), col_names = F, col_types = cols_spec) %>% 
  mutate(ix = str_extract(X4,"(?<=^)([[:print:]]+)(?=\\()")) %>% 
  left_join(., transcripts, by=c("ix")) %>% 
  mutate(X4 = new_name) %>% 
  select(-new_name, -ix)

write_tsv(bed_processed, paste0(OUTPUT_PREFIX,".bed"), col_names = F)
print("Ending GTF postprocessing...")
