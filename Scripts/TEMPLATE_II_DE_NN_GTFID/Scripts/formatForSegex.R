# convert DiffExp_v2... files to files appropriate for uploading to SEGEX database
# produce files for tpm/rpkm and edger/deseq
# input:
# inputFile <- args[1]           # input name
# outputFile <- args[2]          # output name
# epsilon <- as.numeric(args[3]) # epsilon for deseq_ratio calculation
# col_suffix <- args[4]          # suffix which will append to each column
# normaliz <- args[5]            # tpm/rpkm

library(dplyr)
library(readr)
library(tidyr)
library(stringr)

# reducing name for SEGEX
replace_name <- function(name) {
  name %>% 
    str_replace(.,"Gt\\(ROSA\\)26Sor","Gt_ROSA_26Sor") %>% 
    str_replace(.,"^lnc_","nc_") %>%
    str_replace(.,"^ncRNA_","nc_") %>% 
    str_replace(.,"intra-as","intra") %>% 
    str_replace(.,"_chr","_c") %>% 
    str_replace(., "@", "_") %>% 
    str_replace(., "\\(\\+\\)", "p") %>% 
    str_replace(., "\\(-\\)", "m") %>% 
    str_replace(., "random", "") %>% 
    str_replace(., "_JH[0-9]+","_JH") %>% 
    str_replace(.,"LOC100861978_chr4_JH_m", "LOC100861978_chr4m") %>% 
    ifelse(str_detect(., "^nc_"), str_to_lower(.), .) 
}

## MAIN
args <- commandArgs(trailingOnly = TRUE)

inputFile <- args[1]
outputFile <- args[2]
epsilon <- as.numeric(args[3])
col_suffix <- args[4]
normaliz <- args[5]

# by default if not tpm than rpkm normalization
normalization <- ifelse(normaliz == "tpm", normaliz, "rpkm")

data <- read_delim(inputFile, delim="\t", col_names = T)

# separate calculation of the deseq_ratio from the baseMean calculations
deseq_ratio <- data %>% 
  select(starts_with("baseMean_")) %>% 
  mutate(res=(.[[2]]+epsilon)/(.[[1]]+epsilon)) %>% 
  select(deseq_ratio=res) %>% 
  replace(is.na(.), 1)

# print(glimpse(data))
# main table  
res <- data %>% 
  mutate(EdgeR_foldChange = ifelse(is.na(EdgeR_foldChange),1,EdgeR_foldChange),
         DESeq2_foldChange = ifelse(is.na(DESeq2_foldChange),1,DESeq2_foldChange),
         edger_foldChange = ifelse(EdgeR_foldChange < 1, -1*(1/EdgeR_foldChange), EdgeR_foldChange),
         deseq_foldChange = ifelse(DESeq2_foldChange < 1, -1*(1/DESeq2_foldChange), DESeq2_foldChange),
         deseq_padj = ifelse(is.na(DESeq2_padj_FDR), 1, DESeq2_padj_FDR),
         edger_padj= ifelse(is.na(EdgeR_padj_FDR), 1, EdgeR_padj_FDR),
         deseq_pval=ifelse(is.na(DESeq2_pvalue), 1, DESeq2_pvalue),
         edger_pval=ifelse(is.na(EdgeR_pvalue), 1, EdgeR_pvalue),
         id = replace_name(id)) %>% 
  select(id,
         edger_ratio = EdgeR_foldChange,
         deseq_foldChange,
         edger_foldChange,
         contains(paste0(normalization,"_mean")),
         deseq_padj,
         edger_padj,
         deseq_pval,
         edger_pval) %>% 
  cbind(deseq_ratio, .) %>% 
  replace(is.na(.), 0) %>% 
  rename_all(~paste0(., ".", col_suffix))

# because is.infinite() does not work with data frames, we replace all Inf values inplace
res[res == -Inf] <- 1
res[res == Inf] <- 1

# for deseq (rpkm | tpm normalization)
output_filename <- paste0(outputFile, ifelse(normalization == "tpm", "_TPM_DESeq.txt", "_DESeq.txt"))
deseq_output <- res %>% 
  select(one_of(paste0("id.",col_suffix)), matches(paste0("deseq|",normalization))) %>% 
  rename_all(~str_replace(., "deseq_", "")) %>% 
  format(., digits = 8) %>%
  write_tsv(output_filename)

# for edger (rpkm | tpm normalization)
output_filename <- paste0(outputFile, ifelse(normalization == "tpm", "_TPM_EdgeR.txt", "_EdgeR.txt"))
edger_output <- res %>% 
  select(one_of(paste0("id.",col_suffix)), matches(paste0("edger|",normalization))) %>% 
  rename_all(~str_replace(., "edger_", "")) %>% 
  format(., digits = 8) %>%
  write_tsv(output_filename)

