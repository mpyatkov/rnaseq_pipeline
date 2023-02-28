## input params
## path - path to directory with all steps
## output_xls name - name of output file

library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(openxlsx)

args <- commandArgs(T)
path <- args[1]
output_xls <- args[2]

## just for test
# path <- "/projectnb/wax-dk/max/20221110/STAT5_KO_76K/Scripts/"
# output_xls <- "/projectnb/wax-dk/max/20221110/STAT5_KO_76K_summary.xlsx"

samples_path <- paste0(path, "00_Setup_Pipeline/Sample_Labels.txt")

## return group, rpkm, tpm, vector of samples
samples_new <- read_delim(samples_path, delim = ";", col_names = T, comment = "#", trim_ws = TRUE) %>% 
  select(group = 1, condition = 2, sample_id = 3) %>% 
  rowwise() %>% 
  mutate(name = paste0(condition,"_",sample_id),
         rpkm = paste0("rpkm_mean_", condition),
         tpm = paste0("tpm_mean_", condition),
         vec = list(list(name, rpkm, tpm))) %>% 
  group_by(group) %>%
  summarise(rpkm = first(rpkm),tpm = first(tpm), samples = list(name)) %>%
  ungroup() %>%
  select(group, rpkm, tpm, samples)

## input: samples_new
##        DiffExp_v2* data.frame (our usual table from step 9)
## process:
##        for each row from samples_new tried to detect samples inside the file
##        if samples detected then extract column with counts and mean_rpkm and
##        mean_tpm, combine all columns together and return as output
process_file <- function(group_vec_names, file_df){
  
  pmap_dfc(group_vec_names, function(group, rpkm, tpm, samples){
    
    tmp <- file_df %>% select(id, one_of(samples))
    
    if (ncol(tmp) > 1) { # if we have another columns (not only id)
      
      ## attach rpkm and tpm columns to counts
      tmp <- bind_cols(tmp, file_df %>% select(one_of(c(rpkm,tpm))))
      
      ## add group name for each column
      tmp <- tmp %>% rename_at(vars(-id), ~paste0(group,"_",.))
    }
    
    tmp
  }) %>% select(id = `id...1`, !contains("...")) # removing duplicated ids
  
}

## scan files in path directory, extract only DiffExp_v2 + ExonCollapsed
all_files <- list.files(path = path, recursive = T, pattern = "^DiffExp_v2") %>% 
  tibble(file = .) %>% 
  filter(str_detect(file, "ExonCollapsed")) %>% 
  mutate(file = paste0(path,file)) %>% pull(file)


## working "for loop" version (deprecated)
# process_file <- function(group_vec_names, file_df){
# 
# acc <- file_df %>% select(id, one_of("fake_col"))
# groups_acc <- c()
# 
# for(i in 1:nrow(group_vec_names)) {
#   
#   group <- group_vec_names[i, "group"]
#   rpkm <- group_vec_names[i, "rpkm"] %>% pull(rpkm)
#   tpm <- group_vec_names[i, "tpm"] %>% pull(tpm)
#   vsamples <- group_vec_names[i, "samples"] %>% unlist()
#   
#   print(group_vec_names)
#   
#   if (group %in% groups_acc) {next}
#   
#   
#   tmp <- file_df %>% select(id, one_of(vsamples))
#   
#   if (ncol(tmp) > 1) { # if not only id column
#     tmp <- bind_cols(tmp, file_df %>% select(one_of( c(rpkm, tpm) )))
#     tmp <- tmp %>% rename_at(vars(-id), ~paste0(group,"_",.))
#   }
#   
#   acc <- bind_cols(acc, tmp)
#   groups_acc <- c(groups_acc, group)
# }
# 
# acc %>% select(id = `id...1`, !contains("...")) 
# }

## scan all DiffExp files one by one and extract information about counts 
## and rpkm/tpm_mean 
binded_samples <- map(all_files, function(f){
  f_df <- read_tsv(f, col_names = T) 
  process_file(samples_new, f_df)
}) %>% purrr::reduce(., dplyr::left_join) 

# calculates colsums, excludes id-column
calculate_colsums <- function(df) {
  zz_colsums <- df %>% 
    select(-id) %>% 
    as.matrix() %>% 
    colSums() %>% round(., 2) %>% as.vector()
  
  zz_colsums <- tibble(header = colnames(df)[-1], sums = zz_colsums) %>% 
    pivot_wider(names_from = header, values_from = sums) %>% 
    mutate(id = "total") %>% 
    relocate(id, .before = everything())
  
  zz_colsums
}

## calculate colsums separatelly, we are going add them as separate top line 
## in resulted xlsx file
binded_samples_header <- calculate_colsums(binded_samples)

## create xlsx workbook
wb <- createWorkbook()
addWorksheet(wb, sheetName = "summary")
## colNames false for header because we do not need this information inside xlsx
writeData(wb, sheet = "summary", binded_samples_header, startRow = 1, startCol = 1, colNames = FALSE)
writeData(wb, sheet = "summary", binded_samples, startRow = 3, startCol = 1)

## apply bold style for rpkm/tpm columns
rpkm_tpm_cols <- which(str_detect(colnames(binded_samples), "rpkm|tpm") == TRUE)
boldStyle <- createStyle(textDecoration = "Bold") ## create bold style
addStyle(wb, "summary", boldStyle, 3, rpkm_tpm_cols)

## apply number style for column 'total'
comma_style_format <- createStyle(numFmt = "#,##0") # create thousands format
addStyle(wb, "summary", comma_style_format, rows = 1, cols = 2:ncol(binded_samples))

saveWorkbook(wb, output_xls, overwrite = T)

