library(dplyr)
library(stringr)
library(readr)

args <- commandArgs(T)

# # initial index from SEGEX
SEGEX_EXPT_IX <- as.double(args[1])

# sort index to order which feature will be first numerated
# 0- skip renaming,  1 - FGB only, 2 - EC only, 3 - both, FGB - first, EC - second
SEGEX_FEATURE_SORT_IX <- as.double(args[2])

# for test only
# SEGEX_EXPT_IX <- as.double(839)
# SEGEX_FEATURE_SORT_IX <- 1

if (SEGEX_FEATURE_SORT_IX == 0 || SEGEX_FEATURE_SORT_IX > 3){
    warning("Execution is not required. SEGEX_FEATURE_SORT_IX equal to 0 or > 3. SEGEX files will not be renamed.")
    quit(save = "no", status = 0, runLast = F)
}

# extract all information about files and add indexes
all_files <- tibble(full_path = list.files(pattern = "^[[:digit:]+][[:print:]]+\\.txt$", full.names = T, recursive = T)) %>% 
    mutate(file_name = basename(full_path),
           dir_name = dirname(full_path),
           feature = str_extract(file_name, "FullGeneBody|ExonCollapsed"),
           sample_1_name = str_extract(file_name,"(?<=_vs_)([[:print:]]+)(?=_DiffExp)"),
           #sample_2_name = str_extract(file_name,"(?<=^[[:digit:]+]_)([[:print:]]+)(?=_vs_)"),
           sample_2_name = str_extract(file_name,"(?<=[^([:digit:]+)])([:print:]+)(?=_vs)"),
           ix = as.double(str_extract(file_name,"(^[:digit:]+)(?=_)")),
           ix = case_when( SEGEX_FEATURE_SORT_IX == 3 & feature == "ExonCollapsed" ~ ix+100,
                           SEGEX_FEATURE_SORT_IX == 2 & feature == "FullGeneBody" ~ 0,
                           SEGEX_FEATURE_SORT_IX == 1 & feature == "ExonCollapsed" ~ 0,
                           TRUE ~ ix)) %>% 
    filter(ix != 0) %>% 
    arrange(ix) %>% 
    mutate(expt_number = SEGEX_EXPT_IX+row_number()-1) %>% 
    rowwise() %>% 
    mutate(new_file_name = paste0(expt_number, "_", file_name),
           new_full_path = paste0(dir_name,"/", new_file_name))

# output table as requested
output_table <- all_files %>% 
    mutate(platform_number = 4,
           platform = "4_Mouse_RNA-Seq_76K-RefSeqLncRNA-genes",
           exp_desc = str_replace(new_file_name,"_forSEGEXUpload",""),
           exp_desc = str_replace(exp_desc,"\\.txt",""),
           additional_info = paste0("76K_",feature,"_TPM_EdgeR_FeatureCounts_pval1:padj_pval2:pval"),
           study_in_segex = "TBD",
           upload_file = str_replace(new_file_name,"\\.txt","")) %>% 
    select(expt_number, 
           platform_number,
           sample_1_name,
           sample_2_name,
           exp_desc,
           upload_file,
           additional_info,
           platform,
           study_in_segex) %>% 
    write_tsv("Upload_SEGEX_table.tsv", col_names = T)

all_files_renamed <- file.rename(from=all_files$full_path, to=all_files$new_full_path)
sprintf("All required files renamed: %s", all(all_files_renamed))



