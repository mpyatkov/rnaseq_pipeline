# up and down genes summarizing for each comparison group of 09abc directories
# ex: find ../09a_DE_* -iname "*de_ge*counts*" | xargs -n1 -I{} cat {} | grep -v "Output_File" >> ./input_file.csv

library(dplyr)
library(stringr)
library(readr)
library(tidyr)

args <- commandArgs(T)

input_file <- args[1]
comparisons <- args[2]
samples_file <- args[3]
output_name <- args[4]
dataset_label <- args[5]

# read table with DE_gene_counts
up.down.genes <- read_delim(input_file, delim = " ", col_names = F, trim_ws = T) %>% 
  rowwise() %>% 
  mutate(comparison=as.numeric(str_extract(X1, "(\\d)+")),
         feature=str_remove(X1, "[:alpha:]+_[:alpha:]+_[:alnum:]+_[:alpha:]+_"),
         package=ifelse(str_detect(str_to_lower(X2), "deseq"),"DESeq","edgeR"),
         # Condition_1_GM_Numbers=str_extract_all(X2,"G[:alnum:]+_[:alnum:]+", simplify = T)[2],
         # Condition_2_GM_Numbers=str_extract_all(X2,"G[:alnum:]+_[:alnum:]+", simplify = T)[1],
         Dataset_label=dataset_label) %>% 
  select(Dataset_label, comparison, feature, package, Up_genes = X4, Down_genes = X5)

# comparison df
cmp <- read_delim(comparisons, delim = ";", col_names = T, trim_ws = T) %>% 
  select(comparison = Comparison_Number, Condition_1, Condition_2)

# combine samples by dataset_id and samples_id
get_samples_line <- function(smpl){
  tibble(s=smpl) %>% 
    separate(s, c("dataset", "sample"),sep = "_") %>% 
    arrange(sample) %>% 
    group_by(dataset) %>% 
    summarise(samples = str_c(sample, collapse = "")) %>% 
    ungroup() %>%
    unite(res,dataset:samples) %>% pull() %>% paste0(., collapse = "_")
}

# samples df
samples <- read_delim(samples_file, delim=";", col_names = T, trim_ws = T) %>% 
  select(Group, Condition_Name, Sample_ID) %>% 
  group_by(Group) %>% 
  mutate(Samples=get_samples_line(Sample_ID)) %>% 
  ungroup() %>% 
  select(-Sample_ID) %>% 
  distinct()

# combine comparison and samples
cmp_samples <- left_join(cmp, samples, by=c("Condition_1"="Group")) %>% 
  left_join(., samples, by=c("Condition_2"="Group")) %>% 
  select(comparison, Condition_1, Condition_2, 
         Condition_Name_1 = Condition_Name.x, 
         Condition_Name_2 = Condition_Name.y,
         Condition_1_GM_Numbers = Samples.x,
         Condition_2_GM_Numbers = Samples.y)

# final table
left_join(up.down.genes,cmp_samples) %>% 
  select(Dataset_label, comparison,feature,package,Condition_1,Condition_2,
         Condition_Name_1,Condition_Name_2,
         Condition_1_GM_Numbers,Condition_2_GM_Numbers,
         Up_genes,Down_genes) %>% 
  gather(variable, value, -(Dataset_label:Condition_2_GM_Numbers)) %>% 
  unite(combine,feature,variable, sep = ".") %>% 
  spread(combine, value) %>% 
  arrange(comparison) %>% 
#  write_excel_csv2(., output_name) # separator ";"
  write_csv(., output_name) # separator ","


