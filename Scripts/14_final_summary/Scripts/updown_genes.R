# up and down genes summarizing for each comparison group of 09abc directories
# ex: find ../09a_DE_* -iname "*de_ge*counts*" | xargs -n1 -I{} cat {} | grep -v "Output_File" >> ./input_file.csv

library(dplyr)
library(stringr)
library(readr)
library(tidyr)

args <- commandArgs(T)

input_file <- args[1]
comparisons <- args[2]
samples <- args[3]
output_name <- args[4]

# read table with DE_gene_counts
up.down.genes <- read_delim(input_file, delim = " ", col_names = F, trim_ws = T) %>% 
  mutate(comparison=as.numeric(str_extract(X1, "(\\d)+"))) %>% 
  mutate(feature=str_remove(X1, "[:alpha:]+_[:alpha:]+_[:alnum:]+_[:alpha:]+_")) %>% 
  mutate(package=ifelse(str_detect(str_to_lower(X2), "deseq"),"DESeq","edgeR")) %>% 
  select(comparison, feature, package, Up_genes = X4, Down_genes = X5)

# read comparison df
cmp <- read_delim(comparisons, delim = ";", col_names = T, trim_ws = T) %>% 
  select(comparison = Comparison_Number, Condition_1, Condition_2)

# read samples df
samples <- read_delim(samples, delim=";", col_names = T, trim_ws = T) %>% 
  select(Group, Condition_Name) %>% 
  distinct()

# join comparison and samples
cmp_samples <- left_join(cmp, samples, by=c("Condition_1"="Group")) %>% 
  left_join(., samples, by=c("Condition_2"="Group")) %>% 
  select(comparison, Condition_1, Condition_2, Condition_Name_1 = Condition_Name.x, Condition_Name_2 = Condition_Name.y)

# final table
left_join(up.down.genes,cmp_samples) %>% 
  select(comparison,feature,package,Condition_1,Condition_2,Condition_Name_1,Condition_Name_2,Up_genes,Down_genes) %>% 
  gather(variable, value, -(comparison:Condition_Name_2)) %>% 
  unite(combine,feature,variable, sep = ".") %>% 
  spread(combine, value) %>% 
  arrange(comparison) %>% 
  write_csv(., output_name)
