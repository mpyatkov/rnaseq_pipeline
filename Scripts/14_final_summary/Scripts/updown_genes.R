# up and down genes summarizing for each comparison group of 09abc directories
# ex: find ../09a_DE_* -iname "*de_ge*counts*" | xargs -n1 -I{} cat {} | grep -v "Output_File" >> ./input_file.csv

library(dplyr)
library(stringr)
library(readr)
library(tidyr)

args <- commandArgs(T)

comparisons <- args[1]
samples_file <- args[2]
output_name <- args[3]
dataset_label <- args[4]
## CUSTOM_FC <- as.double(args[5])
## CUSTOM_FDR <- as.double(args[6])

CUSTOM_FC <- eval(parse(text=args[5]))
CUSTOM_FDR <- eval(parse(text=args[6]))

files <- list.files(pattern = "forSEGEX")

updown_one <- function(f, FCparam, FDRparam, Dataset_label) {
  comparison <- as.numeric(str_extract(f, "(\\d)+"))
  feature <- str_extract(f, "ExonCollapsed|FullGeneBody|ExonOnly|ExonicOnly|IntronOnly|IntronicOnly")
  package <- ifelse(str_detect(str_to_lower(f), "deseq"),"DESeq","edgeR")
  df <- read_tsv(f, col_names = T, trim_ws = T, comment = "#") %>% 
    select(FC = 3, FDR = 6) %>% 
    filter(abs(FC) > FCparam, FDR < FDRparam) %>% 
    mutate(Up_genes = ifelse(FC>FCparam,1,0),
           Down_genes= ifelse(FC<FCparam,1,0)) %>% 
    summarise(Up_genes = sum(Up_genes), Down_genes = sum(Down_genes))
  cbind(Dataset_label, FCparam, FDRparam, comparison, feature, package,df)
}

## up.down.genes <- dplyr::bind_rows(lapply(files, function(fn){updown_one(fn, CUSTOM_FC, CUSTOM_FDR, dataset_label)}))

## get list of data.frames with different FC and FDR
up.down.genes <- lapply(1:length(CUSTOM_FC), function(ix) {
    
    dplyr::bind_rows(lapply(files, function(fn){updown_one(fn, CUSTOM_FC[ix], CUSTOM_FDR[ix], dataset_label)}))
})

## combine all tables together
up.down.genes <- dplyr::bind_rows(up.down.genes)


# read table with DE_gene_counts
## up.down.genes <- read_delim(input_file, delim = " ", col_names = F, trim_ws = T) %>% 
##   rowwise() %>% 
##   mutate(comparison=as.numeric(str_extract(X1, "(\\d)+")),
##          feature=str_remove(X1, "[:alpha:]+_[:alpha:]+_[:alnum:]+_[:alpha:]+_"),
##          package=ifelse(str_detect(str_to_lower(X2), "deseq"),"DESeq","edgeR"),
##          # Condition_1_GM_Numbers=str_extract_all(X2,"G[:alnum:]+_[:alnum:]+", simplify = T)[2],
##          # Condition_2_GM_Numbers=str_extract_all(X2,"G[:alnum:]+_[:alnum:]+", simplify = T)[1],
##          Dataset_label=dataset_label) %>% 
##   select(Dataset_label, comparison, feature, package, Up_genes = X4, Down_genes = X5)

# comparison df
cmp <- read_delim(comparisons, delim = ";", col_names = T, trim_ws = T, comment = "#") %>% 
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
samples <- read_delim(samples_file, delim=";", col_names = T, trim_ws = T, comment="#") %>% 
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
  select(Dataset_label, FC = FCparam, FDR = FDRparam, comparison,feature,package,Condition_1,Condition_2,
         Condition_Name_1,Condition_Name_2,
         Condition_1_GM_Numbers,Condition_2_GM_Numbers,
         Up_genes,Down_genes) %>% 
  gather(variable, value, -(Dataset_label:Condition_2_GM_Numbers)) %>% 
  unite(combine,feature,variable, sep = ".") %>% 
  spread(combine, value) %>% 
  arrange(comparison) %>% 
#  write_excel_csv2(., output_name) # separator ";"
  write_csv(., output_name) # separator ","


