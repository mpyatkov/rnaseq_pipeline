remotes::install_cran("argparser", upgrade = "never", quiet = T)
library(argparser)

library(tidyverse)

# remotes::install_cran("httr2", upgrade = "never")
# library(httr2)

#remotes::install_cran("rvest", upgrade = "never")
library(rvest)

p <- arg_parser('DAVID GO enrichment download')
p <- add_argument(p,'--input_path', default="../", help="Up/Down files path")
p <- add_argument(p,'--sample_labels', default="../00_Setup_Pipeline/Sample_labels.txt", help="Sample_Labels.txt file")
p <- add_argument(p,'--comparisons', default="../00_Setup_Pipeline/Comparisons.txt", help="Comparisons.txt file")
argv <- parse_args(p)
print(argv)

## TODO:
## take UP/DOWN files path from rnaseq pipeline
## output file should contain number of pcg/nc/total

david_get_results <- function(genes, outname) {
  params_list <- list(idType = "OFFICIAL_GENE_SYMBOL",
                      uploadType = "list",
                      Mode = "paste",
                      ids = paste0(genes, collapse = "\n"),
                      uploadSpecies = 10090,
  
                      pasteBox = paste0(genes, collapse = "\n"),
                      Identifier = "OFFICIAL_GENE_SYMBOL",
                      myLists = "1",
                      speciesSelect = "Mus musculus")
  
  s <- session("https://davidbioinformatics.nih.gov/tools.jsp")
  #s$response$cookies$value
  #SESSIONID = s$response$cookies$value
  
  s <- html_form(s)[[1]] %>%                                      ## choose first form
    html_form_set(!!!params_list) %>%                             ## put params into the form
    session_submit(s, form = .,submit = "B52") %>%                ## submit form and return session
    session_follow_link(i = 7) %>%                                ## click on the 7th link on the page
    session_jump_to("term2term.jsp?annot=28,36,44&currentList=0") ## add only FAT terms for analysis
  
  #full to session jump: "https://david.ncifcrf.gov/term2term.jsp?annot=28,36,44&currentList=0"
  ##Sys.sleep(5)
  
  download_link <- s %>% 
    read_html() %>%                       ## load html from session
    html_elements(css = "a") %>%          ## extract all links
    html_attr("href") %>%                 ## extract all 'href' attrs
    keep(~str_detect(.x,"t2t")) %>%       ## find only links to download file (starts with "t2t")
    str_c("https://david.ncifcrf.gov/",.) ## concatenate url and link to file
  
  ## example of download link: https://david.ncifcrf.gov/data/download/t2t_4057F39AA8F01718198056989.txt
  
  tryCatch(
    download.file(url = download_link, destfile = outname), 
    error = function(e) print(paste('NO COMMON PATHWAYS: ', outname))
  ) 
}

### test
## argv$input_path <- "/projectnb/wax-dk/max/G221_RNASEQ_XE96/"

## TODO: check keep expression below, because I am not sure
## that one term will work with for calculating enrichment score

lf <- list.files(path=argv$input_path, pattern = "^Up|^Down", recursive = T, full.names = T) %>%
  discard(~str_detect(.x, "DESeq|SEGEX_Upload_Files|14_final_summary")) %>% 
  keep(\(f){tmp <- read_tsv(f,col_names = T,show_col_types = F); nrow(tmp)>1})

walk(lf, \(fname){
  
  print(fname)
  
  compar_num <- str_extract(fname,"(?<=09d_DE_)([[:alnum:]]+)(?=_Ref)") %>%
    str_pad(., width = 2, side = "left", pad="0")
  
  output_fname <- basename(fname) %>% tools::file_path_sans_ext()  %>%  paste0(".txt") %>% paste0(compar_num,"_",.)
  print(output_fname)
  
  genes <- read_tsv(fname,col_names = T, show_col_types = F) %>% 
    select(id = 1) %>% 
    filter(!str_detect(id,"nc_")) %>% 
    pull(id)
  
  paste0("Genes: ",length(genes)) %>% print()
  
  david_get_results(genes, output_fname)
  Sys.sleep(5)
})

## save html for debug purpose
## library(xml2)
## s %>% read_html() %>% write_xml(file="test1.html")


## debug only
## big set of genes
# library(readxl)
# genes <- readxl::read_xlsx("/projectnb2/wax-dk/max/DW_DAVID_ES/G224_Cytoplasmic_ChromatinBound_up1.xlsx", col_names = T) %>% 
#   filter(`03_Down_FullGeneBody_G224_ChromatinBoundMale_vs_ChromatinBoundFemale_total2721_pcg1632_lncrna1089` == 1) %>% 
#   pull(id)

## small set of genes for DAVID web
# Fmo3
# Sult3a1
# Sult2a1
# Sult2a5
# Sult3a2
# Sult2a6
# Sult2a2
# Cyp2b13
# Sult2a4
# Cyp3a41b
# Slc22a26
# Cyp2b9
# Cyp2g1
# Cyp3a41a
# Atp6v0d2
# Prom2
# Hao2
# Ntrk2

## for script
# genes <- c("Fmo3","Sult3a1","Sult2a1","Sult2a5","Sult3a2","Sult2a6","Sult2a2","Cyp2b13","Sult2a4","Cyp3a41b","Slc22a26","Cyp2b9","Cyp2g1","Cyp3a41a","Atp6v0d2","Prom2","Hao2","Ntrk2")
