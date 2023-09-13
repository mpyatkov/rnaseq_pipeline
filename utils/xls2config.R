#Rscript xls2config config.xlsx
## Should extract configuration files from xls file
## each sheet - config file
## sheet name - output name for config file

library(readr)
library(readxl)
library(purrr)

args = commandArgs(trailingOnly=TRUE)

#args[1] <- "~/tmp/mnt/G201/Scripts/00_Setup_Pipeline/G201.xlsx"

sheets <- excel_sheets(args[1])

walk(sheets, function(sheet_name){
    print(paste0("Processing config: ", sheet_name))
    df <- read_excel(args[1],col_names = T, sheet = sheet_name)
    if (grepl("index",sheet_name)) {
        write_csv(df,file = sheet_name, col_names = T)
    } else {
        write_delim(df, file = sheet_name, delim = ";", col_names = T)    
    }
})

