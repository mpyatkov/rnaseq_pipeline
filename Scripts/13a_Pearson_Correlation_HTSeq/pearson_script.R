######################################################################################
# Kritika Karri, 08.04.2017
# This script takes the differential expression file (.txt) for Feature Count method starting with the following prefix DiffExp_v2_GeneBody from Output_DiffExp_1a_HTSeq_GeneBody folder for all the DE comparisons.
# The result fo this script is Pearson correlation heatmaps and matrix (.csv) file for every sample replicate and also a merge sample file with correlation values for all samples and their replicates.

##################################################################################

args <- commandArgs(T)
DATASET_LABEL <- args[1]

#---------------------------------------------------------------------------------
wd <- getwd()
if(!is.null(wd))
   setwd(wd)

require(stringr)
require(reshape2)
require(ggplot2)
require(MASS)
library(tools)
require(data.table)
library(ggfortify)
library(tidyverse)
require(dplyr)
library(miscTools)
library(caret)
library(Rtsne)
library(ggrepel)
############################## FUNCTIONS ######################################

delete.na <- function(DF, n=0) {
   DF[rowSums(is.na(DF)) <= n,]
}

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
   cormat[upper.tri(cormat)] <- NA
   return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
   cormat[lower.tri(cormat)]<- NA
   return(cormat)
}


# Analyze_rpkm function takes DiffExp_v2_Genebody (.txt) file from Differential Expression and then calculate pearson correlation and outputs matrix and heatmap
analyze_rpkm <- function(filename)
{
   filename
   data= data.frame(read.table(filename, header = T ))
   
   data_grep <- data[, grepl( "rpkm_" , names( data ) ) ]
   data_filter <- data_grep[data_grep[,grepl( "rpkm_mean" , names(data_grep ) )]>1,]
   data_filter = as.matrix(delete.na(data_filter))
   final <- cor(data_filter)
   final1 <-as.matrix(final)
   file1 <- file_path_sans_ext(filename)
   file_path_sans_ext(file1)
   
   # Write the correlation matrix
   f1 <- paste0("Pearson_Filtered_RPKM_",file1,".csv")
   write.table(final,file=f1) # keeps the rownames
   read.table(f1,header=TRUE,row.names=1) # says first column are rownames
   
   # Create the pearson plots   
   m= round(min(final),3)
   dm = m -0.01
   M= round(max(final),3)
   melted_cormat <- melt(final)
   upper_tri <- get_upper_tri(final)
   melted_cormat <- melt(upper_tri, na.rm = TRUE)
   rpkm_plot <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradientn(colours = rainbow(4),limits=c(dm,M), space= "Lab",
                           #scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           name="Pearson\nCorrelation filtered at RPKM >1") +
      theme_minimal()+ 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 12, hjust = 1))+
      coord_fixed()
   
   f <- paste0("Pearson_Filtered_RPKM",file1,".pdf")
   
   # Save the correlation plot
   ggsave(f,plot= rpkm_plot, device = "pdf",path= wd, width = 30, height = 30, units = "cm")
}

# Analyze_all function calculates the pearson correlation for all genes from DE files 
analyze_all <- function(filename)
{
   
   data= data.frame(read.table(filename, header = T ))
   
   data_grep <- data[, grepl( "rpkm_" , names( data ) ) ]
   final <- cor(data_grep)
   final1 <-as.matrix(final)
   file1 <- file_path_sans_ext(filename)
   file_path_sans_ext(file1)
   
   # Write the correlation matrix
   f1 <- paste0("Pearson_All_",file1,".csv")
   write.table(final,file=f1) # keeps the rownames
   read.table(f1,header=TRUE,row.names=1) # says first column are rownames
   
   # Create the pearson plots   
   m= round(min(final),3)
   dm = m- 0.01
   M= round(max(final),3)
   melted_cormat <- melt(final)
   upper_tri <- get_upper_tri(final)
   melted_cormat <- melt(upper_tri, na.rm = TRUE)
   rpkm_plot <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradientn(colours = rainbow(4),limits=c(dm,M), space= "Lab",
                           
                           #  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           #                      midpoint = 0, limit = c(-1,1), space = "Lab", 
                           name="Pearson\nCorrelation-All genes") +
      theme_minimal()+ 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 12, hjust = 1))+
      coord_fixed()
   
   f <- paste0("Pearson_All_",file1,".pdf")
   
   # Save the correlation plot
   ggsave(f,plot= rpkm_plot, device = "pdf",path= wd, width = 30, height = 30, units = "cm")
}

plot_pca <- function(df_pca, df, title, out_names) {
   # PCA plot
   pca <- autoplot(prcomp(df_pca), data= df, colour= "label", label = F, label.size = 6)
   pca <- pca + 
      ggtitle(paste(title, dim(df_pca)[[2]],sep=" "))+
      geom_label_repel(aes(label = rownames(df_pca), color = label), size = 6)+
      scale_colour_brewer(type="qual", palette = 2)+
      theme_bw()+
      theme(text = element_text(size=14),
            legend.title = element_text(size=15),
            legend.text=element_text(size=15))
   
   ggsave(paste0(out_names,".pdf"),plot= pca, device = "pdf",path= wd, width = 50, height =30, units = "cm")   

}

# create filter conditions for PCA plots
get_conditions <- function(filenames_prefix) {
   listFC <- NULL
   listFDR <- NULL
   cond<- NULL  #cond is significant genes defined as FC >2 and FDR < 0.05
   cond1<- NULL #cond1 is non-significant genes defined as 1.2 < FC < 1/|1.2|  and FDR >0.1
   for(i in 1:length(filenames_prefix))
   {   
      listFC[i] <- paste(filenames_prefix[i],"_edgeRlogFC", sep="")
      listFDR[i] <- paste(filenames_prefix[i],"_edgeRFDR",sep="")
      cond[i] <-  paste("((","abs(",listFC[i],")",">1)","&","(", listFDR[i],"<0.05))", sep=" ")
      cond1[i] <-  paste("((",listFC[i],")","<0.263034","&(","(",listFC[i],")",">","-0.263)","&","(", listFDR[i],">0.1))", sep=" ")
   }
   
   list_data_cond <- (paste(cond, collapse = "|") )                     # merge the list into string
   list_data_cond1 <- (paste(cond1, collapse = "|") )                     # merge the list into string
   
   list(signif=list_data_cond, nonsignif=list_data_cond1)
}

# get out_names for PCA files
get_outnames <- function(filenames) {
    dataset = cbind.data.frame(lapply(filenames, fread, header=TRUE, sep="\t"))
    dataset = dataset[unique(names(dataset))]
    # return out_names
    dataset %>% 
           select(matches("_edgeRlogFC")) %>% 
           colnames(.) %>% 
           str_remove("_edgeRlogFC")
}

# read datasets and return back files for PCA, tSNE and correlation 
read_dataset <- function(filenames, filter_condition = NULL) {
   dataset = cbind.data.frame(lapply(filenames, fread, header=TRUE, sep="\t"))
   dataset = dataset[unique(names(dataset))]
   
   # TODO: need to refactor and remove
   dataset_grep <- dataset[, grepl( "rpkm_" , names( dataset ) ) ]
   
   # choose only required columns
   dataset_grep_filter <- dataset %>%
      select(matches("_edgeRFDR|_edgeRlogFC|rpkm")) %>% 
      select(-contains("_mean_")) 
   
   if (!is.null(filter_condition)) {
      # filter conditions
      dataset_grep_filter <- filter(dataset_grep_filter, eval(parse(text=filter_condition)))
   }
   
   #extract only RPKM columns 
   dataset_grep_filter <- dataset_grep_filter %>% 
      select(matches("rpkm_"))     
   
   df_pca <- t(dataset_grep_filter)
   rownames(df_pca) <- str_remove(rownames(df_pca), "rpkm_")
   
   df <- t(dataset_grep_filter)
   rownames(df) <- str_remove(rownames(df), "rpkm_")
   row <- str_sub(rownames(df), 1, str_length(rownames(df))-1)
   row <- as.factor(t(row))
   df <- insertCol(df, 1, row, "label")
   df <- as.data.frame(df)
   list(df_pca=df_pca, df=df, dg=dataset_grep)
}

# plot PCA for significant, nonsignificant, and all conditions
pca_3plot <- function(filename, dataset_label, number) {
   # get files names without suffixes and prefixes   
   
    out_names <- get_outnames(filename)
   
   # get significant and not-significant genes
   cond <- get_conditions(out_names)
   
   # reassign out_names 
   out_names <- ifelse(length(out_names) > 1, "Merge", out_names)
   
   # dataset without condition
   ds <- read_dataset(filename)
   out_fname <- paste0(number, "_PCA_All_",out_names)
   header <- paste0(dataset_label, ", ", out_fname, ", ", "All genes (without filter), Genes:")
   plot_pca(ds$df_pca, ds$df, header, out_fname)
   
   # dataset with significant genes
   ds <- read_dataset(filename, cond$signif)
   out_fname <- paste0(number, "_PCA_Significant_",out_names)
   header <- paste0(dataset_label, ", ", out_fname, ", ", "Significant genes (|FC|>2 and FDR<0.05), Genes:")
   plot_pca(ds$df_pca, ds$df, header, out_fname)
   
   # dataset with non-significant genes
   ds <- read_dataset(filename, cond$nonsignif)
   out_fname <- paste0(number, "_PCA_Non-significant_", out_names)
   header <- paste0(dataset_label, ", ",  out_fname, ", ", "Non-significant genes (1.2 < |FC| < 1/|1.2|  and FDR >0.1, RPKM >1), Genes:")
   plot_pca(ds$df_pca, ds$df, header, out_fname)
}

####################################### FUNCTION ENDS ###################################

# list all txt files from the current directory
# use the pattern argument to define a common pattern  for import files with regex. Here: .txt
list.filenames.HT <- list()
list.filenames.HT<-list.files(pattern=".txt$")
list.filenames.HT

if(!is.na(list.filenames.HT[1])){
   
   dataset <- read_dataset(list.filenames.HT)
   dataset_grep <- dataset$dg
   

   ## tSNE   

   # if files > 2 than perplexity = 2 else perplexity = 1
   perplexity = ifelse(length(as.vector(list.filenames.HT)) > 1, 2, 1)	

   set.seed(42)
   tsne_model_1 <- Rtsne(as.matrix(dataset$df_pca), check_duplicates=FALSE, pca=TRUE, perplexity=perplexity, theta=0.5, dims=3, set.seed=TRUE)
   
   ## getting the two dimension matrix
   
   d_tsne_1 <- as.data.frame(tsne_model_1$Y)  
   
   tsneplot <- ggplot(d_tsne_1, aes(x=V1, y=V2, z=V3), colour="green") +  
      geom_point(size=1) + 
      guides(colour=guide_legend(override.aes=list(size=3))) +
      xlab("") + ylab("") +
      ggtitle("t-SNE for all genes ") +
      theme_light(base_size=20) +
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank()) +
      scale_colour_brewer(palette = "Set2")
   
   
   tsneplot <- tsneplot + geom_label_repel(aes(label = rownames(dataset$df_pca)),
                                           box.padding   = 0.35, 
                                           point.padding = 0.5,
                                           segment.color = 'grey50') +
      theme_classic()
   ggsave("tSNE_All_Merge.pdf",plot= tsneplot, device = "pdf",path= wd, width = 80, height =40, units = "cm")   
   
   
   ## correlation
   
   final_dataset <- round(cor(dataset_grep),3)
   duplicated.columns <- duplicated(t(final_dataset))
   
   duplicated.rows <- duplicated((final_dataset)) 
   new.matrix <- final_dataset[!duplicated.rows,!duplicated.columns] 
   
   write.table(new.matrix,file= "Pearson_All_Merge_File.csv") # keeps the rownames
   read.table("Pearson_All_Merge_File.csv",header=TRUE,row.names=1) # says first column are rownames
   
   # Create the pearson plots   
   m= round(min(new.matrix),3)
   dm= m-0.01
   M= round(max(new.matrix),3)
   melted_cormat <- melt(new.matrix)
   upper_tri <- get_upper_tri(new.matrix)
   melted_cormat <- melt(upper_tri, na.rm = TRUE)
   rpkm_plot <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+
      #scale_fill_continuous(low= "red",high="green",limits=c(-1,1), space="Lab",
      scale_fill_gradientn(colours = rainbow(4),limits=c(dm,M), space = "Lab",
                           name="Pearson\nCorrelation-All genes\n") +
      theme_minimal()+ 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 12, hjust = 1))
   coord_fixed()
   rpkm_plot
   # f <- paste0("Pearson_All",file1,".pdf")
   
   # Save the correlation plot
   ggsave("Pearson_All_Merge_File.pdf",plot= rpkm_plot, device = "pdf",path= wd, width = 30, height = 30, units = "cm")
   
   ######################## for correlation ######################
   dataset_filter <- dataset_grep[dataset_grep[,grepl( "rpkm_mean" , names(dataset_grep ) )]>1,]
   dataset_filter = as.matrix(delete.na(dataset_filter))
   final_dataset_filter <- round(cor(dataset_filter),3)
   duplicated.columns.filter <- duplicated(t(final_dataset_filter))
   
   duplicated.rows.filter <- duplicated((final_dataset_filter)) 
   new.matrix.filter <- final_dataset_filter[!duplicated.rows.filter,!duplicated.columns.filter] 
   
   write.table(new.matrix.filter,file= "Pearson_Filtered_Merge_File.csv") # keeps the rownames
   read.table("Pearson_Filtered_Merge_File.csv",header=TRUE,row.names=1) # says first column are rownames
   
   # Create the pearson plots   
   m= round(min(new.matrix.filter),3)
   dm = m-0.01
   M= round(max(new.matrix.filter),3)
   melted_cormat_filter <- melt(new.matrix.filter)
   upper_tri_filter <- get_upper_tri(new.matrix.filter)
   melted_cormat_filter <- melt(upper_tri_filter, na.rm = TRUE)
   rpkm_plot_filter <- ggplot(data = melted_cormat_filter, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+
      #scale_fill_continuous(low= "red",high="green",limits=c(-1,1), space="Lab",
      scale_fill_gradientn(colours = rainbow(4),limits=c(dm,M), space = "Lab",
                           name="Pearson\nCorrelation-RPKM >1 genes\n") +
      theme_minimal()+ 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 12, hjust = 1))
   coord_fixed()
   rpkm_plot_filter
   
   # Save the correlation plot
   ggsave("Pearson_Filtered_Merge_File.pdf",plot= rpkm_plot_filter, device = "pdf",path= wd, width = 30, height = 30, units = "cm")
   

   # create 3 pca plot for all genes
   pca_3plot(list.filenames.HT, DATASET_LABEL, 0)
   
   for(i in 1:length(list.filenames.HT)){
      print(paste0("-->",list.filenames.HT[i]))
      number <- str_extract(list.filenames.HT[i], "(\\d)+")   
      # create pca plot for each group
      pca_3plot(list.filenames.HT[i], DATASET_LABEL, number)
      
      analyze_all(list.filenames.HT[i])
      analyze_rpkm(list.filenames.HT[i])

   }
   
} else{
   print("No  files in the folder !!")
   
}




