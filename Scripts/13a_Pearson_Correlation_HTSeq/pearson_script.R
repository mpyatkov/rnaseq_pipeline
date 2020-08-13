######################################################################################
# Kritika Karri, 08.04.2017
# This script takes the differential expression file (.txt) for Feature Count method starting with the following prefix DiffExp_v2_GeneBody from Output_DiffExp_1a_HTSeq_GeneBody folder for all the DE comparisons.
# The result fo this script is Pearson correlation heatmaps and matrix (.csv) file for every sample replicate and also a merge sample file with correlation values for all samples and their replicates.

##################################################################################
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
#library(autoplotly)
#library(DataCombine)
library(tidyverse)
require(dplyr)
library(miscTools)
library("ggfortify")
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
 
 
####################################### FUNCTION ENDS ###################################
 
# list all txt files from the current directory
# use the pattern argument to define a common pattern  for import files with regex. Here: .txt
list.filenames.HT <- list()
list.filenames.HT<-list.files(pattern=".txt$")
list.filenames.HT
list.filenames.HT_prefix <- str_remove(list.filenames.HT, "DiffExp_v2_GeneBody_")
list.filenames.HT_prefix <- str_remove(list.filenames.HT_prefix, "_HTSeq.txt")
list.filenames.HT_prefix

if(!is.na(list.filenames.HT[1])){
   
   dataset = cbind.data.frame(lapply(list.filenames.HT, fread, header=TRUE, sep="\t"))
   dataset_grep <- dataset[, grepl( "rpkm_" , names( dataset ) ) ]
   dataset_grep_pc <- dataset_grep[, -grep("_mean_", colnames(dataset_grep))]
   df <- t(dataset_grep_pc)
   rownames(df) <- str_remove(rownames(df), "rpkm_")
   df1 <- t(dataset_grep_pc)
   row <- str_sub(rownames(df), 1, str_length(rownames(df))-1)
   row <- as.factor(t(row))
   df <- insertCol(df, 1, row, "label")
   df <- as.data.frame(df)
   rownames(df1) <- str_remove(rownames(df1), "rpkm_")
   # PCA plot
   pca <- autoplot(prcomp(df1), data= df, colour= "label", label = TRUE, label.size = 8)
   pca <- pca + ggtitle(paste("All genes (without filter), Genes:", length(rownames(dataset_grep_pc)),sep=" "))
   # correlation
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
   
   # save PCA plot Raw
   ggsave("PCA_All_Merge_File.pdf",plot= pca, device = "pdf",path= wd, width = 60, height =30, units = "cm")
   
   
  ####################################### filtered analysis#########################333
  #create condition statements using the prefix where logFC >1 and FDR < 0.05 for same sample and or condition 
  #across other  samples.
 
  listFC <- NULL
  listFDR <- NULL
  cond<- NULL  #cond is significant genes defined as FC >2 and FDR < 0.05
  cond1<- NULL #cond1 is non-significant genes defined as 1.2 < FC < 1/|1.2|  and FDR >0.1
  for(i in 1:length(list.filenames.HT_prefix))
  {   
    listFC[i] <- paste(list.filenames.HT_prefix[i],"_edgeRlogFC", sep="")
    listFDR[i] <- paste(list.filenames.HT_prefix[i],"_edgeRFDR",sep="")
     cond[i] <-  paste("((","abs(",listFC[i],")",">1)","&","(", listFDR[i],"<0.05))", sep=" ")
     cond1[i] <-  paste("((",listFC[i],")","<0.263034","&(","(",listFC[i],")",">","-0.263)","&","(", listFDR[i],">0.1))", sep=" ")
  }
    
  list_data_cond <- (paste(cond, collapse = "|") )                     # merge the list into string
  list_data_cond1 <- (paste(cond1, collapse = "|") )                     # merge the list into string
  
  
  ###################################################################################### 
  # filtering of the samples.
   
   myVectorOfStrings <- c("_edgeRFDR", "_edgeRlogFC","rpkm")
   matchExpression <- paste(myVectorOfStrings, collapse = "|")  ## to get multiple columns from edgeR_FC, edgeR_FDR, and rpkm
   dataset = cbind.data.frame(lapply(list.filenames.HT, fread, header=TRUE, sep="\t"))
   dataset_grep_filter <- dataset[, grepl( matchExpression,names(dataset))]
   dataset_grep_filter <- dataset_grep_filter[, -grep("_mean_", colnames(dataset_grep_filter))]
   
   ### evaluate the  command that filters basedon the condition created above |logFC|>1, or |FC|>2 AND FDR<0.05
   FC_1 <- eval(parse(text=paste("filter(dataset_grep_filter,", list_data_cond,")", sep="")))  
   
   dataset_grep_FC <- FC_1[,grepl("rpkm_", names(FC_1))]     #extract only RPKM columns
   df_FC <- t(dataset_grep_FC)                               #transpose the matrix 
   rownames(df_FC) <- str_remove(rownames(df_FC), "rpkm_")   #remove the rpkm from colnames for plot
   df2 <- t(dataset_grep_FC)                                 # this is a copy of trasposed matrix for label creation
   row_FC <- str_sub(rownames(df2), 1, str_length(rownames(df2))-1) #remove the numbers for replicates in labels
   row_FC <- as.factor(t(row_FC))                              
   df2 <- insertCol(df2, 1, row_FC, "label")                 #insert label colname for the df2
   df2 <- as.data.frame(df2)
     
   # PCA plot for significant genes
   pca_FC <- autoplot(prcomp(df_FC), data= df2, colour= "label", label = TRUE, label.size = 8)
   pca_FC <- pca_FC+ ggtitle(paste("significant genes (|FC|>2 and FDR<0.05), Genes:", length(rownames(dataset_grep_FC)),sep=" "))
   
   ############ evaluate non-signifcant PCGs ############################
   
   FC_2 <- eval(parse(text=paste("filter(dataset_grep_filter,", list_data_cond1,")", sep=""))) #non siginifcant genes are 0.243< log2FC < -0.243 and FDR >0.1 
   dataset_grep_FC2 <- FC_2[,grepl("rpkm_", names(FC_2))]  
   
   #write.csv(dataset_grep_FC2, "dataset_grep_FC2.txt")
   
   dataset_grep_FC2.1  <- dataset_grep_FC2[which(apply(dataset_grep_FC2,1,function(x) max(x) > 1)),]
   
   df_FC2 <- t(dataset_grep_FC2.1)
   dim(df_FC2)
   rownames(df_FC2) <- str_remove(rownames(df_FC2), "rpkm_")
   df3 <- t(dataset_grep_FC2.1)
   row_FC3 <- str_sub(rownames(df3), 1, str_length(rownames(df3))-1)
   row_FC3 <- as.factor(t(row_FC3))
   df3 <- insertCol(df3, 1, row_FC3, "label")
   df3 <- as.data.frame(df3)
   
   # PCA plot
   pca_FC2 <- autoplot(prcomp(df_FC2), data= df3, colour= "label", label = TRUE, label.size = 8)
   pca_FC2 <- pca_FC2+ ggtitle(paste("Non-significant genes (1.2 < |FC| < 1/|1.2|  and FDR >0.1, RPKM >1), Genes:", length(rownames(dataset_grep_FC2.1)),sep=" "))
   
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
  
  # Save PCA plot
  ggsave("PCA_Significant_File.pdf",plot= pca_FC, device = "pdf",path= wd, width = 60, height =30, units = "cm")
  ggsave("PCA_Non-Significant_File.pdf",plot= pca_FC2, device = "pdf",path= wd, width = 60, height =30, units = "cm")
  
 
 for(i in 1:length(list.filenames.HT)){
 
   analyze_all(list.filenames.HT[i])
   analyze_rpkm(list.filenames.HT[i])
   
   
 }
 
 } else{
   print("No  files in the folder !!")
      
 }
 
 
 
 
