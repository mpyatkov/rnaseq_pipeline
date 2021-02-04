######################################################################################
# Kritika Karri, 08.04.2017
# This script takes the differential expression file (.txt) for Feature Count method starting with the following prefix DiffExp_v2_GeneBody from Output_DiffExp_1a_HTSeq_GeneBody folder for all the DE comparisons.
# The result fo this script is Pearson correlation heatmaps and matrix (.csv) file for every sample replicate and also a merge sample file with correlation values for all samples and their replicates.

##################################################################################

args <- commandArgs(T)
# DATASET_LABEL <- args[1]

EXTARGS=list(dataset_label=args[1],
             detitle=args[2])

#---------------------------------------------------------------------------------
wd <- getwd()
if(!is.null(wd)) {setwd(wd)}

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
library(tidyr)
library(RColorBrewer)
library(tictoc)
############################## FUNCTIONS ######################################

# save correlation plot and table
plot_cor <- function(df, title, method, out_name)
{
   WRITE_FILE <- TRUE

   if (dim(df)[1] <= 1) {
      message(paste0("WARNING: ", out_name, " does not have any significant/non-significant genes, table is empty, correlation plot will not be created"))

      # create empty ggplot
      rpkm_plot<-ggplot() + ggtitle(title) + theme_bw()+
      theme_minimal()+ 
      theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                       size = 13, hjust = 1),
            axis.text.y = element_text(size=13),
            text = element_text(size=12),
            legend.title = element_text(size=14),
            legend.text=element_text(size=14))

    WRITE_FILE <- FALSE

   } else {
   
   # create correlation matrix
   td <- df %>% 
      drop_na() %>% 
      as.matrix() %>% 
      cor(., method = method) 
   
   # remove duplicated rows and columns
   # duplicated.columns <- duplicated(t(td))
   # duplicated.rows <- duplicated(td)
   # td <- td[!duplicated.rows,!duplicated.columns]
   
   # get only lower triangle of matrix
   td[upper.tri(td)] <- NA   
   
   # create pairwise table from matrix [Var1, Var2, value]
   td_melt <- expand.grid(dimnames(td)) %>% 
      cbind(., value = as.vector(td)) %>% 
      drop_na() 
   
   # limits for legend
   mn <- round(min(td, na.rm = T), 3) - 0.01
   mx <- round(max(td, na.rm = T), 3)

   # adaptive font size
   font_size <- function(x) {
   	     f <- floor(-0.34*x+11.36)-1
	     ifelse(f<2,2,f)
   }
   
   # plot correlation
   rpkm_plot <- ggplot(data = td_melt, aes(Var1, Var2, fill = value))+
      ggtitle(title)+
      geom_tile(aes(fill = value)) +
      geom_text(aes(label = round(value, 2)), size = font_size(dim(td)[1]))+
      scale_fill_gradientn(colours = rainbow(4),limits=c(mn,mx), space= "Lab",
                           name=paste0(method,"\ncorrelation")) +
      theme_minimal()+ 
      theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                       size = 13, hjust = 1),
            axis.text.y = element_text(size=13),
            text = element_text(size=12),
            legend.title = element_text(size=14),
            legend.text=element_text(size=14))+
      coord_fixed()
  }
   # Save the correlation plot
   
   ggsave(paste0(out_name,".pdf"),plot= rpkm_plot, device = "pdf",path= wd, width = 40, height =30, units = "cm")
   if (WRITE_FILE) {
      write.table(td, file= paste0(out_name,".csv")) # keeps the rownames
   }
}

plot_pca <- function(df_pca, title, out_names) {
   if (dim(df_pca)[2] <= 2) {
      message(paste0("WARNING: ", out_names, " does not have any significant/non-significant genes, table is empty, PCA will not be created"))
      
      # create empty ggplot
      pca<-ggplot() + ggtitle(title) + theme_bw()+
         theme(text = element_text(size=14),
               legend.title = element_text(size=15),
               legend.text=element_text(size=15))
   } else {
      
      pca <- autoplot(prcomp(subset(df_pca, select = -c(group))), data= df_pca, label = F)+
         ggtitle(title)+
         geom_label_repel(aes(label = rownames(df_pca), colour = group), size = 6)+
         # scale_colour_brewer(type="qual", palette = 2)+
         theme_bw()+
         theme(text = element_text(size=14),
               legend.title = element_text(size=15),
               legend.text=element_text(size=15))
   }
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
   
   list(signif=list(name="signif", cond=list_data_cond), nonsignif=list(name="nonsignif", cond=list_data_cond1))
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
read_dataset <- function(filenames, group_names, filter_condition = NULL, meanonly = F) {
   dataset <- cbind.data.frame(lapply(filenames, fread, header=TRUE, sep="\t"))
   dataset <- dataset[unique(names(dataset))]
   
   # choose only required columns (dataset for pca)
   dataset_grep_filter <- dataset %>%
      select(matches("_edgeRFDR|_edgeRlogFC|rpkm")) 
   # %>% 
   #    select(-contains("_mean_")) 
   # 
   
   # filter conditions
   if (!is.null(filter_condition)) {
      dataset_grep_filter <- filter(dataset_grep_filter, eval(parse(text=filter_condition$cond)))
      if (filter_condition$name == "nonsignif") {
         # additional filtration for rpkm for non-significant
         dataset_grep_filter <- dataset_grep_filter[which(apply(dataset_grep_filter, 1, function(x) max(x) > 1)),]
      }
   }
   
   # dataset for correlation
   dataset_grep_filter <- dataset_grep_filter %>% 
      select(matches("rpkm_")) 
   
   # select columns with "_mean_" only
   if (meanonly == T) {
      for_cor <- dataset_grep_filter %>% 
         select(contains("_mean_")) 
   } else {
      for_cor <- dataset_grep_filter %>% 
         select(-contains("_mean_"))
   }
   
   for_pca <- dataset_grep_filter %>%
      select(-contains("_mean_"))
   
   #extract only RPKM columns 
   df_pca <- t(for_pca)
   rownames(df_pca) <- str_remove(rownames(df_pca), "rpkm_")
   
   # associate rownames(condname_sample) and groups by id(condname_sample)
   groups <- tibble(id = rownames(df_pca)) %>% 
      left_join(., group_names, by=c("id")) %>% 
      select(-id) %>% pull(group)

   # add group column to df_pca
   df_pca <- as.data.frame(df_pca)
   df_pca$group <- groups
   
   list(df_pca=df_pca, df_cor=for_cor)
}

# plot PCA for significant, nonsignificant, and all conditions
pca_3plot <- function(filename, group_names, extargs, number, meanonly = F) {
   
   # get files names without suffixes and prefixes   
   out_names <- get_outnames(filename)
   
   # get significant and not-significant genes
   cond <- get_conditions(out_names)
   
   # reassign out_names 
   out_names <- ifelse(length(out_names) > 1, "Merge", out_names)
   
   # dataset without condition
   ds <- read_dataset(filename, group_names, meanonly = meanonly)
   out_fname <- paste0(number, "_PCA_All_",out_names)
   header <- paste0(extargs$dataset_label, ", ", out_fname, ", ", "\nAll genes (without filter), ",extargs$detitle,"_Genes: ", dim(ds$df_pca)[[2]])
   plot_pca(ds$df_pca, header, out_fname)
   
   # dataset with significant genes
   ds <- read_dataset(filename, group_names, filter_condition = cond$signif, meanonly = meanonly)
   out_fname <- paste0(number, "_PCA_Significant_",out_names)
   header <- paste0(extargs$dataset_label, ", ", out_fname, ", ", "\nSignificant genes (|FC|>2 and FDR<0.05), ",extargs$detitle,"_Genes: ", dim(ds$df_pca)[[2]])
   plot_pca(ds$df_pca, header, out_fname)
   
   # dataset with non-significant genes
   ds <- read_dataset(filename, group_names, filter_condition = cond$nonsignif, meanonly = meanonly)
   out_fname <- paste0(number, "_PCA_Non-significant_", out_names)
   header <- paste0(extargs$dataset_label, ", ",  out_fname, ", ", "\nNon-significant genes (1.2 < |FC| < 1/|1.2|  and FDR >0.1, RPKM >1), ",extargs$detitle,"_Genes: ", dim(ds$df_pca)[[2]])
   plot_pca(ds$df_pca, header, out_fname)
}

# plot PCA for significant, nonsignificant, and all conditions
cor_3plot <- function(filename, group_names, method, extargs, number, meanonly = F) {
   
   # get files names without suffixes and prefixes   
   out_names <- get_outnames(filename)
   
   # get significant and not-significant genes
   cond <- get_conditions(out_names)
   
   # reassign out_names 
   out_names <- ifelse(length(out_names) > 1, "Merge", out_names)
   mname <- ifelse(method == "spearman", "Spearman", "Pearson")
   
   # double starting zero of outname in case of mean-only
   double_number <- function(outname, meanonly) {
      if (meanonly) {
         paste0("0", outname)
      }
      else {
         outname
      }
   }
   
   # dataset without condition
   ds <- read_dataset(filename, group_names, meanonly = meanonly)
   out_fname <- paste0(number, "_", mname, "_All_", out_names)
   out_fname <- double_number(out_fname, meanonly)
   header <- paste0(extargs$dataset_label, ", ", out_fname, ", ", "\nAll genes (without filter), ",extargs$detitle,"_Genes: ", dim(ds$df_cor)[[1]])
   plot_cor(ds$df_cor, header, method, out_fname)

   # dataset with significant genes
   ds <- read_dataset(filename, group_names, filter_condition = cond$signif, meanonly = meanonly)
   out_fname <- paste0(number, "_", mname, "_Significant_", out_names)
   out_fname <- double_number(out_fname, meanonly)
   header <- paste0(extargs$dataset_label, ", ", out_fname, ", ", "\nSignificant genes (|FC|>2 and FDR<0.05), ",extargs$detitle,"_Genes: ", dim(ds$df_cor)[[1]])
   plot_cor(ds$df_cor, header, method, out_fname)
   
   # dataset with non-significant genes
   ds <- read_dataset(filename, group_names, filter_condition = cond$nonsignif, meanonly = meanonly)
   out_fname <- paste0(number, "_", mname, "_Non-significant_", out_names)
   out_fname <- double_number(out_fname, meanonly)
   header <- paste0(extargs$dataset_label, ", ",  out_fname, ", ", "\nNon-significant genes (1.2 < |FC| < 1/|1.2|  and FDR >0.1, RPKM >1), ",extargs$detitle,"_Genes: ", dim(ds$df_cor)[[1]])
   plot_cor(ds$df_cor, header, method, out_fname)
}

####################################### FUNCTION ENDS ###################################

# list all txt files from the current directory
# use the pattern argument to define a common pattern  for import files with regex. Here: .txt
list.filenames.HT <- list()
list.filenames.HT<-list.files(pattern=".txt$")

list.filenames.HT

if(!is.na(list.filenames.HT[1])){
   
   group_names <- read_tsv("./SAMPLE_CONDNAME.tsv", col_names = F) %>% 
      select(id = X3, group = X2)
   
   dataset <- read_dataset(list.filenames.HT, group_names)
   #cor_3plot(list.filenames.HT, "pearson", EXTARGS, 0, meanonly = F)

   ## tSNE   
  
   # if files > 2 than perplexity = 2 else perplexity = 1
   perplexity = ifelse(length(as.vector(list.filenames.HT)) > 1, 2, 1)	

   set.seed(42)
   #tsne_model_1 <- Rtsne(as.matrix(dataset$df_pca), check_duplicates=FALSE, pca=TRUE, perplexity=perplexity, theta=0.5, dims=3, set.seed=TRUE)
   tsne_model_1 <- Rtsne(as.matrix(subset(dataset$df_pca, select = -c(group))), check_duplicates=FALSE, pca=TRUE, perplexity=perplexity, theta=0.5, dims=3, set.seed=TRUE)
   
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
   
   tic("total")

   tic("create 3 pca plot for all genes")
   # create 3 pca plot for all genes
   pca_3plot(list.filenames.HT, group_names, EXTARGS, 0)
   toc()

   tic("create 3 Pearson correlation plot for all genes")
   # create 3 Pearson correlation plot for all genes
   cor_3plot(list.filenames.HT,group_names, "pearson", EXTARGS, 0)
   toc()
   
   tic("create 3 Spearman correlation plot for all genes")
   # create 3 Spearman correlation plot for all genes
   cor_3plot(list.filenames.HT,group_names, "spearman", EXTARGS, 0)
   toc()
   
   tic("create 3 Pearson correlation plot for all genes (meanonly)")
   # create 3 Pearson correlation plot for all genes (meanonly)
   cor_3plot(list.filenames.HT,group_names, "pearson", EXTARGS, 0, meanonly = T)
   toc()
   
   tic("create 3 Spearman correlation plot for all genes (meanonly)")
   # create 3 Spearman correlation plot for all genes (meanonly)
   cor_3plot(list.filenames.HT,group_names, "spearman", EXTARGS, 0, meanonly = T)
   toc()

   for(i in 1:length(list.filenames.HT)){
      print(paste0("-->",list.filenames.HT[i]))
      number <- str_extract(list.filenames.HT[i], "(\\d)+")
      
      # create pca plot for each group
      pca_3plot(list.filenames.HT[i], group_names, EXTARGS, number)
      
      # create correlation plots
      cor_3plot(list.filenames.HT[i], group_names, "pearson", EXTARGS, number)
      cor_3plot(list.filenames.HT[i], group_names, "spearman", EXTARGS, number)
   }
   toc()
   print("DONE")
} else{
   print("No  files in the folder !!")
   
}

