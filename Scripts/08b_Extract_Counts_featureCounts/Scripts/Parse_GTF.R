##################################################################################
#Andy Rampersaud
#10.07.15
#This R script is used to parse a GTF file for use with featureCounts
#This script works within the Extract_Counts.qsub script
#Briefly, we want to generate 2 GTF files:
#1) GTF file for gene symbols where counting will assign reads to all overlapping features
#2) GTF file for gene symbols where counting will assign reads to only 1 overlapping feature
#----------------------------------------------------------------------------------
#Input:
#<GTF_File> 	: full list of gene symbols
#<Subset_List>	: list of gene symbols where counting will assign reads to all overlapping features
#----------------------------------------------------------------------------------
#Usage: 
#Rscript Parse_GTF.R <GTF_File> <Subset_List>
#Example command to run script:
#Rscript Parse_GTF.R RefSeq_GeneBody.gtf GeneSym_List.txt
##################################################################################
#---------------------------------------------------------------------------
#Pass arguments from sh into R
args <- commandArgs(trailingOnly = TRUE)
GTF_File <- args[1]
Subset_List <- args[2]
#---------------------------------------------------------------------------
if(length(args) != 2){
print("Need 2 arguments:")
print("Usage: Rscript Parse_GTF.R <GTF_File> <Subset_List>")
}
#---------------------------------------------------------------------------
print("Print arguments:")
print("-----------------")
print("GTF_File:")
paste(GTF_File,sep="")
print("Subset_List:")
paste(Subset_List,sep="")
print("-----------------")
#---------------------------------------------------------------------------
#Need to set the dir so that this script can work in any folder:
dir <- getwd()
setwd(dir)
#---------------------------------------------------------------------------
paste("Loading ",GTF_File,sep="")
#No header line
GTF_data <- read.table(paste(dir,"/",GTF_File, sep=""), as.is = TRUE, header = FALSE, sep = "\t", fill = TRUE)
#Add colnames:
colnames(GTF_data) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
#dim(GTF_data)
#[1] 665744      9
#The attribute column has the gene symbol and accession number
#Need to parse this column to access that information
#Split by semi-colon:
Split_attribute <- data.frame(do.call('rbind', strsplit(as.character(GTF_data$"attribute"),';',fixed=TRUE)))
#Just use the first 2 columns:
Split_attribute <- Split_attribute[, 1:2]
#colnames(Split_attribute) <- c("gene_id", "transcript_id", "gene_name")
colnames(Split_attribute) <- c("gene_id", "transcript_id")
#Need to remove strings:
Split_attribute$"gene_id" <- gsub("gene_id ","",Split_attribute$"gene_id")
Split_attribute$"transcript_id" <- gsub("transcript_id ","",Split_attribute$"transcript_id")
#Split_attribute$"gene_name" <- gsub("gene_name ","",Split_attribute$"gene_name")
#Combine these columns to the data:
GTF_data <- cbind(GTF_data, Split_attribute)
#Sort by gene symbol
GTF_data <- GTF_data[order(GTF_data$"gene_id"),]
#---------------------------------------------------------------------------
paste("Loading ",Subset_List,sep="")
Subset_data <- read.table(paste(dir,"/",Subset_List, sep=""), header = TRUE)
#dim(Subset_data)
#[1] 1740    1
#Sort by gene symbol
#Since there's only 1 column need: drop = FALSE
Subset_data <- Subset_data[order(Subset_data$"GeneSym_List"), , drop = FALSE]
#---------------------------------------------------------------------------
#Need to compare gene symbol lists:
paste("Number of gene symbols in common:")
common_count <- length(intersect(GTF_data$"gene_id",Subset_data$"GeneSym_List"))
paste (common_count)
#[1] 1713
paste("Percentage overlap:")
total_count <- length(Subset_data$"GeneSym_List")
Percent <- signif((common_count/total_count) * 100 , digits = 2)
paste(Percent, "%", sep="")
#Use setdiff to find the unique GeneSym:
paste("How many GeneSym in the GTF list and not in Subset list:")
paste("(Expect high number)")
length(setdiff(GTF_data$"gene_id",Subset_data$"GeneSym_List"))
#[1] 22484
paste("How many GeneSym in the Subset list and not in GTF list:")
paste("(Due to updated RefSeq gene annotations)")
length(setdiff(Subset_data$"GeneSym_List", GTF_data$"gene_id"))
#[1] 27
#---------------------------------------------------------------------------
#Need to generate GTF file for common_count gene symbols
#Count these gene symbols: assign reads to all overlapping features
common_GeneSym <- intersect(GTF_data$"gene_id",Subset_data$"GeneSym_List")
#Get the corresponding rows from GTF_data
common_GeneSym_data <- GTF_data[GTF_data$"gene_id" %in% common_GeneSym, ]
#dim(common_GeneSym_data)
#[1] 11703    12
#View(common_GeneSym_data)
#Drop columns
#drops <- c("gene_id","transcript_id","gene_name")
drops <- c("gene_id","transcript_id")
common_GeneSym_data <- common_GeneSym_data[,!(names(common_GeneSym_data) %in% drops)]
#Sort by chr and start
common_GeneSym_data <- common_GeneSym_data[order(common_GeneSym_data$"seqname",common_GeneSym_data$"start"),]
#View(common_GeneSym_data)
#Write the common_GeneSym_data to a file:
write.table(common_GeneSym_data, file=paste("GeneSym_assign_all_features.gtf", sep =""), quote=F, sep="\t", col.names=FALSE, row.names=FALSE)
paste("Check out GeneSym_assign_all_features.gtf!",sep="")
#---------------------------------------------------------------------------
#Need to generate GTF file for GeneSym in the GTF list and not in Subset list
#Count these gene symbols: where assign reads to only 1 overlapping feature
unique_GeneSym <- setdiff(GTF_data$"gene_id",Subset_data$"GeneSym_List")
#Get the corresponding rows from GTF_data
unique_GeneSym_data <- GTF_data[GTF_data$"gene_id" %in% unique_GeneSym, ]
#dim(unique_GeneSym_data)
#[1] 654041     12
#View(unique_GeneSym_data)
#Drop columns
#drops <- c("gene_id","transcript_id","gene_name")
drops <- c("gene_id","transcript_id")
unique_GeneSym_data <- unique_GeneSym_data[,!(names(unique_GeneSym_data) %in% drops)]
#Sort by chr and start
unique_GeneSym_data <- unique_GeneSym_data[order(unique_GeneSym_data$"seqname",unique_GeneSym_data$"start"),]
#View(unique_GeneSym_data)
#Write the unique_GeneSym_data to a file:
write.table(unique_GeneSym_data, file=paste("GeneSym_assign_only1_feature.gtf", sep =""), quote=F, sep="\t", col.names=FALSE, row.names=FALSE)
paste("GeneSym_assign_only1_feature.gtf!",sep="")
#---------------------------------------------------------------------------
#Just to check:
#length(common_GeneSym) + length(unique_GeneSym) should equal length(unique(GTF_data$"gene_id"))
#length(common_GeneSym) + length(unique_GeneSym)
#[1] 24197
#length(unique(GTF_data$"gene_id"))
#[1] 24197
#Works!
#---------------------------------------------------------------------------
##################################################################################
