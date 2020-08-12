##############################################################################################
#Andy Rampersaud
#03.12.2016
#This R script is used to compare differentially expressed (DE) gene lists
#The output is a Venn diagram showing the count of common and unique genes
#This scripts works in conjuction with Diff_Genes.R (part of the DiffExp_* job of the RNA-Seq pipeline)
#Output of Diff_Genes.R is used as input for this script
#Tisha adapted the script for lncRNA
#----------------------------------------------------------------------------------
#Notes:
#<File1_Data>		= Output of Diff_Genes.R (Down/Up_Genes_DESeq_LncRNA_ExonCollapsed.txt)
#<File2_Data> 		= Output of Diff_Genes.R (Down/Up_Genes_DESeq_LncRNA_GeneBody.txt)
#<File3_Data> 		= Output of Diff_Genes.R (Down/Up_Genes_DESeq_LncRNA_Exonic_Only.txt)
#<File4_Data>	 	= Output of Diff_Genes.R (Down/Up_Genes_DESeq_LncRNA_Intronic_Only.txt)
#----------------------------------------------------------------------------------
#Usage: 
#Rscript Venn_Count_Method.R <File1_Data> <File2_Data> <File3_Data> <File4_Data>
#Example commands to run script:
#--------------------------------------------------------------------------------
#Rscript Venn_Count_Method.R Down_Genes_DESeq_LncRNA_ExonCollapsed.txt Down_Genes_DESeq_LncRNA_GeneBody.txt Down_Genes_DESeq_LncRNA_Exonic_Only.txt Down_Genes_DESeq_LncRNA_Intronic_Only.txt
#--------------------------------------------------------------------------------
##############################################################################################
#---------------------------------------------------------------------------------------------
#Load packages:
#install.packages("VennDiagram")
library(VennDiagram)
#---------------------------------------------------------------------------------------------
#Pass arguments from sh into R
args <- commandArgs(trailingOnly = TRUE)
# args <- c("Down_Genes_EdgeR_LncRNA_ExonCollapsed_Male_G83_M1M2_Female_G83_M3M4_featureCounts.txt",
#           "Down_Genes_EdgeR_LncRNA_Exonic_Only_Male_G83_M1M2_Female_G83_M3M4_featureCounts.txt",
#           "Down_Genes_EdgeR_LncRNA_GeneBody_Male_G83_M1M2_Female_G83_M3M4_featureCounts.txt",
#           "Down_Genes_EdgeR_LncRNA_Intronic_Only_Male_G83_M1M2_Female_G83_M3M4_featureCounts.txt")

#---------------------------------------------------------------------------------------------
if(length(args) != 4){
print("Need 4 arguments:")
print("Usage: Rscript Venn_Count_Method.R <File1_Data> <File2_Data> <File3_Data> <File4_Data>")
}
#---------------------------------------------------------------------------------------------
print("Print arguments:")
print("-----------------")
for(i in seq(1,length(args))){
  print(paste("File", i, ":"))
  print(args[i])
}
print("-----------------")
#---------------------------------------------------------------------------------------------
#Need to set the dir so that this script can work in any folder:
dir <- getwd()
setwd(dir)
#---------------------------------------------------------------------------------------------
FileDataList <- list()
FileDataLabelList <- list()
for(i in seq(1,length(args))){
  #Read in current File
  File_Data <- read.table(file=args[i], header=TRUE, as.is=TRUE, fill = TRUE)
  #Just rename 1st column:
  colnames(File_Data)[1] <- "Gene.Symbol"
  #Useful to get the file name
  File_Name <- gsub(".txt","", args[i])
  #View(File1_Data)
  #Need a label with only:
  #<Direction><R_Package><feature_counted><Count_Program>
  #Need a series  of ifelse(grepl) commands:
  Direction <- ifelse(grepl("Down",File_Name),"Down","Up")
  # print(Direction)
  R_Package <- ifelse(grepl("DESeq",File_Name),"DESeq","EdgeR")
  # print(R_Package)
  if(grepl("ExonCollapsed", File_Name, ignore.case = FALSE)) {
    feature_counted <- "ExonCollapsed"
  } else if(grepl("GeneBody", File_Name, ignore.case = FALSE)) {
    feature_counted <- "GeneBody"
  } else if(grepl("Exonic_Only", File_Name, ignore.case = FALSE)) {
    feature_counted <- "Exonic_Only"
  } else if(grepl("Intronic_Only", File_Name, ignore.case = FALSE)) {
    feature_counted <- "Intronic_Only"
  }
  #End of if statement
  # print(feature_counted)
  Count_Program <- ifelse(grepl("featureCounts",File_Name),"featureCounts","HTSeq")
  # print(Count_Program)
  #Once I get each part; paste them together
  File_Label <- paste(Direction, R_Package,feature_counted,  Count_Program, sep=".")
  print(File_Label)
  FileDataLabelList[[i]] <- File_Label
  
  #order data frame by gene symbol
  File_Data <- File_Data[order(File_Data$"Gene.Symbol"),]
  FileDataList[[i]] <- File_Data
}


#---------------------------------------------------------------------------------------------
#Now we can compare gene lists
#Use the VennDiagram package
#Creating a venn diagram
#---------------------------------------------------------------------------------------------
GeneSym_ListObj <- lapply(FileDataList, function(i){ return(i$Gene.Symbol) })
#Need to add names to the list object for the Venn diagram:
names(GeneSym_ListObj) <- FileDataLabelList
#Need to color Venn diagrams based in direction of gene regulation
#(red=upregulation) or  (blue=downregulation)
#Need grepl to return true/false and need ignore.case = FALSE
if(sum(grepl("Down", FileDataLabelList, ignore.case = FALSE)) == length(FileDataLabelList)) { #all 4 files says down
  File_Color_list <- topo.colors(length(FileDataLabelList)) #First few colors are blue, then green then yellow
} else if(sum(grepl("Up", FileDataLabelList, ignore.case = FALSE)) == length(FileDataLabelList)) { #all 4 files says up
  File_Color_list <- heat.colors(length(FileDataLabelList)) #Red colors
}
#End of if statement

#Create the Venn diagram:
Venn_Diagram <- venn.diagram(
x = GeneSym_ListObj,
filename=NULL,
lwd = 4,
fill = File_Color_list,
alpha = 0.75,
label.col ="black",
cex = 2.0,
cat.default.pos = "text",
fontfamily = "serif",
fontface = "bold",
cat.col = c("black", "black", "black", "black"),
cat.cex = 1.0,
cat.fontfamily = "serif",
cat.fontface = "bold",
#cat.dist = c(0.06, 0.06, 0.06, 0.6),
cat.pos = 0,
#Adjust so labels are visible:
#List (length = 1/2/3/4 based on set number) of Vectors of length 2 indicating horizontal and vertical justification for each category name
#The 1st number controls left/right position of label
#The 2nd number controls up/down position of label
cat.just=list(c(0.5,0), c(0.5,0), c(0.5,1), c(0.5,-2)),
main = "Gene Symbol Comparison",
main.cex = 2,
sub.cex = 2
);
# grid.draw(Venn_Diagram)
#Options explained:
#http://www.inside-r.org/packages/cran/VennDiagram/docs/venn.diagram
#Create the image file
png(file=paste("Venn_", paste(FileDataLabelList, sep="", collapse="."), ".png", sep=""))
grid.draw(Venn_Diagram)
dev.off()
#---------------------------------------------------------------------------------------------
print("Removing VennDiagram*.log files")
system("rm VennDiagram*.log")
#---------------------------------------------------------------------------------------------
#Also want a table of gene counts for the common and unique gene symbols
paste("Check out Venn diagram!",sep="")
##############################################################################################
