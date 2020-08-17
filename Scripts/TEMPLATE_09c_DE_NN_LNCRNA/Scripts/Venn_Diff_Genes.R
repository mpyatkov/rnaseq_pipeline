##############################################################################################
#Andy Rampersaud
#03.12.2016
#This R script is used to compare differentially expressed (DE) gene lists
#The output is a Venn diagram showing the count of common and unique genes
#This scripts works in conjuction with Diff_Genes.R (part of the DiffExp_* job of the RNA-Seq pipeline)
#Output of Diff_Genes.R is used as input for this script

#Notes:
#<File1_Data>		= Output of Diff_Genes.R (Down_Genes_DESeq.txt or Up_Genes_DESeq.txt)
#<File2_Data> 		= Output of Diff_Genes.R (Down_Genes_EdgeR.txt or Up_Genes_EdgeR.txt)
#<Subtitle>		= Plot subtitle to indicate counting method (GeneBody, Exonic_Only, Intronic_Only) (Use the COL_SUFFIX variable from the DiffExp_* job)
#<DiffExp_Index> 	= Label needed for the Count_Table (Use the DiffExp_Index variable from the DiffExp_* job)

#Usage: 
#Rscript Venn_Diff_Genes.R <File1_Data> <File2_Data> <Subtitle> <DiffExp_Index>
#Example commands to run script:

#File1=Down_Genes_DESeq_GeneBody_Male_G83_M1M2_Female_G83_M3M4_HTSeq.txt
#File2=Down_Genes_EdgeR_GeneBody_Male_G83_M1M2_Female_G83_M3M4_HTSeq.txt
#Subtitle=GeneBody
#DiffExp_Index=DiffExp_1a
#Rscript Venn_Diff_Genes.R ${File1} ${File2} ${Subtitle} ${DiffExp_Index}

#File1=Up_Genes_DESeq_GeneBody_Male_G83_M1M2_Female_G83_M3M4_HTSeq.txt
#File2=Up_Genes_EdgeR_GeneBody_Male_G83_M1M2_Female_G83_M3M4_HTSeq.txt
#Subtitle=GeneBody
#DiffExp_Index=DiffExp_1a
#Rscript Venn_Diff_Genes.R ${File1} ${File2} ${Subtitle} ${DiffExp_Index}

#Load packages:
#install.packages("VennDiagram")
library(VennDiagram)

#Pass arguments from sh into R
args <- commandArgs(trailingOnly = TRUE)
File1 <- args[1]
File2 <- args[2]
Subtitle <- args[3]
DiffExp_Index <- args[4]

# File1 <- "Down_Genes_LncRNA_ExonCollapsed_featureCounts.txt"
# File2 <- "Down_Genes_LncRNA_ExonCollapsed_featureCounts.txt"
# Subtitle <- "LncRNA_ExonCollapsed"
# DiffExp_Index <- "DiffExp_1d"

if(length(args) != 4){
print("Need 4 arguments:")
print("Usage: Rscript Venn_Diff_Genes.R <File1_Data> <File2_Data> <Subtitle> <DiffExp_Index>")
}

print("Print arguments:")
print("-----------------")
print("File1:")
paste(File1,sep="")
print("File2:")
paste(File2,sep="")
print("Subtitle:")
paste(Subtitle,sep="")
print("DiffExp_Index:")
paste(DiffExp_Index,sep="")
print("-----------------")

#Need to set the dir so that this script can work in any folder:
dir <- getwd()
setwd(dir)

#Read in File1
File1_Data <- read.table(file=File1, header=TRUE, as.is=TRUE, fill = TRUE)

#Just rename 1st column:
colnames(File1_Data)[1] <- "Gene.Symbol"
#Useful to get the file name
File1_Name <- gsub(".txt","", File1)
#View(File1_Data)
#Need a label with only:
#<Direction><R_Package><feature_counted><Count_Program>
#Need a series  of ifelse(grepl) commands:
Direction <- ifelse(grepl("Down",File1_Name),"Down","Up")
#print(Direction)
R_Package <- ifelse(grepl("DESeq",File1_Name),"DESeq","EdgeR")
#print(R_Package)
if(grepl("Exonic_Only", File1_Name, ignore.case = FALSE)) {
feature_counted <- "Exonic_Only"
}
#End of if statement
if(grepl("GeneBody", File1_Name, ignore.case = FALSE)) {
feature_counted <- "GeneBody"
}
#End of if statement
if(grepl("Intronic_Only", File1_Name, ignore.case = FALSE)) {
feature_counted <- "Intronic_Only"
}
#End of if statement
if(grepl("ExonCollapsed", File1_Name, ignore.case = FALSE)) {
  feature_counted <- "ExonCollapsed"
}
#End of if statement
#print(feature_counted)
Count_Program <- ifelse(grepl("featureCounts",File1_Name),"featureCounts","HTSeq")
#print(Count_Program)
#Once I get each part; paste them together
File1_Label <- paste(Direction, R_Package,feature_counted,  Count_Program, sep=".")
print(File1_Label)

#Read in File2
File2_Data <- read.table(file=File2, header=TRUE, as.is=TRUE, fill = TRUE)

# interrupt without error (do not draw venn diagram)
if (nrow(File2_Data) == 0 && nrow(File1_Data) == 0) {
       print(paste0("WARNING: Both input files is zero length! Do nothing!"))
       quit(save = "no", status = 0, runLast = F)
}

#Just rename 1st column:
colnames(File2_Data)[1] <- "Gene.Symbol"
#Useful to get the file name
File2_Name <- gsub(".txt","", File2)
#View(File2_Data)
#Need a label with only:
#<Direction><R_Package><feature_counted><Count_Program>
#Need a series  of ifelse(grepl) commands:
Direction <- ifelse(grepl("Down",File2_Name),"Down","Up")
#print(Direction)
R_Package <- ifelse(grepl("DESeq",File2_Name),"DESeq","EdgeR")
#print(R_Package)
if(grepl("Exonic_Only", File2_Name, ignore.case = FALSE)) {
feature_counted <- "Exonic_Only"
}
#End of if statement
if(grepl("GeneBody", File2_Name, ignore.case = FALSE)) {
feature_counted <- "GeneBody"
}
#End of if statement
if(grepl("Intronic_Only", File2_Name, ignore.case = FALSE)) {
feature_counted <- "Intronic_Only"
}
#End of if statement
if(grepl("ExonCollapsed", File2_Name, ignore.case = FALSE)) {
  feature_counted <- "ExonCollapsed"
}
#End of if statement
#print(feature_counted)
Count_Program <- ifelse(grepl("featureCounts",File2_Name),"featureCounts","HTSeq")
#print(Count_Program)
#Once I get each part; paste them together
File2_Label <- paste(Direction, R_Package,feature_counted,  Count_Program, sep=".")
print(File2_Label)
#---------------------------------------------------------------------------------------------
#Sort data_frames by gene symbol
File1_Data <- File1_Data[order(File1_Data$"Gene.Symbol"),]
File2_Data <- File2_Data[order(File2_Data$"Gene.Symbol"),]
#---------------------------------------------------------------------------------------------
#Now we can compare gene lists
#Use the VennDiagram package
#Creating a venn diagram
#---------------------------------------------------------------------------------------------
GeneSym_ListObj <- list(File1_Data$"Gene.Symbol", File2_Data$"Gene.Symbol")
#Need to color Venn diagrams based in direction of gene regulation
#(red=upregulation) or  (blue=downregulation)
#Need grepl to return true/false and need ignore.case = FALSE
if(grepl("Down", File1, ignore.case = FALSE)) {
File1_Color <- "cornflowerblue"
}
#End of if statement
if(grepl("Down", File2, ignore.case = FALSE)) {
File2_Color <- "cadetblue1"
}
#End of if statement
if(grepl("Up", File1, ignore.case = FALSE)) {
File1_Color <- "orangered"
}
#End of if statement
if(grepl("Up", File2, ignore.case = FALSE)) {
File2_Color <- "red"
}
#End of if statement
#Need to add names to the list object for the Venn diagram:
names(GeneSym_ListObj) = c(File1_Label, File2_Label)
#Create the Venn diagram:
Venn_Diagram <- venn.diagram(
x = GeneSym_ListObj,
filename=NULL,
lwd = 4,
fill = c(File1_Color, File2_Color),
alpha = 0.75,
label.col ="black",
cex = 2.0,
cat.default.pos = "text",
fontfamily = "serif",
fontface = "bold",
cat.col = c("black", "black"),
cat.cex = 1.0,
cat.fontfamily = "serif",
cat.fontface = "bold",
cat.dist = c(0.06, 0.06),
cat.pos = 0,
#Adjust so labels are visible:
#List (length = 1/2/3/4 based on set number) of Vectors of length 2 indicating horizontal and vertical justification for each category name
#cat.just=list(c(1,6), c(0.2,1) ),
cat.just=list(c(0.5,3), c(0.5,3) ),
main = "Gene Symbol Comparison",
main.cex = 2,
sub = paste("(", Subtitle, ")",sep=""),
sub.cex = 2
);
#Options explained:
#http://www.inside-r.org/packages/cran/VennDiagram/docs/venn.diagram
#Create the image file
png(file=paste("Venn_", File1_Label, ".", File2_Label, ".png", sep=""))
grid.draw(Venn_Diagram)
dev.off()
#---------------------------------------------------------------------------------------------
#Also want a table of gene counts for the common and unique gene symbols
Common_Count <- length(intersect(File1_Data$"Gene.Symbol", File2_Data$"Gene.Symbol"))
#material that is in the first named set, that is not in the second named set
File1_Unique_Count <- length(setdiff(File1_Data$"Gene.Symbol", File2_Data$"Gene.Symbol"))
File2_Unique_Count <- length(setdiff(File2_Data$"Gene.Symbol", File1_Data$"Gene.Symbol"))
Union_Count <- length(union(File1_Data$"Gene.Symbol", File2_Data$"Gene.Symbol"))
Count_Table <- NULL
sub = paste("(", Subtitle, ")",sep="")
#Need cbind to organize numbers into a table:
Count_Table <- cbind(sub, Common_Count, File1_Unique_Count, File2_Unique_Count, Union_Count)
#Add headers:
colnames(Count_Table) <- c(DiffExp_Index, "Common.Gene.Symbols", paste(File1_Label,".Unique", sep=""), paste(File2_Label, ".Unique", sep=""), "Union.Gene.Symbols")
#Write table to file:
write.table(Count_Table, file=paste("Count.Table.", File1_Label, ".", File2_Label, ".txt", sep=""), quote=F, sep="\t", col.names=TRUE, row.names=F)
#---------------------------------------------------------------------------------------------
print("Removing VennDiagram*.log files")
system("rm VennDiagram*.log")
#---------------------------------------------------------------------------------------------
paste("Check out Venn diagram and Count.Table!",sep="")
#---------------------------------------------------------------------------------------------
##############################################################################################
