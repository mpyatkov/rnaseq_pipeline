#############################
# Andy Rampersaud, 09.21.2015
#############################
#Pass arguments from sh into R
args <- commandArgs(trailingOnly = TRUE)
peak1_total <- as.numeric(args[1])
peak2_total <- as.numeric(args[2])
peak1_common <- as.numeric(args[3])
peak2_common <- as.numeric(args[4])
peak1_unique <- as.numeric(args[5])
peak2_unique <- as.numeric(args[6])
peak1_name <- args[7]
peak2_name <- args[8]
Union_name <- args[9]
Merged_Common_name <- args[10]
##############################################################################
#Function used within functions:
#Calculate peak widths:
#Make an R function to calculate the width of each peak and append to dataframe
peak_width <- function(dataframe){
DHS_width <- dataframe[,3] - dataframe[,2]
dataframe[,dim(dataframe)[2]+1] <- DHS_width
colnames(dataframe)[dim(dataframe)[2]] <- "Peak_Width"
return(dataframe)
}#end of peak_width function
##############################################################################
#Need to set the dir so that this script can work in any folder:
dir <- getwd()
#wdir <- setwd(dir)
#Use variable for Output dir:
Output_Dir <- paste(dir,"/",peak1_name, "_",peak2_name,"_Output", sep ="") 
wdir <- setwd(Output_Dir)
##############################################################################
#Need Peak_Width_Summary function for peak1 and peak2 sites
Peak_Width_Summary <- function(peak_set, peak_name){
Peak_Data <- read.table(paste(peak_set,".bed", sep =""), sep = "\t", header = FALSE)
# print("The dim of Peak_Data:")
# print(dim(Peak_Data))
#Call the function:
Peak_Data <- peak_width(Peak_Data)
# print("The dim of Peak_Data (after peak_width function):")
# print(dim(Peak_Data))
#Make Summary Table:
Summary_table <- rbind(summary(Peak_Data$"Peak_Width"))
#Need to parse Summary_table so that output text file is formatted correctly
Summary_table_out <- rbind(Summary_table[1:6])
peakSet <- c(paste(peak_name," (",dim(Peak_Data)[1],")",sep =""))
Summary_table_out2 <- cbind("Peak Set"=peakSet, Summary_table_out)
colnames(Summary_table_out2)[2:ncol(Summary_table_out2)] <- colnames(Summary_table)
#Need to set row.names = FALSE to avoid the "[1,]" in the second row
write.table(Summary_table_out2, file = (paste(Output_Dir,"/", peak_name,"_Stats.txt",sep="")), quote=FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
#print("Check out Stats.txt!")
}#End of Peak_Width_Summary function
##############################################################################
#Call Peak_Width_Summary function
Peak_Width_Summary("peak1", peak1_name)
Peak_Width_Summary("peak2", peak2_name)
Peak_Width_Summary("Union_Peak_Set", Union_name)
#Only want the Peak_Width_Summary("Merged_Common_Peak_Set")
#when counts of common sites is non-zero
if(peak1_common != 0 & peak2_common != 0){Peak_Width_Summary("Merged_Common_Peak_Set", Merged_Common_name)}#End of if statement
##############################################################################
#Need Peak_Width_Summary function for common and unique sites
Peak_Width_Summary <- function(peak_name, label){
Peak_Data <- read.table(paste(peak_name,label, ".bed", sep =""), sep = "\t", header = FALSE)
# print("The dim of Peak_Data:")
# print(dim(Peak_Data))
#Call the function:
Peak_Data <- peak_width(Peak_Data)
# print("The dim of Peak_Data (after peak_width function):")
# print(dim(Peak_Data))
#Make Summary Table:
Summary_table <- rbind(summary(Peak_Data$"Peak_Width"))
#Need to parse Summary_table so that output text file is formatted correctly
Summary_table_out <- rbind(Summary_table[1:6])
peakSet <- c(paste(peak_name,label," (",dim(Peak_Data)[1],")",sep =""))
Summary_table_out2 <- cbind("Peak Set"=peakSet, Summary_table_out)
colnames(Summary_table_out2)[2:ncol(Summary_table_out2)] <- colnames(Summary_table)
#Need to set row.names = FALSE to avoid the "[1,]" in the second row
write.table(Summary_table_out2, file = (paste(peak_name,label,"_Stats_2.txt",sep="")), quote=FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
#print("Check out Stats.txt!")
}#End of Peak_Width_Summary function
##############################################################################
#Call Peak_Width_Summary function:
if(peak1_common != 0){Peak_Width_Summary(peak1_name, "_common")}#End of if statement
if(peak1_unique != 0){Peak_Width_Summary(peak1_name, "_unique")}#End of if statement
if(peak2_common != 0){Peak_Width_Summary(peak2_name, "_common")}#End of if statement
if(peak2_unique != 0){Peak_Width_Summary(peak2_name, "_unique")}#End of if statement
#Find unique and common percents
peak1_percent_unique <- signif(peak1_unique/peak1_total*100, digits=3)
peak2_percent_unique <- signif(peak2_unique/peak2_total*100, digits=3)
peak1_percent_common <- signif(peak1_common/peak1_total*100, digits=3)
peak2_percent_common <- signif(peak2_common/peak2_total*100, digits=3)
#Print out Overlap_Summary.txt
output_table <- matrix(c(peak1_unique, peak1_percent_unique, peak1_common, peak1_percent_common, peak2_unique, peak2_percent_unique, peak2_common, peak2_percent_common), ncol=4,byrow=TRUE)
colnames(output_table) <- c("Unique_Sites","Unique_Sites(%)","Common_Sites","Common_Sites(%)")
rownames(output_table) <- c(paste(peak1_name," (",peak1_total,")", sep =""),paste(peak2_name," (",peak2_total,")", sep =""))
output_table <- as.table(output_table)
write.table(output_table, file = paste(peak1_name,"_", peak2_name, "_Overlap_Summary.txt", sep =""), quote=FALSE, sep = "\t",col.names=NA)
#Peak Proximity:
##############################################################################
#Write a function to get Peak Proximity statistics
Peak_Proximity <- function(peak_name, peak_total){
Peak_Dis <- read.table(paste(peak_name, "_distance.txt", sep =""), sep = "\t", header = FALSE)
colnames(Peak_Dis) <- c("Distance")
Peak_Dis=as.data.frame(Peak_Dis)
#Get the counts:
Bin_10bp <- length(Peak_Dis[Peak_Dis >= 0 & Peak_Dis < 10])
Bin_100bp <- length(Peak_Dis[Peak_Dis >= 10 & Peak_Dis < 100])
Bin_1KB <- length(Peak_Dis[Peak_Dis >= 100 & Peak_Dis < 1000])
Bin_10KB <- length(Peak_Dis[Peak_Dis >= 1000 & Peak_Dis < 10000])
Bin_100KB <- length(Peak_Dis[Peak_Dis >= 10000 & Peak_Dis < 100000])
Bin_1MB <- length(Peak_Dis[Peak_Dis >= 100000 & Peak_Dis < 1000000])
Bin_gt1MB <- length(Peak_Dis[Peak_Dis >= 1000000])
#Get the percents:
Bin_10bp_percent <- signif(Bin_10bp/peak_total*100, digits=3)
Bin_100bp_percent <- signif(Bin_100bp/peak_total*100, digits=3)
Bin_1KB_percent <- signif(Bin_1KB/peak_total*100, digits=3)
Bin_10KB_percent <- signif(Bin_10KB/peak_total*100, digits=3)
Bin_100KB_percent <- signif(Bin_100KB/peak_total*100, digits=3)
Bin_1MB_percent <- signif(Bin_1MB/peak_total*100, digits=3)
Bin_gt1MB_percent <- signif(Bin_gt1MB/peak_total*100, digits=3)
#Create output:
output_table <- matrix(c(Bin_10bp, Bin_100bp, Bin_1KB, Bin_10KB, Bin_100KB, Bin_1MB, Bin_gt1MB, Bin_10bp_percent, Bin_100bp_percent, Bin_1KB_percent, Bin_10KB_percent, Bin_100KB_percent, Bin_1MB_percent,Bin_gt1MB_percent), ncol=7,byrow=TRUE)
colnames(output_table) <- c("Bin 0 to 10bp","Bin 10 to 100bp","Bin 100 to 1KB","Bin 1KB to 10KB", "Bin 10KB to 100KB", "Bin 100KB to 1MB", "Bin gt 1MB")
rownames(output_table) <- c(paste(peak_name," (",peak_total,") counts", sep =""),paste(peak_name," (",peak_total,") percentages", sep =""))
output_table <- as.table(output_table)
write.table(output_table, file = paste(peak_name,"_","distance_Summary.txt", sep =""), quote=FALSE, sep = "\t",col.names=NA)
}#End of Peak_Proximity function
#Call Peak_Proximity function (no need to assign function output to an object)
Peak_Proximity(peak1_name, peak1_total)
Peak_Proximity(peak2_name, peak2_total)
##############################################################################
