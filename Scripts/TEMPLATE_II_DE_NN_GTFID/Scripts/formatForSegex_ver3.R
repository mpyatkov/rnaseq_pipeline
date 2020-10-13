#Tisha Melia
#Dec 19, 2013: Initial creation of the file
#Jan 14, 2014: Adding pseudocount in fold change calculation
#Jan 15, 2014: Making pseudocount as a command line argument
#May 8 2015: Select for rpkm instead of base means.
#Format output files from the diff exp pipeline for segex
#-----------------------------------------------------------------------------------
#Here's how you use the script

#a) connect to scc4.bu.edu
#b) transfer formatForSegex_ver3.R and the output file of your differential expression (from now on, this will be called the input file)
#c) format files by running the script:
#The script receives three input variables: your input file, output file names, and your pseudocount:

#run the script by typing the following
#Rscript formatForSegex_ver3.R <your input file name here> <your output file name here> <pseudocount>

#Example:
#Rscript formatForSegex_ver2.R input.txt output.txt 1
#In this case, 1)input.txt is the input file name, 2) output.txt is the output file name you want, and 3)pseudocount is set to 1
#d) Check your output file for correctness
#-----------------------------------------------------------------------------------
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 4){
	print("ERROR: you need to supply the 1)input file and 2)output file names 3)pseudocount 4)col_suffix")
	print("Rscript formatForSegex.R <input file> <output file> <pseudocount> <col_suffix>")
} else {
	inputFile <- args[1]
	outputFile <- args[2]
	epsilon <- as.numeric(args[3])
	col_suffix <- args[4]
	
# 	inputFile <- "/dos/waxman-lab/routines/Lab_RNA-Seq_pipeline/temp_count/DiffExp_v2_LncRNA_ExonCollapsed_Female_Male.txt"
# 	outputFile <- "Male_G83_M1M2_vs_Female_G83_M3M4_ExonCollapsed_forSEGEXUpload"
# 	epsilon <- 1
# 	col_suffix <- "LncRNA_ExonCollapsed"
	
	
# 	inputFile <- "/dos/waxman-lab/routines/Lab_RNA-Seq_pipeline/temp_count2/DiffExp_v2_LncRNA_Intronic_Only_Female_Male.txt"
# 	outputFile <- "Male_G83_M1M2_vs_Female_G83_M3M4_Intronic_Only_forSEGEXUpload"
# 	epsilon <- 1
# 	col_suffix <- "LncRNA_Intronic_Only"
	
	#Read in data:
	data <- read.table(file=inputFile, sep="\t", header=TRUE, as.is=TRUE)
	deseq_ratio <- (data[,grep("baseMean_", colnames(data))[2]] + epsilon) / 
                       (data[,grep("baseMean_", colnames(data))[1]] + epsilon)
	deseq_ratio[is.na(deseq_ratio)] <- 1 #genes that are not counted would have ratio as na
	deseq_foldChange <- data$DESeq2_foldChange
	deseq_foldChange[is.na(deseq_foldChange)] <- 1 #genes that are not counted would have fold change as na
	flip.idx <- deseq_foldChange < 1
	deseq_foldChange[flip.idx] <- (1/deseq_foldChange[flip.idx]) * -1
	
	edger_ratio <- 2^data$EdgeR_log2FoldChange
	edger_ratio[is.na(edger_ratio)] <- 1 #genes that are not counted would have ratio as na
	edger_foldChange <- edger_ratio
	edger_foldChange[is.na(edger_foldChange)] <- 1 #genes that are not counted would have fold change as na
	flip.idx <- edger_foldChange < 1 
	edger_foldChange[flip.idx] <- (1/edger_foldChange[flip.idx]) * -1
	output <- cbind("id"=gsub("_chr", "_c", gsub("^ncRNA", "nc", data$id)), #had to shorten the gene name b/c segex limits gene name to 18 chrs
		  "deseq_ratio"=deseq_ratio,
		  "edger_ratio"=edger_ratio,
		  "deseq_foldChange"=deseq_foldChange,
		  "edger_foldChange"=edger_foldChange,
		  data[,grep("rpkm_mean", colnames(data))],
		  "deseq_padj"=data$DESeq2_padj_FDR,
		  "edger_padj"=data$EdgeR_padj_FDR,
		  "deseq_pval"=data$DESeq2_pvalue,
		  "edger_pval"=data$EdgeR_pvalue)
	output <- as.data.frame(output, stringsAsFactors=FALSE)
	for(i in seq(8,ncol(output))){
	output[,i] <- as.numeric(output[,i] )
	change.idx <- is.na(output[,i]) | is.infinite(output[,i])
	output[change.idx,i] <- 1
	}
	output$id <- as.character(output$id)
	output$id[output$id == "Gt(ROSA)26Sor"] <- "Gt_ROSA_26Sor"	
	
	#Need to replace "NA" with zero:
	output[is.na(output)] <- 0
	##Need to sort by gene symbol: already done in the differentialAnalysisDESeq.R
	#output <- output[order(output$id),] 
	#Need to add a suffix to the column names:
	names(output) <- paste(names(output), col_suffix, sep = ".")
	#Write the output file:
	deseqOutput <- output[,grep("id|deseq|rpkm", names(output))]
	names(deseqOutput) <- gsub("deseq_", "", names(deseqOutput))
	
	edgerOutput <- output[,grep("id|edger|rpkm", names(output))]
	names(edgerOutput) <- gsub("edger_", "", names(edgerOutput))
	
	write.table(format(deseqOutput, digits=8, scientific=FALSE), file=paste(outputFile, "_DESeq.txt", sep=""), 
	            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
	
	write.table(format(edgerOutput, digits=8, scientific=FALSE), file=paste(outputFile, "_EdgeR.txt", sep=""), 
	            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
#End of else statement              
}
#-----------------------------------------------------------------------------------
####################################################################################
