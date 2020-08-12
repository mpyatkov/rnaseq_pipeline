#Tisha Melia
#Dec 19, 2013: Initial creation of the file
#Jan 14, 2014: Adding pseudocount in fold change calculation
#Jan 15, 2014: Making pseudocount as a command line argument
#May 8 2015: Select for tpm instead of base means.
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
	
	#Read in data:
	data <- read.table(file=inputFile, sep="\t", header=TRUE, as.is=TRUE)
	deseq_ratio <- (data[,grep("baseMean_", colnames(data))[2]] + epsilon) / 
                       (data[,grep("baseMean_", colnames(data))[1]] + epsilon)
	deseq_foldChange <- data$DESeq_foldChange
	flip.idx <- deseq_foldChange < 1
	deseq_foldChange[flip.idx] <- (1/deseq_foldChange[flip.idx]) * -1
	
	edger_ratio <- 2^data$edgeR_logFC
	edger_foldChange <- edger_ratio
	flip.idx <- edger_foldChange < 1 
	edger_foldChange[flip.idx] <- (1/edger_foldChange[flip.idx]) * -1
	output <- cbind("id"=data$id,
		  "deseq_ratio"=deseq_ratio,
		  "edger_ratio"=edger_ratio,
		  "deseq_foldChange"=deseq_foldChange,
		  "edger_foldChange"=edger_foldChange,
		  data[,grep("tpm_mean", colnames(data))],
		  "deseq_padj"=data$DESeq_padj,
		  "edger_padj"=data$edgeR_FDR,
		  "deseq_pval"=data$DESeq_pvalue,
		  "edger_pval"=data$edgeR_PValue)
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
	deseqOutput <- output[,grep("id|deseq|tpm", names(output))]
	names(deseqOutput) <- gsub("deseq_", "", names(deseqOutput))
	
	edgerOutput <- output[,grep("id|edger|tpm", names(output))]
	names(edgerOutput) <- gsub("edger_", "", names(edgerOutput))
	
	write.table(format(deseqOutput, digits=8, scientific=FALSE), file=paste(outputFile, "_TPM_DESeq.txt", sep=""), 
	            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
	
	write.table(format(edgerOutput, digits=8, scientific=FALSE), file=paste(outputFile, "_TPM_EdgeR.txt", sep=""), 
	            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
#End of else statement              
}
#-----------------------------------------------------------------------------------
####################################################################################
