#----------------------------------------------------------------------------------
#Tisha Melia
#Dec 2012
#April 18, 2013: Tisha added normalized counts to the output
#March 20, 2015: Tisha made some edits to allow for "Intron_Only" features (present in the Intron_Only_Regions.gtf created by Andy Rampersaud)
#April 06, 2015: Andy changed "Intron_Only" to "Intronic_Only" in the Intron_Only_Regions.gtf file, this script was updated accordingly
#April 28, 2015: Andy added system argument "OUTPUT_PREFIX" to label different types of counting:
#May 8, 2015: Tisha added FPKM calculation
#June 30, 2015: Tisha added EdgeR differential expression
#October 15, 2015: Tisha updated DESeq to DESeq2 and capibilities to handle unreplicated samples
#May 25, 2016: Tisha modified the script to adapt for lncRNA differential expression
#----------------------------------------------------------------------------------
rm(list=ls())
library(DESeq2)
library(edgeR)
library(GenomicRanges)
library(readr)
library(dplyr)

#######################################################################
#1. SETTING VARS
#######################################################################


args <- commandArgs(T)
if(length(args) < 8){
	print("diff_analysis failed: not enough inputs")
	quit()
} 
 
CONDITION_1=args[1]
CONDITION_2=args[2]
NUM_REP_CONDITION1=as.numeric(args[3])
NUM_REP_CONDITION2=as.numeric(args[4])
ANNOTATION_FILE=args[5]
READ_COUNT_FOLDER=args[6]
LIB_TYPE="not-used"
OUTPUT_PREFIX = args[7]
GENE_LENGTHS_FILE = args[8]

#######################################################################
#2. READING INPUT FILES
#######################################################################
#read in gene lengts file
geneLength <- read.table(file=paste(READ_COUNT_FOLDER, "/", GENE_LENGTHS_FILE, sep=""), 
                         as.is=TRUE, header=TRUE, colClasses=c("character", "numeric"))
names(geneLength) <- c("GeneSymbol", "length")
###############

curpath<-getwd()

setwd(READ_COUNT_FOLDER)

# get files with counts ex. Male_G186_M1
f_cond_files <- function(p){
  # r <- list.files(pattern = paste0(p,"_[[:alnum:]*][_{1}]"))
  r <- list.files(pattern = paste0("^",p,"_[[:alnum:]]*_[[:alnum:]]*$"))
  r <- r[!grepl("summary", r)]
  r
}

# short version because we parse summary files in DiffExp.qsub
readMapping <- function(fname){
    as.numeric(read.table(file=fname, colClasses="numeric"))
}

readSampleCounts <- function(cond_name) {
  input_dir<-getwd()
  setwd(cond_name)

  cond_files<-f_cond_files(cond_name)
  # sum_files<-f_sum_files(cond_name)
  sum_files<-paste0(cond_files,"_summary")
  
  print(cond_files)
  print(sum_files)

  # read sample counts
  df <-  cbind.data.frame(lapply(cond_files, read_tsv, col_names = F))
  gene.names <- df[,1]

  # only columns with counts (even columns)
  df <- df[,!(seq_len(ncol(df)) %% 2)]

  colnames(df) <- cond_files
  rownames(df) <- gene.names

  # read summary files
  sum_output <- sapply(sum_files, readMapping)
  sum_output <- as.data.frame(t(sum_output))
  colnames(sum_output) <- cond_files
  print(sum_output)

  setwd(input_dir)
  list(data=df, gene.names = gene.names, mapping=sum_output)
}

cond1 <- readSampleCounts(CONDITION_1)
condition1Count <- cond1$data
numMappedCondition1List <- cond1$mapping

cond2 <- readSampleCounts(CONDITION_2)
condition2Count <- cond2$data
numMappedCondition2List <- cond2$mapping

countTable <- cbind.data.frame(condition1Count, condition2Count, stringsAsFactors=FALSE)
rownames(countTable) <- rownames(condition1Count)

setwd(curpath)

#######################################################################
#3. FORMAT THE INPUT FILES
#######################################################################
#create the meta data
designMat = data.frame(
    row.names = colnames( countTable ),
    condition = factor(c(rep(CONDITION_1, NUM_REP_CONDITION1), rep(CONDITION_2, NUM_REP_CONDITION2)),
                          levels=c(CONDITION_1, CONDITION_2)),
    libType = rep(LIB_TYPE, ncol(countTable ) ) )


###############################################################
#4. DIFFERENTIAL EXPRESSION CALCULATION
###############################################################
deseqMat <- DESeqDataSetFromMatrix(countData = countTable,
                                   colData = designMat,
                                   design = ~ condition)
suppressWarnings(deseqOutput <- DESeq(deseqMat))
baseMeanPerCondition <- sapply( levels(designMat$condition), 
                                function(cond){
                                    currData <- counts(deseqOutput,normalized=TRUE)[,deseqOutput$condition == cond]
                                    if(is.null(ncol(currData))){
                                        meansPerRow=currData
                                    }else {
                                        meansPerRow=rowMeans(currData)
                                    }
                                    return(meansPerRow)
                                })
colnames(baseMeanPerCondition) <- paste("baseMean", colnames(baseMeanPerCondition), sep="_")

temp <- as.data.frame(results(deseqOutput))
#change fold change = NA to 0, as how edgeR handled it
temp$log2FoldChange[is.na(temp$log2FoldChange)] <- 0
result <- cbind(id= rownames(baseMeanPerCondition), 
                baseMean= temp$baseMean, 
                baseMeanPerCondition,
                foldChange= 2^temp$log2FoldChange,
                temp[,c("log2FoldChange", "pvalue", "padj")])
#DJW request padj would be renamed to padj_FDR
names(result)[ncol(result)] <- "padj_FDR"
names(result) <- paste(c(rep("",4), rep("DESeq2_",4)), names(result), sep="")


dge <- DGEList(counts=countTable, group=factor(designMat$condition)) 
dge <- calcNormFactors(dge)
if(NUM_REP_CONDITION1 == 1 & NUM_REP_CONDITION2 == 1){
    dge <- estimateGLMCommonDisp(dge, method="deviance", robust="TRUE", subset=NULL) 
    resultEdger=exactTest(dge)
} else {
    dge <- estimateCommonDisp(dge) 
    dge <- estimateTagwiseDisp(dge) 
    resultEdger <- exactTest(dge, pair=c(CONDITION_1, CONDITION_2)) 
}
sortedResultEdger <- topTags(resultEdger, n=nrow(resultEdger))
mergedResult <- cbind("id"=rownames(sortedResultEdger), 
                      foldChange= 2^sortedResultEdger$table$logFC,
                      sortedResultEdger$table[,c("logFC", "PValue", "FDR")])
#rename column to mimic DESeq2 column names
names(mergedResult) <- c("id", "foldChange", "log2FoldChange", "pvalue", "padj_FDR")
names(mergedResult) <- paste("EdgeR", names(mergedResult), sep="_")
mergedResult <- merge(result, mergedResult, by.x="id", by.y="EdgeR_id", sort=FALSE)

mergedResult <- cbind("edgeR_logFC_copy"=mergedResult$EdgeR_log2FoldChange, mergedResult)
mergedResult <- cbind("edgeR_FDR_copy"=mergedResult$EdgeR_padj_FDR, mergedResult)


new_edgeRFC <-(paste(CONDITION_1,"_", CONDITION_2,"_edgeRlogFC",sep=""))  # vairable with prefix for logFC
new_edgeRFDR <-(paste(CONDITION_1,"_", CONDITION_2,"_edgeRFDR",sep=""))    # variable with prefix for FDR

colnames(mergedResult)[colnames(mergedResult) == "edgeR_logFC_copy"] <- new_edgeRFC
colnames(mergedResult)[colnames(mergedResult) == "edgeR_FDR_copy"] <- new_edgeRFDR

#### RPKM and TPM old

temp<-as.data.frame(rownames(condition1Count))
names(temp) <- "GeneSymbol"
# geneLength <- merge(temp, geneLength, by.x="GeneSymbol", by.y="GeneSymbol", sort=FALSE)[]
geneLength <- left_join(temp, geneLength, by=c("GeneSymbol"))

rpkm <- function(x, totalLength, mappedReadsNum) {
  return(x/(unlist(mappedReadsNum) * 1e-09 * totalLength))
}

rpkmCondition1 <- sapply( seq(1, ncol(condition1Count)) , function(i){
  return( rpkm(condition1Count[,i], geneLength$length, numMappedCondition1List[i]) )
})

rpkmCondition1 <- cbind(rpkmCondition1, apply(rpkmCondition1, 1, mean))
colnames(rpkmCondition1) <- c(paste("rpkm", colnames(condition1Count), sep="_"),
                              paste("rpkm_mean",CONDITION_1, sep="_"))

sum1 <- sapply(seq(1, ncol(rpkmCondition1)), function(i){ return(sum(rpkmCondition1[,i]))})

tpmCondition1 <- sapply( seq(1, ncol(condition1Count)) , function(i){
  return(((rpkm(condition1Count[,i], geneLength$length, numMappedCondition1List[i]))/sum1[i])*10^6)
})

tpmCondition1 <- cbind(tpmCondition1, apply(tpmCondition1, 1, mean))
colnames(tpmCondition1) <- c(paste("tpm", colnames(condition1Count), sep="_"),
                             paste("tpm_mean",CONDITION_1, sep="_"))

rpkmCondition2 <- sapply( seq(1, ncol(condition2Count)) , function(i){
  return( rpkm(condition2Count[,i], geneLength$length, numMappedCondition2List[i]) )
})

rpkmCondition2 <- cbind(rpkmCondition2, apply(rpkmCondition2, 1, mean))
colnames(rpkmCondition2) <- c(paste("rpkm", colnames(condition2Count), sep="_"),
                              paste("rpkm_mean",CONDITION_2, sep="_"))

sum2 <- sapply(seq(1, ncol(rpkmCondition2)), function(i){ return(sum(rpkmCondition2[,i]))})

tpmCondition2 <- sapply( seq(1, ncol(condition2Count)) , function(i){
  return(((rpkm(condition2Count[,i], geneLength$length, numMappedCondition2List[i]))/sum2[i])*10^6)
})

tpmCondition2 <- cbind(tpmCondition2, apply(tpmCondition2, 1, mean))
colnames(tpmCondition2) <- c(paste("tpm", colnames(condition2Count), sep="_"),
                              paste("tpm_mean",CONDITION_2, sep="_"))

########### NEW 
# #resort geneLength 
# geneLength <- merge(temp, geneLength, by.x="GeneSymbol", by.y="GeneSymbol", sort=FALSE)[]
# # sum(geneLength$GeneSymbol != rownames(condition1Count)) #should be 0

# edger_rpkm <- edgeR::rpkm(dge$counts, geneLength$length)
# edger_tpm <- t(t(edger_rpkm)/colSums(edger_rpkm))*(10^6)


# add_mean_and_rename <- function(m, names, prefix, cond_name) {
#   m1 <- m[, names]
#   m1 <- cbind(m1, apply(m1,1,mean))
#   colnames(m1) <- c(paste(prefix, names, sep="_"), paste(prefix,"mean", cond_name, sep="_"))
#   m1
# }

# edger_rpkm1 <- add_mean_and_rename(edger_rpkm, colnames(condition1Count), "rpkm", CONDITION_1)
# edger_rpkm2 <- add_mean_and_rename(edger_rpkm, colnames(condition2Count), "rpkm", CONDITION_2)
# edger_tpm1 <- add_mean_and_rename(edger_tpm, colnames(condition1Count), "tpm", CONDITION_1)
# edger_tpm2 <- add_mean_and_rename(edger_tpm, colnames(condition2Count), "tpm", CONDITION_2)

# #make sure the gene length is ordered the same way as the count
# temp <- as.data.frame(rownames(condition1Count))
# names(temp) <- "GeneSymbol"

###############################################################
#4. PREPPING OUTPUT FILES
###############################################################
#get the count, rpkm and diff. exp result
mergedResult <- as.data.frame(mergedResult)
countNorm <- counts(deseqOutput, normalized=TRUE )
colnames(countNorm) <- paste("norm_", colnames(countNorm), sep="")
# for edger output rpkm, tpm
# mergedResult <- cbind(countTable, edger_rpkm1, edger_tpm1, edger_rpkm2, edger_tpm2 , countNorm, mergedResult) #row order is OK
# for old version

mergedResult <- cbind(countTable, rpkmCondition1, tpmCondition1, rpkmCondition2, tpmCondition2 ,countNorm, mergedResult)
mergedResult <- subset(mergedResult, select=-c(id))

# rownames to columns
mergedResult <- cbind(id = rownames(mergedResult), mergedResult)
rownames(mergedResult) <- NULL

#mergedResult <- mergedResult[,c("id", names(mergedResult)[-match("id", names(mergedResult))])]

# TODO: fix this to some more general: remove all XLOC for LncRNA output
library(dplyr)
mergedResult<-mergedResult%>%filter(!grepl("XLOC", id))

# drop all rows which contains NA ()
# mergedResult <- na.omit(mergedResult)

#Print output file location:
#print(paste("output file is in: ", READ_COUNT_FOLDER, "/",OUTPUT_PREFIX,"_",CONDITION_1,"_", CONDITION_2, ".txt", sep=""))
#Create the output file:
write.table(mergedResult,file=paste(READ_COUNT_FOLDER, "/",OUTPUT_PREFIX,"_", CONDITION_1, "_", CONDITION_2, ".txt", sep=""), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
#----------------------------------------------------------------------------------
