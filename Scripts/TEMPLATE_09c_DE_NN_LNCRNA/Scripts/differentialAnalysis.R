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

#######################################################################
#0. FUNCTIONS
#######################################################################


rpkm <- function(x, totalLength, mappedReadsNum) {
  return(x/(unlist(mappedReadsNum) * 1e-09 * totalLength))
}

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

# CONDITION_1="Female"
# CONDITION_2="Male"
# NUM_REP_CONDITION1=2
# NUM_REP_CONDITION2=2
# ANNOTATION_FILE="/dos/waxman-lab/ncRNA/analysis/output/ncRNA_output_filtered_final_gene.txt"
# # READ_COUNT_FOLDER="/dos/waxman-lab/routines/Lab_RNA-Seq_pipeline/temp_count"
# LIB_TYPE="not-used"
# # OUTPUT_PREFIX = "DiffExp_v2_LncRNA_ExonCollapsed"
# # GENE_LENGTHS_FILE = "GTF_Files/lengths/ncRNA_exon_for_counting_lengths.txt"
# 
# READ_COUNT_FOLDER="/dos/waxman-lab/routines/Lab_RNA-Seq_pipeline/temp_count2"
# OUTPUT_PREFIX = "DiffExp_v2_LncRNA_Intronic_Only"
# GENE_LENGTHS_FILE = "GTF_Files/lengths/intronic_only_gene_models_ncRNA_for_counting_lengths.txt"


print("Arguments for differentialAnalysisDESeq.R:")
print(CONDITION_1)
print(CONDITION_2)
print(NUM_REP_CONDITION1)
print(NUM_REP_CONDITION2)
print(ANNOTATION_FILE)
print(READ_COUNT_FOLDER)
print(OUTPUT_PREFIX)
print(GENE_LENGTHS_FILE)


#######################################################################
#2. READING INPUT FILES
#######################################################################
#read in gene lengts file
geneLength <- read.table(file=paste(READ_COUNT_FOLDER, "/", GENE_LENGTHS_FILE, sep=""), 
                         as.is=TRUE, header=TRUE, colClasses=c("character", "numeric"))
names(geneLength) <- c("GeneSymbol", "length")


numMappedCondition1List<- NULL
#loading input files
condition1Count <- NULL
for(i in 1:NUM_REP_CONDITION1){
	temp <- read.delim2(file=paste(READ_COUNT_FOLDER, "/", CONDITION_1, (i-1), ".out", sep=""), header=FALSE,
                        as.is=TRUE, colClasses=c("character", "numeric"), col.names=c("gene", "count"))

  numFragmentFile <- read.table(file=paste(READ_COUNT_FOLDER, "/",CONDITION_1, (i-1), ".summary", sep=""), header=TRUE, sep="\t", colClasses=c("character", "numeric"))
  #number of fragment is in the second column, It's unfortunately hard coded b/c the name of the column changes according to the bam file name
  numFragment <- sum(numFragmentFile[numFragmentFile$Status %in% c("Assigned", "Unassigned_NoFeatures"), 2])

  if(i == 1){
		condition1Count <- temp$count
		numMappedCondition1List <- numFragment
	}else {
		condition1Count <- cbind(condition1Count, temp$count)
		numMappedCondition1List <- cbind(numMappedCondition1List, numFragment)
	}
	
}
#Need to convert to data frame to add rownames:
#(Instances when there's only 1 replicate)
condition1Count <- data.frame(condition1Count)
rownames(condition1Count) <- temp$gene 
colnames(condition1Count) <- paste(rep(CONDITION_1, NUM_REP_CONDITION1),seq(1,NUM_REP_CONDITION1), sep="")
colnames(numMappedCondition1List) <- paste(rep(CONDITION_1, NUM_REP_CONDITION1),seq(1,NUM_REP_CONDITION1), sep="")


numMappedCondition2List<- NULL
condition2Count <- NULL
for(i in 1:NUM_REP_CONDITION2){
	temp <- read.delim2(file=paste(READ_COUNT_FOLDER, "/", CONDITION_2, (i-1), ".out", sep=""), header=FALSE,
                        as.is=TRUE, colClasses=c("character", "numeric"), col.names=c("gene", "count"))

	numFragmentFile <- read.table(file=paste(READ_COUNT_FOLDER, "/",CONDITION_2, (i-1), ".summary", sep=""), header=TRUE, sep="\t", colClasses=c("character", "numeric"))
	#number of fragment is in the second column, It's unfortunately hard coded b/c the name of the column changes according to the bam file name
	numFragment <- sum(numFragmentFile[numFragmentFile$Status %in% c("Assigned", "Unassigned_NoFeatures"), 2])
	
	if(i == 1){
		condition2Count <- temp$count
		numMappedCondition2List <- numFragment
	}else {
		condition2Count <- cbind(condition2Count, temp$count)
		numMappedCondition2List <- cbind(numMappedCondition2List, numFragment)
	}
}
#Need to convert to data frame to add rownames:
#(Instances when there's only 1 replicate)
condition2Count <- data.frame(condition2Count)
rownames(condition2Count) <- temp$gene 
colnames(condition2Count) <- paste(rep(CONDITION_2, NUM_REP_CONDITION2),seq(1,NUM_REP_CONDITION2), sep="")
colnames(numMappedCondition2List) <- paste(rep(CONDITION_2, NUM_REP_CONDITION2),seq(1,NUM_REP_CONDITION2), sep="")


countTable <- cbind(condition1Count, condition2Count)
countTable <- as.data.frame(countTable, stringsAsFactors=FALSE)
rownames(countTable) <- rownames(condition1Count)

#make sure the gene length is ordered the same way as the count
temp <- as.data.frame(rownames(condition1Count))
names(temp) <- "GeneSymbol"
#resort geneLength 
geneLength <- merge(temp, geneLength, by.x="GeneSymbol", by.y="GeneSymbol", sort=FALSE)[]
# sum(geneLength$GeneSymbol != rownames(condition1Count)) #should be 0
rpkmCondition1 <- sapply( seq(1, ncol(condition1Count)) , function(i){
  return( rpkm(condition1Count[,i], geneLength$length, numMappedCondition1List[i]) )
})
rpkmCondition1 <- cbind(rpkmCondition1, apply(rpkmCondition1, 1, mean))
colnames(rpkmCondition1) <- c(paste("rpkm", colnames(condition1Count), sep="_"), 
                              paste("rpkm_mean",CONDITION_1, sep="_"))


sum1 <- sapply(seq(1, ncol(rpkmCondition1)), function(i){ return(sum(rpkmCondition1[,i]))})
print(paste("sum1:",sum1,sep=""))
tpmCondition1 <- sapply( seq(1, ncol(condition1Count)) , function(i){
  return(((rpkm(condition1Count[,i], geneLength$length, numMappedCondition1List[i]))/sum1[i])*10^6)
})



rpkmCondition2 <- sapply( seq(1, ncol(condition2Count)) , function(i){
  return( rpkm(condition2Count[,i], geneLength$length, numMappedCondition2List[i]) )
})
rpkmCondition2 <- cbind(rpkmCondition2, apply(rpkmCondition2, 1, mean))
colnames(rpkmCondition2) <- c(paste("rpkm", colnames(condition2Count), sep="_"),
                              paste("rpkm_mean",CONDITION_2, sep="_"))


sum2 <- sapply(seq(1, ncol(rpkmCondition2)), function(i){ return(sum(rpkmCondition2[,i]))})
print(paste("sum2:",sum2,sep=""))
tpmCondition2 <- sapply( seq(1, ncol(condition2Count)) , function(i){
  return(((rpkm(condition2Count[,i], geneLength$length, numMappedCondition2List[i]))/sum2[i])*10^6)
})

tpmCondition2 <- cbind(tpmCondition2, apply(tpmCondition2, 1, mean))
colnames(tpmCondition2) <- c(paste("tpm", colnames(condition2Count), sep="_"),
                              paste("tpm_mean",CONDITION_2, sep="_"))



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




###############################################################
#4. PREPPING OUTPUT FILES
###############################################################
#get the count, rpkm and diff. exp result
mergedResult <- as.data.frame(mergedResult)
countNorm <- counts(deseqOutput, normalized=TRUE )
colnames(countNorm) <- paste("norm_", colnames(countNorm), sep="")
mergedResult <- cbind(countTable, rpkmCondition1,tpmCondition1, rpkmCondition2, tpmCondition2 ,countNorm, mergedResult) #row order is OK

#read in lncRNA annotation
lncRNAAnnotation <- read.table(ANNOTATION_FILE, sep="\t", as.is=TRUE, header=TRUE);dim(lncRNAAnnotation)#15558    23
#add lncRNA annotation to the mergedResult
mergedResultLncRNA <- merge(mergedResult, lncRNAAnnotation, by.x="id", by.y="ncRNAId", all.y=TRUE, sort=FALSE);dim(mergedResultLncRNA)#15558    47
#subset columns to the columsn that I need and reorder them
mergedResultLncRNA <- mergedResultLncRNA[,c("id", "hits_oldLincs_accession", "type", "chr", "start", "end", "putativeStrand", "numExon",
                                            "hits_noncode_accession", "hits_smallRNA_accession", "hits_pseudogene_accession", "longestOrfLength", 
                                            #exclude column named id
                                            names(mergedResult)[-match("id", names(mergedResult))])]





#sort result
mergedResultLncRNA <- mergedResultLncRNA[order(mergedResultLncRNA$chr, mergedResultLncRNA$start),]
lncRNAChr1toChr9 <- mergedResultLncRNA[mergedResultLncRNA$chr %in% paste("chr", seq(1,9), sep=""),]
lncRNAChr10toChr19 <- mergedResultLncRNA[mergedResultLncRNA$chr %in% paste("chr", seq(10,19), sep=""),]
lncRNAChrXtoChrY <- mergedResultLncRNA[mergedResultLncRNA$chr %in% c("chrX", "chrY"),]
unique(lncRNAChr1toChr9$chr) #[1] "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9"
unique(lncRNAChr10toChr19$chr) # [1] "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19"
unique(lncRNAChrXtoChrY$chr) #[1] "chrX" "chrY"
lncRNA <- rbind(lncRNAChr1toChr9, lncRNAChr10toChr19, lncRNAChrXtoChrY)

# sum( abs(mergedResultLncRNA$DESeq_log2FoldChange)>1 & mergedResultLncRNA$DESeq_padj < 0.05, na.rm=TRUE)#0
# sum( abs(mergedResultLncRNA$edgeR_logFC)>=1 & mergedResultLncRNA$edgeR_FDR <= 0.05)#125
# fivenum(mergedResultLncRNA$edgeR_logFC)
# fivenum(as.data.frame(results(deseqOutput))$log2FoldChange)

#Print output file location:
print(paste("output file is in: ", READ_COUNT_FOLDER, "/",OUTPUT_PREFIX,"_",CONDITION_1,"_", CONDITION_2, ".txt", sep=""))
#Create the output file:
write.table(mergedResultLncRNA,file=paste(READ_COUNT_FOLDER, "/",OUTPUT_PREFIX,"_", CONDITION_1, "_", CONDITION_2, ".txt", sep=""), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
#----------------------------------------------------------------------------------
