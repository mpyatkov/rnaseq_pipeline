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
#----------------------------------------------------------------------------------
rm(list=ls())
library(DESeq2)
library(edgeR)
library(GenomicRanges)

#######################################################################
#0. FUNCTIONS
#######################################################################

read.GTF <- function(file, feature=c("exon", "CDS", "intron", "utr"), 
                     attributes=c("gene_id", "transcript_id", "gene_name", 
                                  "exon_number", "tss_id")){

  cat("load GTF file ... \n")
  GTF <- read.table(file, sep="\t", header=FALSE, stringsAsFactors=FALSE, 
                    colClasses=c(rep("character", 3), rep("numeric", 2), rep("character", 4)))
  cat("parse attributes ... \n")
  attrs <- GTF[,9]
  attrs <- strsplit(attrs, split=" |;")
  for(i in 1:length(attrs)){
    ai <- match(attributes, attrs[[i]])
    if(sum(is.na(ai))==0){
      break
    }
  }
  Attrs <- do.call("cbind", lapply(attrs, function(x)x[ai+1]))
  Attrs <- t(Attrs)
  colnames(Attrs) <- attributes
  gtf.idx <- GTF[,7] == "."
  GTF[gtf.idx,7] = "*"
  GR <- GRanges(seqnames=GTF[,1], IRanges(start=GTF[,4], end=GTF[,5]), strand=GTF[,7], source=GTF[,2], feature=GTF[,3], score=GTF[,6], frame=GTF[,8], Attrs)
  



  if(feature[1]=="all"){
    return(GR)
  }else{ #split by feature
    #exon
    if(length(grep("exon", feature, ignore.case=TRUE))>0){
      cat("extract exon ...\n")
      GRexon <- GR[GR@elementMetadata$feature=="exon"]
      #GRL <- list(exon=sort(GRexon))
      GRL <- list(exon=GRexon)
    }
    if(length(grep("Intronic_Only", feature, ignore.case=TRUE))>0){
      cat("extract Intron only ...\n")
      GRintronOnly <- GR[GR@elementMetadata$feature=="Intronic_Only"]
      #GRL <- list(exon=sort(GRexon))
      GRL <- list(intronOnly=GRintronOnly)
    }
    #CDS, stop_codon will be included
    if(length(grep("CDS", feature, ignore.case=TRUE))>0){
      cat("extract CDS ...\n")
      cds <- GR[GR@elementMetadata$feature %in% c("CDS", "stop_codon")]
      cds <- GRanges(seqnames=paste(seqnames(cds), cds@elementMetadata$transcript_id, sep="|"), ranges=ranges(cds), strand=strand(cds))
      GRcds <- reduce(cds)
      chrid <- do.call("rbind", strsplit(as.character(seqnames(GRcds)), split="\\|"))
      GRcds <- GRanges(seqnames=chrid[,1], ranges=ranges(GRcds), strand=strand(GRcds), transcript_id=chrid[,2])
      GRL <- c(GRL, CDS=GRcds)
    }
    if(length(grep("utr|intron", feature)  )>0){
      #intron
      #define gaps
      cat("extract intron/utr ...\n")
      GRexon1 <- GRanges(seqnames=paste(seqnames(GRexon), GRexon@elementMetadata$transcript_id, sep="|"), ranges=ranges(GRexon), strand=strand(GRexon))
      GRgaps <- gaps(GRexon1)
      GRgaps <- GRgaps[start(GRgaps)!=1]
      #define intron, gaps
      chrid <- do.call("rbind", strsplit(as.character(seqnames(GRgaps)), split="\\|"))
      intron <- GRanges(seqnames=chrid[,1], ranges=ranges(GRgaps), strand=strand(GRgaps), transcript_id=chrid[,2])
      if(length(grep("intron", feature, ignore.case=TRUE))>0){
        GRL <- c(GRL, intron=intron)
      }
      if(length(grep("utr", feature, ignore.case=TRUE))>0){
        GRLcds <- split(GRcds, GRcds@elementMetadata$transcript_id)
        GRLexon <- split(GRexon, GRexon@elementMetadata$transcript_id)
        #filter
        GRLexon <- GRLexon[names(GRLexon) %in% names(GRLcds)]
        if(sum(names(GRLcds)!=names(GRLexon))>0)stop("exon and CDS not match")
        exons <- start(GRLexon)
        exone <- end(GRLexon)
        cdss <- start(GRLcds)
        cdse <- end(GRLcds)
        tids <- names(GRLexon)
        strands <- strand(GR)[match(tids, GR@elementMetadata$transcript_id)]
        seqs <- seqnames(GR)[match(tids, GR@elementMetadata$transcript_id)]
        
        u1 <- GRanges(seqs, IRanges(min(exons), min(cdss)-1), strand=strands, transcript_id=tids)
        u1 <- u1[start(u1)!=(end(u1)+1)]
        u2 <- GRanges(seqs, IRanges(max(cdse)+1, max(exone)), strand=strands, transcript_id=tids)
        u2 <- u2[(start(u2)-1)!=end(u2)]
        
        UTR1 <- diffGR(u1, intron)
        UTR2 <- diffGR(u2, intron)
        utr5 <- UTR1[strand(UTR1)=="+"]
        utr5 <- sort(c(utr5, UTR2[strand(UTR2)=="-"]))
        utr3 <- UTR1[strand(UTR1)=="-"]
        utr3 <- sort(c(utr3, UTR2[strand(UTR2)=="+"]))
        GRL <- c(GRL, utr3=utr3, utr5=utr5)
      }
      #GRL <- c(GRL, attributes=list(unique(Attrs)))
    }
    return(GRL)
  } 
}

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
GTF_FILE=args[5]
READ_COUNT_FOLDER=args[6]
LIB_TYPE="not-used"
OUTPUT_PREFIX = args[7]
GENE_LENGTHS_FILE = args[8]

# CONDITION_1="TD206_83L_total_tumor"
# CONDITION_2="TD206_83L_F480"
# NUM_REP_CONDITION1=1
# NUM_REP_CONDITION2=1
# GTF_FILE="GTF_Files/RefSeq_GeneBody.gtf"
# READ_COUNT_FOLDER="count"
# LIB_TYPE="not-used"
# OUTPUT_PREFIX = "DESeq_v2_RefSeq_GeneBody"
# GENE_LENGTHS_FILE = "lengths/Exon_Regions_Lengths.txt"

# CONDITION_1="HypoxMale"
# CONDITION_2="GH30min"
# NUM_REP_CONDITION1=2
# NUM_REP_CONDITION2=2
# GTF_FILE="GTF_Files/RefSeq_GeneBody.gtf"
# READ_COUNT_FOLDER="count"
# LIB_TYPE="not-used"
# OUTPUT_PREFIX = "DESeq_v2_RefSeq_GeneBody"
# GENE_LENGTHS_FILE = "lengths/Exon_Regions_Lengths.txt"

print("Arguments for differentialAnalysisDESeq.R:")
print(CONDITION_1)
print(CONDITION_2)
print(NUM_REP_CONDITION1)
print(NUM_REP_CONDITION2)
print(GTF_FILE)
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
	#If using HTseq (last 5 rows are special counters)
	#get rid of the last 5 rows
	#temp <- temp[1:(nrow(temp)-5),]
	#The DiffExp.qsub already removed these rows
	#If using featureCounts from DiffExp.qsub (need to keep all rows)

  numMappedReads <- read.table(file=paste(READ_COUNT_FOLDER, "/",CONDITION_1, (i-1), "_num_mapped_reads.txt", sep=""), colClasses="numeric")
	
  if(i == 1){
		condition1Count <- temp$count
		numMappedCondition1List <- numMappedReads
	}else {
		condition1Count <- cbind(condition1Count, temp$count)
		numMappedCondition1List <- cbind(numMappedCondition1List, numMappedReads)
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
	#If using HTseq (last 5 rows are special counters)
	#get rid of the last 5 rows
	#temp <- temp[1:(nrow(temp)-5),]
	#The DiffExp.qsub already removed these rows
	#If using featureCounts from DiffExp.qsub (need to keep all rows)
	numMappedReads <- read.table(file=paste(READ_COUNT_FOLDER, "/",CONDITION_2, (i-1), "_num_mapped_reads.txt", sep=""), colClasses="numeric")
	
	if(i == 1){
		condition2Count <- temp$count
		numMappedCondition2List <- numMappedReads
	}else {
		condition2Count <- cbind(condition2Count, temp$count)
		numMappedCondition2List <- cbind(numMappedCondition2List, numMappedReads)
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
names(result) <- paste(c(rep("",4), rep("DESeq_",4)), names(result), sep="")


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
mergedResult <- cbind("id"=rownames(sortedResultEdger), sortedResultEdger$table[,c("logFC", "PValue", "FDR")])
names(mergedResult) <- paste("edgeR", names(mergedResult), sep="_")

mergedResult <- merge(result, mergedResult, by.x="id", by.y="edgeR_id", sort=FALSE)

###############################################################
#4. PREPPING OUTPUT FILES
###############################################################
gtf_file <- read.GTF(file=GTF_FILE, feature=c("all"))
lookupTable <- cbind("transcript_id"=as.character(elementMetadata(gtf_file)$transcript_id), 
                     "gene_name"=as.character(elementMetadata(gtf_file)$gene_name))
lookupTable <- as.data.frame(lookupTable, stringsAsFactors = FALSE)
lookupTable <- unique(lookupTable)
lookupTableSplitted <- split(lookupTable, lookupTable$gene_name)

mergedResult <- as.data.frame(mergedResult)
countNorm <- counts(deseqOutput, normalized=TRUE )
colnames(countNorm) <- paste("norm_", colnames(countNorm), sep="")
mergedResult <- cbind(countTable, rpkmCondition1,tpmCondition1, rpkmCondition2, tpmCondition2, countNorm, mergedResult) #row order is OK
mergedResult <- mergedResult[order(mergedResult$id), ]
lookupTableSplitted <- lookupTableSplitted[order(names(lookupTableSplitted))]
accessionList <- sapply(seq(1,length(lookupTableSplitted)), function(i){
  return(paste(lookupTableSplitted[[i]]$transcript_id, sep="", collapse=";"))
})

mergedResult <- cbind("edgeR_logFC_copy"=mergedResult$edgeR_logFC, mergedResult)
mergedResult <- cbind("edgeR_FDR_copy"=mergedResult$edgeR_FDR, mergedResult)
new_edgeRFC <-(paste(CONDITION_1,"_", CONDITION_2,"_edgeRlogFC",sep=""))  # vairable with prefix for logFC
new_edgeRFDR <-(paste(CONDITION_1,"_", CONDITION_2,"_edgeRFDR",sep=""))    # variable with prefix for FDR
colnames(mergedResult)[colnames(mergedResult) == "edgeR_logFC_copy"] <- new_edgeRFC
colnames(mergedResult)[colnames(mergedResult) == "edgeR_FDR_copy"] <- new_edgeRFDR


mergedResult <- cbind("accession"=accessionList, mergedResult)#row order is OK
mergedResult <- cbind("id"=mergedResult$id, mergedResult[,-match("id", names(mergedResult))]) 


# sum( abs(mergedResult$DESeq_log2FoldChange)>1 & mergedResult$DESeq_padj < 0.05)#0
# sum( abs(mergedResult$edgeR_logFC)>=1 & mergedResult$edgeR_FDR <= 0.05)#125
# fivenum(mergedResult$edgeR_logFC)
# fivenum(as.data.frame(results(deseqOutput))$log2FoldChange)

#Print output file location:
print(paste("output file is in: ", READ_COUNT_FOLDER, "/",OUTPUT_PREFIX,"_",CONDITION_1,"_", CONDITION_2, ".txt", sep=""))
#Create the output file:
write.table(mergedResult,file=paste(READ_COUNT_FOLDER, "/",OUTPUT_PREFIX,"_", CONDITION_1, "_", CONDITION_2, ".txt", sep=""), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
#----------------------------------------------------------------------------------
