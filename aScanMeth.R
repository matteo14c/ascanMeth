# parameters
args <- commandArgs(trailingOnly = TRUE)

# initial values
h1_val <- NULL
h2_val <- NULL
n_val <- NULL
c_val <- NULL

if (length(args) > 0) {
	i <- 1
	while (i <= length(args)) {
		if (args[i] == "--h1") {
			h1_val <- args[i+1]
			i <- i + 1
		} else if (args[i] == "--h2") {
			h2_val <- args[i+1]
			i <- i + 1
		} else if (args[i] == "--n") {
			n_val <- as.numeric(args[i+1]) # Convert to numeric
			i <- i + 1
		} else if (args[i] == "--c") {
			c_val <- as.numeric(args[i+1]) # Convert to numeric
			i <- i + 1
		} else {
			warning(paste("Unknown argument:", args[i]))
		}
    	i <- i + 1
	}
}

header<-c("chr","start","end","mod","score","strand","start","end","color","coverage","ratio","nmod","canonical","nother","ndelete","nfail","ndiff","nNo")

# Input files. Could be parameters
f1<-h1_val
f2<-h2_val

# Coverage threshold
covThr<-ifelse(is.null(c_val),5,c_val)
numThr<-ifelse(is.null(n_val),20,n_val)
# 

# Add header to bedMethyl. Just to make it more readable and say what's what
header<-c("chr","start","end","mod","score","strand","start","end","color","coverage","ratio","nmod","canonical","nother","ndelete","nfail","ndiff","nNo");

# Read input files and attach header;
meth1<-read.table(f1,header=FALSE);
meth2<-read.table(f2,header=FALSE);

colnames(meth1)=header;
colnames(meth2)=header;


# Compute  common CpGs. Since the file are sorted and in the same order, they can be juxtaposed later
commonPositionsM1<- meth1$chr %in% meth2$chr & meth1$start %in% meth2$start & meth1$end %in% meth2$end
commonPositionsM2<- meth2$chr %in% meth1$chr & meth2$start %in% meth1$start & meth2$end %in% meth1$end



# Merge the 2 files at common CpG sites. 
# NOTE: What should we do with the additional data/non shared CpGs?
# See/check if there in an excess of the number of CpGs in mom VS dad?
# We could use the same logic used for differential methylation and implement
# some kind of test to check significant increase in CpGs
commonCpG<-data.frame(chr=meth1$chr[commonPositionsM1], start=meth1$start[commonPositionsM1], end=meth1$end[commonPositionsM1], mod=meth1$mod[commonPositionsM1], covM1=meth1$coverage[commonPositionsM1], nmodM1=meth1$nmod[commonPositionsM1], covM2=meth2$coverage[commonPositionsM2], nmodM2=meth2$nmod[commonPositionsM2])


# Make a copy of the dataset
# Being lazy, since filters applied later on
#commonCpGunf<-commonCpG;

# Filter coverage
commonCpG<-commonCpG[commonCpG$covM1 >=covThr & commonCpG$covM2 >=covThr ,]

# Compute equal methyl ratio (expected, this is m in the formula used by aScan)
commonCpG$equalRatioMod=(commonCpG$nmodM1+commonCpG$nmodM2)/2
# Compute value for H1: c1*log(c1/m)
commonCpG$ascan1= ifelse(commonCpG$nmodM1==0, 0, commonCpG$nmodM1 * log (commonCpG$nmodM1/commonCpG$equalRatioMod ) )
# compute value for H2: c2*log(c2/m)
commonCpG$ascan2= ifelse(commonCpG$nmodM2==0, 0, commonCpG$nmodM2 * log (commonCpG$nmodM2/commonCpG$equalRatioMod ) )
# compute aScanStatistics: 2*(H1+H2);
commonCpG$ascanStat= 2*(commonCpG$ascan1 + commonCpG$ascan2)
# log meth ratio: used for plots
commonCpG$logMethRatio= log2(commonCpG$nmodM1/commonCpG$nmodM2)


# From here onward we could check whatever!
# If stats are computed throughout the table, then, easy to iterate N-rows at a time and compute the p-values.

# Compute p-values and FDR.
# Note: chunKsize sets the number of consecutive, non overlapped CpGs to test.
# Questions: should we remove from testing unmethylated CpGs? I fathom not

computeStatsChunks<-function(commonCpG,chunkSize=20,mod="m")
{
	# keep only the modification I need
	commonCpG<-commonCpG[commonCpG$mod==mod,]
	# make vectors into matrices, number of rows is chunkSize
	# I will keed only a multiple of chunkSize values
	# the last N will be discarded. If chunkSize=20, up to 19 will be discarded.
	keep<-nrow(commonCpG) - nrow(commonCpG) %% chunkSize
	# adjust to multiple of chunKsize
	commonCpG<-commonCpG[1:keep,]

	# use R matrix to make matrices of chunkSize rows
	fakeMatStat<-matrix(commonCpG$ascanStat,nrow=chunkSize)
	fakeMatStart<-matrix(commonCpG$start,nrow=chunkSize)
	fakeMatEnd<-matrix(commonCpG$end,nrow=chunkSize)
	fakeMatChr<-matrix(commonCpG$chr,nrow=chunkSize)
	fakeMatMod1<-matrix(commonCpG$nmodM1,nrow=chunkSize)
	fakeMatMod2<-matrix(commonCpG$nmodM2,nrow=chunkSize)
	fakeMatCoV1<-matrix(commonCpG$covM1,nrow=chunkSize)
        fakeMatCoV2<-matrix(commonCpG$covM2,nrow=chunkSize)

	
	# Sum by columns do the chunsize consecutive elements are selected
	sumStat<-colSums(fakeMatStat)
	sumM1<-colSums(fakeMatMod1)
	sumM2<-colSums(fakeMatMod2)
	covM1<-colSums(fakeMatCoV1)
	covM2<-colSums(fakeMatCoV2)


	
	# assemble results table
	results<-data.frame(chr=fakeMatChr[1,] , start=fakeMatStart[1,],end=fakeMatEnd[chunkSize,],M1=sumM1,C1=covM1,M2=sumM2,C2=covM2,stat=sumStat)
	
	# compute pvalues
	results$pvalue<-pchisq(results$stat, chunkSize, lower.tail = FALSE)
	#fdr
	results$FDR<-p.adjust(results$pvalue)
	results$logFC <- log2(results$M1/results$M2)
	
	return(results)
}

# results for m: m indicates 5mC, 5-methylcytosine
results_m <- computeStatsChunks(commonCpG,chunkSize=numThr)

# results for h: h indicates 5hmC, 5-hydroxymethylcytosine
results_h <- computeStatsChunks(commonCpG,mod="h",chunkSize=numThr)


# I selected a region with a very significant corrected p-values, actually 3 regions.
# These will be plotted

#plot(commonCpG$start[commonCpG$start >=193406526 & commonCpG$end <=193419432 & commonCpG$mod=="m"], 
#commonCpG$logMethRatio[commonCpG$start >=193406526 & commonCpG$end <=193419432 & commonCpG$mod=="m"],
#xlab="ch4_T2T", ylab="log2(meth1/meth2)",pch=20)

#abline(h=1,col="red",lty=20)
#abline(h=-1,col="blue",lty=20)

# Significant hypermethylation of H1 is observed, data saved here
#singnR<-commonCpG[commonCpG$start >=193408987 & commonCpG$end <=193416629 & commonCpG$mod=="m" ,c("start","covM1","nmodM1","covM2","nmodM2")]

# write results
write.table(results_m,file="results_M_CpG.csv",row.names=F);
write.table(results_h,file="results_H.CpG.csv",row.names=F);

# compute results with smaller chunkSize
#results_m_10 <- computeStatsChunks(commonCpG,chunkSize=10)

######################################################################
# different CpG density!
# Here I tag the 2 haplotypes: h1 and h2;
# Merge and order the data;
# Filter: 1 modification, sufficient coverage;
# Take N CpGs at a time, and test if they are equally distributed (N/2)
# between h1 and h2

countCpG<-function(meth1,meth2,chunkSize=40,mod="m",cov=5)
{

	# Tagging
	meth1$haplotype="h1"
	meth2$haplotype="h2"
	# Merging.
	merged=rbind(meth1,meth2)
	merged_m=merged[merged$mod==mod,]
	# Sorting and filtering
	merged_m=merged_m[order(merged_m$start),]
	merged_m=merged_m[merged_m$coverage>=cov,]
	keep=nrow(merged_m) -nrow(merged_m) %% chunkSize
	# R magic to divide data in matrices with chunkSize rows
	M_chunk_m=matrix(merged_m$haplotype[1:keep],nrow=chunkSize)
	M_chunk_start=matrix(merged_m$start[1:keep],nrow=chunkSize)
	M_chunk_end=matrix(merged_m$end[1:keep],nrow=chunkSize)
	M_chunk_chr=matrix(merged_m$chr[1:keep],nrow=chunkSize)

	# Contingency table of h1 and h2
	countF<-function(v){return(c(sum(v=="h1"),sum(v=="h2")))}
	contingency=matrix(apply(M_chunk_m,2,countF),ncol=2,byrow=T)
	colnames(contingency)=c("h1","h2")
	# Assemble results
	results=data.frame(chr=M_chunk_chr[1,], start=M_chunk_start[1,], end=M_chunk_end[nrow(M_chunk_end),], h1=contingency[,"h1"],h2=contingency[,"h2"])
	
	# Compute stats
	mHalf=chunkSize/2
	results$h1stat= ifelse(results$h1==0, 0, results$h1 * log (results$h1/mHalf))
	results$h2stat= ifelse(results$h2==0, 0, results$h2 * log (results$h2/mHalf))
	results$aScanstat=2*(results$h1stat+results$h2stat)
	# Compute p-value
	results$pvalue<-pchisq(results$aScanstat, 2, lower.tail = FALSE)
        # Compute FDR
        results$FDR<-p.adjust(results$pvalue)
	
	return(results)
}

resCpGdens<-countCpG(meth1,meth2,chunkSize=numThr);
write.table(resCpGdens,"results_D.CpG.csv",row.names=FALSE)




