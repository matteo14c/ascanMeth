#
library(data.table)
# parameters
args <- commandArgs(trailingOnly = TRUE)

# initial values
h1_val <- NULL
h2_val <- NULL
#n_val <- NULL
c_val <- NULL


#parse parameters
if (length(args) > 0) {
	i <- 1
	while (i <= length(args)) {
		if (args[i] == "--h1") {
			h1_val <- args[i+1]
			i <- i + 1
		} else if (args[i] == "--h2") {
			h2_val <- args[i+1]
			i <- i + 1
		#} else if (args[i] == "--n") {
		#	n_val <- as.numeric(args[i+1]) # Convert to numeric
		#	i <- i + 1
		} else if (args[i] == "--c") {
			c_val <- as.numeric(args[i+1]) # Convert to numeric
			i <- i + 1
		} else {
			warning(paste("Unknown argument:", args[i]))
		}
    	i <- i + 1
	}
}

header<-c("chr","start","end",
          "mod","score","strand","start","end",
          "color","coverage","ratio","nmod",
          "canonical","nother","ndelete","nfail","ndiff","nNo",
          "chr_bed","start_bed","end_bed","suffix",
          "name_CpG","len_CpG","num_CpG","numGC","pCpG","rCpG","ratioCpG","included")

# Input files. 
f1<-"hp1.with_cpg_islands_all.tsv" #h1_val
f2<-"hp2.with_cpg_islands_all.tsv"# h2_val

# Coverage threshold
covThr<-ifelse(is.null(c_val),10,c_val)
# 


# Read input files and attach header;
meth1<-read.table(f1,header=FALSE);
meth2<-read.table(f2,header=FALSE);
colnames(meth1)=header;
colnames(meth2)=header;
meth1$name_CpG = as.character(meth1$name_CpG)
meth2$name_CpG = as.character(meth2$name_CpG)


filterAndJoin <-function(meth1,meth2,covThr=10,mod="m")
{
  
  meth1 <- meth1[meth1$mod==mod & meth1$coverage>=covThr,]
  meth2 <- meth2[meth2$mod==mod & meth2$coverage>=covThr,]
  # compute common position
  
  commonPositionsM1<- meth1$chr %in% meth2$chr & 
    meth1$start %in% meth2$start & 
    meth1$end %in% meth2$end
  
  commonPositionsM2<- meth2$chr %in% meth1$chr & 
    meth2$start %in% meth1$start & 
    meth2$end %in% meth1$end
  

  # change meth1 names for CpG
  
  selName <- meth1$name_CpG!="."
  meth1$name_CpG[selName]=
    paste(meth1$chr_bed[selName],meth1$start_bed[selName],meth1$end_bed[selName],
          meth1$len_CpG[selName],meth1$num_CpG[selName],sep="_")
  
  # merge
  commonCpG<-data.frame(chr=meth1$chr[commonPositionsM1], 
                        start=meth1$start[commonPositionsM1], 
                        end=meth1$end[commonPositionsM1], 
                        mod=meth1$mod[commonPositionsM1], 
                        covM1=meth1$coverage[commonPositionsM1], 
                        nmodM1=meth1$nmod[commonPositionsM1], 
                        covM2=meth2$coverage[commonPositionsM2], 
                        nmodM2=meth2$nmod[commonPositionsM2],
                        CpGName=meth1$name_CpG[commonPositionsM1]
                        )
  return(commonCpG)
}


makeChunks <- function(commonCpG, targetReads = 400) {
  # skip non-island block
  if (length(unique(commonCpG$CpGName)) == 1 && unique(commonCpG$CpGName) == ".") {
    commonCpG$CpGName <- paste("nonCpG_island",commonCpG$chr,sep="_")
    commonCpG$chunkSizeUsed <- NA_integer_
    return(commonCpG)
  }
  
  # compute representative coverage
  totCoverage <- quantile(commonCpG$covM1 + commonCpG$covM2, probs = 0.1, na.rm = TRUE)
  chunkSize <- max(1L, ceiling(targetReads / totCoverage))
  
  # assign chunks within this island
  localChunk <- ceiling(seq_len(nrow(commonCpG)) / chunkSize)
  
  commonCpG$chunkSizeUsed <- chunkSize
  commonCpG$CpGName <- paste0(unique(commonCpG$CpGName), "_", localChunk)
  
  return(commonCpG)
}

refineNonCpG <- function(chunked_commonCpG, chunkSize = 6) {
  if (!startsWith(chunked_commonCpG$CpGName[1], "nonCpG_")) return(chunked_commonCpG)
  
  grp <- ceiling(seq_len(nrow(chunked_commonCpG)) / chunkSize)
  
  chunked_commonCpG$CpGName <- paste0(chunked_commonCpG$CpGName[1], "_", grp)
  chunked_commonCpG$chunkSizeUsed <- chunkSize
  
  return(chunked_commonCpG)
}

computeStats <- function(chunked_commonCpG) {
  chunkSize <- length(chunked_commonCpG$end)
  
  sumM1 <- sum(chunked_commonCpG$nmodM1)
  sumM2 <- sum(chunked_commonCpG$nmodM2)
  
  covM1 <- sum(chunked_commonCpG$covM1)
  covM2 <- sum(chunked_commonCpG$covM2)
  
  equalRatio <- (sumM1 + sumM2) / 2
  
  ascan1 <- ifelse(sumM1 == 0, 0, sumM1 * log(sumM1 / equalRatio))
  ascan2 <- ifelse(sumM2 == 0, 0, sumM2 * log(sumM2 / equalRatio))
  ascanStat <- 2 * (ascan1 + ascan2)
  
  methPc1 <- (sumM1+0.5) / (covM1)
  methPc2 <- (sumM2+0.5) / (covM2)
  
  logMethRatio <- log2(methPc1 / methPc2)
  
  results <- list(
    chr    = chunked_commonCpG$chr[1],
    start  = chunked_commonCpG$start[1],
    end    = chunked_commonCpG$end[chunkSize],
    cpG    = chunked_commonCpG$CpGName[1],
    chunk  = chunkSize,
    M1     = sumM1,
    C1     = covM1,
    M2     = sumM2,
    C2     = covM2,
    methPc1 = methPc1,
    methPc2 = methPc2,
    stat   = ascanStat,
    logFC  = logMethRatio
  )
  
  results$pvalue <- pchisq(results$stat, df = chunkSize, lower.tail = FALSE)
  
  return(results)
}

computeResults  <- function(inputL)
{
  splitL=split(inputL,inputL$CpGName)
  list_res=lapply(splitL, computeStats)
  return(list_res)
}


# filter and join
commonCpG <-filterAndJoin(meth1,meth2)
# make chunks on islands
chunked_commonCpG <- lapply(split(commonCpG, commonCpG$CpGName), makeChunks)
# compute 3rdQ, I need to unlist
vals <- unlist(lapply(chunked_commonCpG, function(x) unique(x$chunkSizeUsed)))
chunk_to_use  <- quantile(vals,0.75,na.rm = TRUE)
# make chunks on non islands
chunked_commonCpG <- lapply(chunked_commonCpG,refineNonCpG,chunkSize=chunk_to_use)
# compute results
#chunked_commonCpG <- do.call(rbind, chunked_commonCpG)
#rownames(chunked_commonCpG)  <-  NULL
chunked_commonCpG_results  <- lapply(chunked_commonCpG, computeResults)

flat_list <- unlist(chunked_commonCpG_results, recursive = FALSE)
datatable <- rbindlist(flat_list,use.names = TRUE, fill = TRUE)
datatable$fdr  <- p.adjust(datatable$pvalue)
datatable$significant   <- ifelse(datatable$fdr<=0.01 &
                                    abs(datatable$logFC)>=2 &
                                    datatable$chunk!=1 &
                                    (datatable$methPc1 >= 0.6 | datatable$methPc2>=0.6),
                                        "significant",
                                        "not_significant")
datatable    <-  datatable[order(datatable$start),]
fwrite(datatable, "final_output.csv",sep = "\t")
#####
#Add different cpg density
  

