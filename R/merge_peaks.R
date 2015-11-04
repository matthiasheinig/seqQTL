#' Merge MACS peaks calls from files in a directory
#'
#' @param dir path to the directory that holds the peak call files
#' @param FDR threshold for the false discovery rate
#' @param pattern regular expression to identify the peak files
#' @return a GRanges object with the merged peak regions
#' @export
merge.peaks <- function(dir="results/current/macs", FDR = 5, pattern="peaks.xls$") { 
  require(GenomicRanges)
  
  name = basename(base)
  peak.caller = unlist(strsplit(name, "-"))[1]
  tissue = unlist(strsplit(name, "-"))[2]
  modification = unlist(strsplit(name, "-"))[3]

  print(paste("read MACS directory:",dir))

  files = list.files(path=dir, pattern=pattern, full.names = F)

  for (i in 1:length(files)) {

    cat(files[i], "(FDR=",FDR,"%)\n")
    macs.GR <- read.table(file.path(dir, files[i]), header=T, sep="\t", stringsAsFactors=F)
    macs.GR <- macs2GRange(macs.GR)
    runValue(strand(macs.GR)) = "+"

    # filter regions based on macs FDR
    macs.GR.fdr = subset(macs.GR, elementMetadata(macs.GR)$fdr < FDR & as.character(seqnames(macs.GR)) != "*")
    
    cat(name, "add=", length(macs.GR.fdr),"\n")
    
    if (i>1) { macs.GRlst <- c(macs.GRlst, GRangesList(macs.GR.fdr)) }
    if (i==1)  { macs.GRlst <- GRangesList(macs.GR.fdr) }
   
  }

  macs.GRlst.union = macs.GRlst[[1]]
  for (i in 2:length(macs.GRlst)) {
    macs.GRlst.union = union(macs.GRlst.union, macs.GRlst[[i]])
  }
  

  #obtain overlap (YES/NO) of single samples with union 
  OL <- do.call( "rbind" , lapply ( macs.GRlst , function (x)  countOverlaps(macs.GRlst.union,x))  )
  elementMetadata(macs.GRlst.union)$countOL = colSums(OL)
 
  #keep peaks present in more than one sample (countOverlaps > 1) 
  macs.GRlst.union.VAL = subset(macs.GRlst.union, elementMetadata(macs.GRlst.union)$countOL > 1) 

  # definiton of ID: chr:start..end
  elementMetadata(macs.GRlst.union.VAL)$ID = paste(seqnames(macs.GRlst.union.VAL),":",
                                                   start(macs.GRlst.union.VAL),"..",end(macs.GRlst.union.VAL),sep="")

  fname = paste(name, ".gff", sep="")
  ## GR2gff(macs.GRlst.union.VAL, file=file.path(outdir, fname))
  ## print(paste("write file:",file.path(outdir, fname)))
  
  return(macs.GRlst.union.VAL)
}


#' Transform a MACS output peakfile into a GRanges object 
#'
#' @param peaks a data.frame with output from MACS
#' @param chroms character vector with chromosome names to be used. Default is
#'        to use all chromosome names from the MACS file
#' @return a GRanges object with the MACS peak calls
macs2GRange <- function(peaks, chroms = NULL) {
  myrange=GRanges(peaks[,1],
    IRanges(peaks[,2], peaks[,3],
            names=paste(peaks[,1],(peaks[,5]+peaks[,2]),sep=":")),
    count=peaks[,6], 
    score=peaks[,7], 
    FE=peaks[,8], 
    fdr=peaks[,9], 
    summit=(peaks[,5]+peaks[,2]) )
  if (is.null(chroms)) chroms = seqlevels(myrange)
  myrange=myrange[which(seqnames(myrange) %in% chroms)]
  return(myrange)
}
