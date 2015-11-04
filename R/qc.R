#' Summary of SPP quality metrics from files in a directory
#'
#' @param dir path to directory that holds the SPP quality results
#' @param pattern regular expression to identify the SPP files
#' @return a table of SPP quality metrics
#' @export
spp.qc.summary <- function(dir, pattern) { 
  qc.files = list.files(dir, pattern="_qc.txt$", full=T)
  qc = NULL
  for (qc.file in qc.files) {
    qc = rbind(qc, read.table(qc.file))
  }
  ## set the colnames (they are not included in the files)
  colnames(qc) = c("Filename", "numReads", "estFragLen", "corr_estFragLen",
            "phantomPeak", "corr_phantomPeak", "argmin_corr", "min_corr",
            "NSC", "RSC", "QualityTag")
  return(qc)
}

#' Plot SPP quality metrics
#'
#' @param qc.table data.frame created with \code{\link{spp.qc.summary}}
#' @export
plot.spp.qc <- function(qc.table) {
  require(ggplot2)
  require(reshape)
  require(grid)

  tab$estFragLen = as.numeric(sapply(strsplit(tab$estFragLen, ","), "[", 1))

  a = qplot(x=sample, y=numReads, geom="bar", stat="identity", data=tab, fill=numReads > 1e7) + coord_flip() + geom_hline(aes(yintercept=1e7, colour="red")) + scale_y_log10()

  b = qplot(x=sample, y=estFragLen, geom="bar", stat="identity", data=tab, ylab="SPP strand shift") + coord_flip()

  c = qplot(x=sample, y=bioanalyzer_shift_size, geom="bar", stat="identity", data=tab, ylab="Bioanalyzer strand shift") + coord_flip()
  
  d = qplot(x=sample, y=NSC, geom="bar", stat="identity", data=tab, fill=NSC > 1.05) + coord_flip() + geom_hline(aes(yintercept=1.05, colour="red")) 

  e = qplot(x=sample, y=RSC, geom="bar", stat="identity", data=tab, fill=RSC > 0.8) + coord_flip() + geom_hline(aes(yintercept=0.8, colour="red"))

  grid.newpage()
  pushViewport(viewport(layout=grid.layout(1, 5, widths=unit(c(2.2, 0.7, 0.7, 0.7, 0.7), "null"))))
  vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
  
  no.yaxis = theme(axis.line=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())
  no.xaxis.label = theme(axis.title.x=element_blank(), legend.position="none")
  no.legend = theme(legend.position="none")
  
  print(a + no.legend, vp=vplayout(1,1))
  print(b + no.yaxis, vp=vplayout(1,2))
  print(c + no.yaxis, vp=vplayout(1,3))
  print(d + no.yaxis + no.legend, vp=vplayout(1,4))
  print(e + no.yaxis + no.legend, vp=vplayout(1,5))
  
}
