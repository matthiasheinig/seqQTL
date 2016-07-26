#' Helper function to extract text using regular expressions
#'
#' @param pattern regular expression with grouping
#' @param string character vector from which to extract the text
#' @param perl logical indicating whether regular expression is perl style
#' @return character vector with the extracted text
extract <- function (pattern, string, perl = TRUE) {
  r <- paste(".*", pattern, ".*", sep = "")
  matched <- grep(r, string, perl = perl)
  result <- rep(NA, length(string))
  result[matched] <- sub(r, "\\1", string[matched], perl = perl)
  return(result)
}


#' Read gene annotations in gff format
#'
#' @param filename gff filename
#' @param gffAttrNames character vector with attribute names that should be
#'        parsed from the gtf file. Usually these are "ID", "Parent"
#'        and the like. By default no attributes are parsed.
#' @return a GRanges object with all entries from the gff file with metadata
#'         columns for source, type and the attributes parsed
#' @export
gff2GR <- function(filename, gffAttrNames=NULL) {
  # read gff into genomic ranges
  # require(GenomicRanges)
  regions = scan(filename, what=list(character(), character(), character(), numeric(), numeric(), character(), character(), character(), character()), comment.char="#", sep="\t")

  strand = regions[[7]]
  strand[strand == "."] = "*"
  strand[strand == "1"] = "+"
  strand[strand == "-1"] = "-"

  gr = GRanges(seqnames=regions[[1]],
    ranges=IRanges(start=regions[[4]], end=regions[[5]]),
    strand=strand)

  src = regions[[2]]
  type = regions[[3]]
  score = regions[[6]]
  df = DataFrame(src, type, score)
  
  if (!is.null(gffAttrNames)) {
    df = cbind(df, DataFrame(sapply(gffAttrNames, function(n)
      extract(paste(n , "=(.+?)(;|$)", sep=""), regions[[9]]))))
  }
  elementMetadata(gr) = df
  
  return(gr)
}


#' Export GRanges to gff format
#'
#' @param regions GRanges object to be saved
#' @param filename name of the gff file
#' @param feature.type character for the feature type field
#' @param src character for the source type field
#' @param score numeric or character for the score field
#' @param phase numeric or character for the phase field
#' @export
GR2gff <- function(regions, filename, feature.type="experimental_feature", src="GenomicRanges", score=".", phase=".") {
  # require(GenomicRanges)

  strnd = as.character(strand(regions))
  strnd[strnd == "*"] = "."
  
  tab = data.frame(as.character(seqnames(regions)), src, feature.type, as.numeric(start(regions)), as.numeric(end(regions)), score, strnd, phase, makeGffAttributes(as.data.frame(elementMetadata(regions))), stringsAsFactors=F)

  write.table(tab, file=filename, sep="\t", quote=F, row.names=F, col.names=F)
}


#' Helper function to generate the gtf attributes field
#'
#' @param df data.frame with attribute information
#' @param cols character vector with column names of the attributes to be used.
#'        With the default value of NULL all columns will be used.
makeGtfAttributes <- function(df, cols=NULL) {
  if (is.null(cols))
    cols = colnames(df)
  # make sure that gene_id and transcript_id are the first two columns
  mandatory = c("gene_id", "transcript_id")
  o = match(c(mandatory, setdiff(cols, mandatory)), cols)
  if (any(is.na(o[1:length(mandatory)]))) {
    warning("mandatory gtf attributes gene_id or transcript_id missing")
    o = o[!is.na(o)]
  }
  cols = cols[o]
  return(paste(apply(sapply(cols, function(s) {
           content = df[,s]
           if (is.character(content) | is.factor(content)) {
             content = paste('"', content, '"', sep="")
           }
           paste(gsub(".", "_", s, fixed=T), content, sep=" ")
         }), 1, paste, collapse="; "), ";", sep=""))
}


#' Read gene annotations in gtf format
#'
#' @param filename gtf filename
#' @param gtfAttrNames character vector with attribute names that should be
#'        parsed from the gtf file. Usually these are "gene_id", "gene_name",
#'        "transcript_id" and the like. By default no attributes are parsed.
#' @return a GRanges object with all entries from the gtf file with metadata
#'         columns for source, type and the attributes parsed
#' @export
gtf2GR <- function(filename, gtfAttrNames=NULL) {
  # read gff into genomic ranges
  # require(GenomicRanges)
  regions = scan(filename, what=list(character(), character(), character(), numeric(), numeric(), character(), character(), character(), character()), comment.char="#", sep="\t")

  strand = regions[[7]]
  strand[strand == "."] = "*"

  gr = GRanges(seqnames=regions[[1]],
    ranges=IRanges(start=regions[[4]], end=regions[[5]]),
    strand=strand)

  src = regions[[2]]
  type = regions[[3]]

  df = data.frame(src, type, stringsAsFactors=F)
  if (!is.null(gtfAttrNames)) {
    df = data.frame(df, sapply(gtfAttrNames, function(n)
      extract(paste(n , " (.+?)(;|$)", sep=""), regions[[9]])),
      stringsAsFactors=F)
  }
  elementMetadata(gr) = df
  return(gr)
}

#' Export GRanges to gtf format
#'
#' @param regions GRanges object to be saved
#' @param filename name of the gtf file
#' @param feature.type character for the feature type field
#' @param src character for the source type field
#' @param score numeric or character for the score field
#' @param phase numeric or character for the phase field
#' @param attributes character vector of columns to be saved
#' @param ... additional arguments passed to \code{\link{write.table}}
#' @export
GR2gtf <- function(regions, filename, feature.type="experimental_feature", src="GenomicRanges", score=".", phase=".", attributes=NULL, ...) {
  # require(GenomicRanges)

  strnd = as.character(strand(regions))
  strnd[strnd == "*"] = "."
  
  tab = data.frame(as.character(seqnames(regions)), src, feature.type, as.numeric(start(regions)), as.numeric(end(regions)), score, strnd, phase, makeGtfAttributes(as.data.frame(elementMetadata(regions)), cols=attributes), stringsAsFactors=F)

  write.table(tab, file=filename, sep="\t", quote=F, row.names=F, col.names=F, ...)
}


#' Read genomic ranges from bed files (zero based coordinates)
#'
#' @param filename name of the bed file
#' @param nfields number of fields in the bed file to be parsed. Default is to
#'        automatically detect the number of fields.
#' @return GRanges object with the coordinates from the bed file
#' @export
bed2GR <- function(filename, nfields=NULL) {
  # read gff into genomic ranges
  # require(GenomicRanges)
  what = list(character(), numeric(), numeric(), character(), numeric(), character())
  if (is.null(nfields)) {
    hdr = read.csv(filename, sep="\t", header=F, nrows=3)
    nfields = ncol(hdr)
  }
  if (nfields > length(what)) {
    for (i in (length(what) + 1):nfields) {
      what[[i]] = character()
    }
  }
  regions = scan(filename, what=what[1:nfields])

  if (nfields >= 6) {
    strand = regions[[6]]
    strand[strand == "."] = "*"
  } else {
    strand = "*"
  }

  ## bed files are zero based so add one to the start
  gr = GRanges(seqnames=regions[[1]],
    ranges=IRanges(start=regions[[2]] + 1, end=regions[[3]]),
    strand=strand)

  if (nfields >= 4) {
    names(gr) = regions[[4]]
    if (nfields >= 5) {
      m = DataFrame(score=as.numeric(regions[[5]]))
      if (nfields >= 7) {
        for (i in 7:nfields) {
          m = DataFrame(m, regions[[i]])
          colnames(m)[ncol(m)] = paste("field", i, sep="")
        }
      }
      elementMetadata(gr) = m
    }
  }
  
  return(gr)
}


#' Helper function to generate the gff attributes field
#'
#' @param df data.frame with attribute information
#' @param cols character vector with column names of the attributes to be used.
#'        With the default value of NULL all columns will be used.
makeGffAttributes <- function(df, cols=NULL) {
  if (ncol(df) == 0)
    return(rep("", nrow(df)))
  if (is.null(cols))
    cols = colnames(df)
  return(apply(sapply(cols, function(s) paste(gsub(".", "_", s, fixed=T), df[,s], sep="=")), 1, paste, collapse=";"))
}



#' Count reads from a bam file in genomic ranges
#'
#' @param bam.file filename of the bam file. The bam file must be sorted and
#'        indexed
#' @param granges GRanges object in which to count the reads
#' @param min.mapq minimal mapping quality. Default is to ignore mapping quality
#' @param read.width the length of a read to be used for overlaps. Default is to #'        use just the start position of the read.
#' @return count vector of the same length as the GRanges object
#' @export
countBamInGRanges <- function(bam.file, granges, min.mapq=NULL, read.width=1) {
  # require(GenomicRanges)
  # require(Rsamtools)
  
  rds.counts <- numeric(length(granges))
  seq.names <- as.character(unique(seqnames(granges)))
  seq.names.in.bam <- names(scanBamHeader(bam.file)[[1]]$targets)
  for (seq.name in seq.names) {
    if (seq.name %in% seq.names.in.bam) {
      granges.subset <- granges[seqnames(granges)==seq.name]
      strand(granges.subset) <- "*"
      rds <- scanBam(bam.file,param=ScanBamParam(what=c("pos","mapq"),which=range(granges.subset)))
      if (!is.null(min.mapq)) {
        mapq.test <- rds[[1]]$mapq >= min.mapq & !is.na(rds[[1]]$mapq)
      } else {
        mapq.test = rep(T, length(rds[[1]]$mapq))
      }
      if (sum(mapq.test) > 0) {
        rds.ranges <- GRanges(seq.name,IRanges(start=rds[[1]]$pos[mapq.test],width=read.width))
        rds.counts.seq.name <- countOverlaps(granges.subset,rds.ranges)
        rds.counts[as.logical(seqnames(granges)==seq.name)] <- rds.counts.seq.name
      } else {
        rds.counts[as.logical(seqnames(granges)==seq.name)] <- 0
      }
    } else {
      rds.counts[as.logical(seqnames(granges)==seq.name)] <- 0
    }
  }
  rds.counts
}

getBins <- function(chr=NULL, n=NULL, bin.size=NULL, genome=Rnorvegicus, offset=0) {
  stopifnot(!all(c(is.null(n), is.null(bin.size)), "specify either bin size or number of bins"))
  if (is.null(chr)) {
    chr = seqnames(genome)
  }
  if (!is.null(n)) {
    bin.size = floor((seqlengths(genome)[chr] - offset) / n)
    names(bin.size) = chr
    n = rep(n, length(chr))
    names(n) = chr
  } else {
    n = floor((seqlengths(genome)[chr] - offset) / bin.size)
    names(n) = chr
    bin.size = rep(bin.size, length(chr))
    names(bin.size) = chr
  }
  
  g = GRanges()
  for (ch in chr) {
    g = c(g, GRanges(seqnames=ch, IRanges(start=0:(n[ch] - 1) * bin.size[ch] + 1 + offset, width=bin.size)))
  }
  return(g)
}


#' Coverage profile from a bam file across genomic ranges
#'
#' @param bam.file filename of the bam file. The bam file must be sorted and
#'        indexed
#' @param granges GRanges object in which to count the reads. Note: all ranges
#'        must be of the same width
#' @param min.mapq minimal mapping quality. Default is to ignore mapping quality
#' @param reads.collapsed logical indicating whether duplicated reads were
#'        removed and counts were added as _x to the readnames (this is the case
#'        for some of the tools from the Rajewski lab)
#' @param width length of the sequencing reads, by default this will be extracted
#'        from the bam file
#' @return count matrix with rows for each element of the GRanges object
#' @export
coverageBamInGRanges <- function(bam.file, granges, min.mapq=NULL, reads.collapsed=FALSE, width=NULL) {
  # require(GenomicRanges)
  # require(Rsamtools)

  # first check that all granges have the same width
  w = width(granges[1])
  stopifnot(all(width(granges) == w))
  
  seq.names <- as.character(unique(seqnames(granges)))
  seq.names.in.bam <- names(scanBamHeader(bam.file)[[1]]$targets)

  grange.coverage = matrix(0, nrow=length(granges), ncol=w)
  for (seq.name in seq.names) {
    if (seq.name %in% seq.names.in.bam) {
      granges.subset <- granges[seqnames(granges)==seq.name]
      strand(granges.subset) <- "*"
      what = c("pos", "mapq", "qwidth")
      if (reads.collapsed) {
        what = c(what, "qname")
      }
      rds <- scanBam(bam.file,param=ScanBamParam(what=what, which=range(granges.subset)))
      if (!is.null(min.mapq)) {
        mapq.test <- rds[[1]]$mapq >= min.mapq & !is.na(rds[[1]]$mapq)
      } else {
        mapq.test = rep(T, length(rds[[1]]$mapq))
      }
      if (sum(mapq.test) > 0) {
        if (is.null(width)) {
          width = rds[[1]]$qwidth[mapq.test]
        }
        rds.ranges <- GRanges(seq.name, IRanges(start=rds[[1]]$pos[mapq.test], width=width))
        if (reads.collapsed) {
          multiply = as.numeric(sapply(strsplit(rds[[1]]$qname[mapq.test], "_x"), "[", 2))
          select = unlist(lapply(1:length(multiply), function(i) rep(i, multiply[i])))
          rds.ranges = rds.ranges[select]
        }
        # set the seqlength, so the coverage Rle gets the right length
        len = seqlengths(granges)[seq.name]
        if (is.na(len)) {
          len = max(c(end(granges), end(rds.ranges)))
        }
        seqlengths(rds.ranges)[seq.name] = len

        coverage.seq.name <- coverage(rds.ranges)[[1]]
        v = Views(coverage.seq.name, start=start(granges.subset), end=end(granges.subset))
        cvg = t(sapply(v, as.numeric))
        grange.coverage[as.logical(seqnames(granges)==seq.name),] <- cvg
      }
    } 
  }
  # reverse the ones on the minus strand
  minus = as.logical(strand(granges) == "-")
  if (any(minus)) {
    grange.coverage[minus,] = t(apply(grange.coverage[minus,], 1, rev))
  }
  return(grange.coverage)
}


#' Get a count matrix from bam files that match a certain filename pattern
#'
#' @param granges GRanges object in which to count reads
#' @param dir name of the directory that contains the bam files
#' @param pattern regulator expression to identify bam files in the directory
#' @return a count matrix with rows for each element of the GRanges object and
#'         columns for each bam file
#' @export
get.count.matrix <- function(granges, dir, pattern) {

  bam.files = list.files(dir, pattern, full=T)
  counts = matrix(nrow=length(granges), ncol=length(bam.files))
  colnames(counts) = basename(bam.files)

  for (bam.file in bam.files) {
    counts[,bam.file] = countBamInGRanges(bam.file, granges)
  }

  return(counts)
}


#' Get a count matrix from htseq-count files that match a filename pattern
#'
#' @param dir name of the directory that contains the bam files
#' @param pattern regulator expression to identify bam files in the directory
#' @return a count matrix with rows for each element of the GRanges object and
#'         columns for each bam file
#' @export
get.htseq.count.matrix <- function(dir, pattern) {
  files = list.files(dir, pattern, full=T)
  expr = NULL
  for (f in files) {
    x = read.csv(f, sep="\t", header=F)
    colnames(x) = c("ID", basename(f))
    if (is.null(expr)) {
      expr = x
    } else {
      expr = merge(expr, x, all.x=T, all.y=T)
    }
  }
  expr[is.na(expr)] = 0
  rownames(expr) = expr[,"ID"]
  expr = as.matrix(expr[,-match("ID", colnames(expr))])
  return(expr)
}


#' Normalize to reads per kilobase per million sequences
#'
#' @param annotation GRanges object with exons that have a "gene_id" field in
#'        the meta data
#' @param counts read count matrix per gene. The rownames have to be the same as
#'        those used in the gene_id field.
#' @return matrix of RPKM values
#' @export
get.rpkm <- function(annotation, counts) {
  bygene = split(annotation, values(annotation)[,"gene_id"])
  len = sum(width(reduce(bygene)))
  genes = intersect(names(bygene), rownames(counts))
  counts = counts[genes,]
  total = colSums(counts)
  len = len[genes]
  rpkm = counts / rep(total, each=nrow(counts)) / rep(len, ncol(counts)) * 1e9 
  return(rpkm)
}


