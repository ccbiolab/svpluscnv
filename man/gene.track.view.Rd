% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene.track.view.r
\name{gene.track.view}
\alias{gene.track.view}
\title{Gene track visualization}
\usage{
gene.track.view(
  chrom = NULL,
  start = NULL,
  stop = NULL,
  symbol = NULL,
  upstr = NULL,
  dnstr = NULL,
  genome.v = "hg19",
  cex.text = 0.6,
  addtext = TRUE,
  plot = TRUE,
  summary = TRUE,
  ...
)
}
\arguments{
\item{chrom}{(character) Chromosome (e.g. chr9)}

\item{start}{(numeric) Genomic coordinate from specified chromosome to start plotting}

\item{stop}{(numeric) Genomic coordinate from specified chromosome to stop plotting}

\item{symbol}{(character) Gene acceoted hgnc symbol to retrieve coordinates and area plotting ()}

\item{upstr}{(numeric) Distance upstream specified gene to extend the area plotted}

\item{dnstr}{(numeric) Distance downstream specified gene to extend the area plotted}

\item{genome.v}{(character) Reference genome version to draw chromosome limits and centromeres (hg19 or hg38)}

\item{cex.text}{(numeric) The magnification to be used for transcript RefSeq text added}

\item{addtext}{(logic) Whether to include transcript RefSeq ids in the plot}

\item{plot}{(logic) Whether to generate plot in open device}

\item{summary}{(logic) Whether to produce a data.table output with transcript information}

\item{...}{Additional graphical parameters}
}
\value{
A data.frame with gene isoform annotations and/or plot into open device
}
\description{
Creates a track visualization of a genomic region defined by gene boundaries or custom provided
}
\examples{

# obtain the coordinates of a desired genomic regionbased on a known gene locus 
refSeqGene <- gene.symbol.info(refseq_hg19,"PTPRD")
chrom <- refSeqGene$chrom
start <- refSeqGene$start - 150000;
stop <- refSeqGene$stop + 50000;

gene.track.view(symbol="PTPRD", genome.v="hg19")
}
\keyword{CNV,}
\keyword{segmentation}
