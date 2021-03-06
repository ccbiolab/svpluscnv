% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate.input.data.r
\name{validate.svc}
\alias{validate.svc}
\title{Initialization of SVC data}
\usage{
validate.svc(sv.df)
}
\arguments{
\item{sv.df}{(data.frame) structural variant table including the following fields: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass}
}
\value{
an instance of the class 'svcnvio' containing SV data derived from SVC data type;  A unique id (uid) column is also added
}
\description{
This function validates and reformats the SV (structural variant) calls input. It is used internaly by 'svpluscnv' functions that require this type of data.
A few formatting rules are enforced:
1) The input must obtain 8 columns in the following order(sample ID, chromosome of origin, strand of origin, position of origin,, chromosome of destination, strand of destination, position of destination, SV class)
2) SV classes accepted: DEL(deletion), DUP(duplication), INS(insertion), TRA(translocation), INV(inversion) and BND(break end)
3) Any variant in which chromosome of origin and destination differ are encoded as TRA (translocation)
4) pos1 < pos2 is enforced for all variants in which chromosome of origin and destination are the same
5) The class BND can be used to operate with complex events as long as both break ends are the same chromosome
}
\examples{

validate.svc(svdat_lung_ccle)
}
\keyword{SV,}
\keyword{structural}
\keyword{variants}
