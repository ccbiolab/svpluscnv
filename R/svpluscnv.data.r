#' Data class refSeqDat
#' 
#' Class to store refseq data from UCSC containing exon level info for known transcripts
#' 
#' @param data (data.table): transcript information
#' @param exonStarts (list): every transcript exonic start position
#' @param exonStarts (list): every transcript exonic end position
#' @param genome.v (character): the genome version encoding transcript data
#' @return an instance of the class 'refSeqDat' containing transcript exonic coordinates
#' @export

refSeqDat <- setClass("refSeqDat", representation(
    data  = "data.table",
    exonStarts = "list",
    exonEnds= "list",
    genome.v="character"
))

setMethod("show","refSeqDat",function(object){
    writeLines(paste("An object of class refSeqDat from svpluscnv with ",nrow(object@data),"transcipts from",object@genome.v,"genome version"))
})


#'
#' Return coordinates of an specified gene
#' 
#' @param object (refSeqDat) An object of class refSeqDat containing gene transcript mapping. svpluscnv includes two selfloaded objects: refseq_hg19 & refseq_hg38
#' @param symbol (character) a valid HGNC gene symbol included in the refseq object
#' @export
#' @docType methods
#' @return A list containing chr, start, end coordinates
#' @rdname gene.symbol.info-methods

setGeneric("gene.symbol.info", function(object, symbol) standardGeneric("gene.symbol.info"))

#' @rdname gene.symbol.info-methods
setMethod("gene.symbol.info", "refSeqDat", function(object, symbol){
    DT <- object@data[which(object@data$name2 == symbol)]
    return(list(
        chrom = unique(DT$chrom),
        start = min(DT$txStart),
        stop = max(DT$txEnd)
    ))
})


utils::globalVariables(c("refseq_hg19", "refseq_hg38"))

#' Reference transcript and exon annotations for hg19 
#' 
#' refSeq annotations for hg19 version from UCSC (http://genome.ucsc.edu/cgi-bin/hgTables)
#' 
#' @name refseq_hg19
#' @docType data
#' @keywords genes, transcripts, exons
#' 
"refseq_hg19"


#' Reference transcript and exon annotations for hg38
#' 
#' refSeq annotations for hg38 version from UCSC (http://genome.ucsc.edu/cgi-bin/hgTables)
#' 
#' @name refseq_hg38
#' @docType data
#' @keywords genes, transcripts, exons
#' 
"refseq_hg38"


utils::globalVariables(c("segdat_lung_ccle", "svdat_lung_ccle","cnv_blacklist_regions","nbl_segdat","nbl_svdat"))

#' Lung CCLE CNV data
#' 
#' CCLE CNV segmentation data from LUNG tissue cell lines (DepMap): https://depmap.org/portal/download/
#' @name segdat_lung_ccle
#' @docType data
#' @keywords CNV segmentation
"segdat_lung_ccle"

#' Lung CCLE SVC data
#' 
#' CCLE translocation data from LUNG tissue cell lines (DepMap): https://depmap.org/portal/download/
#' @name svdat_lung_ccle
#' @docType data
#' @keywords SVs
"svdat_lung_ccle"

#' Low coverage regions
#' 
#' @name cnv_blacklist_regions
#' @docType data
#' @keywords CNV segmentation
"cnv_blacklist_regions"

#' TARGET Neuroblastoma CNV
#' 
#' TARGET CNV segmentation: https://target-data.nci.nih.gov/
#' @name nbl_segdat
#' @docType data
#' @keywords CNV segmentation, SVs
"nbl_segdat"

#' TARGET Neuroblastoma SVC
#' 
#' TARGET CGI structural variants: https://target-data.nci.nih.gov/
#' 
#' @name nbl_svdat
#' @docType data
#' @keywords  SVs
"nbl_svdat"



