#' Data class cnvmat
#' 
#' Class to store breakpoint annotations
#' 
#' @param cnvmat (data.frame): matrix containing average CNV per gene (rows) for each sample (columns)
#' @param genesgr (S4): a GenomicRanges object with genomic feature annotations such as gene coordinates
#' @param cnv (S4) an object of class svcnvio containing data type 'cnv' validated by validate.cnv
#' @param param (list):
#' @return an instance of the class 'genecnv' containing gene level copy number info
#' @export

genecnv <- setClass("genecnv",
                        representation(
                            cnvmat  = "matrix",
                            genesgr = "GRanges",
                            cnv = "svcnvio",
                            param = "list"
                        ))


setMethod("show","genecnv",function(object){
    writeLines(paste("An object of class genecnv from svpluscnv containing gene level CNV data
                \nNumber of samples=",ncol(object@cnvmat),
                "\nAltered genes=",nrow(object@cnvmat)))
})


#' Gene-level CNV
#' 
#' Obtains a gene-level copy number matrix from a segmentation profile.
#'  
#' @param cnv (S4) an object of class svcnvio containing data type 'cnv' initialized by validate.cnv
#' @param genome.v (hg19 or hg38) reference genome version to draw chromosome limits and centromeres
#' @param genesgr (S4) a GenomicRanges object containing gene annotations (if not NULL overides genome.v). It must containg 'strand' and a metadata field 'gene_id' with unique values. Seqnames are expected in the format (chr1, chr2, ...) 
#' @param chrlist (character) list of chromosomes to include chr1, chr2, etc...
#' @param fill.gaps (logical) whether to fill the gaps in the segmentation file using gap neighbour segmean average as log ratio
#' @param verbose (logical) 
#' @return an instance of the class 'genecnv' containing gene level copy number info
#' @keywords CNV, segmentation, genes
#' @export
#' @examples
#' 
#' ## validate input data.frames
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' gene.cnv(cnv)

gene.cnv <- function(cnv, 
                     genome.v="hg19",
                     genesgr=NULL,
                     chrlist=NULL, 
                     fill.gaps=FALSE,
                     verbose=TRUE){

stopifnot(cnv@type == "cnv")
cnv2<-cnv
if(fill.gaps) cnv2 <- segment.gap(cnv2, chrlist=chrlist, verbose=verbose)
cnvdat <- cnv2@data
    
if(is.null(chrlist)) chrlist <- unique(cnvdat$chrom)
chrlist <- chr.sort(chrlist)




if(!is.null(genesgr)){
    if(anyDuplicated(genesgr@elementMetadata$gene_id) > 0) stop("The genesgr provided object contains duplicated gene_id values")
}else{
    genesgr <- get.genesgr(genome.v=genome.v)
}

cnvdat_gr <- with(cnvdat, GRanges(chrom, IRanges::IRanges(start=start, end=end)))

hits <- GenomicAlignments::findOverlaps(genesgr,cnvdat_gr)

overlaps_all <- pintersect(genesgr[queryHits(hits),], cnvdat_gr[subjectHits(hits),])
width_overlap <- width(overlaps_all)

df <- data.table(cnvdat[subjectHits(hits),c("sample","segmean")],genesgr[queryHits(hits)]@elementMetadata$gene_id,width_overlap)
colnames(df) <- c("sample","segmean","gene_id","width")

a <- sapply(unique(cnvdat$sample), 
          function(i) df[sample == i, .(CN=mean(segmean)), by = "gene_id"],
          simplify=FALSE)

newfunc <- function(dfi) { 
    cn<- dfi$CN
    names(cn) <- dfi$gene_id
    return(cn)
    }

b<- lapply(a, function(x) newfunc(x)[genesgr@elementMetadata$gene_id] )
cnvmat <- do.call(cbind,b)

out <- genecnv(
    cnvmat=cnvmat,
    genesgr=genesgr,
    cnv=cnv,
    param=list(genome.v=genome.v,
               chrlist=chrlist, 
               fill.gaps=fill.gaps,
               verbose=verbose
    )
)
return(out)

}




#' Amplifications and deletions
#' 
#' Retrieve amplification and deletion events from a 'genecnv.obj' generated by 'gene.cnv' function
#' 
#' @param genecnv.obj (genecnv) an instance of the class 'genecnv' containing gene level copy number info
#' @param logr.cut (numeric) the log-ratio cutoff above which genes are considered amplified (e.g 2 = 8 copies for amplification and 0.5 copies for deep deletions, in diploid regions)
#' @return (list) A list of lists including amplified.list, amplified.rank, deepdel.list and deepdel.rank
#' @keywords CNV, segmentation, genes
#' @export
#' @examples
#' 
#' ## validate input data.frames
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' genecnv.obj <- gene.cnv(cnv)
#' 
#' geneampdel <- amp.del(genecnv.obj, logr.cut = 2)
#' lapply(geneampdel,head)

amp.del <- function(genecnv.obj, logr.cut=2){
    
    amp_list <- apply(genecnv.obj@cnvmat, 1, function(x) names(which(x >= 2)))
    amp_list <- amp_list[which(unlist(lapply(amp_list,length)) > 0)]
    amp_rank <- sort(unlist(lapply(amp_list,length)),decreasing=TRUE)

    del_list <- apply(genecnv.obj@cnvmat, 1, function(x) names(which(x <= -2)))
    del_list <- del_list[which(unlist(lapply(del_list,length)) > 0)]
    del_rank <- sort(unlist(lapply(del_list,length)),decreasing=TRUE)
    
    return(list(amplified.list = amp_list,
                amplified.rank = amp_rank,
                deepdel.list = del_list,
                deepdel.rank = del_rank))
}

