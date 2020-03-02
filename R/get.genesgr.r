#' Genes GRanges 
#' 
#' Retrieves a GRanges object containinng gene annotations for an specified genome version 
#' 
#' @param genome.v (hg19 or GRCh37 and hg38 or GRCh38) reference genome version to retrieve gene annotations 
#' @param chrlist (character)  
#' @return a GRanges class object from the specified human genome version 
#' @keywords CNV, segmentation, genes
#' @export
#' @examples
#' 
#' get.genesgr(genome.v = "hg19",chrlist=NULL)
#' 

get.genesgr<- function(genome.v="hg19",chrlist=NULL){

    if(genome.v %in% c("hg19","GRCh37")){
        genesgr = GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene, columns="gene_id")
    }else if(genome.v %in% c("hg38","GRCh38")){
        genesgr = GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg38.knownGene, columns="gene_id")
    }else{stop("Unspecified, or non available genome")}
    
    if(is.null(chrlist)){ 
        chrlist <- paste("chr",c(1:22,"X","Y"),sep="")
    }

    err <- capture.output(
        genesgr@elementMetadata$gene_id <- mapIds(org.Hs.eg.db, genesgr@elementMetadata$gene_id,  'SYMBOL','ENTREZID'),
        type="message")
    
    genesgr <- genesgr[which(!is.na(genesgr$gene_id))]
    genesgr <- genesgr[which(lapply(genesgr@elementMetadata$gene_id,length) > 0)]
    genesgr <- genesgr[which(as.character(genesgr@seqnames) %in% chrlist)]

    return(genesgr)
}
