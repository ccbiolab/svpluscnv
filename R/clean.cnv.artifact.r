#' CNV segmentation gap filling
#' 
#' Fills the gaps in a segmentation data.frame. Chromosome limits are defined for the complete segmentation dataset then segments fill the missing terminal regions. 
#' The CN log-ratio of the added segments is set to the average of the closest neighbours in each sample.
#' 
#' @param cnv (S4) an object of class svcnvio containing data type 'cnv' initialized by validate.cnv
#' @param minsize (numeric) the minimum gap size required to fill the gap
#' @param chrlist (character) list of chromosomes to include chr1, chr2, etc...
#' @param verbose (logical) whether to return internal messages
#' @return a data.frame containing CNV data
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' ## validate input data.frames
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' cnv2 <- segment.gap(cnv)
#' cnv2

segment.gap <- function(cnv,
                        minsize=5000,
                        chrlist=NULL,
                        verbose=FALSE){
    
    stopifnot(cnv@type == "cnv")
    cnvdat <- cnv@data
    
    chrlims <- chromosome.limit.coords(cnv)
    if(is.null(chrlist)) chrlist <- chrlims$chrom
    
    chrlims_df<- data.frame(chrlims)
    rownames(chrlims_df) <- chrlims_df$chrom
    
    cnvdat_df <- data.frame(cnvdat)
    
    if(verbose){
        message("Filling gaps is the segmentation data.frame")
        pb <- txtProgressBar(style=3)
        cc <-0
        tot <- nrow(cnvdat_df)
    }
    newsegments<-list()
    if(cnvdat_df[1,"start"] > chrlims_df[cnvdat_df[1,"chrom"],"begin"]){
        newsegments[["1"]] <- data.frame(cnvdat_df[1,c("sample","chrom")],chrlims_df[cnvdat_df[1,"chrom"],"begin"],cnvdat_df[1,"start"]-1,0,cnvdat_df[1,"segmean"])
    }
    
    for(i in 2:nrow(cnvdat_df)){
        if(cnvdat_df[i,"chrom"] == cnvdat_df[i-1,"chrom"] ){
            if( cnvdat_df[i,"start"] - cnvdat_df[i-1,"end"] > minsize){
                newsegments[[as.character(i)]] <- data.frame(cnvdat_df[i,c("sample","chrom")],cnvdat_df[i-1,"end"]+1,cnvdat_df[i,"start"]-1,0,mean(cnvdat_df[c(i,i-1),"segmean"]) )
            }
        }else{
            if(cnvdat_df[i,"start"] > chrlims_df[cnvdat_df[i,"chrom"],"begin"]){ 
                newsegments[[as.character(i)]] <- data.frame(cnvdat_df[i,c("sample","chrom")],chrlims_df[cnvdat_df[i,"chrom"],"begin"],cnvdat_df[i,"start"]-1,0,cnvdat_df[i,"segmean"])
            }
            if(cnvdat_df[i-1,"end"] < chrlims_df[cnvdat_df[i-1,"chrom"],"end"]){
                newsegments[[as.character(i)]] <- data.frame(cnvdat_df[i-1,c("sample","chrom")],cnvdat_df[i-1,"end"]+1,chrlims_df[cnvdat_df[i-1,"chrom"],"end"],0,cnvdat_df[i,"segmean"])
            }
        }
        if(verbose) cc <- cc+1
        if(verbose) setTxtProgressBar(pb, cc/tot)
    }
    if(cnvdat_df[i,"end"] < chrlims_df[cnvdat_df[i,"chrom"],"end"]) newsegments[[as.character(i)]] <- data.frame(cnvdat_df[i,c("sample","chrom")],cnvdat_df[i,"end"]+1,chrlims_df[cnvdat_df[i,"chrom"],"end"],0,cnvdat_df[i-1,"segmean"])
    if(verbose) close(pb)
    
    newsegments <- lapply(newsegments, setNames, colnames(cnvdat_df)[1:6])
    
    segout <- rbind(cnvdat_df[,1:6], do.call(rbind,newsegments))
    out <- validate.cnv(segout)
    
    return(out)
}


#' CNV artifact detection and filtering
#' 
#' Detects identical or near-identical CNV segments across multiple samples susceptible of representing common variants or technical artifacts. Then those segments CNV log-ratio is replaced by the flanking segments average 
#' 
#' @param cnv (S4) an object of class svcnvio containing data type 'cnv' validated by validate.cnv
#' @param n.reps (numeric) number of samples with identical segment to consider artifact
#' @param cnv.size (numeric) only smaller segments will be modified in the cnv data.frame
#' @param pc.overlap (numeric) minimun percentage overlap for a pair of segments to be consider identical 
#' @param fill.gaps (logical) whether to fill gaps from the segmentaed file after filtering artifacts
#' @param minsize (numeric) the minimum gap size required to fill the gap. Only used if 'fill.gaps=TRUE'
#' @param verbose  (logical) whether to print internal messages
#' @return a data.frame containing CNV data
#' @keywords CNV, segmentation, filter
#' @export
#' @examples
#' 
#' ## validate input data.frame
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' cnvcl <- clean.cnv.artifact(cnv)
#' cnvcl

clean.cnv.artifact<- function(cnv,
                              n.reps=4,
                              cnv.size=2000000,
                              pc.overlap=0.99,
                              fill.gaps=TRUE,
                              minsize=5000,
                              verbose=TRUE){

    stopifnot(cnv@type == "cnv")
    cnvdat <- cnv@data
    
all_artifacts_l <-list()

cnvdat_short <- cnvdat[which(cnvdat$end - cnvdat$start < cnv.size),]

for(chr in unique(cnvdat$chrom)){

  if(verbose) cat("\r",chr)

  segchr <- cnvdat_short[which(cnvdat_short$chrom == chr),]
  segchr.gr <- with(segchr, GRanges('chrom', IRanges(start=start, end=end)))
  hits = GenomicAlignments::findOverlaps(segchr.gr,segchr.gr)
  overlaps <- pintersect(segchr.gr[queryHits(hits)], segchr.gr[subjectHits(hits)])
  
  percentOverlapA <- width(overlaps) / width(segchr.gr[queryHits(hits)])
  percentOverlapB <- width(overlaps) / width(segchr.gr[subjectHits(hits)])
  hits_p <- as.data.frame(hits[intersect(which(percentOverlapA >= pc.overlap),which(percentOverlapB >= pc.overlap)),])
  reps <- aggregate(subjectHits~queryHits,hits_p,paste,simplify=FALSE)
  reps_list <- reps$subjectHits
  names(reps_list) <- reps$queryHits
  reps_list_collapse <- lapply(lapply(reps_list,sort),paste,collapse=" ")
  groups_a <- table(unlist(reps_list_collapse))
  all_artifacts <- as.numeric(unlist(strsplit(names(which(groups_a > n.reps))," ")))
  all_artifacts_l[[chr]] <- segchr[all_artifacts,]
}

all_artifacts <- do.call(rbind,unname(all_artifacts_l))
toremove <- unite(all_artifacts, "newcol", c("sample","chrom","start","end"), remove=FALSE,sep=":")$newcol
allsegids <- unite(cnvdat, "newcol", c("sample","chrom","start","end"), remove=FALSE,sep=":")$newcol
cnvdat_clean <- svcnvio(data = cnvdat[which(!allsegids %in% toremove),],type = "cnv")

if(fill.gaps){ 
    segclean_fill <-  segment.gap(cnvdat_clean, minsize=minsize, verbose=verbose)
    return(segclean_fill)
  }else{
    return(cnvdat_clean)
  }

}



