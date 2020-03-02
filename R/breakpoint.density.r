#' Data class breaks
#' 
#' Class to store breakpoint annotations in association with genomic features (e.g. gene loci)
#' 
#' @param breaks (data.table): the breakpoint info containing data.table, this will be occupied by the CNV segmentation data in the case of cnv.break.annot or SV for sv.break.annot. Unique random string rownames are added to the returned breaks data.frame.
#' @param burden (numeric): a vector containing the total number of breakpoints in each sample 
#' @param param (list): a list of parametres provided 
#' @return an instance of the class 'breaks' containing breakpoint and breakpoint burden information
#' @export
breaks <- setClass("breaks", representation(
                            breaks  = "data.table",
                            burden = "numeric",
                            param = "list"
                        ))


setMethod("show","breaks",function(object){
    writeLines(paste("An object of class breaks from svpluscnv containing",object@param$datatype,"breakpoints:
                \nNumber of samples=",length(object@burden),
                "\nTotal number of breakpoints =",nrow(object@breaks)))
})


#' Identify CNV breakpoints
#' 
#' Identify CNV breakpoints filtered by the change in copy number log-ratio between contiguous segments
#' 
#' @param cnv (S4) an object of class svcnvio containing data type 'cnv' initialized by validate.cnv
#' @param fc.pct (numeric) copy number change between 2 consecutive segments: i.e (default) cutoff = 0.2 represents a fold change of 0.8 or 1.2
#' @param break.width (numeric) the maximum distance between a segment end and the subsequent segment start positions beyond which breakpoints are discarded
#' @param min.cnv.size (numeric) The minimun segment size (in base pairs) to include in the analysis 
#' @param min.num.probes (numeric) The minimun number of probes per segment to include in the analysis 
#' @param chrlist (character) list of chromosomes to include chr1, chr2, etc...
#' @param low.cov (data.frame) a data.frame (chr, start, end) indicating low coverage regions to exclude from the analysis
#' @param clean.brk (numeric) identical breakpoints across multiple samples tend to be artifacts; remove breaks > N 
#' @param verbose (logical) whether to return  
#' @return an instance of the class 'breaks' containing breakpoint and breakpoint burden information
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' # initialized CNV data
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' cnv.breaks(cnv)
#' 


cnv.breaks <- function(cnv,
                       fc.pct = 0.2,
                       break.width = 10000,
                       min.cnv.size = NULL,
                       min.num.probes = NULL,
                       chrlist = NULL,
                       low.cov = NULL, 
                       clean.brk = NULL,
                       verbose = TRUE){
    

stopifnot(cnv@type == "cnv")
cnvdat <- cnv@data
    
if(is.null(chrlist)) chrlist <- unique(cnvdat$chrom)
chrlist <- chr.sort(chrlist)
    
brk.burden <- rep(0,length(unique(cnvdat$sample)))
names(brk.burden) <- unique(cnvdat$sample)
    
if(!is.null(min.cnv.size)) cnvdat <- cnvdat[which(cnvdat$end - cnvdat$start >= min.cnv.size),]
if(!is.null(min.num.probes)) cnvdat <- cnvdat[which(cnvdat$probes  >= min.num.probes),]
    
lastrow <- nrow(cnvdat)
pos <- round(apply(cbind(cnvdat[2:(lastrow),"start"], cnvdat[1:(lastrow-1),"end"]),1,mean))
chrom <- cnvdat[2:(lastrow),"chrom"]
sample <- cnvdat[2:(lastrow),"sample"]
width <- cnvdat[2:(lastrow),"start"] - cnvdat[1:(lastrow-1),"end"]
FC <-  (2^cnvdat[1:(lastrow-1),"segmean"]) / (2^cnvdat[2:lastrow,"segmean"])
uid <- paste("brk_",createRandomString(nrow(cnvdat)-1,8),sep="")
breakpoints <- data.table(sample,chrom,pos,width,FC,uid)
colnames(breakpoints) <- c("sample","chrom","pos","width","FC","uid")

break_idx <- c(which( log2(FC) >= log2(1+fc.pct)),which( log2(FC) < log2(1 - fc.pct)))
    
samechr <- which(apply(cbind(cnvdat[1:(lastrow-1),"chrom"],cnvdat[2:(lastrow),"chrom"]),1,anyDuplicated) == 2)

samesample <-  which(apply(cbind(cnvdat[1:(lastrow-1),"sample"],cnvdat[2:(lastrow),"sample"]),1,anyDuplicated) == 2)

if(is.null(break.width)) break.width <- Inf
brwidthin <- which(width < break.width)
    
breakpoints <- breakpoints[Reduce(intersect, list(break_idx,samechr,samesample,brwidthin)),]
    

if(!is.null(low.cov)){
    message("Filtering breakpoints in low coverage regiomns")
    colnames(low.cov) <- c("chrom","start","end")
    low_cov_GR = with(low.cov, GRanges(chrom, IRanges(start=start, end=end)))
    breakpoints_GR = with(breakpoints, GRanges(chrom, IRanges(start=start, end=end)))
    overlapgr <- GenomicAlignments::findOverlaps(breakpoints_GR,low_cov_GR,ignore.strand=TRUE)
    breakpoints <- breakpoints[setdiff(1:nrow(breakpoints),queryHits(overlapgr)),]
}
    
if(!is.null(clean.brk)){
    breakids <- unite(breakpoints[,c(2:4)],"newcol")$newcol
    breakids.freq <- sort(table(breakids),decreasing=TRUE)
    breakpoints <- breakpoints[which(breakids %in% names(which(breakids.freq < clean.brk))),]
}
    
brk.burden.sub <- table(breakpoints$sample)
brk.burden[names(brk.burden.sub)] <- brk.burden.sub
    
return(breaks(breaks=breakpoints,
            burden=brk.burden,
            param=list(
                datatype=cnv@type,
                fc.pct = fc.pct,
                min.cnv.size = min.cnv.size,
                min.num.probes=min.num.probes,
                low.cov=low.cov, 
                clean.brk=clean.brk
            )
            )
       )
}



#' Identify SVC breakpoints
#' 
#' Transform structural varian (SVC) data.frame into a 'breaks' object 
#' 
#' @param svc (S4) an object of class svcnvio containing data type 'svc' initialized by validate.svc
#' @param low.cov (data.table) a data.table (chrom, start, end) indicating low coverage regions to exclude from the analysis
#' @return an instance of the class 'breaks' containing breakpoint and breakpoint burden information
#' @keywords Structural variants
#' @export
#' @examples
#' 
#' ## Obtain breakpoints from SV calls data
#' svc <- validate.svc(svdat_lung_ccle)
#' 
#' svc.breaks(svc)



svc.breaks <- function(svc, low.cov=NULL){
    
stopifnot(svc@type == "svc")
svcdat <- svc@data
    
brk.burden <- rep(0,length(unique(svcdat$sample)))
names(brk.burden) <- unique(svcdat$sample)


uid<- paste("brk_",createRandomString(nrow(svcdat)*2,8),sep="")
svcdat.breaks <- data.table(c(svcdat$sample,svcdat$sample),
                           c(svcdat$chrom1,svcdat$chrom2),
                           c(svcdat$pos1,svcdat$pos2),
                           c(svcdat$strand1,svcdat$strand2),
                           c(svcdat$svclass,svcdat$svclass),
                           c(svcdat$uid,svcdat$uid),
                           uid)

colnames(svcdat.breaks) <- c("sample","chrom","pos","strand","svclass","svcuid","uid")
if(!is.null(low.cov)){
    low.cov.df <- data.table(low.cov[,1:3])
    colnames(low.cov.df) <- c("chrom","start","end")
    
    svc_ranges <- with(svcdat.breaks, GRanges(chrom, IRanges(start=pos, end=pos)))
    low.cov_ranges <- with(low.cov.df, GRanges(chrom, IRanges(start=start, end=end)))
    
    low.cov_ranges = GenomicAlignments::findOverlaps(svc_ranges,low.cov_ranges)
    
    svcdat.breaks <- svcdat.breaks[which(!svcdat.breaks$id %in% queryHits(low.cov_ranges)),]
}else{
    svcdat.breaks <- svcdat.breaks
}

brk.burden.sub <- table(svcdat.breaks$sample)
brk.burden[names(brk.burden.sub)] <- brk.burden.sub
    

return(breaks(breaks=svcdat.breaks,
            burden=brk.burden,
            param=list(
                datatype=svc@type,
                low.cov=low.cov
            )
        )
    )
    
}




#' Breakpoint density map
#' 
#' Generating a genomic map based on a defined bin size and sliding window and counts the number of breakpoints mapped onto each bin. This function is used internally by svpluscnv::shattered.regions and svpluscnv::shattered.regions.cnv
#' 
#' @param brk (breaks) An instance of the class 'breaks' obtained from CNV segmentation data (svpluscnv::cnv.breaks) or Structural Variant calls (svpluscnv::svc.breaks). 
#' @param chr.lim (data.frame) 3 column table (chrom, begin, end) indicating the chromosome most distal coordinates with coverage. Also returned by the function svpluscnv::chromosome.limit.coords.
#' @param genome.v (hg19 or hg38) reference genome version to draw chromosome limits and centromeres
#' @param window.size (numeric) size in megabases of the genmome bin onto which breakpoints will be mapped 
#' @param slide.size (numeric) size in megabases of the sliding genomic window; if slide.size < window.size the genomic bins will overlap
#' @param verbose (logical) whether to return internal messages
#' @return a matrix of samples (rows) and genomic bins (cols) qith the number of breakpoints mapped in heach cell
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' # initialize CNV data
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' # obtain CNV breakpoints
#' brk <- cnv.breaks(cnv)
#' 
#' break.density(brk)


break.density <- function(brk, 
                          chr.lim=NULL, 
                          genome.v = "hg19",
                          window.size = 10, 
                         slide.size=2,
                         verbose=TRUE){
if(is.null(chr.lim)){
    chr.lim<- d3gb.chr.lim(genome.v=genome.v)
}else{
    stopifnot(ncol(chr.lim) == 3)   
}
    
    chr.begin <- chr.lim$begin
    chr.end <- chr.lim$end
    names(chr.begin) <- names(chr.end) <- chr.lim$chrom
    
  # make sure both chr.lim and breaks have same chromosome names 
  seqnames <- intersect(chr.lim$chrom,brk@breaks$chr)
  stopifnot(length(seqnames) > 0) 
  
  # a template vector to save breakpoint counts 
  templatevector <- brk@burden
  templatevector[]<-0
  
  WS <- window.size * 1e+6
  SS <- slide.size * 1e+6
  offset <- window.size/slide.size
  
    chrlist <- chr.sort(chr.lim$chrom)
  
  # count breaks for each chromosome for each fragment
  fragment <- list()
  for(chr in  chrlist){

    if(verbose) cat("\r",chr)

    chr_breaks <- brk@breaks[which(brk@breaks$chrom == chr),]
    frag <- seq(chr.begin[chr],chr.end[chr]+SS,SS)
    
    for(i in (1+offset):length(frag)){
      start <- frag[i - offset]
      stop <- frag[i]
      fragment[[paste(chr,start,stop)]] <- templatevector
      break.position <- chr_breaks$pos
      res_bp <- table(chr_breaks[intersect(which(break.position > start),which(break.position < stop)),"sample"])
      fragment[[paste(chr,start,stop)]][names(res_bp)] <- res_bp
    }
  }
  if(verbose) cat("Done!\n")

  return( do.call(cbind,fragment))
  
}




#' Breakpoint matching
#' 
#' Match common breakpoints from two different datasets or data types based on their co-localization in the genome. 
#' 
#' @param brk1 (S4) an object of class breaks as returned by `svc.breaks` and `cnv.breaks`
#' @param brk2 (S4) an object of class breaks as returned by `svc.breaks` and `cnv.breaks` to compare against brk1
#' @param maxgap (numeric) distance (base pairs) limit for nreakpoints to be consider colocalized 
#' @param plot (logical) whether to plot into open device
#' @param verbose (logical) whether to return internal messages
#' @return an object containing co-localizing breakpoints from two input 'breaks'  
#' @keywords CNV, SV, genomic breakpoints
#' @export
#' @examples
#' 
#' # initialize CNV and SVC data
#' cnv <- validate.cnv(segdat_lung_ccle)
#' svc <- validate.svc(svdat_lung_ccle)
#' 
#' ## Obtain breakpoints from CNV and SVC
#' brk1 <- cnv.breaks(cnv)
#' brk2 <- svc.breaks(svc)
#' 
#' common.brk <- match.breaks(brk1, brk2)
#' 



match.breaks <- function(brk1, 
                         brk2, 
                         maxgap=100000,
                         verbose=FALSE,
                         plot=TRUE){
    
    common_samples <- intersect(names(brk1@burden),names(brk2@burden))
    stopifnot(length(common_samples) > 0) 
    
    brk1_match <- brk2_match <- res <- list()
    for(id in common_samples){
        
        brk1_i <- brk1@breaks[which(brk1@breaks$sample == id),]
        brk_ranges1 <- with(brk1_i, GRanges(chrom, IRanges(start=pos, end=pos)))
        
        brk2_i <- brk2@breaks[which(brk2@breaks$sample == id),]
        brk_ranges2 <- with(brk2_i, GRanges(chrom, IRanges(start=pos, end=pos)))
        
        
        options(warn=-1)
        seg_seg = GenomicAlignments::findOverlaps(brk_ranges1, brk_ranges2, maxgap=maxgap)
        options(warn=0)
        
        brk_match1 <- sort(unique(queryHits(seg_seg)))
        brk_match2 <- sort(unique(subjectHits(seg_seg)))
        
        res[[id]] <- data.table(id,length(brk_match1), nrow(brk1_i), length(brk_match2), nrow(brk2_i))
        colnames(res[[id]]) <- c("sample","matched.brk1", "total.brk1", "matched.brk2", "total.brk2")
        
        brk1_match[[id]] <- brk1_i[brk_match1,]
        brk2_match[[id]] <- brk2_i[brk_match2,]
    }

    restab <- do.call(rbind,res)
    
    if(plot){
      def.par <- par(no.readonly = TRUE)
      par(mfrow=c(2,1))
      restab <- restab[order(restab$total.brk2)]
      m2 <- sprintf("%.1f",100*mean(na.omit(restab$matched.brk2/restab$total.brk2))) 
      barplot(rbind(restab$matched.brk2, restab$total.brk2 - restab$matched.brk2),
              border=NA,las=2,xlab="",horiz=FALSE,cex.main=.7,cex.names=.4, 
              names=restab$sample,ylab="#samples" )
      legend("top",paste(brk2@param$datatype," breaks matched by ",
                         brk1@param$datatype,
                         " breaks\n","Average = ",m2,"%",sep=""),bty='n')
      grid(ny=NULL,nx=NA)

            restab <- restab[order(restab$total.brk1)]
      m2 <- sprintf("%.1f",100*mean(na.omit(restab$matched.brk1/restab$total.brk1))) 
      barplot(rbind(restab$matched.brk1, restab$total.brk1 - restab$matched.brk1),
              border=NA,las=2,xlab="",horiz=FALSE,cex.main=.7,cex.names=.4, 
              names=restab$sample,ylab="#samples")
      legend("top",paste(brk1@param$datatype,
                         " breaks matched by ",brk2@param$datatype,
                         " breaks\n","Average = ",m2,"%",sep=""),bty='n')
      grid(ny=NULL,nx=NA)
      par(def.par)
    }
    
    return(list(
        brk1_match = do.call(rbind,brk1_match),
        brk2_match = do.call(rbind,brk2_match),
        restab= restab))
}

