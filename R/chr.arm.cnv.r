#' Chromosome arm mean CNV
#'
#' Obtains a matrix with the weighted average CN per chromosome arm 
#' @param cnv (S4) an object of class svcnvio containing data type 'cnv' validated by validate.cnv
#' @param genome.v (character) (hg19 or hg38) reference genome version to draw chromosome limits and centromeres
#' @param verbose (logical) whether to return internal messages
#' @return a matrix of chromosome arms (rows) versus samples (cols) with average segment logRs per cell
#' @keywords CNV, segmentation, chromosome arm
#' @export 
#' @examples
#' 
#' # initialize CNV data
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' arm_mat <- chr.arm.cnv(cnv, genome.v="hg19")
#' dim(arm_mat)


chr.arm.cnv <- function(cnv,
                    genome.v="hg19",
                    verbose=FALSE){
 
    stopifnot(cnv@type == "cnv")
    cnvdat <- cnv@data
    
    if(genome.v %in% c("GRCh37","hg19")){ 
        bands <- GRCh37.bands
    }else if(genome.v %in% c("GRCh38","hg38")){ 
        bands <- GRCh38.bands
    }else{stop("Genome version not provided")}
  
    centromeres_start <- bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"start"]
    centromeres_end <- bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"end"]
    names(centromeres_start) <-  names(centromeres_end) <- paste("chr",bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"chr"],sep="")
  
    chr.lim <- chromosome.limit.coords(cnv)
    chrarms <- rbind(cbind(chr.lim$begin,centromeres_start[chr.lim$chrom]),cbind(centromeres_end[chr.lim$chrom],chr.lim$end))
    chrarms <- data.table(rownames(chrarms),chrarms,c(paste(chr.lim$chrom,"p",sep=""), paste(chr.lim$chrom,"q",sep="")))
    colnames(chrarms) <- c("chrom","start","end","arm")
  
    chrarms <- chrarms[which(chrarms$end -chrarms$start > 0),]

    chrarmsGR <- with(chrarms,GRanges(chrom, IRanges(start=start, end=end)))

    cnvdat_gr <- with(cnvdat, GRanges(chrom, IRanges(start=start, end=end)))
    hits <- GenomicAlignments::findOverlaps(chrarmsGR,cnvdat_gr)
  
    armcnvmat <- matrix(ncol=length(unique(cnvdat$sample)), nrow=nrow(chrarms) )
    colnames(armcnvmat) <- unique(cnvdat$sample)
    rownames(armcnvmat) <- chrarms$arm

    for(i in unique(queryHits(hits))){ 
        arm <- chrarms[i,"arm"][[1]]
    
        if(verbose) cat("\r",arm)
    
        armdf <- cnvdat[subjectHits(hits)[which(queryHits(hits) == i)],]
        armdf[which(armdf$start < chrarms[i,"start"]),"start"] <- chrarms[i,"start"]
        armdf[which(armdf$end > chrarms[i,"end"]),"end"] <- chrarms[i,"end"]

        arm.width <- armdf$end - armdf$start
        armdf <- data.table(armdf,arm.width)
        armlength <- aggregate(arm.width~sample,armdf,sum)[,2]
        names(armlength) <- aggregate(arm.width~sample,armdf,sum)[,1]
        part <- armdf$segmean * armdf$arm.width / armlength[armdf$sample]
    
        armdf <- data.table(armdf,arm.width,part,armlength[armdf$sample])
    
        meanArmSegment <- aggregate(part~sample,armdf,sum)
    
        num <-  as.numeric(meanArmSegment[,2])
        names(num) <- as.character(meanArmSegment[,1])
        armcnvmat[arm,names(num)] <- num
    }
    return(armcnvmat)
}



