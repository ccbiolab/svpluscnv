#' Data class svcnvio
#' 
#' Class to store CNV segmentation data
#' 
#' @param data (data.table): cnv or svc data.table to be validated by 'validate.cnv' or 'validate.svc' respectivelly
#' @param type (character): the data type  "cnv" or "svc" defined by "validate.cnv" or "validate.svc" respectivelly
#' @seealso Additional data format information in the man pages of validate.cnv and validate.svc
#' @return an instance of the class 'svcnvio' containing SV data derived from CNV or SVC data types;  A unique id (uid) column is also added
#' @export

svcnvio <- setClass("svcnvio", representation(
    data  = "data.table",
    type = "character"
))

setMethod("show","svcnvio",function(object){
    writeLines(paste("An object of class svcnvio from svpluscnv storing",object@type,"data from",length(unique(object@data$sample)),"samples"))
})

#' Initialization of SVC data
#' 
#' This function validates and reformats the SV (structural variant) calls input. It is used internaly by 'svpluscnv' functions that require this type of data.
#' A few formatting rules are enforced:
#' 1) The input must obtain 8 columns in the following order(sample ID, chromosome of origin, strand of origin, position of origin,, chromosome of destination, strand of destination, position of destination, SV class)
#' 2) SV classes accepted: DEL(deletion), DUP(duplication), INS(insertion), TRA(translocation), INV(inversion) and BND(break end)
#' 3) Any variant in which chromosome of origin and destination differ are encoded as TRA (translocation)
#' 4) pos1 < pos2 is enforced for all variants in which chromosome of origin and destination are the same
#' 5) The class BND can be used to operate with complex events as long as both break ends are the same chromosome
#' 
#' @param sv.df (data.frame) structural variant table including the following fields: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass
#' @return an instance of the class 'svcnvio' containing SV data derived from SVC data type;  A unique id (uid) column is also added
#' @keywords SV, structural variants
#' @export
#' @examples
#' 
#' validate.svc(svdat_lung_ccle)


validate.svc <- function(sv.df){
    
    stopifnot(ncol(sv.df) >= 8)
    uid <-  paste("svc_",createRandomString(nrow(sv.df),8),sep="")
    svc <- data.table(remove.factors(sv.df[,1:8]),uid)
    
    colnames(svc) <- c("sample","chrom1","pos1","strand1","chrom2","pos2","strand2","svclass","uid")
    if(length(grep("chr",svc[1]$chrom1)) == 0) svc$chrom1 <- paste("chr",svc$chrom1,sep="")
    if(length(grep("chr",svc[1]$chrom2)) == 0) svc$chrom2 <- paste("chr",svc$chrom2,sep="")
    
    stopifnot(is.numeric(svc$pos1))
    stopifnot(is.numeric(svc$pos2))
    stopifnot(is.character(svc$chrom1))
    stopifnot(is.character(svc$chrom2))
    stopifnot(is.character(svc$sample))
    
    svc[grep("INV",svc$svclass)]$svclass <- "INV"
    svc[grep("DUP",svc$svclass)]$svclass <- "DUP"
    
    extrachr <- which(unlist(lapply(apply(svc[,c("chrom1","chrom2")],1,unique),length)) == 2) 
    svc[extrachr]$svclass <- "TRA"
    
    wrong_class <- setdiff(unique(svc$svclass),c("DEL","DUP","TRA","INV","INS","BND"))
    try(if(length(wrong_class) > 0) message(paste("SV classes not accepted:", paste(wrong_class,collapse=","), "will be set as BND") ))
    svc[which(!svc$svclass %in% c("DEL","DUP","TRA","INV","INS","BND"))]$svclass <- "BND"
    
    # ensure that pos1 is upstream pos2
    intrachr <- which(unlist(lapply(apply(svc[,c("chrom1","chrom2")],1,unique),length)) == 1) 
    intrachr_rev <- intersect(which(svc$pos2 -svc$pos1 < 0),intrachr)
    
    
    if(length(intrachr_rev) > 0){
        svcrev <- svc[intrachr_rev,c(1,2,6,7,5,3,4,8,9)]
        colnames(svcrev) <- c("sample","chrom1","pos1","strand1","chrom2","pos2","strand2","svclass","uid")
        svc <- rbind(svcrev,svc[setdiff(1:nrow(svc),intrachr_rev)])
    }
    
    stopifnot(nrow(svc) > 0)
    
    return(svcnvio(
        data=svc,
        type="svc"
    ))
    
}

#' Chromosome ordering
#' 
#' A function to order a list of chromosomes 
#' 
#' @param chrlist (character): a vector containing chromosome names (chr1, chr2...chrX,chrY  ) 
#' @return a character vector of sorted chromosomes
#' @keywords CNV, segmentation, genes
#' @export
#' @examples
#' 
#' chrlist <- paste("chr",c("X","Y",sample(1:22)),sep="")
#' chr_sorted <- chr.sort(chrlist)


chr.sort <- function(chrlist){ 
    chrunique <- sort(gsub("chr","",unique(chrlist)))
    chrsort <- paste("chr",chrunique[suppressWarnings(order(as.numeric(chrunique) ))],sep="")
    return(chrsort)
}


#' Initialization of CNV data
#' 
#' This function validates and reformats the CNV segmentation data type containing copy number log-ratios. It is used internaly by 'svpluscnv' functions that require this type of data.
#' 
#' @param cnv.df (data.frame) segmentation data with at least 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @return an instance of the class 'svcnvio' containing segmentation data derived from CNV data type;  A unique id (uid) column is also added
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' validate.cnv(segdat_lung_ccle)


validate.cnv <- function(cnv.df){
    
    stopifnot(ncol(cnv.df) >= 6)
    uid <-  paste("cnv_",createRandomString(nrow(cnv.df),8),sep="")
    cnvdat <- data.table(cnv.df[,1:6],uid)
    
    colnames(cnvdat) <- c("sample","chrom","start","end","probes","segmean","uid")
    if(length(grep("chr",cnvdat[1,2])) == 0) cnvdat[,"chrom"] <- paste("chr",cnvdat$chrom,sep="")
    stopifnot(is.numeric(cnvdat$start))
    stopifnot(is.numeric(cnvdat$end))
    stopifnot(is.numeric(cnvdat$segmean))
    stopifnot(is.character(cnvdat$sample))
    stopifnot(is.character(cnvdat$chrom))
    
    chrlist <- chr.sort(unique(cnvdat$chrom))
        
    cnvdat <- cnvdat[order(cnvdat$start),]
    cnvdat <- cnvdat[order(match(cnvdat$chrom, chrlist)),]
    cnvdat <- cnvdat[order(cnvdat$sample),]
    
    stopifnot(nrow(cnvdat) > 0)
    
    return(svcnvio(
        data=cnvdat,
        type="cnv"
    ))
    
}


#' Chromosome limit map
#' 
#' Obtain chromosome start and end positions based on mapped regions from CNV segmentation data
#' 
#' @param cnv (S4) an object of class svcnvio containing data type 'cnv' initialized by validate.cnv
#' @keywords CNV, segmentation, mapping
#' @return data.table indicating start and end mapped positions of each chromosome
#' @export
#' @examples
#' 
#' ## validate input data.frame
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' chr.lim <- chromosome.limit.coords(cnv)

chromosome.limit.coords <- function(cnv){
    
    stopifnot(cnv@type == "cnv")
    cnvdat <- cnv@data
    
    chrlist <- chr.sort(unique(cnvdat$chrom))
    chrmin <- chrmax <- list()
    for(chr in chrlist){
        if(chr %in% cnvdat$chrom){
            chrmin[[chr]] <- min(cnvdat[which(cnvdat$chrom == chr)]$start) 
            chrmax[[chr]] <- max(cnvdat[which(cnvdat$chrom == chr)]$end) 
        }
    }
    begin <- unlist(chrmin)
    end <- unlist(chrmax)
    chr.lim <- data.table(chrlist,begin,end)
    colnames(chr.lim) <- c("chrom","begin","end")
    return(chr.lim)
}




