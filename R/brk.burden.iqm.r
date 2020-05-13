#' Data class break.iqm
#' 
#' Class to store breakpoint annotations in association with genomic features (e.g. gene loci)
#' 
#' @param summary (data.table): the frequency of gains and losses in each defined genomic bin
#' @param brk.mat (numeric): a matrix of genomic bins versus samples
#' @param chrlimits (data.frame): a table containing the chromosome limit coordinates and global genomic coordinates
#' @param plot (graphical): a recorded plot object
#' @param param (list): a list of parametres provided 
#' @return an instance of the class 'cnvfreq' 
#' @export

break.iqm <- setClass("break.iqm", representation(
    summary  = "data.table",
    brk.mat = "matrix",
    chrlimits = "data.table",
    plot = "recordedplot",
    param = "list"
))


setMethod("show","break.iqm",function(object){
    writeLines(paste("An object of class break.iqm from svpluscnv containing the following stats:
                \nNumber of samples=",nrow(object@brk.mat)))
})


#' Evaluates the breakpoint burden based on a instance 'breaks' produced by svpluscnv::scv_breaks or svpluscnv::cnv_breaks. 
#' Breakpoint densities are calculated for each chromosome arm and the inter quantile mean (svpluscnv::IQM) of al chromosome arms is reported for each sample.
#' A Graphical output is generated indicating every sample's arm burden ordered by their IQM. 
#' 
#' @param brk (breaks) An instance of the class 'breaks' obtained from CNV segmentation data (svpluscnv::cnv.breaks) or Structural Variant calls (svpluscnv::svc.breaks).
#' @param sample.col (character) A vector of valid colors. Names must match sample column from 'brk'. If null a gradiant color based on breakpoint burden IQM will be used. 
#' @param chr.lim (data.frame) 3 column table (chrom, begin, end) indicating the chromosome most distal coordinates with coverage. Also returned by the function svpluscnv::chromosome.limit.coords.
#' @param genome.v (hg19 or hg38) reference genome version to draw chromosome limits and centromeres
#' @param min.arm.size (numeric) minimum size in base pairs for a chromosome arm to be included in the analysis. Size will be calculated based on the 'genome.v' centromere location (excluding centromere bands). Chromosome start and en locations can be provided in 'chr.lim'.
#' @param bp.unit (numeric) The genomic size unit in base pairs to report brekpoint densities. This parameter is also used for the y axis of the plot. 
#' @param verbose (logical) whether to return internal messages
#' @return an instance of the class 'cnvfreq' and optionally a plot into open device
#' @keywords structural variants, mutational burden, chromosomal instability
#' @export
#' @examples
#' 
#' # initialize CNV data
#' svc <- validate.svc(nbl_svdat)
#' 
#' # obtain CNV breakpoints
#' brk <- cnv.breaks(cnv)
#' 
#' brk.burden.iqm(brk)


brk.burden.iqm <- function(brk,
                           sample.col = NULL,
                           min.arm.size = 2e7,
                           bp.unit=1e7,
                           genome.v="hg19",
                           chr.lim= NULL,
                           plot=TRUE,
                           verbose=TRUE){
    
stopifnot(isS4(brk))

if(genome.v %in% c("GRCh37","hg19")){ 
    bands <- remove.factors(GRCh37.bands)
}else if(genome.v %in% c("GRCh38","hg38")){
    bands <- remove.factors(GRCh38.bands)
}else{stop("Genome version not provided")}
    
centromeres_start <- bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"start"]
centromeres_end <- bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"end"]
chromosome_start <- sapply(as.character(unique(bands$chr)), function(i) min( bands$start[which(bands$chr == i)] ))
chromosome_end <- sapply(as.character(unique(bands$chr)), function(i) max( bands$end[which(bands$chr == i)] ))
names(chromosome_start) <-  names(chromosome_end) <-  names(centromeres_start) <-  names(centromeres_end) <- paste("chr",bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"chr"],sep="")
    
if(!is.null(chr.lim)){
    centromeres_start <- centromeres_start[chr.lim$chrom]
    centromeres_end <- centromeres_end[chr.lim$chrom]
    chromosome_start <- chromosome_start[chr.lim$chrom]
    chromosome_end <- chromosome_end[chr.lim$chrom]
    chromosome_end[] <- chr.lim$end
    chromosome_start[] <- chr.lim$begin
}else{
    chr.lim<- data.table(names(centromeres_start),centromeres_start,chromosome_end)
    colnames(chr.lim) <- c("chrom","begin","end")
}
    
mapped_p <- names(which(centromeres_start -chromosome_start > min.arm.size))
mapped_q <- names(which(chromosome_end -centromeres_end > min.arm.size))
p.arm.df <- data.frame(mapped_p,chromosome_start[mapped_p],centromeres_start[mapped_p])
q.arm.df <- data.frame(mapped_q,centromeres_end[mapped_q],chromosome_end[mapped_q])
colnames(p.arm.df) <- colnames(q.arm.df) <- c("chrom","start","end")

p.arm.gr <- with(p.arm.df,GRanges(chrom, IRanges(start=start, end=end)))
q.arm.gr <- with(q.arm.df,GRanges(chrom, IRanges(start=start, end=end)))
    
breaks.gr <- with(brk@breaks, GRanges(chrom, IRanges(start=pos,end=pos)))

p.hits <- GenomicAlignments::findOverlaps(breaks.gr,p.arm.gr)
q.hits <- GenomicAlignments::findOverlaps(breaks.gr,q.arm.gr)

p.armname <- paste(mapped_p,"p",sep="")
q.armname <- paste(mapped_q,"q",sep="")
arm.size <- c(p.arm.df$end -p.arm.df$start, q.arm.df$end -q.arm.df$start)
names(arm.size) <- c(p.armname,q.armname)

template <- rep(0, length(c(p.armname,q.armname)))
names(template) <- c(p.armname,q.armname)
arm.brk.dens <- sapply(unique(brk@breaks$sample), function(i) template, simplify=FALSE)

p.hits.info <- data.table(brk@breaks$sample[queryHits(p.hits)],p.armname[subjectHits(p.hits)])
q.hits.info <- data.table(brk@breaks$sample[queryHits(q.hits)],q.armname[subjectHits(q.hits)])

total.brk <- list()
for(sample.id  in names(arm.brk.dens)){
    input <- c(table(p.hits.info$V2[which(p.hits.info$V1 == sample.id)]),
               table(q.hits.info$V2[which(q.hits.info$V1 == sample.id)]))
    arm.brk.dens[[sample.id]][names(input)] <- input*bp.unit/arm.size[names(input)]
    total.brk[[sample.id]] <- sum(input)
}

arm.brk.iqm <- log10(1+sort(unlist(lapply(arm.brk.dens,IQM))))

if(is.null(sample.col)){
    sample.col <- rep("green",length(unique(brk@breaks$sample)))
    names(sample.col) <- unique(brk@breaks$sample)
    sample.col.tmp <- map2color(arm.brk.iqm,pal <- colorRampPalette(c("darkgreen","orange","red"))(256))
    names(sample.col.tmp) <- names(arm.brk.iqm)
    sample.col[names(sample.col.tmp)] <- sample.col.tmp
    }

datavector <- log10(1+unlist(lapply(arm.brk.dens[names(arm.brk.iqm)],sort)))
datacolor <- unlist(sapply(names(arm.brk.iqm), function(i) rep(sample.col[i], length(template)),simplify=FALSE))
names(datacolor) <- names(datavector)

# plot 

npoints <- length(template)
plot(datavector,pch=20,xaxt='n',yaxt='n',col="white",xlab="",ylab='',
     xaxt='n',bty='n',xlim=c(100,length(datavector)-100))
altcol<-"grey95"
for(i in 1:length(arm.brk.dens)){ 
    rect((i-1)*npoints,-10,i*npoints,50,col=altcol,border=NA)
    if(altcol == "grey95"){ altcol <- "grey85"
    }else{altcol <- "grey95"}
}
abline(h=seq(-2,6,0.5),lty=1,lwd=.2,col="black")

points(datavector,pch=20,cex=0.3,col=datacolor)
axis(2,labels=sprintf("%.2f",10^(seq(-2,4,0.5))-1),at=seq(-2,4,0.5),las=3,family="Courier",font=1,line=0,cex.axis=1.2,las=1)

mtext(paste("log10(1+breaks/",bp.unit,")",sep=""),side=2,line=4,cex=1.3)
lines(seq(npoints/2,length(datavector),length(datavector)/length(arm.brk.iqm)),log2(1+arm.brk.iqm) )
# save plot
p <- recordPlot()

nbreaks <- table(brk@breaks$sample)[names(arm.brk.iqm)]
nbreaks.map <- unlist(total.brk)[names(arm.brk.iqm)]
brk.dens <- (nbreaks.map*bp.unit/sum(arm.size))[names(arm.brk.iqm)]

summary <- data.table(names(arm.brk.iqm), 
                      arm.brk.iqm, 
                      sample.col[names(arm.brk.iqm)], 
                      as.numeric(nbreaks), 
                      nbreaks.map, 
                      brk.dens )
colnames(summary) <- c("sample","brk.iqm","color","total breaks","nbreaks mapped","overal density")


return(break.iqm(
    summary = summary,
    brk.mat = do.call(rbind,arm.brk.dens),
    chrlimits = chr.lim,
    plot=p,
    param = list(
        min.arm.size= min.arm.size,
        bp.unit=bp.unit,
        genome.v= genome.v,
        verbose= verbose
    )
)
)
}



