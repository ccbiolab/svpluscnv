#' Data class cnvfreq
#' 
#' Class to store breakpoint annotations in association with genomic features (e.g. gene loci)
#' 
#' @param freqsum (data.table): the frequency of gains and losses in each defined genomic bin
#' @param chrlimits (data.frame): a table containing the chromosome limit coordinates and global genomic coordinates
#' @param bin.mat (numeric): a matrix of genomic bins versus samples
#' @param plot (graphical): a recorded plot object
#' @param param (list): a list of parametres provided 
#' @return an instance of the class 'cnvfreq' 
#' @export

cnvfreq <- setClass("cnvfreq", representation(
    freqsum  = "data.table",
    chrlimits = "data.frame",
    bin.mat = "matrix",
    plot = "recordedplot",
    param = "list"
))


setMethod("show","cnvfreq",function(object){
    writeLines(paste("An object of class cnvfreq from svpluscnv containing the following stats:
                \nNumber of samples=",ncol(object@bin.mat),
                "\nNumber of genomic bins =",nrow(object@bin.mat)))
})


#' CNV frequency map
#' 
#' Creates a map of CNVs using genome binning and plots CNV frequency across the genome. This function optionally returns text, graphical or both outputs.
#' 
#' @param cnv (S4) an object of class svcnvio containing data type 'cnv' initialized by validate.cnv
#' @param fc.pct (numeric) percentage CNV gain/loss for a segment to be considered changed (i.e. 0.2 = 20 percent change 0.8 < segmean && segmean > 1.2)
#' @param genome.v (character) (hg19 or h38) reference genome version to draw chromosome limits and centromeres
#' @param ploidy (logical) whether to apply ploidy correction; the function med.segmean will be used to obtain each sample's ploidy logR then this value substracted to each sample's logR values
#' @param g.bin (numeric) size in megabases of the genmome bin to compute break density 
#' @param sampleids (character) vector containing list of samples to include in plot. if set to NULL, all samples in the input will be used
#' @param cex.axis,cex.lab,label.line (numeric) plot parameters
#' @param verbose (logical) whether to return internal messages
#' @return an instance of the class 'cnvfreq' and optionally a plot into open device
#' @keywords CNV, segmentation, plot
#' @export
#' @examples
#' 
#' ## validate input data.frame
#' cnv <- validate.cnv(nbl_segdat)
#' 
#' cnv.freq(cnv, genome.v = "hg19")

cnv.freq <- function(cnv,
                     fc.pct= 0.2,
                     genome.v= "hg19",
                     ploidy=FALSE,
                     g.bin= 1,
                     sampleids=NULL,
                     cex.axis= 1,
                     cex.lab= 1,
                     label.line= -1.2,
                     verbose=TRUE){
  
stopifnot(cnv@type == "cnv")
cnvdat <- cnv@data
    
if(!is.null(sampleids)) cnvdat <- cnvdat[which(cnvdat$sample %in% sampleids),]
  
if(ploidy){
    ploidy_val <- med.segmean(cnv)
    cnvdat$segmean <- cnvdat$segmean - ploidy_val[cnvdat$sample]
}

stopifnot(genome.v %in% c("hg19","hg38","GRCh37","GRCh38"))
if(genome.v %in% c("hg19","GRCh37")){ bands <- GRCh37.bands
}else if(genome.v %in% c("hg38","GRCh38")){ bands <- GRCh38.bands}

centromeres <- bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"start"]
names(centromeres) <- paste("chr",bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"chr"],sep="")

# define chromosome mapped limits and the global genome coordinates for each chromosome start
chrlimits <-   chromosome.limit.coords(cnv)
offset <- c(0,vapply(seq_len(nrow(chrlimits)-1), 
                     function(i) sum(chrlimits[seq_len(i),"end"]) + i*g.bin,1))
chrlabelpos <- offset + chrlimits$end/2
chrlimits <- data.frame(offset,as.data.frame(chrlimits),chrlabelpos)
rownames(chrlimits) <- chrlimits$chrom

g.bin.mb <- g.bin*1e6
  
if(verbose) message("Generating binned genome map ")

chrbins <- list()

for(chr in rownames(chrlimits)){
    seqpos <- seq(chrlimits[chr,"begin"],chrlimits[chr,"end"]+g.bin.mb,g.bin.mb)
    ranges <-  t( vapply(seq(2,length(seqpos)), function(i) c(seqpos[i-1],seqpos[i]),double(2)) )
    chrcol<- rep(chr,length(seqpos)-1)
    segcol_del <- segcol_gain <- rep("grey",length(chrcol))
    segcol_del[which(ranges[,2] <= centromeres[chr])] <- "lightblue"
    segcol_del[which(ranges[,2] > centromeres[chr])] <- "blue"
    segcol_gain[which(ranges[,2] <= centromeres[chr])] <- "salmon"
    segcol_gain[which(ranges[,2] > centromeres[chr])] <- "red"
    chrbins[[chr]] <- data.table(chrcol,ranges,segcol_del,segcol_gain)
}
  
  chrbins.df <- do.call(rbind,unname(chrbins) )
  chrbins.df<- data.table(chrbins.df,unite(chrbins.df[,c(1,2,3)],paste)$paste)
  colnames(chrbins.df) <- c("chr","start","end","segcol_del","segcol_gain","binid")


  if(verbose) message("Calculating mean segmean per genomic bin")
  # find overlaps between bins and cnv segments
  binsGR <- with(chrbins.df, GRanges(chr, IRanges(start=start, end=end)))
  segGR <- with(cnvdat, GRanges(chrom, IRanges(start=start, end=end)))
  hits <-GenomicAlignments::findOverlaps(binsGR,segGR)
  
  outmat <- matrix(ncol=length(unique(cnvdat$sample)),nrow=nrow(chrbins.df))
  colnames(outmat) <- unique(cnvdat$sample)
  rownames(outmat) <- chrbins.df$binid

  for(i in seq_len(nrow(chrbins.df)) ){
    segtmp<- cnvdat[subjectHits(hits)[which(queryHits(hits) == i)],]
    if(nrow(segtmp)>0){
      a <- aggregate(segmean~sample,segtmp, sum)  
      outmat[i,a$sample]<- a$segmean
    }else{
      outmat[i,a$sample]<- NA
    }
  }

  if(verbose) message("Calculating gain/loss frequencies per genomic bin")
  outmat[which(is.na(outmat),arr.ind=TRUE)] <- 0
  
  outmat_gain<-outmat_loss<-outmat
  outmat_gain[]<-outmat_loss[]<-0
  nsamples <- ncol(outmat_gain)
  
  outmat_gain[which(outmat > log2(1+fc.pct), arr.ind=TRUE)] <-  1
  outmat_loss[which(outmat < log2(1-fc.pct), arr.ind=TRUE)] <-  1
  freq.gains <- apply(outmat_gain,1,sum)/nsamples
  freq.loss <- apply(outmat_loss,1,sum)/nsamples
  
    plot.end<- chrlimits$offset[nrow(chrlimits)]+chrlimits$end[nrow(chrlimits)]
    bin.loc <- chrlimits[chrbins.df[names(freq.gains),on="binid"]$chr,"offset"] + chrbins.df[names(freq.gains),,on="binid"]$start

    if(verbose) message("Plotting ...")
    altcols <- rep(c(rgb(0.1,0.1,0.1,alpha=0.1),rgb(0.8,0.8,0.8,alpha=0.1)),12)
    altcols2<- rep(c(rgb(0.1,0.1,0.1,alpha=1),rgb(0.4,0.4,0.4,alpha=1)),12)
  
    plot(x=NULL,y=NULL,xlim=c(0,plot.end),ylim=c(-1,1),bty='n',xaxt='n',yaxt='n',xlab="",ylab="")
    for(i in seq_len(length(chrlimits$offset)) ) rect( chrlimits$offset[i],-1,chrlimits$offset[i]+chrlimits$end[i],1, col=altcols[i],border=NA )
    points(bin.loc,freq.gains,type='h',col=chrbins.df$segcol_gain)
    points(bin.loc,-freq.loss,type='h',col=chrbins.df$segcol_del)
    lines(c(0,plot.end),c(0,0),col="lightgrey")
    lines(c(0,plot.end),c(0.5,0.5),col="lightgrey",lty=3)
    lines(c(0,plot.end),c(-0.5,-0.5),col="lightgrey",lty=3)
    mtext(gsub("chr","",rownames(chrlimits))[seq(1,nrow(chrlimits),2)],side=1,at=chrlimits$chrlabelpos[seq(1,nrow(chrlimits),2)],las=1,col=altcols2[seq(1,nrow(chrlimits),2)],line=label.line,cex=cex.lab)
    mtext(gsub("chr","",rownames(chrlimits))[seq(2,nrow(chrlimits),2)],side=3,at=chrlimits$chrlabelpos[seq(2,nrow(chrlimits),2)],las=1,col=altcols2[seq(2,nrow(chrlimits),2)],line=label.line,cex=cex.lab)
    mtext("Frequency",side=4,line=1)
    mtext("#samples",side=2,line=1)
    axis(4,c(100,50,0,50,100),at=c(-1,-0.5,0,0.5,1),las=1,pos=plot.end, cex.axis=cex.axis)
    axis(2,c(nsamples,round(nsamples/2),0,round(nsamples/2),nsamples),at=c(-1,-0.5,0,0.5,1),las=1, pos=0, cex.axis=cex.axis)
    p <- recordPlot()


summary <- data.table(chrbins.df[,c("chr","start","end")],bin.loc,freq.gains,freq.loss)

return(cnvfreq(
            freqsum = summary,
            bin.mat = outmat,
            chrlimits = chrlimits,
            plot=p,
            param = list(
                fc.pct= fc.pct,
                genome.v= genome.v,
                g.bin= g.bin,
                sampleids=sampleids,
                cex.axis= cex.axis,
                cex.lab= cex.lab,
                label.line= label.line   
                )
            )
       )
}

