#' SV integrated visualization
#' 
#' Integrated visualization of SVC and CNV data for defined genomic locations. CNV and SVC data is overlayed into a sample-based track visualization map.
#' 
#' @param cnv (S4) an object of class svcnvio containing data type 'cnv' initialized by validate.cnv
#' @param svc (S4) an object of class svcnvio containing data type 'svc' initialized by validate.svc
#' @param chrom (character) chromosome (e.g chr9)
#' @param start (numeric) genomic coordinate from specified chromosome to start plotting
#' @param stop (numeric) genomic coordinate from specified chromosome to stop plotting
#' @param sampleids (character) a vector containing a list of sample ids represented in svc and/or cnv objects to be plotted
#' @param cnvlim (numeric) limits for color coding of background CNV log-ratios. Use to modify the CNV color contrast at different levels.
#' @param addlegend (character) One of 'sv' (show SV type legend), 'cnv' (show CNV background color legend) or 'both'.
#' @param cex.legend (numeric) The cex values for each legend
#' @param interval (numeric) The axis interval in base pairs
#' @param addtext (character) a vector indicating what SV types should include text labels indicating brakpoint partners genomic locations. The added labels are point breakpoint locations outside the plot area. (e.g. c("TRA","INV") )
#' @param cex.text (numeric) The magnification to be used for SV text info added
#' @param plot (logic) whether to produce a graphical output
#' @param summary (logic) whether the function shoud return CNV segment 'segbrk' and SV 'svbrk' breakpoints tabular output
#' @param ... additional plot parameters from graphics plot function 
#' @return a data.frame with CNV and SVN breakpoint annotations and/or plot into open device
#' @keywords structural variant, CNV, segmentation
#' @export
#' @examples
#' 
#' ## validate input data.frames
#' cnv <- validate.cnv(segdat_lung_ccle)
#' svc <- validate.svc(svdat_lung_ccle)
#'
#' # obtain the coordinates of a desired genomic regionbased on a known gene locus 
#' refSeqGene <- gene.symbol.info(refseq_hg19,"PTPRD")
#' start <- refSeqGene$start - 150000;
#' stop <- refSeqGene$stop+ 50000;
#' chrom <- refSeqGene$chrom
#' 
#' sv.model.view(cnv, svc, chrom, start, stop)
#' 



sv.model.view <- function(cnv, svc, chrom, start, stop, 
                          sampleids=NULL,
                          cnvlim=c(-2,2), 
                          addlegend='both',
                          cex.legend=1,
                          interval=NULL,
                          addtext=NULL,
                          cex.text=.8,
                          plot=TRUE,
                          summary=TRUE,
                          ...){
    

 stopifnot(!is.null(chrom) && !is.null(start) && !is.null(stop))

    stopifnot(cnv@type == "cnv")
    cnvdat <- cnv@data
    
    stopifnot(svc@type == "svc")
    svcdat <- svc@data
    
    if(!is.null(sampleids)){
        missing.samples <- setdiff(sampleids,c(svcdat$sample,cnvdat$sample))
        if(length(missing.samples) == length(unique(sampleids))){
            stop("None of the samples provided were found in 'sv' and 'cnv' input data!")
        }else if(length(missing.samples) > 0){
            warning(paste("The following samples provided are not found in 'sv' and 'cnv' input data:", paste(missing.samples,collapse=" "),sep=" "))
        }
        svcdat<-svcdat[which(svcdat$sample %in% intersect(sampleids,svcdat$sample)),]
        cnvdat<-cnvdat[which(cnvdat$sample %in% intersect(sampleids,cnvdat$sample)),]
    }

    genegr <- with(data.frame(chrom,start,stop), GRanges(chrom, IRanges(start=start, end=stop))) 
    
    # Find samples with SV breaks within defined genomic region
    sv1gr = with(svcdat, GRanges(chrom1, IRanges(start=pos1, end=pos1))) 
    sv2gr = with(svcdat, GRanges(chrom2, IRanges(start=pos2, end=pos2))) 
    
    sv_hits1 = GenomicAlignments::findOverlaps(sv1gr,genegr)
    sv_hits2 = GenomicAlignments::findOverlaps(sv2gr,genegr)
    svtab <- svcdat[sort(unique(c(queryHits(sv_hits1),queryHits(sv_hits2)))),]
    svBreakSamples <- unique(svtab$sample)
    if(length(svBreakSamples) == 0) warning("Thre is no SV breakpoints in the defined genomic region")
        
    # obtain SVs for plotting with different colors for each svclass
    svcolormap = setNames(c("blue", "red", "orange", "black", "green","grey20"), 
                   c("DEL", "DUP", "INV", "TRA", "INS", "BND"))
    svcolor <- svcolormap[svtab$svclass]
    svtab_plot <- data.table(svtab,svcolor)
    svtab_plot_seg <- svtab_plot[which(svtab_plot$svclass != "TRA")]
    svtab_plot_tra <- svtab_plot[which(svtab_plot$svclass == "TRA")]
    
    # Find samples with CNV segment breaks within defined genomic region
    seg1br  = with(cnvdat, GRanges(chrom, IRanges(start=start, end=start))) 
    seg2br  = with(cnvdat, GRanges(chrom, IRanges(start=end, end=end))) 
    seg_hits1 = GenomicAlignments::findOverlaps(seg1br,genegr)
    seg_hits2 = GenomicAlignments::findOverlaps(seg2br,genegr)
    segBreakSamples <- unique(cnvdat[sort(unique(c(queryHits(seg_hits1),queryHits(seg_hits2))))]$sample)
    if(length(segBreakSamples) == 0) warning("Thre is no CNV segment breakpoints in the defined genomic region")    
    segbrk <- cnvdat[sort(unique(c(queryHits(seg_hits1),queryHits(seg_hits2))))]
    
    if(plot==TRUE){
        # Find overlap between all CNV segments and the defined genomic region for plotting

        seggr <- with(cnvdat, GRanges(chrom, IRanges(start=start, end=end))) 
        hits_seg = GenomicAlignments::findOverlaps(seggr,genegr)
        seg_plot <- cnvdat[queryHits(hits_seg)]
        segcolor <- map2color(seg_plot$segmean,
                  pal=colorRampPalette(c("lightblue","white","salmon"))(256),
                  limits=cnvlim)
        seg_plot <- data.table(seg_plot,segcolor)
    
        if(!is.null(sampleids)){
            sample_order <- 1:length(sampleids)
            names(sample_order) <- sampleids
        }else{
            sample_order <- 1:length(unique(c(svBreakSamples,segBreakSamples)))
            names(sample_order) <- unique(c(svBreakSamples,segBreakSamples))
        }
    
        if(!is.null(addlegend)){
            plot_ylim <- length(sample_order)*10/100+length(sample_order)
            legend_ypos <- plot_ylim - length(sample_order)*3/100 
            if(length(sample_order) < 10) plot_ylim <- length(sample_order) +1
        }else{
            plot_ylim <- length(sample_order)
        }
        
        plot(x=NULL,y=NULL,xlim=range(c(start,stop)),ylim=range(c(0,plot_ylim)),
             xaxt='n',yaxt='n',xlab='',ylab='',bty='n', ...)
    
        mtext(side=2,at=sample_order-0.5,text=names(sample_order),las=2,line = 0.5, ...)
    
        for(sid in names(sample_order)){
            ypos <- sample_order[sid]
            polygon(rbind(
                c(start-1e7,ypos+0.02),
                c(start-1e7,ypos-0.98),
                c(stop+1e7,ypos-0.98),
                c(stop+1e7,ypos+0.02)),
                col=rep(c("grey80","grey80"),length(sample_order))[ypos],border=NA)
        }
        
        for(sid in names(sample_order)){
            seg_sample_plot <- seg_plot[which(seg_plot$sample == sid),]
            ypos <- sample_order[sid]
            for(i in 1:nrow(seg_sample_plot)){
                polygon(rbind(
                    c(seg_sample_plot[i]$start,ypos),
                    c(seg_sample_plot[i]$start,ypos-1),
                    c(seg_sample_plot[i]$end,ypos-1),
                    c(seg_sample_plot[i]$end,ypos)
                ),col=seg_sample_plot[i]$segcolor,border=NA)
            }
        }
    
    
        for(sid in unique(svtab_plot_tra$sample)){
            svtab_plot_tra_i <- svtab_plot_tra[which(svtab_plot_tra$sample == sid),]
            ypos <- sample_order[sid]
            addrnorm <- rep(c(0,0.3,-0.3,0.1,-0.1,0.2,-0.2),nrow(svtab_plot_tra_i))
            for(i in 1:nrow(svtab_plot_tra_i)){
                if(svtab_plot_tra_i[i]$chrom2 == chrom){ 
                    points(svtab_plot_tra_i[i]$pos2,ypos-0.5+addrnorm[i],pch=10)
                    lines(c(svtab_plot_tra_i[i]$pos2,svtab_plot_tra_i[i]$pos2),c(ypos,ypos-1),lwd=1,lty=3)
                    if("TRA" %in% addtext){
                        text(svtab_plot_tra_i[i]$pos2,ypos-0.5+addrnorm[i],
                             paste("  ",svtab_plot_tra_i[i]$chrom1,":",svtab_plot_tra_i[i]$pos1,sep=""),
                             pos=4,offset=0,cex=cex.text)
                    }
                }            
                if(svtab_plot_tra_i[i,"chrom1"] == chrom){
                    points(svtab_plot_tra_i[i]$pos1,ypos-0.5+addrnorm[i],pch=10)
                    lines(c(svtab_plot_tra_i[i]$pos1,svtab_plot_tra_i[i]$pos1),c(ypos,ypos-1),lwd=1,lty=3)
                    if("TRA" %in% addtext) {
                        text(svtab_plot_tra_i[i]$pos1,ypos-0.5+addrnorm[i],
                             paste("  ",svtab_plot_tra_i[i]$chrom2,":",svtab_plot_tra_i[i]$pos2,sep=""),
                             pos=4,offset=0,cex=cex.text)
                    }
                }
            }
        }
        
        for(sid in unique(svtab_plot_seg$sample)){
            svtab_plot_seg_i <- svtab_plot_seg[which(svtab_plot_seg$sample == sid)]
            ypos <- sample_order[sid]
            addrnorm <- rep(c(0,0.2,-0.2,0.1,-0.1,0.3,-0.3),nrow(svtab_plot_seg_i))
            for(i in 1:nrow(svtab_plot_seg_i)){
                polygon(rbind(
                    c(svtab_plot_seg_i[i]$pos1,ypos-0.4-addrnorm[i]),
                    c(svtab_plot_seg_i[i]$pos1,ypos-0.6-addrnorm[i]),
                    c(svtab_plot_seg_i[i]$pos2,ypos-0.6-addrnorm[i]),
                    c(svtab_plot_seg_i[i]$pos2,ypos-0.4-addrnorm[i])
                ),col=NA,border=svtab_plot_seg_i[i]$svcolor)
                
                if(svtab_plot_seg_i[i]$svclass %in% addtext){
                    if(svtab_plot_seg_i[i]$pos1 < start){
                        text(start,ypos-0.5-addrnorm[i],
                             paste("<-",svtab_plot_seg_i[i]$pos1,sep=""),
                             pos=4,offset=0,cex=cex.text)
                    }
                    if(svtab_plot_seg_i[i]$pos2 > stop){
                        text(stop,ypos-0.5-addrnorm[i],
                             paste(svtab_plot_seg_i[i]$pos2,"->",sep=""),
                             pos=2,offset=0,cex=cex.text)
                    }
                }
            }
        }
        
        if(is.null(interval)) interval <- round((stop - start)/5000) * 1000
        xlabs <- seq(floor(start/10000)*10000, ceiling(stop/10000)*10000,interval)
        axis(1, at = xlabs,labels=TRUE, lwd.ticks=1.5, pos=0,...)

        if(is.null(cex.legend)) cex.legend <- 1
        
        if(addlegend %in% c("sv","both")) {
            fillx <- c("white", "white", "white", "white", "white",NA)
            borderx <- c("blue", "red","orange","green","grey20",NA)
            pchx <- c(NA,NA,NA,NA,NA,10)
            names(fillx) <- names(borderx) <- names(pchx) <- c("DEL", "DUP", "INV","INS","BND", "TRA")
            svclassin <- sort(unique(svtab_plot$svclass))
            legend(x= start, y =legend_ypos+0.2, legend = svclassin, bg=NA,
                   bty = "n", fill = fillx[svclassin], border=borderx[svclassin], 
                   pch = pchx[svclassin], horiz = TRUE, x.intersp=0.2, cex = cex.legend)
        }
        if(addlegend %in% c("cnv","both")) {
            colkey(colorRampPalette(c("lightblue","white","salmon"))(256),clim = c(-4,4),side=3,add=TRUE,side.clab=1,length=0.5,shift=0.2,lwd.ticks = 2, dist = -0.12)
        }
    }
    if(summary){
        return(list(svbrk=svcdat,segbrk=segbrk))
    }
}



