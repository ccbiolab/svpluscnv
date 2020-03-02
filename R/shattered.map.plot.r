#' Shattered regions genomic map
#' 
#' Plots a genome wide map of shattered region frequencies
#' 
#' @param chromo.regs.obj (chromo.regs) An object of class chromo.regs 
#' @param conf (character) either 'hc' for high confidence objects or else all included
#' @param genome.v (character)  reference genome version to draw chromosome limits and centromeres either hg19 or hg38 accepted
#' @param freq.cut the value to draw an horizontal line; use 'freq.p.test' to obtain a threshold for statisticaly significant hot spots 
#' @param add.legend the position of the legend in the plot; if null, no legend will be draw
#' @return a plot into open device
#' @keywords chromosome shattering, genome map
#' @export
#' @examples
#' 
#' 
#' ## validate input data.frames
#' cnv <- validate.cnv(segdat_lung_ccle)
#' svc <- validate.svc(svdat_lung_ccle)
#' 
#' ## obtain shattered regions
#' chromo.regs.obj <- shattered.regions(cnv,svc)
#' 
#' shattered.map.plot(chromo.regs.obj)

shattered.map.plot <- function(chromo.regs.obj,
                          conf="hc",
                          genome.v = "hg19",
                          freq.cut=NULL,
                          add.legend="top"){


    if(genome.v %in% c("hg19","GRCh37")){ bands <- GRCh37.bands
    }else if(genome.v %in% c("hg38","GRCh38")){ bands <- GRCh38.bands}
    
    centromeres <- bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"start"]
    names(centromeres) <- paste("chr",bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"chr"],sep="")

    chrlengths <- vapply(unique(bands$chr), function(i) max(bands$end[which(bands$chr == i)]), 1)
    names(chrlengths) <- paste("chr",unique(bands$chr),sep="")

    chrlist <- unique(do.call(rbind,strsplit(colnames(chromo.regs.obj@high.density.regions)," "))[,1])
    stopifnot( length(which(!chrlist %in% names(chrlengths))) == 0 )
    
    if(conf == "hc") {
        highDensitiRegionsFreq <- apply(chromo.regs.obj@high.density.regions.hc,2,sum)
    }else{
        highDensitiRegionsFreq <- apply(chromo.regs.obj@high.density.regions,2,sum)
    }

p_chrcols  <- rep(c("salmon4","salmon4"),12)
q_chrcols  <- rep(c("salmon","salmon"),12)
names(p_chrcols) <- names(q_chrcols) <- chrlist
chrom <- do.call(rbind,strsplit(names(highDensitiRegionsFreq)," "))[,1]
coloresBarplot  <- rep("white",length(chrom))

parm <- which(as.numeric(do.call(rbind,strsplit(names(highDensitiRegionsFreq)," "))[,3]) - centromeres[chrom] > 0)
qarm <- which(as.numeric(do.call(rbind,strsplit(names(highDensitiRegionsFreq)," "))[,2]) - centromeres[chrom] < 0)
coloresBarplot[parm] <- p_chrcols[names(parm)]
coloresBarplot[qarm] <- q_chrcols[names(qarm)]


axislab <- chrstarts<-  chrend <- chrlengths[chrlist]
tab <- data.table(do.call(rbind,strsplit(names(highDensitiRegionsFreq)," ")),names(highDensitiRegionsFreq))
colnames(tab) <- c("chrom","start","end","regid")
tab$start <-as.numeric(tab$start)
tab$end <-as.numeric(tab$end)

for(i in unique(tab$chrom)) chrend[i] <- max(tab[which(tab[,1] == i),3])
for(i in 0:(length(chrend)-1) ) axislab[i+1] <- chrend[i+1]/2 + sum( chrend[0:i])
for(i in 0:(length(chrend)-1) ) chrstarts[i+1] <- sum(chrend[0:i])
data <- cbind( (tab$end + tab$start) / 2 + chrstarts[tab$chrom], highDensitiRegionsFreq)


altcols <- rep(c(rgb(0.1,0.1,0.1,alpha=0.1),rgb(0.8,0.8,0.8,alpha=0.1)),12)
altcols2<- rep(c(rgb(0.1,0.1,0.1,alpha=1),rgb(0.4,0.4,0.4,alpha=1)),12)
ctrmr <- chrstarts+centromeres[names(chrstarts)]

plot(data[,1:2],type='h',col=coloresBarplot,xaxt='n',lwd=1.5,ylim=c(0, max(data[,2])+5),
     las=1,bty='n',yaxt='n',family="Arial",ylab="",xlab="")
for(i in 1:length(chrstarts) ) rect( chrstarts[i],0,chrstarts[i]+chrlengths[i],1000, col=altcols[i],border=NA )
mtext(gsub("chr","",names(axislab)),side=1,at=axislab,las=1,col=altcols2,cex=c(rep(1,17),rep(0.8,5),1) )
if(!is.null(freq.cut)) lines(c(0,chrstarts["chrX"]+chrlengths["chrX"]),c(freq.cut,freq.cut),lty=3,col="black")    
axis(2,las=1,pos= 0, cex=1.2)
axis(4,las=1,pos= max(data[,1])+10000, cex=1.2, at=axTicks(2), labels=sprintf("%.2f",axTicks(2)/dim(chromo.regs.obj@high.density.regions)[1]) )
mtext("Frequency",side=4,line=1.5)
mtext("#samples",side=2,line=1)
if(!is.null(add.legend)) legend(add.legend,c("short (p) arm","long (q) arm"),border=c("salmon","salmon4"),fill=c("salmon","salmon4"),bty='n',ncol=2)
}

