#' Percent genome change calculation
#' 
#' Calculates the percentage of genome changed using CNV segmentation profiles. Genome change is defined based on the fold change CNV log-ratio between a sampele and a reference. 
#' 
#' @param cnv (S4) an object of class svcnvio containing data type 'cnv' initialized by validate.cnv
#' @param fc.pct (numeric) percentage CNV gain/loss for a segment to be considered changed (e.g. 0.2 = 20 percent change 0.8 < segmean && segmean > 1.2)
#' @param discard.sex (logical) whether sex chromosomes should be included
#' @return (numeric) vector containing percent genome changed values (0-1)
#' @seealso Additional data format information in the man pages of validate.cnv
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' ## validate input CNV data.frames
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' pct_changed <- pct.genome.changed(cnv)
#' head(pct_changed)

pct.genome.changed <- function(cnv, 
                               fc.pct=0.2, 
                               discard.sex=TRUE){

cnvdat <- cnv@data
if(discard.sex == TRUE) cnvdat <- cnvdat[which(!cnvdat$chrom %in% c("chrX","chrY")),]

width <- cnvdat$end - cnvdat$start
segmean <- cnvdat$segmean
sample <- cnvdat$sample
df <- data.table(sample,width,segmean)
idx_changed <- c(which(df$segmean < log2(1-fc.pct)),which(df$segmean >= log2(1+fc.pct)))
idx_normal <- setdiff(1:nrow(df),idx_changed)
df_normal <- df[idx_normal,]
df_changed <-  df[idx_changed,]
  
length_changed_df <- aggregate(width~sample ,df_changed,sum)
length_normal_df <- aggregate(width~sample ,df_normal,sum)
  
nochange <- setdiff(length_normal_df$sample,length_changed_df$sample)
fullchange <- setdiff(length_changed_df$sample,length_normal_df$sample)
nochange_x <- rep(0,length(nochange))
names(nochange_x) <- nochange
fullchange_x <- rep(0,length(fullchange))
names(fullchange_x) <- fullchange
  
length_changed <- c(length_changed_df[,2],nochange_x)
names(length_changed)<- c(length_changed_df[,1],nochange)
  
length_normal <- c(length_normal_df[,2],fullchange_x)
names(length_normal)<- c(length_normal_df[,1],fullchange)
  
pct.change<- length_changed/apply(cbind(length_normal[names(length_changed)],length_changed),1,sum)
  
return(pct.change)
}


