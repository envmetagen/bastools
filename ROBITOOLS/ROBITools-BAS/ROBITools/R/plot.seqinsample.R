#' @include 02_class_metabarcoding.data.R
NULL

#' Plot sequence abundance in samples
#' 
#' Plots relative abundances of a set of sequences in all samples (log10 transformed)
#' 
#' 
#' @param x          a \code{\link{metabarcoding.data}} object
#' @param seqset     a vetcor with sequences names
#' @param seqtype    a string indicating what type of sequences are displayed
#' @param controls   a vector indicating the negative controls names in the x object.
#'                   Default is \code{NULL}
#' 
#' @return returns a plot with the log10 transformed relative porportion of 
#'         selected MOTUs in each samples. If the number of samples is > 96,
#'         then the plot is displayed in 4 panels
#' 
#' @examples
#' 
#' data(termes)
#' 
#' seqset = rownames(termes$motus)[which(termes$motus$genus_name=="Anoplotermes")]
#' plot.seqinsample(termes, seqset, "Anoplotermes")
#' 
#' controls = rownames(termes)[grep("r", rownames(termes))]
#' seqset = rownames(termes$motus)[which(termes$motus$best_identity<0.7)]
#' plot.seqinsample(termes, seqset, "Not assigned", controls)
#' 
#' @seealso \code{\linkS4class{taxonomy.obitools}}, and method \code{\link{taxonmicank}}
#'
#' @author Lucie Zinger
#' @keywords metabarcoding
#' 
#' @export
#' 

plot.seqinsample = function(x, seqset, seqtype, controls=NULL){
 
  require(vegan)
    
  x.freq = vegan::decostand(x$reads,"total",1)
  
  if(!is.null(controls)){
    controls.ind = match(controls, rownames(x.freq))
  }
  
  if(nrow(x.freq)>96){
    x.freq.parse = seq(0,round(nrow(x$samples), digit=0),
                       round(nrow(x$samples)/4, digit=0))
    
    layout(matrix(c(1,2,3,1,4,5),3,2), height=c(0.3,1,1))
    par(oma=c(1,1,1,0), mar=c(3,3,1,1))
    
    #legend
    breaks = seq(log10(1e-4),log10(1), length.out=100)
    plot(breaks, rep(1,100), col=topo.colors(100), pch=15, cex=2, ylim=c(0,1.5),
         xaxt="n", yaxt="n", bty='n')
    text(breaks[seq(1,100,10)], rep(0.7,length(seq(1,100,10))), 
         round(10^breaks[seq(1,100,10)],4))
    mtext("Seqence frequencies:", side=3, line=0, cex=0.8)
    
    #plot
    for(i in 1:(length(x.freq.parse)-1)) {
      range = (x.freq.parse[i]+1):(x.freq.parse[i]+round(nrow(x$samples)/4, digit=0))
      mat = x.freq[range,seqset]
      image(log10(mat),col = topo.colors(100), xaxt="n", yaxt="n", breaks=c(breaks,0))
      
      if(!is.null(controls)){
        if(length(na.omit(match(controls.ind, range)))!=0){
          abline(v=seq(0,1,l=round(nrow(x$samples)/4, digit=0))[match(controls.ind, range)],col="red", lty=3)
        }}
      
      axis(side=1,at=seq(0,1,l=round(nrow(x$samples)/4,digit=0)),
           labels=rownames(x$samples)[range],
           las=2, cex.axis=0.3)
    }
    mtext(side=2, paste(seqtype, "n = ", length(seqset)), outer=T, cex=0.7, font=3)
    mtext(side=1, "Samples", cex=0.7, outer=T)
  
  } else {
    layout(matrix(c(1,2,1,2),2,2), height=c(0.3,1))
    par(oma=c(1,1,1,0), mar=c(3,3,1,1))
    
    #legend
    breaks = seq(log10(1e-4),log10(1), length.out=100)
    plot(breaks, rep(1,100), col=topo.colors(100), pch=15, cex=2, ylim=c(0,1.5),
         xaxt="n", yaxt="n", bty='n')
    text(breaks[seq(1,100,10)], rep(0.7,length(seq(1,100,10))), 
         round(10^breaks[seq(1,100,10)],4))
    mtext("Seqence frequencies:", side=3, line=0, cex=0.8)
    
    image(log10(x.freq[,seqset]),col = topo.colors(100), xaxt="n", yaxt="n", breaks=c(breaks,0))
    
    if(!is.null(controls)){
      abline(v=seq(0,1,l=round(nrow(x$samples), digit=0))[controls.ind],col="red", lty=3)
      }
    axis(side=1,at=seq(0,1,l=round(nrow(x$samples),digit=0)),
         labels=rownames(x$samples),
         las=2, cex.axis=0.3)
  mtext(side=2, paste(seqtype, "n = ", length(seqset)), outer=T, cex=0.7, font=3)
  mtext(side=1, "Samples", cex=0.7, outer=T)
  }
}

