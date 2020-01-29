#' @include 02_class_metabarcoding.data.R
NULL

#' Plot PCR plates
#' 
#' Plots samples localization in PCR plates, and points out problematic samples if provided.
#' 
#' @param x a \code{\link{metabarcoding.data}} object
#' @param samples a character vector containing names of problematic samples. Default is \code{NULL}
#' @param different a boolean indicating whether different tags where used in forward and reverse to identify samples. Default is \code{TRUE}
#' @param ... arguments ot be passed to methods, such as graphical parameters
#' 
#' @return \code{\link{plot.PCRplate}} returns a plot displaying no more than 4 PCR plates, with problematic sample localization
#' 
#' @examples
#' \dontshow{# switch the working directory to the data package directory}
#' \dontshow{setwd(system.file("extdata", package="ROBITools"))}
#' 
#' data(termes)
#' 
#' # reading the termes_ngsfilt.txt file
#' termes.ngs=import.ngsfilter.data('termes_ngsfilt.txt', platewell="position")
#' 
#' # including ngsfilter data into termes data
#' attr(termes, "samples") = termes.ngs[rownames(termes),]
#' 
#' #plot PCR plate plan
#' col = rep("green", nrow(termes))
#' col[grep("r", rownames(termes))] = "red"
#' plot.PCRplate(termes, col=col)
#' 
#' #highlighting location of samples with low identification score
#' 
#' #low quality taxonomic assignements identification
#' library(plotrix)
#' weighted.hist(termes$motus$best_identity, colSums(termes$reads), breaks = 20, ylab = "Nb reads", xlab = "Ecotag scores", xaxis=F)
#' axis(1, labels = T)
#' lowqual.seq = rownames(termes$motus)[termes$motus$best_identity < 0.7]
#' 
#' #identification and localization (in PCR plate) of samples with high proportions of low quality taxonomic assignements
#' termes.freq= normalize(termes, MARGIN=1)$reads
#' hist(log10(rowSums(termes.freq[,lowqual.seq]) + 1e-05), breaks = 20, xlab = "Prop low quality reads")
#' lowqual.sample = rownames(termes)[log10(rowSums(termes.freq[, lowqual.seq]) + 1e-05) > -0.5]
#' 
#' plot.PCRplate(termes, lowqual.sample, col=col)
#'   
#' @seealso \code{\link{import.metabarcoding.data}}
#'
#' @author Lucie Zinger
#' @keywords DNA metabarcoding
#' @export
#' 
plot.PCRplate = function(x, samples=NULL, col="cyan2", different=T, ...) {
  
  if(length(grep("xPlate", colnames(x$samples)))==0 | 
       length(grep("yPlate", colnames(x$samples)))==0) {
    stop("samples/controls position in PCR plates (xPlate and yPlate) are not defined")
  }
  
  if(length(grep("tagF", colnames(x$samples)))==0 | 
       length(grep("tagR", colnames(x$samples)))==0) {
    stop("tags (tagF and tagR) are not defined")
  }
  
  nplate = max(x$samples$nbPlate)
  
  if(nplate>4) {
    stop("Cannot plot more than 4 plates")
  }
   
    plot(x$samples$xPlate, -x$samples$yPlate, pch=19, xaxt="n", yaxt="n", col=col,
         xlim=c(-5,17), ylab="y plate", xlab= "x plate", ylim=c(-4.5*8-5,0), ...)
    if(different==T) {
      text(-3, -unique(x$samples$yPlate[order(x$samples$yPlate)]), unique(x$samples$tagF[order(x$samples$yPlate)]), cex=0.5)
      text(unique(x$samples$xPlate[order(x$samples$xPlate)]), -5, unique(x$samples$tagR[order(x$samples$xPlate)]), cex=0.5, srt=90)
      }
      abline(h=-seq(8.5,8*nplate+0.5,8), lty=2, col="grey")
      segments(c(0,13), rep(min(-x$samples$yPlate),2), c(0,13), c(0,0), lty=2, col="grey")
  
    #plot problematic samples
    if(!is.null(samples)) {
      points(x$samples[samples,"xPlate"], -x$samples[samples,"yPlate"], pch="x")
    }  
}
