#'@include ROBIBarcodes.R
#'@include logo.R
NULL

#' Draw a scatter plot of the reverse mismatches as a function of forward mismatches.
#' 
#' The \code{mismatchplot} function draws a scatter plot of the number of mismatches
#' observed in an ecoPCR result for the reverse primer as a function of the mismatches
#' for the reverse primer. Each point for a pair (forward_mismatch,reverse_mismatch) is
#' drawn as a circle having a surface proportional to the aboundance of this pair in the
#' results. If a grouping factor is specified, then the circle is replaced by a pie chart.
#' 
#' @param ecopcr an ecoPCR result data.frame as returned by the \code{\link{read.ecopcr.result}}
#'               function. 
#'               
#' @param group  a factor decribing classes amongst the amplicon described in the ecoPCR
#'               result
#'               
#' @param col    a vector describing the colored used for the bubble or the pie charts
#' 
#' @param legend a character vector describing the legend for each modality of the
#'               grouping factor. By default the factor levels are used for the legend
#'               
#' @param legend.cex the expension factor for the legend text
#' 
#' @param inset  the distance to the margin of the legend box (see the \code{\link{legend}} 
#'               documentation)
#'               
#' @param view.legend if set to \code{FALSE} the legend corresponding to the groups is not dispayed.
#' 
#' @param maxerror allows for specifying the maximum of errors to display on the graph.
#' 
#' @examples
#' 
#' # Load the ROBITools library
#' library(ROBITools)
#' 
#' # Load the default taxonomy
#' taxo = default.taxonomy()
#' 
#' # Load the sample ecoPCR data file
#' data(GH.ecopcr)
#' 
#' # Computes classes associated to each taxid 
#' orders = as.factor(taxonatrank(taxo,GH.ecopcr$taxid,'order',name=T))
#' 
#' # Plot the graph
#' mismatchplot(GH.ecopcr,group=orders)
#' 
#' @seealso \code{\link{read.ecopcr.result}}
#' @author Eric Coissac
#' @export
mismatchplot = function(ecopcr,group=NULL,
                        col=NULL,legend=NULL,
                        legend.cex=0.7,inset=c(0.02,0.02),
                        view.legend=TRUE,
                        maxerror=NA) {
  
  maxforward_error = max(ecopcr$forward_mismatch)
  maxreverse_error = max(ecopcr$reverse_mismatch)
  if (is.na(maxerror))
    maxerror=max(maxforward_error,maxreverse_error)

  if (is.null(group))
    group=factor(rep("all",dim(ecopcr)[1]))
  else
    group=as.factor(group)
  
  if (is.null(legend))
    legend = levels(group)
  
  actualheight= maxerror + 1
  actualwidth = maxerror + 1
  
  if (length(levels(group)) > 1 & view.legend)
    actualwidth = actualwidth + 2
  
  whitepaper(actualwidth,actualheight,xmin=-0.5,ymin=-0.5,asp=1)
  
  axis(1,at=0:maxerror,
       labels=0:maxerror)
  
  axis(2,at=0:maxerror,
       labels=0:maxerror)
  
  
  data = aggregate(group,by=list(forward=factor(ecopcr$forward_mismatch,levels=0:maxerror),
                                 reverse=factor(ecopcr$reverse_mismatch,levels=0:maxerror)),
                   table)
  
  data <- data[rowSums(data[,c(-1,-2),drop=FALSE])>0, , drop=FALSE]
  
  if (is.null(col)) 
    col <- c("white", "lightblue", "mistyrose", "lightcyan", 
             "lavender", "cornsilk")
  
  
  value=data[,c(-1,-2),drop=FALSE]
  x = as.integer(data[,1]) - 1
  y = as.integer(data[,2]) - 1
  diam = sqrt(rowSums(value))
  radius = diam / max(diam) / 2
  
  hide = mapply(pie.xy,x,y,
                data=lapply(1:(dim(value)[1]),function(y) value[y,]),
                radius=radius,
                label="",MoreArgs=list(col=col))
  
  
  if ((length(levels(group))) > 1 & view.legend)
    legend('topright',legend=legend,fill=col, cex=legend.cex, inset=inset)
  
}