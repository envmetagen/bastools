#' @include 02_class_metabarcoding.data.R
NULL

#' Read ngsfilter text file
#' 
#' Reads  the text file used for assigning reads to samples with the
#'  \code{ngsfilter} command of the \strong{OBITools} package. 
#' 
#' @param file a string containing the file name for the \code{ngsfilter} command.
#' @param platewell a string corresponding to the tag used for storing the sample location
#'                  in the PCR plate. Should be of the form "nbPlate_Well" (e.g. "01_A02").
#'                  Default is \code{NULL}
#' @return \code{\link{import.ngsfilter.data}} returns a \code{\link{data.frame}} instance
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
#' colnames(termes$samples)
#'   
#' @seealso \code{\link{import.metabarcoding.data}} and \code{\link{read.obitab}} for other methods of data importation 
#'
#' @author Lucie Zinger
#' @keywords DNA metabarcoding
#' @export
#' 
import.ngsfilter.data = function(file, platewell=NULL) {
  raw = read.table(file, sep="\t")
  
  #get samples names
  names = raw[,2]
  
  #form first part of the output table (default ngsfilter text input)
  out = raw[,-c(2,3,ncol(raw))]
  colnames(out) = c("Experiment", "primerF", "primerR")

  #add tags
  out[,c("tagF", "tagR")] = do.call("rbind", strsplit(as.vector(raw[,3]), "\\:"))
  
  #collect nb and names of additionnal information
  max.add = max(unlist(lapply(strsplit(gsub("^F @ ","", raw[, ncol(raw)]), "; "), length)))
  names.add = unique(unlist(lapply(strsplit(unlist(strsplit(gsub("^F @ ","", raw[, ncol(raw)]), "; ")), "="), "[[",1)))
  
  #form table of additionnal info
  form = lapply(strsplit(gsub("^F @ ","", raw[, ncol(raw)]), "; "), strsplit, "=")
  additionnals = as.data.frame(do.call("rbind", lapply(form, function(y) {
    val = rep(NA, , max.add)
    names(val) = names.add
    val[match(unlist(lapply(y, "[[", 1)), names(val))] = gsub(";", "",unlist(lapply(y, "[[", 2)))
    val
  })))
  
  #create PCR plate coordinates
  if(!is.null(platewell)) {
    form = strsplit(as.vector(additionnals[, platewell]), "_")
    nbPlate = as.numeric(gsub("^0", "", unlist(lapply(form, "[[", 1))))
    wellPlate = unlist(lapply(form, "[[", 2))
    xPlate = as.numeric(gsub("[A-Z]", "", wellPlate))
    yPlate = as.numeric(as.factor(gsub("[0-9]*", "", wellPlate))) + 8*nbPlate
    
    additionnals = additionnals[,-grep(platewell, colnames(additionnals))]
    out = data.frame(out, additionnals, nbPlate, wellPlate, xPlate, yPlate)
  }
  else {
    additionnals[,ncol(additionnals)] = gsub(";","", additionnals[,ncol(additionnals)])
    out = data.frame(out, additionnals)
  }

  rownames(out) = names
  return(out)
}
