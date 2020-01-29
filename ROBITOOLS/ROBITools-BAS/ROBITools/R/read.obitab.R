#' @include 02_class_metabarcoding.data.R
NULL


#' Reads a data file produced by the obitab command
#'
#' Read a data file issued from the convertion of a fasta 
#' file to a tabular file by the obitab command
#' 
#' @param file a string containing the file name of the obitab file.
#' @param sep  Column separator in the obitab file. 
#'             The default separator is the tabulation.
#'
#' @return a \code{data.frame} instance containing the obitab file
#'
#' @examples
#' require(ROBITools)
#' 
#' \dontshow{# switch the working directory to the data package directory}
#' \dontshow{setwd(system.file("extdata", package="ROBITools"))}
#' 
#' # read the termes.tab file
#' termes=read.obitab('termes.tab')
#' 
#' # print the dimensions of the data.frame
#' dim(termes)
#'   
#' @seealso \code{\link{import.metabarcoding.data}}
#' @author Eric Coissac
#' @export
#'
read.obitab <-
function(filename,sep='\t') {

   data=read.delim(filename,sep=sep,strip.white=T,check.names =F)
   data
   
}

