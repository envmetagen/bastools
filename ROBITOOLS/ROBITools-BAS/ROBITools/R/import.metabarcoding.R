#' @include read.obitab.R
#' @include 02_class_metabarcoding.data.R
NULL

#' Read a data file produced by the \code{obitab} command 
#' 
#' Read a data file issued from the conversion of a \strong{fasta} 
#' file to a tabular file by the \code{obitab} command of the 
#' \strong{OBITools} package 
#' 
#' @param file a string containing the file name of the obitab file.
#' @param sep  Column separator in the obitab file. 
#'             The default separator is the tabulation.
#' @param sample A regular expression allowing to identify columns 
#'               from the file describing abundances of sequences per sample
#' @param sample.sep Separator between combined sample name.
#' @param attribute Separator used to split between sample 'tag' and sample name.
#' 
#' @return a \code{\link{metabarcoding.data}} instance
#' 
#' @examples
#' require(ROBITools)
#' 
#' \dontshow{# switch the working directory to the data package directory}
#' \dontshow{setwd(system.file("extdata", package="ROBITools"))}
#' 
#' # read the termes.tab file
#' termes=import.metabarcoding.data('termes.tab')
#' 
#' # print the number of samples and motus described in the file
#' dim(termes)
#'   
#' @seealso \code{\link{metabarcoding.data}}
#'
#' @author Eric Coissac
#' @keywords DNA metabarcoding
#' @export
#' 
import.metabarcoding.data = function(file,sep='\t',sample="sample",sample.sep="\\.",attribute=":") {
	
	data=read.obitab(file,sep=sep)
	
	# get the colnames matching the sample pattern
	
	column=colnames(data)
	pat = paste('(^|',sample.sep,')',sample,'[',sample.sep,attribute,']',sep='')
	scol= grep(pat,column)
	
	# reads informations about samples
	
	reads  = data[,scol]
	names  = colnames(reads)
	names  = strsplit(names,split=attribute)
	
			# for sample name just remove the first part of the col names
			# usally "sample:"
	
	sample.names = sapply(names,function(a) paste(a[-1],collapse=attribute))	
	
	reads=t(reads)
	rownames(reads)=sample.names
	
	# sample's data
	
	sample.data = data.frame(t(data.frame(strsplit(sample.names,split=attribute))))
	rownames(sample.data)=sample.names
	colnames(sample.data)=strsplit(names[[1]][1],split=attribute)
	
	
	# motus information 

	motus = data[,-scol]
	
	motus.id = motus$id
	
	rownames(motus)=motus.id
	colnames(reads)=motus.id
	
	
	return(metabarcoding.data(reads,sample.data,motus))
	
}


#pcr = gh[,grep('^sample',colnames(gh))]
#pcr.names = colnames(pcr)
#pcr.names = sub('sample\\.','',pcr.names)
#sequencer = rep('Solexa',length(pcr.names))
#sequencer[grep('454',pcr.names)]='454'
#sequencer=factor(sequencer)
#
#tmp = strsplit(pcr.names,'\\.[A-Z](sol|454)\\.')
#
#sample = sapply(tmp,function(x) x[1])
#locality = factor(sapply(strsplit(sample,'_'),function(x) x[1]))
#sample = factor(sample)
#repeats= factor(sapply(tmp,function(x) x[2]))
#
#tmp = regexpr('[A-Z](454|sol)',pcr.names)
#run=factor(substr(pcr.names,tmp,tmp+attr(tmp,"match.length")-1))
#
#pcr.metadata = data.frame(run,sequencer,locality,sample,repeats)
#
#rownames(pcr.metadata)=pcr.names


