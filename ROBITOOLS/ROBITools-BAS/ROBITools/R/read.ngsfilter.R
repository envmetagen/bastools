#' @include 02_class_metabarcoding.data.R
NULL

#' Read an OBITools ngsfilter file
#' 
#' Reads a ngsfilter file as formatted for the OBITools. For now, needs to be tab delimited till the "F" column. 
#' Any additionnal information needs to be space delimited.
#' 
#' @seealso \code{\link{import.metabarcoding.data}}
#' @author Lucie Zinger
#' @keywords data import
#' @export
#' 

read.ngsfilter <- function(filename, decimal='.', as.is=!stringsAsFactors, stringsAsFactors = default.stringsAsFactors()) {
	
	t<-read.table(file=filename, header=F, sep="\t", as.is=T)
	beg <- t[,1:5]
	
	colnames(beg) <- c('experiment','sample','tags','forward_primer','reverse_primer')
	if (length(unique(beg$sample))==nrow(beg)) 
		rownames(beg) <- beg$sample
	end <- t[,c(2,6)]
	
	#F <- unlist(lapply(end$V6, function(x) strsplit(x,"@")[[1]][1]))
	rawextras <- unlist(lapply(end$V6, function(x) strsplit(x,"@")[[1]][2]))
	
	rawextras <- lapply(rawextras, function(s) strsplit(s, '; ')[[1]])
	rawextras <- lapply(rawextras, function(l) unlist(lapply(l, function(s) sub("^ +","",s))))
	rawextras <- lapply(rawextras, function(l) unlist(lapply(l, function(s) sub(" +$","",s))))
	
	
	rawextras <- lapply(rawextras, function(l) unlist(lapply(l, function(s) strsplit(s,"="))))
	
	
	columnnames <- unique(unlist(lapply(rawextras, function(l) l[seq(1,length(l),2)])))
	
	m <- matrix(nrow=nrow(end), ncol=length(columnnames))
	colnames(m) <- columnnames
	m <- as.data.frame(m)    
	
	
	#print(head(rawextras))
	
	
	tt <- lapply(rawextras, function(l) list(l[seq(1,length(l),2)],l[seq(2,length(l),2)]))
	invisible(lapply(1:length(tt), function(i){m[i,tt[[i]][[1]]] <<- tt[[i]][[2]]}))
	
	invisible(lapply(colnames(m), function(n) m[,n] <<- type.convert(m[,n], dec=decimal, as.is=as.is)))
	
	ngs = cbind(beg, m)
	rownames(ngs) = ngs$sample
	class(ngs)<-c('ngsfilter.data',class(ngs))
	
	return(ngs)
}
