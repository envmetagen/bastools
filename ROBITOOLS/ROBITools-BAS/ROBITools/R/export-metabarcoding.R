#' @include 02_class_metabarcoding.data.R
NULL


# TODO: Add comment
# 
# Author: coissac
###############################################################################

require(utils)

expand.metabarcoding.data=function(data,minread=1) {
	resultonesample=function(sample) {
		mo= data@reads[sample,] >= minread
		s = data@samples[rep(sample,sum(mo)),]
		r = as.numeric(data@reads[sample,mo])
		m = data@motus[mo,]
		
		result = data.frame(s,frequency=r,m, 
				stringsAsFactors =FALSE,
				row.names = NULL)
		
		result
	}
	
	res = lapply(1:data@scount, resultonesample)
	
	do.call(rbind,res)	
}

#setGeneric("utils::write.csv")
write.csv.metabarcoding.data = function(...) {
	Call <- match.call(expand.dots = TRUE)
	if (!is.null(Call[["minread"]])) {
		minread = Call[["minread"]]
		Call = Call[!names(Call)=="minread"]
	}
	else
		minread = 1
	data = eval.parent(Call[[2L]])
	data = expand.metabarcoding.data(data,minread)
	Call[[1L]] <- as.name("write.csv")
	Call[[2L]] <- as.name("data")
	eval(Call)
}

#setGeneric("utils::write.csv2")
write.csv2.metabarcoding.data = function(...) {
	Call <- match.call(expand.dots = TRUE)
	if (!is.null(Call[["minread"]])) {
		minread = Call[["minread"]]
		Call = Call[!names(Call)=="minread"]
	}
	else
		minread = 1
	data = eval.parent(Call[[2L]])
	data = expand.metabarcoding.data(data,minread)
	Call[[1L]] <- as.name("write.csv2")
	Call[[2L]] <- as.name("data")
	eval(Call)
	
}
