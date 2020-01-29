#' @include 02_class_metabarcoding.data.R
NULL


# TODO: Add comment
# 
# Author: coissac
###############################################################################


read.sumatra = function(filename) {
	data = read.table(filename,sep="\t",header=FALSE)
	score = data[,3]
	name.first = mapply(min,as.character(s[,1]),as.character(s[,2]))
	name.second= mapply(max,as.character(s[,1]),as.character(s[,2]))
	sname = as.character(interaction(data[,1],data[,2]))
}
