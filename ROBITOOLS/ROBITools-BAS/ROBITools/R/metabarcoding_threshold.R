#' @include 02_class_metabarcoding.data.R
NULL


#' @export
setGeneric("marginalsum", function(data,MARGIN="sample", na.rm = FALSE) {
			return(standardGeneric("marginalsum"))
		})


#' Computes marginal sums over read counts.
#' 
#' Method \code{marginalsum} computes marginal sums over read counts of
#' a \code{\link{metabarcoding.data}} instance.
#' 
#' @param data The \code{\linkS4class{metabarcoding.data}} instance
#'             on which marginal sums have to be computed.
#' @param MARGIN Indicates if the sums have to be computed across 
#'               samples or motus. 
#'               Allowed values are :
#'               \itemize{
#'                 \item{'sample' or 1} for computing sum across samples
#'                 \item{'motu' or 2} for computing sum across motus
#'                 }
#' @param na.rm  Logical. Should missing values be omitted from the 
#'               calculations?
#'               
#' @return Returns the vector of marginal sums as a \code{numeric} vector
#' 
#' @examples
#' # load termite data set from the ROBITools sample data
#' data(termes)
#' 
#' # Computes marginal sums per sample
#' ssum = marginalsum(termes,MARGIN="sample")
#' 
#' # Computes marginal sums per MOTU
#' msum = marginalsum(termes,MARGIN="motu")
#' 
#' @seealso \code{\linkS4class{metabarcoding.data}}
#'
#' @docType methods
#' @rdname marginalsum-methods
#' @aliases marginalsum-methods,metabarcoding.data
#' @author Aurelie Bonin
#' 
setMethod("marginalsum", "metabarcoding.data", function(data,MARGIN='sample', na.rm = FALSE) {
	
	if (MARGIN == 'sample')
		MARGIN=1
	
	if (MARGIN == 'motu')
		MARGIN=2
	
	readcount = reads(data)
  if (MARGIN==1)    
	  margesum = rowSums(readcount,na.rm=na.rm)
  else
    margesum = colSums(readcount,na.rm=na.rm)
	
	
	return(margesum)
})

rowSums.metabarcoding.data = function (x, na.rm = FALSE, dims = 1L) {
  print("coucou")
}

#' @export
setGeneric("normalize", function(data,MARGIN='sample',as.matrix=FALSE) {
			return(standardGeneric("normalize"))
		})


#' Normalizes read counts by sample or by MOTU.
#' 
#' Method \code{normalize} computes a normalized read aboundancy matrix
#' (relative frequency matrix) of a \code{\link{metabarcoding.data}} instance.
#' Normalization can be done according aboundancies per sample or per MOTU.
#' 
#' @param data The \code{\linkS4class{metabarcoding.data}} instance
#'             on normalisation have to be computed.
#' @param MARGIN Indicates if the sums have to be computed across 
#'               samples or motus. 
#'               Allowed values are :
#'               \itemize{
#'                 \item{'sample' or 1} for computing sum across samples
#'                 \item{'motu' or 2} for computing sum across motus
#'                 }
#' @param as.matrix Logical indicating if the normalized aboundancies
#'               must be returned as a simple \code{matrix} (TRUE) or as a new
#'               instance of the \code{\linkS4class{metabarcoding.data}} class
#'               (FALSE, the default case).
#'               
#' @return Returns a new instance of \code{\linkS4class{metabarcoding.data}}
#'         or a \code{numeric} matrix according to the \code{return.as.matrix}
#'         parameter.
#' 
#' @examples
#' # load termite data set from the ROBITools sample data
#' data(termes)
#' 
#' # Computes normalized aboundancies per sample
#' termes.norm = normalize(termes,MARGIN="sample")
#' 
#' # Computes normalized aboundancies per sample and
#' # stores the result as a new layer into the thermes
#' # structure
#' termes$normalized = normalize(termes,MARGIN="sample",as.matrix=TRUE)
#' 
#' @seealso \code{\linkS4class{metabarcoding.data}}
#'
#' @docType methods
#' @rdname normalize-methods
#' @aliases normalize-methods,metabarcoding.data
#' @author Aurelie Bonin
#' 
setMethod("normalize", "metabarcoding.data", function(data,MARGIN="sample",as.matrix=FALSE) {
	
	if (MARGIN == 'sample')
		MARGIN=1
	
	if (MARGIN == 'motu')
		MARGIN=2
	
	readcount = reads(data)
	margesum  = marginalsum(data,MARGIN,na.rm=TRUE)
	
	readcount = sweep(readcount,MARGIN,margesum, FUN="/")
	
  if (as.matrix)
    newdata=readcount
  else
	  newdata = copy.metabarcoding.data(data,reads=readcount)
	
	return(newdata)
})

#' @export
setGeneric("threshold", function(data,MARGIN="sample",threshold=0.97) {
			return(standardGeneric("threshold"))
		})

#' Compute the cumulative threshold of read aboundances.
#' 
#' The method \code{threshold} of the class \code{\linkS4class{metabarcoding.data}}
#' computes the thresold to be used for conserving just a part of the global
#' signal. This thresold is computed by ranking aboundances by decreasing order.
#' The cululative sums of these ranked abondencies are computed and the aboundance
#' corresponding to the first sum greater than the threshold is returned as result.
#' 
#' @param data The \code{\linkS4class{metabarcoding.data}} instance
#'             on normalisation have to be computed.
#' @param MARGIN Indicates if the sums have to be computed across 
#'               samples or motus. 
#'               Allowed values are :
#'               \itemize{
#'                 \item{'sample' or 1} for computing sum across samples
#'                 \item{'motu' or 2} for computing sum across motus
#'                 }
#' @param threshold a numeric value between 0 and 1 indicating which part of 
#'                  the signal must be conserved. Default value is setup to
#'                  0.97 (97% of the total signal).
#'                  
#' @return a numeric vector containing the limit aboundancy to consider for
#'         each sample or each MOTU according to the value of the \code{MARGIN} 
#'         parameter.
#'                  
#' @examples
#' # load termite data set from the ROBITools sample data
#' data(termes)
#' 
#' # computes threshold value to used for keep 95% of 
#' # the reads per MOTU
#' 
#' t = threshold(termes,MARGIN='motu',threshold=0.95)
#'  
#' @seealso \code{\linkS4class{metabarcoding.data}}, \code{\link{threshold.mask}}
#'
#' @docType methods
#' @rdname threshold-methods
#' @aliases threshold-methods,metabarcoding.data
#' @author Aurelie Bonin
#'
setMethod("threshold", "metabarcoding.data", function(data,MARGIN="sample",threshold=0.97) {
	
	
	onethreshold=function(x,threshold) {
		s = x[order(-x)]
		cs= cumsum(s) 
    total=cs[length(cs)]
    if (total > 0) {
		  cs= cs / total
		  cs = cs > threshold
		  t = s[cs][1]
    }
    else t=0
		
    return(t)
	}
	
	
	if (MARGIN == 'sample')
		MARGIN=1
	
	if (MARGIN == 'motu')
		MARGIN=2
	
	readcount = reads(data)
	
	t = apply(readcount,MARGIN,onethreshold,threshold)
	
	return(t)
})

#' @export
setGeneric("threshold.mask", function(data,MARGIN,threshold=0.97,operator='<') {
			return(standardGeneric("threshold.mask"))
		})

#' Computes a cumulatif thresold mask for filtering read aboundancies.
#' 
#' The method \code{threshold.mask} of the class \code{\linkS4class{metabarcoding.data}}
#' computes a logical matrix of the same size than the read matrix of the data parameter.
#' Each cell of this matrix contains a \code{TRUE} or a \code{FALSE} value according to the
#' relationship existing between the read abondancy and the corresponding theshold as computed
#' by the \code{\link{theshold}} method.
#' 
#' (computed value) = (read aboundancy) operator (threshold value)
#' 
#' for a cell in the result matrix, \code{(read aboundancy)} is extracted from the read layer.
#' \code{operator} is a comparaison operator and \code{(threshold value)} is estimated with the
#' \code{\link{theshold}} method.
#' 
#' @param data The \code{\linkS4class{metabarcoding.data}} instance
#'             on normalisation have to be computed.
#' @param MARGIN Indicates if the sums have to be computed across 
#'               samples or motus. 
#'               Allowed values are :
#'               \itemize{
#'                 \item{'sample' or 1} for computing sum across samples
#'                 \item{'motu' or 2} for computing sum across motus
#'                 }
#' @param threshold a numeric value between 0 and 1 indicating which part of 
#'                  the signal must be conserved. Default value is setup to
#'                  0.97 (97% of the total signal).
#' @param operator is a logical comparison operator.
#' 
#' @return A logical matrix usable for selecting cell in the read aboundancy matrix.
#'                   
#' @seealso \code{\linkS4class{metabarcoding.data}}, \code{\link{threshold.mask}}, \code{\link{threshold}}
#'
#' @docType methods
#' @rdname threshold-mask-methods
#' @aliases threshold.mask-methods,metabarcoding.data
#' @author Aurelie Bonin
#' 
setMethod("threshold.mask", "metabarcoding.data", function(data,MARGIN,threshold=0.97,operator='<') {
	
	
	if (MARGIN == 'sample')
		MARGIN=1
	
	if (MARGIN == 'motu')
		MARGIN=2
	
	readcount = reads(data)
	
	t = threshold(data,MARGIN,threshold)
	mask = apply(readcount,c(2,1)[MARGIN],operator,t)
	
	if (MARGIN==2)
		mask = t(mask)
	
	return(mask)
})


#' @export
setGeneric("const.threshold.mask", function(data,MARGIN,threshold=0.01,operator='<') {
			return(standardGeneric("const.threshold.mask"))
		})

#' Computes a constant thresold mask for filtering read aboundancies.
#' 
#' The method \code{const.threshold.mask} of the class \code{\linkS4class{metabarcoding.data}}
#' computes a logical matrix of the same size than the read matrix of the data parameter.
#' Each cell of this matrix contains a \code{TRUE} or a \code{FALSE} value according to the
#' relationship existing between the read abondancy and the global theshold.
#' 
#' (computed value) = (normalized read aboundancy) operator (threshold value)
#' 
#' for a cell in the result matrix, \code{(normalized read aboundancy)} is extracted from the read layer
#' after normalization.
#' \code{operator} is a comparaison operator and \code{(threshold value)} is estimated with the
#' \code{\link{theshold}} method.
#' 
#' @param data The \code{\linkS4class{metabarcoding.data}} instance
#'             on normalisation have to be computed.
#' @param MARGIN Indicates if the sums have to be computed across 
#'               samples or motus. 
#'               Allowed values are :
#'               \itemize{
#'                 \item{'sample' or 1} for computing sum across samples
#'                 \item{'motu' or 2} for computing sum across motus
#'                 }
#' @param threshold a numeric value between 0 and 1 indicating which part of 
#'                  the signal must be conserved. Default value is setup to
#'                  0.01 (1% of the normalized signal).
#' @param operator is a logical comparison operator.
#' 
#' @return A logical matrix usable for selecting cell in the read aboundancy matrix.
#'                   
#' @seealso \code{\linkS4class{metabarcoding.data}}, \code{\link{threshold.mask}}, \code{\link{normalize}}
#'
#' @docType methods
#' @rdname const-threshold-mask-methods
#' @aliases const.threshold.mask-methods,metabarcoding.data
#' @author Aurelie Bonin
#' 
setMethod("const.threshold.mask", "metabarcoding.data", function(data,MARGIN,threshold=0.01,operator='<') {
  
  
  if (MARGIN == 'sample')
    MARGIN=1
  
  if (MARGIN == 'motu')
    MARGIN=2
  
  readcount = normalize(data,MARGIN,as.matrix=TRUE)	
  
  mask = do.call(operator,list(readcount,threshold))
  
  return(mask)
})

#' @export
setGeneric("threshold.set", function(data,
				MARGIN,
				threshold=0.97,
				operator='<',
				value=0,
				normalize=TRUE,
				mask.fun=threshold.mask) {
			return(standardGeneric("threshold.set"))
		})


setMethod("threshold.set", "metabarcoding.data", function(data,
		MARGIN,
		threshold=0.97,
		operator='<',
		value=0,
		normalize=TRUE,
		mask.fun=threshold.mask) {
	
	
	
	if (MARGIN == 'sample')
		MARGIN=1
	
	if (MARGIN == 'motu')
		MARGIN=2
	
	readcount = reads(data)
	
	if (normalize)
		data = normalize(data,c(2,1)[MARGIN])
	
	mask = mask.fun(data,MARGIN,threshold,operator)
	
	readcount[mask] = value
	
	newdata = copy.metabarcoding.data(data,reads=readcount)
	
	return(newdata)
	
})
