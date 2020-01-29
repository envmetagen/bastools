#' @include 02_class_metabarcoding.data.R
NULL

# TODO: Add comment
# 
# Author: coissac
###############################################################################

#' @export
setGeneric("rarefy", function(x,n,first.pass=0.95,pseudo.count=0,...) {
			return(standardGeneric("rarefy"))
		})

setMethod("rarefy", "ANY", function(x,n,first.pass=0.95,pseudo.count=0,sum=NA) {
				
	if (is.na(sum))
		sum=sum(x)

	if (sum < sum(x))
		stop("sum parameter must be greater or equal to sum(x)")
  
  grey = sum-sum(x)
	
	probs = x + pseudo.count
  
  if (grey > 0)
	  probs = c(probs,grey)
	
	# Just to ensure at least one execution of the loop
	n1 = n * 2
	
	while(n1 > n)
		n1 = rpois(1,n * first.pass)
	
	rep1 = as.vector(rmultinom(1,n1,probs))
	n2  = sum(rep1)
	
	levels = 1:length(probs)
	
	rep2= as.vector(table(factor(sample(levels,
						   	                      n - n2,
							                        replace=TRUE, 
							                        prob = probs),
					                            levels=levels)))
	
	rep1 = (rep1 + rep2)
  
  if (grey > 0)
    rep1 = rep1[-length(rep1)]
  
	return(rep1)
})


setMethod("rarefy", "metabarcoding.data", function(x,n,first.pass=0.95,pseudo.count=0,MARGIN='sample') {
		
	if (MARGIN == 'sample')
		MARGIN=1
	
	if (MARGIN == 'motu')
		MARGIN=2
	
  dreads= dim(x@reads)
  rreads= matrix(0,nrow = dreads[1] , ncol = dreads[2])
	
	if (MARGIN == 1)
      for (i in 1:dreads[1]) {
        rreads[i,]=rarefy(x@reads[i,],
                          n=n,
                          first.pass=first.pass,
                          pseudo.count=pseudo.count)
      }
    
# 		rreads = t(apply(reads,1,rarefy,n=n,
# 					 	 first.pass=first.pass,
# 						 pseudo.count=pseudo.count))
	else
	  for (i in 1:dreads[2]) {
	    rreads[,i]=rarefy(x@reads[,i],
	                      n=n,
	                      first.pass=first.pass,
	                      pseudo.count=pseudo.count)
	  }

# rreads =   as.matrix(apply(reads,2,rarefy,n=n,
# 								   first.pass=first.pass,
# 								   pseudo.count=pseudo.count))

  rreads=as.matrix(rreads)
			
	rownames(rreads) = rownames(x@reads)
	colnames(rreads) = colnames(x@reads)
	
	newdata = copy.metabarcoding.data(x,reads=rreads)
	
	return(newdata)
	
})
	