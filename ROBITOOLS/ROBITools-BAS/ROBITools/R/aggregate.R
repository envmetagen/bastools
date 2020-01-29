#' @include 02_class_metabarcoding.data.R
NULL

# TODO: Add comment
# 
# Author: coissac
###############################################################################

#' @export
aggregate.metabarcoding.data=function(x, by, FUN,..., 
		                              MARGIN='sample',
									  default.layer=NULL,
									  layers=NULL) {	
  
	uniq.value = function(z) {
		
		if (is.null(z) | 
				any(is.na(z)) |
				length(z)==0)
			ans = NA
		else {
      if (all(z==z[1]))
			  ans = z[1]
		  else
			  ans = NA
		}
		if (is.factor(z))
			ans = factor(ans,levels=levels(z))
		
		return(ans)
	}			
	
	
	#
	# Deals with the supplementaty aggregate arguments
	#
			
	if (is.null(default.layer))
		default.layer=uniq.value
	
	
	if (is.null(layers)) {
		layers = as.list(rep(c(default.layer),length(x@layers)))
		names(layers)=layer.names(x)
	}
	else {
		for (n in layer.names(x))
			if (is.null(layers[[n]]))
				layers[[n]]=default.layers
	}
	
	if (MARGIN == 'sample')
		MARGIN=1
	
	if (MARGIN == 'motu')
		MARGIN=2
	
	reads = x@reads
	  
	if (MARGIN==1) {
		# prepare the aggrevation arguments for the read table
		# from the function arguments
	  dotted = list(...)
	  if (length(dotted) > 0)
	    aggr.args = list(reads,by=by,FUN=FUN,...=dotted,simplify=FALSE)
	  else
	    aggr.args = list(reads,by=by,FUN=FUN,simplify=FALSE)
	  		
		# Aggregate the read table
		ragr = do.call(aggregate,aggr.args)
		
		# extrat new ids from the aggregated table
		ncat = length(by)
		ids  = as.character(interaction(ragr[,1:ncat,drop=FALSE]))
		
		# remove the aggregations modalities to rebuild a correct
		# reads table
		ragr = as.matrix(ragr[,-(1:ncat),drop=FALSE])
		dragr= dim(ragr)
		cragr= colnames(ragr)
		ragr = as.numeric(ragr)
		dim(ragr)=dragr
		colnames(ragr)=cragr		
		rownames(ragr)=ids
		
		#
		# Apply the same aggragation to each layer
		#
		
		ln = layer.names(x)
		
    la = vector(mode="list",length(ln))
    names(la)=ln
    
		for (n in ln) {
			f = layers[[n]]
			if (is.factor(x[[n]])){
			  isfact = TRUE
			  lf = levels(x[[n]])
			  df = dim(x[[n]])
			  m = matrix(as.character(x[[n]]))
			  dim(m)=df
			}
			else
			  m = x[[n]]

      aggr.args = list(m,by=by,FUN=f,simplify=FALSE)
			lagr = do.call(aggregate,aggr.args)
			lagr = as.factor.or.matrix(lagr[,-(1:ncat),drop=FALSE])

			if (isfact){
			  df = dim(lagr)
			  lagr = factor(lagr,levels=lf)
			  dim(lagr)=df
			}
			
      rownames(lagr)=ids
			la[[n]]=lagr
		}
		
		# aggragate the sample table according to the same criteria
		#
		# TODO: We have to take special care of factors in the samples
		#       data.frame
		
		sagr = aggregate(samples(x),by,uniq.value,simplify=FALSE)
		
		# move the first columns of the resulting data frame (the aggregations
		# modalities to the last columns of the data.frame
		sagr = sagr[,c((ncat+1):(dim(sagr)[2]),1:ncat),drop=FALSE]
		larg = c(lapply(sagr,unlist),list(stringsAsFactors=FALSE))	
		sagr = do.call(data.frame,larg)
		
		# set samples ids to the ids computed from modalities
		sagr$id=ids
		rownames(sagr)=ids
		
		# build the new metabarcoding data instance
		newdata = copy.metabarcoding.data(x,reads=ragr,samples=sagr)
		
	}
	else {
		# prepare the aggregation arguments for the read table
		# from the function arguments
		# BECARFUL : the reads table is transposed
		#            standard aggregate runs by row and we want
		# 			 aggregation by column
    
    dotted = list(...)
    if (length(dotted) > 0)
		  aggr.args = list(t(reads),by=by,FUN=FUN,...=dotted,simplify=FALSE)
    else
      aggr.args = list(t(reads),by=by,FUN=FUN,simplify=FALSE)
    
		
    # Aggregate the read table
		ragr = do.call(aggregate.data.frame,aggr.args)
		
		# extrat new ids from the aggregated table
		ncat = length(by)
		ids  = as.character(interaction(ragr[,1:ncat,drop=FALSE]))
		
		# remove the aggregations modalities to rebuild a correct
		# reads table
    
		ragr = t(ragr[,-(1:ncat),drop=FALSE])		
    dragr= dim(ragr)
    rragr= rownames(ragr)
    ragr = as.numeric(ragr)
    dim(ragr)=dragr
    colnames(ragr)=ids
    rownames(ragr)=rragr
		
		#
		# Apply the same aggragation to each layer
		#
		
		ln = layer.names(x)
    
    la = vector(mode="list",length(ln))
    names(la)=ln
    
    for (n in ln) {
			f = layers[[n]]
      
      if (is.factor(x[[n]])){
			  isfact = TRUE
        lf = levels(x[[n]])
        df = dim(x[[n]])
        m = matrix(as.character(x[[n]]))
        dim(m)=df
      }
      else
          m = x[[n]]

      aggr.args = list(t(m),by=by,FUN=f,simplify=FALSE)
			lagr = do.call(aggregate,aggr.args)
			lagr = t(as.factor.or.matrix(lagr[,-(1:ncat),drop=FALSE]))
      
      if (isfact){
        df = dim(lagr)
        lagr = factor(lagr,levels=lf)
        dim(lagr)=df
			}

      colnames(lagr)=ids
			la[[n]]=lagr
		}
    
		# aggragate the motus table according to the same criteria
		magr = aggregate(motus(x),by,uniq.value,simplify=FALSE)
		
		# move the first columns of the resulting data frame (the aggregations
		# modalities to the last columns of the data.frame
		magr = magr[,c((ncat+1):(dim(magr)[2]),1:ncat),drop=FALSE]
		larg = c(lapply(magr,unlist),list(stringsAsFactors=FALSE))	
		magr = do.call(data.frame,larg)
		
		# set motus ids to the ids computed from modalities
		magr$id=ids
		rownames(magr)=ids
		
		# build the new metabarcoding data instance
		newdata = copy.metabarcoding.data(x,reads=ragr,motus=magr,layers=la)
	}
	
	return(newdata)			
}

