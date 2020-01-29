#' @include ROBITools.R
#' @include s3objects.R
#' @import ROBITaxonomy 
NULL

require(ROBITaxonomy)

#
# FOR THE DEVELOPPER : we have to check that the code doesn't relies on the
#                      fact that the xx@samples$sample column is not always 
#                      identical to the rownames(xx@samples)

setClassUnion("characterOrNULL",c("character","NULL"))
setClassUnion("matrixOrfactorL",c("matrix","factor"))

#
# We specialize data.frame in two subclasses motus.frame and samples.frame
# for this we add to function insuring the type checking and the cast from
# data.frame
#

is.motus.frame= function(x) any(class(x)=="motus.frame")
is.samples.frame= function(x) any(class(x)=="samples.frame")

as.motus.frame= function(x) {
  if (! is.data.frame(x))
    stop("only cast from data.frame is allowed")
  if (! is.motus.frame(x))
    x = addS3Class(x,"motus.frame")
  
  return(x)
}


as.samples.frame= function(x) {
  if (! is.data.frame(x))
    stop("only cast from data.frame is allowed")
  if (! is.samples.frame(x))
    x = addS3Class(x,"samples.frame")
  return(x)
}

samples.frame=as.samples.frame
motus.frame=as.motus.frame

as.factor.or.matrix = function(x) {
  if (is.matrix(x))
    return(x)
  
  if (is.factor(x)){
    if (length(dim(x))!=2)
      stop('Just factor with two dimensions are allowed')
    return(x)
  }
  
  if (!is.data.frame(x))
    stop('Just matrix, 2D factor and data.frame can be casted')
  
  tps   = sapply(x,class)
  allna = sapply(x, function(y) all(is.na(y)))
  
  if (all(tps==tps[[1]] | allna)) {
    tps = tps[[1]] 
  }
  else
    stop('all the column of the data.frame must have the same type')

  tps = tps[[1]]
  
  x = as.matrix(x)
  dx = dim(x)
  if (tps=='factor')
    x = factor(x)
    dim(x)=dx
  
  return(x)
}

#' DNA metabarcoding experiment description class
#' 
#' A S4 class describing a DNA metabarcoding experiment. It groups
#' three data frames describing samples, motus and occurrences of 
#' MOTUs per sample
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{reads}:}{Matrix of class \code{"numeric"},  
#'                         containing the counts of reads per samples
#'                         \itemize{
#'                           \item{1 samples per line}
#'                           \item{1 sequence per column}
#'                         }
#'                        }
#'                        
#'    \item{\code{samples}:}{Object of class \code{"data.frame"}, describing samples
#'                         \itemize{
#'                           \item{1 samples per line}
#'                           \item{1 property per column}
#'                         }
#'                       }
#'                       
#'    \item{\code{motus}:}{Object of class \code{"data.frame"}, describing MOTUs (sequences)
#'                         \itemize{
#'                           \item{1 MOTU per line}
#'                           \item{1 property per column}
#'                         }
#'                       }
#'                       
#'    \item{\code{layers}:}{Object of class \code{"list"}, containing a set of data layers 
#'                          linking motus and samples. Each element of the list is a matrix
#'                          of the same size than the \code{read} slot with
#'                         \itemize{
#'                           \item{1 samples per line}
#'                           \item{1 sequence per column}
#'                         }
#'                       }
#'        
#'    \item{\code{scount}:}{Object of class \code{"integer"}, containing the count of sample}             
#'        
#'    \item{\code{mcount}:}{Object of class \code{"integer"}, containing the count of MOTUs}             
#'        
#'    \item{\code{sample.margin}:}{Vector of class \code{"numeric"},  describing the total count of 
#'                          sequence per sample. By default this slot is set by applying sum
#'                         to the reads data.frame lines}             
#'
#'    \item{\code{taxonomy}:}{Object of class \code{"taxonomy.obitools"}, linking the DNA metabarcoding
#'                            experiment to a taxonomy}             
#'        
#'    \item{\code{taxid}:}{Vector of class \code{"character"}, list of MOTUs' attributes to manage as taxid}             
#'  }
#'
#' @seealso \code{\link{taxonomy.obitools}},
#' @name metabarcoding.data
#' @rdname metabarcoding-data-class
#' @keywords DNA metabarcoding
#' @author Eric Coissac
#' @exportClass metabarcoding.data

setClass("metabarcoding.data",
		
		
		#
		# Attribute declaration
		#
		
		 representation(reads         = "matrix", 
						        samples       = "data.frame",
  		  		        motus         = "data.frame",
        						layers        = "list",
    				        scount        = "integer",
        						mcount        = "integer",
 				 		        sample.margin = "numeric",
 						        taxonomy      = "obitools.taxonomyOrNULL",
						        taxid         = "characterOrNULL"
						),
				
		#
		# Check object structure 
		#
				
		validity = function(object) { 
			
						## object : nom reserve !
			
			#
			# Check that reads / samples and motus data.frames
			# have compatible sizes
			#
			#   reads line count   = samples line count
			#   reads column count = motus   line count
			
			rsize = dim(object@reads)
			ssize = dim(object@samples)
			msize = dim(object@motus)
			csize = length(object@sample.margin)
			
			if (rsize[1] != ssize[1] & 
				rsize[2] != msize[1] &
				rsize[1] != csize)
				return(FALSE)
			
			
			# if no layer, object is ok
			
			if (length(object@layers)==0)
				return(TRUE)
			
			# otherwise we check the size of each layer as we
			# did for reads
			
			return(! any(sapply(object@layers, 
					   function(l) any(dim(l)!=c(ssize[1],msize[1])))))
			
		}
)



#
#' metabarcoding.data constructor
#' 
#' @docType methods
#' @rdname initialize-methods
#' @aliases initialize-methods,metabarcoding.data
setMethod("initialize",
		  "metabarcoding.data",
			function(.Object, reads,samples,motus,
					          taxonomy=NULL,taxid=NULL,
							  sample.margin=NA,
							  layers=list()) {
        
        rn = rownames(reads)
        cn = colnames(reads)
        
				.Object@reads   <- reads
        
				# .Object@samples <- as.samples.frame(samples)
        .Object@samples <- samples
        row.names(.Object@samples) = rn
        
				#.Object@motus   <- as.motus.frame(motus)
        .Object@motus   <- motus
        row.names(.Object@motus) = cn
        
        
        # Set colnames and rownames to each layers
        layers = lapply(layers, function(x) {colnames(x)=cn 
                                             rownames(x)=rn 
                                             return(x)})
		.Object@layers   <- layers
        				
				# Precompute sample count and motu count
				
		.Object@scount = dim(.Object@samples)[1]
		.Object@mcount = dim(.Object@motus)[1]
				
		.Object@taxonomy = taxonomy
		.Object@taxid = taxid
				
		if (is.null(sample.margin))
			.Object@sample.margin = rowSums(reads)
		else
			.Object@sample.margin = sample.margin
        
        names(.Object@sample.margin) = rn
				
		validObject(.Object) ## valide l'objet
        
		return(.Object)
		})


#
# metabarcoding.data getters
#

#' @export
setGeneric("reads", function(obj) {
			return(standardGeneric("reads"))
		})

#' Extracts the matrix describing MOTUs abondances
#' 
#' Extract the the matrix describing MOTUs abondances (read counts) 
#' from a \code{\link{metabarcoding.data}} instance.
#' 
#' @param   obj a \code{\link{metabarcoding.data}} instance
#' @return  a matrix containing data about reads 
#' 
#' @examples
#' # load termite data set from the ROBITools sample data
#' data(termes)
#' 
#' # Extract the matrix describing MOTUs abondances
#' d = reads(termes)
#' 
#' head(d)
#' 
#' @seealso \code{\link{metabarcoding.data}},
#'          \code{\link{motus}}, \code{\link{samples}}
#'          
#' @docType methods
#' @rdname read-methods
#' @aliases read-methods,metabarcoding.data
#' @author Eric Coissac
#' 
setMethod("reads", "metabarcoding.data", function(obj) {
			return(obj@reads)
		})


# get samples data.frames

#' @export
setGeneric("samples", function(obj) {
			return(standardGeneric("samples"))
		})

#' Extracts the samples description data.frame
#' 
#' Extract the sample description data.frame from a 
#' \code{\link{metabarcoding.data}} instance.
#' 
#' @param   obj a \code{\link{metabarcoding.data}} instance
#' @return  a data.frame containing data about sample 
#' 
#' @examples
#' # load termite data set from the ROBITools sample data
#' data(termes)
#' 
#' # Extract the data frame describing samples
#' d = samples(termes)
#' 
#' head(d)
#' 
#' @seealso \code{\link{metabarcoding.data}},
#'          \code{\link{motus}}, \code{\link{reads}}
#'          
#' @docType methods
#' @rdname samples-methods
#' @aliases samples-methods,metabarcoding.data
#' @author Eric Coissac
#' 
setMethod("samples", "metabarcoding.data", function(obj) {
			return(obj@samples)
		})


#' @export
setGeneric("motus", function(obj) {
			return(standardGeneric("motus"))
		})

#' Extracts the MOTU descriptions \code{data.frame}
#' 
#' Extract the MOTUs description \code{data.frame} from a 
#' \code{\link{metabarcoding.data}} instance.
#' 
#' @param   obj a \code{\link{metabarcoding.data}} instance
#' @return  a data.frame containing data about MOTU 
#' 
#' @examples
#' # load termite data set from the ROBITools sample data
#' data(termes)
#' 
#' # Extract the data.frame describing MOTUs
#' d = motus(termes)
#' 
#' head(d)
#' 
#' @seealso \code{\link{metabarcoding.data}},
#'          \code{\link{reads}}, \code{\link{samples}}
#'          
#' @docType methods
#' @rdname motu-methods
#' @aliases motu-methods,metabarcoding.data
#' 
setMethod("motus", "metabarcoding.data", function(obj) {
			return(obj@motus)
		})


# get sample count 

setGeneric("sample.count", function(obj) {
			return(standardGeneric("sample.count"))
		})

setMethod("sample.count", "metabarcoding.data", function(obj) {
			return(obj@scount)
		})

# get motu count 

setGeneric("motu.count", function(obj) {
			return(standardGeneric("motu.count"))
		})

setMethod("motu.count", "metabarcoding.data", function(obj) {
			return(obj@mcount)
		})

# dim method

setMethod("dim", "metabarcoding.data", function(x) {
			return(c(x@scount,x@mcount))
		})


setMethod('[', "metabarcoding.data", function(x,i=NULL,j=NULL,...,drop=TRUE) {
			
      # special case if samples are not specified (dimension 1)
			if (!hasArg(i))
				i = 1:x@scount
			
			# special case if motus are not specified (dimension 2)
			if (!hasArg(j))
				j = 1:x@mcount

      # special case if the layer attribut is specified
			args = list(...)
			
			if (!is.null(args$layer))
          return(x[[args$layer]][i,j])	
      
      #####################
      #
      # normal case 
      #
			  
			r = x@reads[i,j,drop=FALSE]
      
			if (sum(dim(r) > 1)==2 | ! drop)
			{
        
        # we do the selection on the motus and samples description data.frame
        
				m = x@motus[j,,drop=FALSE]
				s = x@samples[i,,drop=FALSE]
        
        # we do the selection on each layers
        l = lapply(x@layers,function(l) l[i,j,drop=FALSE])
        
				newdata = copy.metabarcoding.data(x, reads=r, samples=s, motus=m, layers=l)
			}
			else
			{
				newdata = as.numeric(x@reads[i,j])
			}
			
			return(newdata)
			
		})

setMethod('[<-', "metabarcoding.data",
          function (x, i, j, ..., value) {
            if (!hasArg(i))
              i = 1:x@scount
            
            if (!hasArg(j))
              j = 1:x@mcount
            
            args = list(...)
            
            if (is.null(args$layer))
              x@reads[i, j]=value
            else
              
              x[[args$layer]][i,j]=value
            
            return(x)
          })

	

#################################################
#
# User interface function to create 
# metabarcoding.data objects
#
#################################################

#'@export
metabarcoding.data = function(reads,samples,motus,
							  taxonomy=NULL,taxid=NULL,
							  sample.margin=NULL,
                layers=list()) {
	rd = new('metabarcoding.data',
			reads=reads,
			samples=samples,
			motus=motus,
			taxonomy=taxonomy,
			taxid=taxid,
			sample.margin=sample.margin,
      layers=layers
	)
	
	return(rd)
}

copy.metabarcoding.data = function(data, 
                                   reads=NULL,
                                   samples=NULL,motus=NULL,
										               taxonomy=NULL,taxid=NULL,
										               sample.margin=NULL,
                                   layers=NULL) {
    
    
    
		if (is.null(reads))
			reads = data@reads
		
		if (is.null(samples))
			samples = data@samples
		
		if (is.null(motus))
			motus = data@motus
		
		if (is.null(taxonomy))
			taxonomy = data@taxonomy
		
		if (is.null(taxid))
			taxid = data@taxid
		
		if (is.null(sample.margin))
			sample.margin = data@sample.margin
		
		if (is.null(layers))
		  layers = data@layers
		
		
		rd = new('metabarcoding.data',
				reads=reads,
				samples=samples,
				motus=motus,
				taxonomy=taxonomy,
				taxid=taxid,
				sample.margin=sample.margin,
        		layers=layers
		)
		
		return(rd)
}

#' @export
setGeneric('rownames')

#' @export
setMethod("rownames", "metabarcoding.data", function(x, do.NULL = TRUE, prefix = "col") {
  return(rownames(x@reads,do.NULL,prefix))
})

#' @export
setGeneric('colnames')

#' @export
setMethod("colnames", "metabarcoding.data", function(x, do.NULL = TRUE, prefix = "col") {
  return(colnames(x@reads,do.NULL,prefix))
})
