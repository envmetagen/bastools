#' @include 02_class_metabarcoding.data.R
#' @import igraph
NULL

require(igraph)

# pos = expand.grid(x,y)

#' Computes the pairwise distance matrix as a data.frame where
#' 
#' @param x a vector for the X coordinates
#' @param y a vector for the Y coordinates
#' @param labels a vector with the sample names
#' 
#' @return a data.frame instance of three columns
#'    - a : The label of the first sample
#'    - b : The label of the second sample
#'    - dist : The euclidian distance beween sample a and b
#' 
#' @examples
#' data(termes)
#' termes.ok = termes[,colSums(termes$reads)>0]
#' pos = expand.grid(1:3 * 10,1:7 * 10)
#' labels = rownames(termes.ok)
#' d = dist.grid(pos[,1],pos[2],labels)
#' 
#' @export
dist.grid = function(x,y,labels=NULL){
  pos = data.frame(x,y)
  
  if (is.null(labels))
    labels = as.character(interaction(pos))
  else
    labels = as.character(labels)

  llabels=length(labels)
  dpos=dist(pos)

  a = rep(labels[1:(llabels-1)],(llabels-1):1)
  b = do.call(c,(lapply(2:llabels, function(i) labels[i:llabels])))

  return(data.frame(a,b,dist=as.vector(dpos)))
}

#' Builds the list of sample groups included in a circle around a central sample
#' 
#' @param dtable a distance table between samples as 
#'               computed by \code{\link{dist.grid}}
#' @param radius the radius of the circle
#' @param center a \code{logical} value indicating if the center of
#'               the group must be included in the group
#'               
#' @return a list of vectors containing the labels of the group members 
#'          
#' @examples
#' data(termes)
#' termes.ok = termes[,colSums(termes$reads)>0]
#' pos = expand.grid(1:3 * 10,1:7 * 10)
#' labels = rownames(termes.ok)
#' d = dist.grid(pos[,1],pos[2],labels)
#' groups = dist.center.group(d,20)
#' 
#' @export
dist.center.group=function(dtable,radius,center=TRUE) {
  
  fgroup = function(c) {
    ig = dtable[(dtable[,1]==c | dtable[,2]==c) & dtable[,3] <= radius,]
    return(union(ig[,1],ig[,2]))
  }
  
  pos = as.character(union(dtable[,1],dtable[,2]))
  
  g = lapply(pos,fgroup)
  names(g) = pos
  
  if (!center)
    g = mapply(setdiff,g,pos)
  
  return(g)
  
}

#' Builds the list of sample groups including samples closest than a define distance
#' 
#' A graph is build by applying the threshold \code{dmax} to the distance matrix
#' A group is a clique max in this graph. Consequently all member pairs of a group
#' are distant by less or equal to \code{dmax}.
#' 
#' @param dtable a distance table between samples as 
#'               computed by \code{\link{dist.grid}}
#' @param dmax the maximum distance between two samples
#'               
#' @return a list of vectors containing the labels of the group members 
#'          
#' @examples
#' data(termes)
#' termes.ok = termes[,colSums(termes$reads)>0]
#' pos = expand.grid(1:3 * 10,1:7 * 10)
#' labels = rownames(termes.ok)
#' d = dist.grid(pos[,1],pos[2],labels)
#' groups = dist.clique.group(d,20)
#' 
#' @export
dist.clique.group=function(dtable,dmax,center=True) {
  gp = igraph::graph.edgelist(as.matrix(dtable[dtable$dist <= dmax,c('a','b')]),directed=FALSE)
  g  = igraph::maximal.cliques(gp)
  return(lapply(g, function(i) igraph::V(gp)$name[i]))
}

#' Computes the univariate M statistics 
#' 
#' @param w the weigth matrix indicating the presence probability of each motu
#'          in each samples. Each line corresponds to a sample and each column
#'          to a MOTU. \code{rownames} of the \code{w} matrix must be the sample
#'          names.  It is nice but not mandatory if the \code{colnames} refer to the MOTU id.
#'          
#' @param groups the list of considered groups as computed by the \code{\link{dist.center.group}}
#'               function
#'
#' @seealso \code{\link{dist.center.group}}
#' @seealso \code{\link{m.weight}}
#' 
#' @examples
#' data(termes)
#' termes.ok = termes[,colSums(termes$reads)>0]
#' pos = expand.grid(1:3 * 10,1:7 * 10)
#' labels = rownames(termes.ok)
#' d = dist.grid(pos[,1],pos[2],labels)
#' groups = dist.center.group(d,20)
#' w = m.weight(termes.ok)
#' m = m.univariate(w,groups)
#' 
#' @references Marcon, E., Puech, F., and Traissac, S. (2012). 
#'             Characterizing the relative spatial structure of point patterns. 
#'             International Journal of Ecology, 2012.
#'             
#' @export
m.univariate = function(w,groups) {
  
  nunivar = function(members,center) {
    g = w[members,]
    
    wn = colSums(g)
    wa = sum(wn)
 
    wn = wn - center
    wa = wa - center
    
    p = wn / wa * center
        
    return(p)
  }

  centers = lapply(names(groups),function(x) w[x,])
  
  Wf = colSums(w)
  Wa = sum(Wf)

  Denom.univar = colSums(w * (sweep(-w,2,Wf,'+') / (Wa - w)))
  Num.univar   = rowSums(mapply(nunivar,groups,centers))
  
  Munivar=Num.univar/Denom.univar
  Munivar[Denom.univar==0]=0
  
  return(Munivar)
}


#' Computes the bivariate M statistics 
#' 
#' The function computes the bivariate M statiscics for a set of target species around a set of
#' focus species.
#' 
#' @param w1 the weigth matrix indicating the presence probability of each motu
#'          used as focus species in each samples. Each line corresponds to a sample and each column
#'          to a MOTU. \code{rownames} of the \code{w} matrix must be the sample
#'          names. It is nice but not mandatory if the \code{colnames} refer to the MOTU id.
#'          
#' @param w2 the weigth matrix indicating the presence probability of each motu
#'          used as target species in each samples. Each line corresponds to a sample and each column
#'          to a MOTU. \code{rownames} of the \code{w} matrix must be the sample
#'          names. It is nice but not mandatory if the \code{colnames} refer to the MOTU id.
#'          if \code{w2} is not set, w1 is also used as target species. in this case the diagonal
#'          of the matrix return contains the univariate M statistic for the diferent species.
#'          
#' @param groups the list of considered groups as computed by the \code{\link{dist.center.group}}
#'               function
#'               
#' @return a matrix of M bivariate statistics with one focus species by row and one target species 
#'         by columns If \code{w2} is not specified the diagonal of the matrix is equal to the univariate
#'         M statistic of the corresponding species.
#'
#' @seealso \code{\link{dist.center.group}}
#' @seealso \code{\link{m.weight}}
#' 
#' @examples
#' data(termes)
#' termes.ok = termes[,colSums(termes$reads)>0]
#' pos = expand.grid(1:3 * 10,1:7 * 10)
#' labels = rownames(termes.ok)
#' d = dist.grid(pos[,1],pos[2],labels)
#' groups = dist.center.group(d,20)
#' w = m.weight(termes.ok)
#' m = m.bivariate(w,groups)
#' 
#' @references Marcon, E., Puech, F., and Traissac, S. (2012). 
#'             Characterizing the relative spatial structure of point patterns. 
#'             International Journal of Ecology, 2012.
#'             
#' @export
m.bivariate = function(w1,w2=NULL,groups) {
  
  nunbivar = function(members,center) {
    g = w2[members,]
    
    wn = colSums(g)
    wa = sum(wn)
    
    if (self){
      mwn  = wn %*% t(rep(1,length(wn)))
      diag(mwn)= wn - center
      wa = wa - center
      wna = mwn/wa
      p = sweep(wna,2,center,'*')
      #p = center %*% wna
    }
    else {
      wna= matrix(wn/wa,nrow=1)
      p = center %*% wna
    }    
    
    return(p)
  }
  
  if (is.null(w2)){
    self = TRUE
    w2=w1
  }
  else {
    self = FALSE
  }
  
  centers = lapply(names(groups),function(x) w[x,])

  Wf = colSums(w1)
  Wn = colSums(w2)
  Wa = sum(Wn)
  
  if (self){      
    Wn = sweep(-w1,2,Wn,'+')
    Wna = Wn/(Wa - w1)
    Denom.bivar = t(w1) %*% Wna
  }
  else {
    Wna= t(Wn/Wa)
    Denom.bivar = Wf %*% Wna    
  }
  
  Num.bivar = matrix(0,nrow=ncol(w1),ncol=ncol(w2))

  ng = length(groups)
  
  for (i in 1:ng) {
    Num.bivar = Num.bivar + nunbivar(groups[[i]],centers[[i]])
  }
  
  Mbivar=Num.bivar/Denom.bivar

  Mbivar[Denom.bivar==0]=0
  
  return(Mbivar)
}

#' Computes a weigth matrix from a \code{\linkS4class{metabarcoding.data}}
#' 
#' The weight can be considered as a propability of presence of a MOTU in a 
#' given sample. This function defines this probability as the fraction of
#' the maximal occurrence frequency over all samples. 
#' 
#' @param data a \code{\linkS4class{metabarcoding.data}} instance
#' 
#' @return a weight matrix usable for M statistics
#' 
#' @examples
#' data(termes)
#' termes.ok = termes[,colSums(termes$reads)>0]
#' w = m.weight(termes.ok)
#' 
#' @export
m.weight = function(data) {
  ndata = normalize(data,MARGIN='sample')
  fmax=apply(ndata$reads,2,max)
  w = sweep(ndata$reads,2,fmax,'/')
  rownames(w)=rownames(ndata)
  colnames(w)=colnames(ndata)
  return(w)
}

#' Simulate null distribion of the M statistics by Monte-Carlo
#' 
#' Computes the null empirical distribution of the M statistics
#' by shuffling MOTUs among location.
#' 
#' @param w the weigth matrix indicating the presence probability of each motu
#'          in each samples. Each line corresponds to a sample and each column
#'          to a MOTU. \code{rownames} of the \code{w} matrix must be the sample
#'          names.
#' @param groups the list of considered groups as computed by the \code{\link{dist.center.group}}
#'               function
#' @param resampling the number of simulation to establish the null distribution
#' 
#' @return a matrix of M score under the null hypothesis of random distribution of MOTUs
#'         with a MOTUs per line and a culumn per simulation
#'
#' @examples       
#' data(termes)
#' termes.ok = termes[,colSums(termes$reads)>0]
#' pos = expand.grid(1:3 * 10,1:7 * 10)
#' labels = rownames(termes.ok)
#' d = dist.grid(pos[,1],pos[2],labels)
#' groups = dist.center.group(d,20)
#' w = m.weight(termes.ok)
#' dnull = dm.univariate(w,groups)
#' 
#' @export              
dm.univariate = function(w,groups,resampling=100) {
  
  shuffle = function(w){
    wr =apply(w,2,function(y) sample(y,length(y),replace=FALSE))
    rownames(wr)=rownames(w)
    return(wr)
  }
  
  msim = function(x) {
    return(m.univariate(shuffle(w),groups))
  }
  
  dnull = mapply(msim,1:resampling)
  
  rownames(dnull) = colnames(w)
  
  return(dnull)
}

#' Test the significance of the M statistics by Monte-Carlo
#' 
#' Computes computes the p.value the M statistics asociated to a MOTU
#' by shuffling MOTUs among location.
#' 
#' @param w the weigth matrix indicating the presence probability of each motu
#'          in each samples. Each line corresponds to a sample and each column
#'          to a MOTU. \code{rownames} of the \code{w} matrix must be the sample
#'          names.
#' @param groups the list of considered groups as computed by the \code{\link{dist.center.group}}
#'               function
#' @param resampling the number of simulation to establish the null distribution
#' 
#' @param alternative a character value in \code{c('two.sided','less','greater')}
#'               - two.sided : the m stat is check against the two side of the empirical
#'                             M distribution
#'               - less : test if the M stat is lesser than the M observed in the the empirical
#'                             M distribution (exlusion hypothesis)
#'               - greater : test if the M stat is greater than the M observed in the the empirical
#'                             M distribution (aggregation hypothesis)
#' 
#' @return a vector of p.value with an attribute \code{m.stat} containing the actual M stat
#'         for each MOTUs
#'
#' @examples       
#' data(termes)
#' termes.ok = termes[,colSums(termes$reads)>0]
#' pos = expand.grid(1:3 * 10,1:7 * 10)
#' labels = rownames(termes.ok)
#' d = dist.grid(pos[,1],pos[2],labels)
#' groups = dist.center.group(d,20)
#' w = m.weight(termes.ok)
#' pval = m.univariate.test(w,groups)
#' 
#' @export              
m.univariate.test = function(w,groups,resampling=100,alternative='two.sided') {
  dnull = dm.univariate(w,groups,resampling)
  m = m.univariate(w,groups)
  pnull = sapply(1:dim(dnull)[1],function(y) 1 - ecdf(dnull[y,])(m[y]))
  
  p.value=NULL
  
  if (alternative=='two.sided') {
    p.value = mapply(min,pnull,1 - pnull)
  }
  
  if (alternative=='less') {
    p.value = pnull
  }

  if (alternative=='greater') {
    p.value = 1 - pnull
  }
  
  # Set p.value to 1 if the MOTU occurres in only one place
  n = colSums(w > 0)
  p.value[n==1]=1
  
  names(p.value) = colnames(w)
  attr(p.value,'m.stat')=m
  
  return(p.value)
}
