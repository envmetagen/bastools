#' @include 02_class_metabarcoding.data.R
NULL

#' Read frequencies krigging
#' 
#' Extrapolates read frequencies from a \code{\link{metabarcoding.data}} object in space for a finer resolution 
#' 
#' @param x a vector or matrix from a row-normalized read table 
#'          \code{\link{metabarcoding.data}} object
#' @param min.coord a vector of length = 2 indicating the minimum values of x and y
#'                  coordinates to be used for the predicted grid
#' @param max.coord a vector of length = 2 indicating the maximum values of x and y
#'                  coordinates to be used for the predicted grid
#' @param grid.grain an integer indicating the resolution (i.e. nb of subpoints) in x and y
#'                  coordinates required for the predicted grid
#' @param coords a dataframe containing the x and y coordinates of the abundances 
#'               from x to be extrapolated.
#' @param otus.table a motus data.frame containing motus informations of x
#' @param cutoff a cutoff below which abundances are set to 0. 
#'               This threshold also determines the value to be added to 0 values for log10 
#'               transformation
#' @param return.metabarcoding.data if \code{TRUE}, returns a \code{\link{metabarcoding.data}} object. Default is \code{FALSE}
#'
#' @return either a dataframe or a S3 object with a structure similar to \code{\link{metabarcoding.data}} object. 
#'         The number of samples corresponds to the predicted points.
#'         The two last columns (if \code{return.metabarcoding.data==F}) or sample data.frame contains x y coordinates of the predicted grid
#'         The all but last two columns (if \code{return.metabarcoding.data==F}) or read matrix contains the predicted log10 transformed relative abundances 
#'         instead of reads counts
#'         If \code{return.metabarcoding.data==F} the motus data.frame contains the motus informations from x
#'
#' @examples
#' 
#' data(termes)
#' #Create dummy spatial coordinates
#' attr(termes, "samples")[c("x", "y")] = expand.grid(1:7,1:3)
#' 
#' #compute frequencies
#' attr(termes, "layers")[["reads.freq"]] = normalize(termes, MARGIN=1)$reads
#' 
#' # Getting extrapolations
#' termes.pred = extrapol.freq(attr(termes, "layers")[["reads.freq"]], min.coord=c(1,1), max.coord=c(7,3), 
#'                             grid.grain=100,termes$samples[,c("x", "y")], termes$motus, cutoff=1e-3)
#'
#' head(termes.pred$reads)
#' @seealso \code{\link{map.extrapol.freq}} as well as \code{sp} and \code{gstat} packages
#' @author Lucie Zinger
#' @export

extrapol.freq = function(x, min.coord, max.coord, grid.grain=100, coords, otus.table, cutoff=1e-3, return.metabarcoding.data = FALSE) {
  require(gstat)
  require(sp)
  
  #predicted grid setting
  new.x = seq(min.coord[1], max.coord[1], length.out = grid.grain)
  new.y = seq(min.coord[2], max.coord[2], length.out = grid.grain)
  grid.p=expand.grid(new.x, new.y)
  colnames(grid.p)=c("x", "y")
  S=sp::SpatialPoints(grid.p); sp::gridded(S)<-TRUE
  m=gstat::vgm(50, "Exp", 100)
  
  #krigging
  preds = apply(x, 2, function(otu) {
    otu[otu<cutoff] = cutoff
    spj=cbind(coords,otu)
    colnames(spj)=c("x", "y", "otu")
    spj.g=gstat::gstat(id="Log10.freq", formula=log10(otu)~1,locations=~x+y,data=spj,model=m)
    gstat::predict.gstat(spj.g, grid.p, quiet=T)$Log10.freq.pred
  })
  
  #formatting the output
  colnames(preds) = rownames(otus.table)
  rownames(preds) = paste("s", 1:nrow(grid.p), sep=".")
  row.names(grid.p) = rownames(preds)
  
  if(return.metabarcoding.data==F) {
    out = data.frame(preds, grid.p)
  } else{ 
    out = metabarcoding.data(preds, grid.p, otus.table)
  }
  return(out)
}


#' Maps of krigged log10-transformed frequencies
#' 
#' Maps the output of extrapol.freq
#' 
#' 
#' @param x an extrapol.freq output
#' @param path the path of the folder to export the map. Default is \code{NULL} and map is printed in Rplot/quartz
#' @param col.names a vector containing the names of the columns to be used for defining the file name. Typically
#'                  the column names containing the taxonomic information and/or sequence/motus id.
#' @param index  an integer indicating column number of the motu/sequence to be plotted.
#' @param cutoff  lower motu frequency accepted to consider motu abundance as different
#'                from 0. Should be the same than the one used in extrapol.freq
#' @param add.points a 3-column data.frame containing factor levels and associated x and y coordinates
#'                   to be added to the map. Typically taxa observed in the field.
#' @param adj a value used for adjusting text position in the map. Default is \code{4}
#'
#' @return a map/png file displaying motus distribution.
#'
#' @examples
#' 
#' data(termes)
#' attr(termes, "samples")[c("x", "y")] = expand.grid(1:7,1:3)
#' 
#' #compute frequencies
#' attr(termes, "layers")[["reads.freq"]] = normalize(termes, MARGIN=1)$reads
#' 
#' # Getting extrapolations
#' termes.pred = extrapol.freq(attr(termes, "layers")[["reads.freq"]], 
#' grid.grain=100,termes$samples[,c("x", "y")], termes$motus, cutoff=1e-3)
#' 
#' #mapping the distribution of the 3 most abundant sequences (caution, mfrow does not work for lattice's levelplot)
#' map.extrapol.freq(termes.pred, path=NULL, col.name=NULL, 1, cutoff=1e-3)
#' map.extrapol.freq(termes.pred, path=NULL, col.name=NULL, 2, cutoff=1e-3)
#' map.extrapol.freq(termes.pred, path=NULL, col.name=NULL, 3, cutoff=1e-3)
#' 
#' #dummy observationnal data
#' termes.obs = data.frame(x=c(2,3,5), y=c(2.7,2,2.6), taxa = rep("Isoptera Apicotermitinae", 3))
#' map.extrapol.freq(termes.pred, path=NULL, col.name=NULL, 3, cutoff=1e-3, add.points=termes.obs)
#' 
#' @seealso \code{\link{extrapol.freq}}, and \code{levelplot} from \code{lattice} package
#' @author Lucie Zinger
#' @export

map.extrapol.freq = function(x, path=NULL, col.name=NULL, index, cutoff=1e-3, add.points=NULL, adj=4) {
  
  require(lattice)
  
  if(!is.null(path)) {
    x.motus = apply(x$motus,2,as.character)
    name = gsub("\\.", "_", paste(gsub(", ", "_", toString(x.motus[index,col.name])), x.motus[index,"id"], sep="_"))
    file.out = paste(path, "/", name, ".png", sep="")
  }
  
  z=x$reads[,index]
  z[abs(z)>abs(log10(cutoff))]=log10(cutoff)
  z[z>0] = 0
  spj=as.data.frame(cbind(x$samples,z))
  colnames(spj)=c("x", "y", "z")
  
  map.out=levelplot(z~x+y, spj, col.regions=topo.colors(100), 
                    at=seq(log10(cutoff),log10(1), by=0.2), 
                    colorkey=list(at=seq(log10(cutoff),log10(1), by=0.2),
                                  labels=list(at=seq(log10(cutoff),log10(1), by=0.2),
                                              labels=round(10^seq(log10(cutoff),log10(1), by=0.2),3))),
                    aspect = "iso", contour=F, main=list(label=x$motus[index, "id"], cex=0.7))         
  
  if(!is.null(path)) {
    png(file=file.out, width=800, height=800)  
    print(map.out)
    if(!is.null(add.points)) {
      n = (max(spj[,"y"])-min(spj["y"]))/length(unique(spj[,"y"]))*adj
      trellis.focus("panel", 1, 1, highlight=FALSE)
      lpoints(add.points[,"x"], add.points[,"y"], cex=0.7, lwd=3, col="red")
      ltext(add.points[,"x"], add.points[,"y"]+n, add.points[,-match(c("x", "y"), colnames(add.points))], col="red", cex=1.5)
      trellis.unfocus()
    }
    dev.off() 
    
  } else {
    print(map.out)
    if(!is.null(add.points)) {
      n = (max(spj[,"y"])-min(spj["y"]))/length(unique(spj[,"y"]))*adj
      trellis.focus("panel", 1, 1, highlight=FALSE)
      lpoints(add.points[,"x"], add.points[,"y"], cex=0.7, lwd=3, col="red")
      ltext(add.points[,"x"], add.points[,"y"]+n, add.points[,-match(c("x", "y"), colnames(add.points))], col="red", cex=1)
      trellis.unfocus()
    }
  }
  
}

                                 
                                 
                                 
                                 