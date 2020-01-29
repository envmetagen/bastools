#' @include 02_class_metabarcoding.data.R
NULL

#11.03.2011
#L.Zinger

#######################
#function anosim.pw
#######################
#computes pairwise anosim computation
#input:
#dat: dissimilarity matrix
#g: factor defining the grouping to test
#permutations: nb of permutation to access anosim statistics
#p.adjust.method: method of correction for multiple-testing
#
#output: a distance-like table containing:
#in the upper triangle: the anosims R values
#in the lower triangle: the adjusted p-values


### start

anosim.pw<-function(dat, g, permutations, p.adjust.method, ...) {
	require(vegan)
	#data.trasformation
	dat<-as.matrix(dat)
	g<-factor(g)
	
	#empty object for result storage
	ano<-matrix(NA, nrow=nlevels(g), ncol=nlevels(g), dimnames=list(levels(g),levels(g)))
	p.val.tmp<-NULL
	#running anosims
	for(i in 1:(nlevels(g)-1)) for(j in (i+1):nlevels(g)){
			tmp<-anosim(as.dist(dat[c(which(g==levels(g)[i]),which(g==levels(g)[j])),
									c(which(g==levels(g)[i]),which(g==levels(g)[j]))]),
					c(rep(levels(g)[i], length(which(g==levels(g)[i]))),
							rep(levels(g)[j], length(which(g==levels(g)[j])))), permutations)
			ano[i,j]<-tmp$statistic
			p.val.tmp<-append(p.val.tmp, tmp$signif)
		}
	
	#p value correction for multiple comparison
	p.val.tmp<-p.adjust(p.val.tmp, p.adjust.method )
	
	#put the corrected p values in the anosim table
	tmp<-NULL
	tmp2<-NULL
	for(i in 1:(nlevels(g)-1)) for(j in (i+1):nlevels(g)){
			tmp<-append(tmp,i)
			tmp2<-append(tmp2,j)
		}
	for(i in 1:length(p.val.tmp)){
		ano[tmp2[i],tmp[i]]<-p.val.tmp[i]}
	
	return(ano)
}

### end




#23 Nov 2012
#L.Zinger
###################
#function MOTUtable
###################
# Generates ready-to-use MOTU tables and basic statistics on samples (i.e. sequencing depth, raw richness, and invsimpson index)
#input:
#x: an obitable output (samples should be indicated as e.g. "sample.A01r" in column names)
#y: the column name by which that data are to be aggregated. Should be e.g. "cluster" or "species_name"
#outputs:
#x.otu: the ready-to-use MOTU table
#x.rawstats: basic statistics on samples

### start

MOTUtable<-function(x, y) {
	
	require(vegan)
	nom<-as.character(substitute(x))
	
	tmp<-x[,c(grep(y, colnames(x)), grep("sample", colnames(x)))]
	tmp2<-t(aggregate(tmp[,-1], by=list(tmp[,1]), sum))
	x.otu<-tmp2[-1,]
	colnames(x.otu)<-paste(y,tmp2[1,], sep=".")
	
	x.rawstats<-data.frame(Nb_ind=rowSums(x.otu), Raw_richness=specnumber(x.otu, MARGIN=1), Raw_eveness=diversity(x.otu, "invsimpson", MARGIN=1) )
	#may have a pb in the rowSums depending on the R version (allows or not non-numeric)
	
	assign(paste(nom, y, sep="."),x.otu,env = .GlobalEnv)
	assign(paste(nom, y, "rawstats", sep="."),x.rawstats,env = .GlobalEnv)
}

### end




#26 Nov 2012
#F.Boyer
###################
#function reads.frequency & filter.threshold
###################
#can be used to filter the table of reads to have the sequences that represents at least 95% of the total reads by sample
#
#e.g. reads.treshold(reads.frequency(metabarcodingS4Obj@reads), 0.95)


filter.threshold <- function(v, threshold) {
	o <- order(v, decreasing=T)
	ind <- which(cumsum(as.matrix(v[o]))>threshold)
	v[-o[seq(min(length(o), 1+length(o)-length(ind)))]] <- 0
	v
}

reads.threshold  <- function (reads, threshold, by.sample=T) {
	res <- apply(reads, MARGIN=ifelse(by.sample, 1, 2), filter.threshold, thr=threshold)
	if (by.sample) res <- t(res)
	data.frame(res)
}

reads.frequency <- function (reads, by.sample=T) {
	res <- apply(reads, MARGIN=ifelse(by.sample, 1, 2), function(v) {v/sum(v)})
	if (by.sample) res <- t(res)
	data.frame(res)
}


#06 Jan 2013
#F.Boyer
###################
#function removeOutliers
###################
#given a contengency table and a distance matrix
#returns the list of samples that should be removed in order to have only 
#distances below thresold
#can't return only one sample
#
#e.g. intraBad <- lapply(levels(sample.desc$sampleName), function(group) {samples<-rownames(sample.desc)[sample.desc$sampleName==group]; removeOutliers(contingencyTable[samples,], thr=0.3, distFun = function(x) vegdist(x, method='bray'))})



#require(vegan)
removeOutliers <- function(m, thr=0.3, distFun = function(x) vegdist(x, method='bray') ) {
	distMat <- as.matrix(distFun(m))
	maxM <- max(distMat)
	theBadGuys =c()
	
	while (maxM>thr) {
		bad <- apply(distMat, MARGIN=1, function(row, maxM) {any(row==maxM)}, maxM=maxM)    
		bad <- names(bad)[bad]
		bad <- apply(distMat[bad,], MARGIN=1, mean)
		badGuy <- names(bad)[bad==max(bad), drop=F][1]
		
		theBadGuys <- c(theBadGuys, badGuy)
		
		stillok <- rownames(distMat) != badGuy
		distMat <- distMat[stillok, stillok, drop=F]
		maxM <- max(distMat)
	}
	
	if (length(theBadGuys) >= (nrow(m)-1)) {
		theBadGuys <- rownames(m)    
	}
	theBadGuys
}


#31.05.2013
#L.Zinger
#getAttrPerS, a function allowing to get the values of a sequence attribute per sample
#(e.g. best_identities, etc...) the output is a list with one dataframe per sample.
#This dataframe contains:
#	 first column (named as attr): the attribute value for each sequence present in the sample
#	 second column (named weight): the corresponding number of reads in the sample

getAttrPerS=function(x,attr){
	#x: a metabarcoding object
	#attr: a character object corresponding to the attribute
	#for which values per sample are needed (should be equal to a colname in x@motus)
	
	if(class(x)[1]!= "metabarcoding.data") {
		stop("x is not a metabarcoding S4 object")
	}
	
	if(is.character(attr)==F) {
		stop("attr is not a character object")
	}
	
	x.motus = motus(x)
	x.reads = reads(x)
	
	otu = apply(x.reads, 1, function(y) x.motus[match(names(y[which(y!=0)]),x.motus$id), grep(attr, colnames(x.motus))])
	reads = apply(x.reads, 1, function(y) y[which(y!=0)])
	
	output = mapply(cbind, otu, reads)
	output = lapply(output, function(y) {
				colnames(y)=c(attr,"weight")
				return(y)
			})
	return(output)
}
### end getAttrPerS

