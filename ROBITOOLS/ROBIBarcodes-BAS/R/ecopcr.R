#'@include ROBIBarcodes.R
NULL

#   column 1 : accession number
#   column 2 : sequence length
#   column 3 : taxonomic id
#   column 4 : rank
#   column 5 : species taxonomic id
#   column 6 : scientific name
#   column 7 : genus taxonomic id
#   column 8 : genus name
#   column 9 : family taxonomic id
#   column 10 : family name
#   column 11 : super kingdom taxonomic id
#   column 12 : super kingdom name
#   column 13 : strand (direct or reverse)
#   column 14 : first oligonucleotide
#   column 15 : number of errors for the first strand
#   column 16 : Tm for hybridization of primer 1 at this site
#   column 17 : second oligonucleotide
#   column 18 : number of errors for the second strand
#   column 19 : Tm for hybridization of primer 1 at this site
#   column 20 : amplification length
#   column 21 : sequence
#   column 22 : definition  

#' Read the result file produced by the ecoPCR program.
#' 
#' @export
read.ecopcr.result = function(file) 
{
  
  split.line = function(line) {
    l = strsplit(line,split=" +\\| +")[[1]]
    l = c(l[1:21],paste(l[-c(1:21)],sep="|"))
    return(l)
  }
  
  if (missing(file) && !missing(text)) {
    file <- textConnection(text)
    on.exit(close(file))
  }
  if (is.character(file)) {
    file <- file(file, "rt")
    on.exit(close(file))
  }
  if (!inherits(file, "connection")) 
    stop("'file' must be a character string or connection")
  if (!isOpen(file, "rt")) {
    open(file, "rt")
    on.exit(close(file))
  }
  
  line = readLines(file,1)
  while (length(grep('^#',line))==1) {
    line = readLines(file,1)
  }
  
  pushBack(line,file)
  
  lines = lapply(readLines(file),split.line)
  nlines = length(lines)
  AC                = sapply(1:nlines,function(x) lines[[x]][1])
  seq_length        = as.integer(sapply(1:nlines,function(x) lines[[x]][2]))
  taxid             = as.integer(sapply(1:nlines,function(x) lines[[x]][3]))
  rank              = as.factor(sapply(1:nlines,function(x) lines[[x]][4]))
  species           = type.convert(sapply(1:nlines,function(x) lines[[x]][5]),na.string="###")
  species_name      = sapply(1:nlines,function(x) lines[[x]][6])
  genus             = type.convert(sapply(1:nlines,function(x) lines[[x]][7]),na.string="###")
  genus_name        = sapply(1:nlines,function(x) lines[[x]][8])
  family            = type.convert(sapply(1:nlines,function(x) lines[[x]][9]),na.string="###")
  family_name       = sapply(1:nlines,function(x) lines[[x]][10])
  superkingdom      = type.convert(sapply(1:nlines,function(x) lines[[x]][11]),na.string="###")
  superkingdom_name = sapply(1:nlines,function(x) lines[[x]][12])
  strand            = as.factor(sapply(1:nlines,function(x) lines[[x]][13]))
  forward_match     = sapply(1:nlines,function(x) lines[[x]][14])
  forward_mismatch  = as.integer(sapply(1:nlines,function(x) lines[[x]][15]))
  forward_tm        = as.double(sapply(1:nlines,function(x) lines[[x]][16]))
  reverse_match     = sapply(1:nlines,function(x) lines[[x]][17])
  reverse_mismatch  = as.integer(sapply(1:nlines,function(x) lines[[x]][18]))
  reverse_tm        = as.double(sapply(1:nlines,function(x) lines[[x]][19]))
  amplicon_length   = as.integer(sapply(1:nlines,function(x) lines[[x]][20]))
  sequence          = sapply(1:nlines,function(x) lines[[x]][21])
  definition        = sapply(1:nlines,function(x) lines[[x]][22])
  
  eco = data.frame(AC,seq_length,taxid,rank,
                   species,species_name,
                   genus,genus_name,
                   family,family_name,
                   superkingdom,superkingdom_name,
                   strand,
                   forward_match,forward_mismatch,forward_tm,
                   reverse_match,reverse_mismatch,reverse_tm,
                   amplicon_length,sequence,definition
  )
  
  return(eco)
}

ecopcr.frequencies = function(matches,group=NULL) {
  compute = function(matches) {
    w = as.matrix(do.call(rbind,strsplit(as.character(matches),'')))
    d = dim(w)
    w=factor(w,levels=c('A','C','G','T'))
    dim(w)=d
    w=t(w)
    freq = mapply(function(x) table(w[x,]),1:d[2])
    freq = freq[c('A','C','G','T'),]
    csum = colSums(freq)
    freq = sweep(freq,2,csum,'/')
    attr(freq,'count')=length(w)
    return(freq)
  }
  if (is.null(group))
    return(compute(matches))
  else {
    lmatches = aggregate(matches,by=list(group=as.factor(group)),as.character)
    w = lmatches$x
    names(w)=lmatches$group
    lf = lapply(w,compute)
    return(lf)
  }
}

#' @export
ecopcr.forward.frequencies = function(ecopcr,group=NULL) {
  return(ecopcr.frequencies(ecopcr$forward_match,group))
}

#' @export
ecopcr.reverse.frequencies = function(ecopcr,group=NULL) {
  return(ecopcr.frequencies(ecopcr$reverse_match,group))
}

#' @export
dna.shanon = function(freq,base=2) {
  shanon = log(4)/log(base) - colSums(-freq *log(freq) / log(base),na.rm=TRUE)
  return(sweep(freq,2,shanon,'*')) 
}


ecopcr.shanon = function(matches,group=NULL,base=2) {
  
  if (is.null(group)) {
    freq = ecopcr.frequencies(matches)
    return(dna.shanon(freq))
  }
  else {
    lf = lapply(ecopcr.frequencies(matches,group),dna.shanon)
    return(lf)
  }
}

#' @export
ecopcr.forward.shanon = function(ecopcr,group=NULL,base=2) {
  return(ecopcr.shanon(ecopcr$forward_match,group,base))
}

#' @export
ecopcr.reverse.shanon = function(ecopcr,group=NULL,base=2) {
  return(ecopcr.shanon(ecopcr$reverse_match,group,base))
}
