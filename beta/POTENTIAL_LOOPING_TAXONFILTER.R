filterlist<-list()
a<-c(0,0.1,0.25,0.5,1,2,3,5,10)
for(i in 1:length(a)){
  filterlist[[i]]<-taxon.filter.solo.df(taxatab,taxonpc=a[i])
}

negativesList<-list()
for(i in 1:length(filterlist)){
  negativesList[[i]]<-suppressMessages(negs.stats(taxatab=filterlist[[i]],ms_ss = ms_ss,real=real,ex_hominidae=F,printnegs = F))
}

dxns.in.negs.List<-list()
for(i in 1:length(negativesList)){
dxns.in.negs<-count.dxns.by.taxon(negativesList[[i]]$Extraction_Negative)
dxns.in.negs.List[[i]]<-sum(dxns.in.negs$n.samples)
}

do.call(c,dxns.in.negs.List)

out<-data.frame(taxonpc=a,dxns=do.call(c,dxns.in.negs.List))


require(tidyverse)

long<-reshape2::melt(taxatab)
long<-long[long$value>0,]

long$sample_type<-ms_ss[match(long$variable,ms_ss$ss_sample_id),"sample_type"]
long$taxon<-as.character(long$taxon)

splitdf<-split(long, f = long$taxon)

#only include dfs with at least one negs and real

##neg sample types
negtypes<-unique(ms_ss$sample_type)[!unique(ms_ss$sample_type) %in% real]

splitdf2<-list()

##dfs with negs
for(i in 1:length(splitdf)){
  if(!TRUE %in% (negtypes %in% splitdf[[i]]$sample_type)) {
    splitdf2[[i]]<-NULL
  } else {
    splitdf2[[i]]<-splitdf[[i]]
  }
}

splitdf2<-splitdf2 %>% discard(is.null)

splitdf3<-list()

##dfs with real
for(i in 1:length(splitdf2)){
  if(!TRUE %in% (real %in% splitdf2[[i]]$sample_type)) {
    splitdf3[[i]]<-NULL
  } else {
    splitdf3[[i]]<-splitdf2[[i]]
  }
}

splitdf3<-splitdf3 %>% discard(is.null)

taxon.filter.solo.df2<-function(taxatab){
  taxon.filter.solo.df(taxatab,taxonpc=2)
}

splitdf4<-lapply(splitdf3,FUN = taxon.filter.solo.df2)

