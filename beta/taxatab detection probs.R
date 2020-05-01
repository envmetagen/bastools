spident=99
topS=2
gpident=97
topG=4
fpident=96
topF=5
abspident=94 
topAbs=7

stepstotake=c("step2","step3","step4","step5","step6","step7","step8","step9","step10","step11","step12","step13","step14")


taxatab<-taxatab.list[[16]]

taxatab<-binarise.taxatab(taxatab,t = T)

taxatab2<-reshape2::melt(taxatab)

#get quants for detection probability for each species
qdata.list<-list()
for(i in 1:nrow(taxatab)){
  taxatab.temp<-taxatab[taxatab$taxon==taxatab$taxon[i],]

  if(sum(taxatab.temp[,-1])>0) {
  
    qdata<-as.data.frame(t(quantile(taxatab.temp[,-1,drop=F],probs = c(0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99))))
    colnames(qdata)<-"quant"
    qdata$qcount<-0
    qdata$quant2<-as.numeric(gsub("%","",rownames(qdata)))
    #get equivalent number of samples
    taxatab.m<-reshape2::melt(taxatab.temp)
    
    for(j in 1:nrow(qdata)){
      qdata$qcount[j]<-nrow(taxatab.m[taxatab.m$value<=qdata[j,1],])
    }
    
    qdata$taxon<-taxatab$taxon[i]
    
    qdata.list[[i]]<-qdata
  } else message("Excluding ", taxatab$taxon[i], " from plot: No Detections")
}

qdata.df<-do.call(rbind,qdata.list)

ggplot(qdata.df,aes(x = qcount,y = quant2,colour=taxon)) + geom_line()

taxatab2$nsamples<-ncol(taxatab[,-1])

ggplot(taxatab2,aes(x = nsamples,y = value,colour=taxon)) + geom_point()
