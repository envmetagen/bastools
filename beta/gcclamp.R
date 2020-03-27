gc.clamp<-function(ecopcroutput_clean,primer.seq,primer.direction){
  
  #split query forward primer
  if(primer.direction=="f"){
    sst = split.primer.into.cols(primer.seq = primer.seq,ecopcroutput = ecopcroutput_clean,primer.direction = "f")
  }
  
  if(primer.direction=="r"){
    sst = split.primer.into.cols(ecopcroutput_clean,primer.seq,primer.direction = "r")
  }
  
  #first, only get the last 5 bases
  a<-length(colnames(sst))
  sst2<-sst[,(a-4):a]
  #first, only get the last 5 bases
  pflist<-pflist[(a-4):a]
  
  #******************************OK to here
  
  
  
  #then those that are C or G
  sst3<-sst2[,pflist %in% c("C","G") ]
  
  
  
  #then get columns with C, G, R, Y, S, K, 
  
  #split subject primer into list
  pflist<-list()
  for(i in 1:nchar(primer.seq)){
    pflist[[i]]<-substr(primer.seq,start = i,stop = i)
  }
  
  #IUPAC RULES IN LIST in terms of C G .......have been using this, but what if primer binding sites have ambiguities. e.g. if a site is S or R?
  IUPAC=list(R=c("A","G"),
             Y=c("C","T"),
             S=c("G","C"),
             W=c("A","T"),
             K=c("G","T"),
             M=c("A","C"),
             B=c("C","G","T"),
             D=c("A","G","T"),
             H=c("A","C","T"),
             V=c("A","C","G"),
             N=c("C","G","T","A"))
  
  #replace subject primer baseswith IUPAC ambiguities
  suppressWarnings(for(i in 1:length(pflist)){
    if(pflist[[i]]=="R") pflist[[i]]<-IUPAC[["R"]]
    if(pflist[[i]]=="Y") pflist[[i]]<-IUPAC[["Y"]]
    if(pflist[[i]]=="S") pflist[[i]]<-IUPAC[["S"]]
    if(pflist[[i]]=="W") pflist[[i]]<-IUPAC[["W"]]
    if(pflist[[i]]=="K") pflist[[i]]<-IUPAC[["K"]]
    if(pflist[[i]]=="M") pflist[[i]]<-IUPAC[["M"]]
    if(pflist[[i]]=="B") pflist[[i]]<-IUPAC[["B"]]
    if(pflist[[i]]=="D") pflist[[i]]<-IUPAC[["D"]]
    if(pflist[[i]]=="H") pflist[[i]]<-IUPAC[["H"]]
    if(pflist[[i]]=="V") pflist[[i]]<-IUPAC[["V"]]
    if(pflist[[i]]=="N") pflist[[i]]<-IUPAC[["N"]]
  }  )

  
  
  
  sst2[] <- lapply(sst2, as.character)

  #ignore non GC bases
  nongc<-c("R","Y","S","W","K","M","B","D","H","V","N","A","T")
  for(i in 1:length(nongc)){
    sst2[sst2==nongc[i]]<-"N"
  }
  
  #check if query in subject primer
  sst3<-list()
  for(i in 1:length(pflist)){
  sst3[[i]]<-sst2[,i] %in% pflist[[i]]
  }

  sst4<-as.data.frame(do.call(cbind,sst3))  


