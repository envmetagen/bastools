

mismatch.table<-function(ecopcroutput_clean,primer.seq=NULL,primer.direction){

  #split query forward primer
  if(primer.direction=="f"){
    sst <- as.data.frame(t(as.data.frame(strsplit(as.character(ecopcroutput_clean$forward_match), ""))))}
  if(primer.direction=="r"){
    sst <- as.data.frame(t(as.data.frame(strsplit(as.character(ecopcroutput_clean$reverse_match), ""))))}
  sst$family_name<-ecopcroutput_clean$family_name
  rownames(sst)<-NULL

  #split subject primer into list
  pflist<-list()
  for(i in 1:nchar(primer.seq)){
    pflist[[i]]<-substr(primer.seq,start = i,stop = i)
  }

  #IUPAC RULES IN LIST
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
  })

  #check if query primer is in subject primer
  #function to apply %in% to one data.frame
  matchdf<-function(df,pflist){
    out<-data.frame(matrix(ncol=length(pflist),nrow = length(df[,1])))
    for(i in 1:length(pflist)){
      out[,i]<-df[,i] %in% pflist[[i]]
    }
    return(out)
  }

  #apply
  out2<-matchdf(sst,pflist)

  #convert true/false to 1/0
    out3<-out2*1
    #swap 0s for 1s

      out3[out3==1]<-2
      out3[out3==0]<-1
      out3[out3==2]<-0

      return(out3)
}
