rm.detections<-function(df,absolute_in_rep,sumreps){
  message("If a detection was identified in more than one PCR replicate, all detections are kept.
          
          If a detection was identified in only one PCR replicate, but it had a read count above 
          absolute_in_rep, it is kept. If a detection was identified in only one PCR replicate, 
          and had a read count below absolute_in_rep, it is put to zero. 
          
          If only one PCR replicate is available,read counts below absolute_in_rep are put to zero.
          
          Setting absolute_in_rep=1000000000, effectively turns off one-replicate option, in which case 
          detections are only kept if they appear in both samples.")
  DNA_samples<-unique(gsub("(.*)-(.*)-(.*)-(.*)-(.*)$", "\\1-\\2-\\3", colnames(df)))
  
  clean<-list()
  for(i in 1:length(DNA_samples)){
    df2<-as.data.frame(df[,grep(DNA_samples[i],colnames(df))])
    
    # test<-data.frame(matrix(ncol = 4,nrow = 10),x=c(4000,100,4000,100,4000,100,0,0,0,0),
    #                  y=c(4000,100,100,4000,0,0,100,4000,0,0),
    #                  z=c(20,10000,20,20,4000,100,0,0,0,0),
    #                  a=c(20,10000,20,20,4000,100,1000,4000,2,1))
    # test<-test[,5:8]
    # rep1<-test[,1]>absolute_in_rep | test[,1]<absolute_in_rep & test[,2]!=0 
    # test[!rep1,1]<-0
    # rep2<-test[,2]>absolute_in_rep | test[,2]<absolute_in_rep & test[,1]!=0 
    # #more than 2 reps
    # rep1<-test[,1]>absolute_in_rep | test[,1]<absolute_in_rep & rowSums(test[,-1]>0)>0 
    # rep2<-test[,2]>absolute_in_rep | test[,2]<absolute_in_rep & rowSums(test[,-2]>0)>0 
    # rep3<-test[,3]>absolute_in_rep | test[,3]<absolute_in_rep & rowSums(test[,-3]>0)>0 
    # rep4<-test[,4]>absolute_in_rep | test[,4]<absolute_in_rep & rowSums(test[,-4]>0)>0 
    
    if(length(colnames(df2))>4) stop("Can only handle 4 PCR replicates")
    
    if(length(colnames(df2))>1){
      #for each rep: keep detection if rep has a detection above threshold, 
      #or if the detection is below the threshold but one of the other reps also has reads
      rep1<-df2[,1]>absolute_in_rep | df2[,1]<absolute_in_rep & df2[,2]!=0 
      rep2<-df2[,2]>absolute_in_rep | df2[,2]<absolute_in_rep & df2[,1]!=0 
      if(length(colnames(df2))>2){
      rep3<-df2[,3]>absolute_in_rep | df2[,3]<absolute_in_rep & rowSums(df2[,-3]>0)>0} 
      if(length(colnames(df2))>3){
      rep4<-df2[,4]>absolute_in_rep | df2[,4]<absolute_in_rep & rowSums(df2[,-4]>0)>0} 
      #put rep 1 to zero if below threshold
      df2[!rep1,1]<-0
      #put rep 2 to zero if below threshold
      df2[!rep2,2]<-0
      #put rep 3 to zero if below threshold
      if(length(colnames(df2))>2){
      df2[!rep3,3]<-0}
      #put rep 4 to zero if below threshold
      if(length(colnames(df2))>3){
      df2[!rep2,4]<-0}
      
      clean[[i]]<-df2
    } else {
      #rep one has detection above threshold
      rep1<-df2>absolute_in_rep 
      #put rep 1 to zero if below threshold
      rownames(df2)<-rownames(df)
      colnames(df2)<-grep(DNA_samples[i],colnames(df),value = T)
      df2[!rep1,1]<-0
      clean[[i]]<-df2
    }
  }
  
  #optional sum replicates
  if(sumreps==T){
    for(i in 1:length(clean)){
      if(length(colnames(clean[[i]]))>1){
        clean[[i]]$sum<-rowSums(clean[[i]])
        colnames(clean[[i]])<-gsub("sum",DNA_samples[i],colnames(clean[[i]]))
        clean[[i]]<-as.data.frame(clean[[i]][,-grep("^(.*)-(.*)-(.*)-(.*)-(.*)$",colnames(clean[[i]]))])
        colnames(clean[[i]])<-DNA_samples[i]
      }
      if(length(colnames(clean[[i]]))<2){
        colnames(clean[[i]])<-gsub("(.*)-(.*)-(.*)-(.*)-(.*)$", "\\1-\\2-\\3", colnames(clean[[i]]))
      }
    }
  }
  
  out<-do.call(cbind,clean)
}
