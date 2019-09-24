#make count function
do.counts<-function(df, step){
  
  #Make empty df
  df_c<-data.frame(matrix(nrow=1,ncol = 7))
  colnames(df_c)<-c("step","n.reads","n.detections","n.pcr.reps","n.DNA_samples","n.biomaterial","n.taxa")
  
  df<-as.data.frame(df)
  
  #remove samples and taxa with no reads
  if(length(colnames(df))>1){
  df2<-df[rowSums(df)!=0,]
  df3<-as.data.frame(df2[,colSums(df2)!=0])
  if(length(colnames(df3))<2){
  colnames(df3)<-colnames(df2)[colSums(df2)>0]
  df2<-df3
  } }
  
  if(length(colnames(df))==1){
    df2<-as.data.frame(df[!df[1]==0,])
    if(length(colnames(df2))==1){  
    df2<-as.data.frame(df2[,!colSums(df2)==0])
    colnames(df2)<-colnames(df)
    }}
  
  if(length(colnames(df))==0){
    df2<-df}

  if(length(colnames(df2))>0){
  #count reads
  count.reads.df<-sum(df2)
  #count detections
  count.detections.df<-length(df2[df2>0])
  #count replicates
  count.pcr.replicates<-length(colnames(df2))
  #count DNA samples
  count.DNA_samples<-length(unique(gsub("(.*)-(.*)-(.*)-(.*)-(.*)$", "\\1-\\2-\\3", colnames(df2))))
  #count biomaterial
  count.bio_samples<-length(unique(gsub("(.*)-(.*)-(.*)-(.*)-(.*)$", "\\1-\\2", colnames(df2))))
  #count.taxa
  count.taxa<-length(rownames(df2))
  #combine
  df_c[1,]<-c(step,count.reads.df,count.detections.df,count.pcr.replicates,count.DNA_samples,
              count.bio_samples,count.taxa)
  }
  
  return(df_c)
}
