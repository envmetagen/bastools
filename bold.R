#' @export
#'
get_BOLD_BAS<-function(groups){
#would be good to apply taxise to ncbi download

#make folder for all
dir.create(path = "./BOLD")
setwd("./BOLD")

#check the level of each group
bb<-c("kingdom","phylum","class","order")
a2<-taxize::classification(groups,db = "col")
cc<-data.frame(ncol(4))
for(i in 1:length(a2)){
  cc[i,1:4]<-bb %in% a2[[i]]$rank[length(a2[[i]]$rank)]
}

#get families, if group is greater than family
a3<-list()
for(i in 1:length(cc$V1)){
  if("TRUE" %in% cc[i,]==TRUE){
    a3[[i]]<-taxize::downstream(groups[i], db = "col", downto = "family")
  }
}

#Loop to download sequences for each group
if(length(a3)>0){
for(i in 1:length(a3)){
  #make folder for each group
  dir.create(path = names(a3[[i]][1]))
  setwd(names(a3[[i]][1]))
    time1<-Sys.time()
    ctn<-vector()
    a4<-a3[[i]][[1]]$childtaxa_name
    ctn<-c(ctn,a4)
    ctn<-ctn[!ctn=="Not assigned"]
    ctn<-ctn[order(... = ctn,na.last = F)]
  #download sequences for each family
  for(j in 1:length(ctn)){
      out <- bold::bold_seq(ctn[j])
      time2<-Sys.time()
 ##Make fasta file for each family, if anything was found
      if(!length(out)==0){
          df1<-data.table::rbindlist(l = out)
          df1$def<-paste0(df1$id,"_",df1$name)
          df2<-df1[,c(5,4)]
          invisible(seqRFLP::dataframe2fas(x = df2,file = paste0(ctn[j],"_BOLD.fasta")))
          message(c("Downloaded ",length(df2$def), " sequences for ", paste0(ctn[j]), " in ",
              round(difftime(time2,time1),digits = 2)," mins from BOLD"))
       }
  }
    setwd("../../BOLD")
}
}


if(length(a3)==0){
  for(i in 1:length(groups)){
  #download data
  time1<-Sys.time()
  groups<-groups[order(... = groups,na.last = F)]
  out <- bold::bold_seq(groups[i])
  time2<-Sys.time()

  ##Output fasta files
  if(!length(out)==0){
  df1<-data.table::rbindlist(l = out)
  df1$def<-paste0(df1$id,"_",df1$name)
  df2<-df1[,c(5,4)]
  invisible(seqRFLP::dataframe2fas(x = df2,file = paste0(groups[i],"_BOLD.fasta")))
  message(c("Downloaded ",length(df2$def), " sequences for ", groups[i], " in ",round(difftime(time2,time1),digits = 2),
            " mins from BOLD"))
  }
  }
}
}

