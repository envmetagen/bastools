google.read.ss<-function(sheetname){
  ss_info<-googlesheets::gs_title(sheetname)
  rowlength<-ss_info$ws$row_extent[match("ss_used_indexes_only",table = ss_info$ws$ws_title)]
  ss_data<-ss_info %>% googlesheets::gs_read(ws = "ss_used_indexes_only",range=paste0("A20:H",rowlength))
  ss_data<-as.data.frame(ss_data[!is.na(ss_data$Sample_Name),])  
}

google.read.master<-function(sheetname){
  ss_info<-googlesheets::gs_title(sheetname)
  ss_data<-ss_info %>% googlesheets::gs_read(ws = "Master_Samplesheet")
  ss_data<-as.data.frame(ss_data[!is.na(ss_data$Sample_Name),])  
}

google.read.master.url<-function(sheeturl,out=NULL){
  library(dplyr)
  url2<-stringr::str_split(sheeturl,"/d/")[[1]][2]
  url2<-stringr::str_split(url2,"/")[[1]][1]
  ss_info<-googlesheets::gs_key(url2)
  ss_data<-ss_info %>% googlesheets::gs_read(ws = "Master_Samplesheet")
  ss_data<-as.data.frame(ss_data[!is.na(ss_data$Sample_Name),])  
  if(!is.null(out)){
  write.table(ss_data,file = paste0(gsub(" ","_",ss_info$sheet_title),".txt"),quote = F,row.names = F,sep = "\t")
  message(paste("file saved as",paste0(gsub(" ","_",ss_info$sheet_title),".txt")))}
  return(ss_data)
}