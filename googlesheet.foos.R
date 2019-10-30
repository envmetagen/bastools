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

#read (and write) a processing sheet 
google.make.experiment.sheet<-function(outDir,sheeturls,experiment_id){
master<-list()
headers<-c("barcode_id","Primer_set","Primer_F","Primer_R","Min_length","Max_length","ss_sample_id","experiment_id")
for(i in 1:length(sheeturls)){
  master[[i]]<-google.read.master.url(sheeturls[i])
  if(length(headers)!=sum(headers %in% colnames(master[[i]]))){
    stop (c("one of the following headers missing: ", paste(headers)))}
  master[[i]]<-master[[i]][,headers]
}

#make a processing sheet
experimentsheet<-as.data.frame(data.table::rbindlist(master))
experimentsheet<-experimentsheet[experimentsheet$experiment_id==experiment_id,]

#write file
write.table(experimentsheet,paste0(outDir,experiment_id,"_experiment_sheet.txt"),sep = "\t",quote = F,row.names = F)
message(paste("file saved as",paste0(outDir,experiment_id,"_experiment_sheet.txt")))

return(experimentsheet)
}