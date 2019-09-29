google.read.ss<-function(sheetname){
  ss_info<-googlesheets::gs_title(sheetname)
  rowlength<-ss_info$ws$row_extent[match("ss_used_indexes_only",table = ss_info$ws$ws_title)]
  ss_data<-ss_info %>% googlesheets::gs_read(ws = "ss_used_indexes_only",range=paste0("A20:H",rowlength))
  ss_data<-ss_data[!is.na(ss_data$Sample_Name),]  
}

google.read.master<-function(sheetname){
  ss_info<-googlesheets::gs_title(sheetname)
  ss_data<-ss_info %>% googlesheets::gs_read(ws = "Master_Samplesheet")
  ss_data<-ss_data[!is.na(ss_data$Sample_Name),]  
}