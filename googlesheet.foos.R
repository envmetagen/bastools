google.read.ss<-function(sheetname){
  ss_info<-googlesheets4::gs_title(sheetname)
  rowlength<-ss_info$ws$row_extent[match("ss_used_indexes_only",table = ss_info$ws$ws_title)]
  ss_data<-ss_info %>% googlesheets::gs_read(ws = "ss_used_indexes_only",range=paste0("A20:H",rowlength))
  ss_data<-as.data.frame(ss_data[!is.na(ss_data$Sample_Name),])  
}

google.read.master<-function(sheetname){
  ss_info<-googlesheets::gs_title(sheetname)
  ss_data<-ss_info %>% googlesheets::gs_read(ws = "Master_Samplesheet")
  ss_data<-as.data.frame(ss_data[!is.na(ss_data$Sample_Name),])  
}

google.read.master.url<-function(sheeturl,out=NULL,ws="Master_Samplesheet"){
  url2<-stringr::str_split(sheeturl,"/d/")[[1]][2]
  url2<-stringr::str_split(url2,"/")[[1]][1]
  ss_info<-googlesheets4::sheets_get(ss = url2)
  if(ws == "ENA_sample_data") {
    ss_data<-googlesheets4::read_sheet(ss = url2,sheet = ws,col_types = "c") 
    colnames(ss_data)<-ss_data[2,]
    ss_data<-ss_data[3:length(ss_data$sample_alias),]
    ss_data<-as.data.frame(ss_data[!is.na(ss_data$sample_alias),])
    } else{
      if(ws == "ENA_library_data") { 
      ss_data<-googlesheets4::read_sheet(ss = url2,sheet = ws,col_types = "c") 
      ss_data<-as.data.frame(ss_data[!is.na(ss_data$sample_alias),])
      } else{
        if(ws == "Library_data") { 
          ss_data<-googlesheets4::read_sheet(ss = url2,sheet = ws,col_types = "c") 
          ss_data<-as.data.frame(ss_data[!is.na(ss_data$run_alias),])
        } else{
          ss_data<-googlesheets4::read_sheet(ss = url2,sheet = ws,col_types = "c") 
          ss_data<-as.data.frame(ss_data[!is.na(ss_data$Sample_Name),])
          }
        }
    }


  if(!is.null(out)){
  write.table(ss_data,file = paste0(gsub(" ","_",ss_info$name),"_",gsub(" ","_",ws),".txt"),quote = F,row.names = F,sep = "\t")
  message(paste("file saved as",paste0(gsub(" ","_",ss_info$name),"_",gsub(" ","_",ws),".txt")))
  }
  return(ss_data)
}

#read mastersheets and make a processing sheet 
google.make.experiment.sheet<-function(outDir,sheeturls,experiment_id){
  master<-list()
  headers<-c("barcode_id","Primer_set","Primer_F","Primer_R","Min_length","Max_length","ss_sample_id","experiment_id")
  for(i in 1:length(sheeturls)){
    master[[i]]<-google.read.master.url(sheeturls[i])
    if(length(headers)!=sum(headers %in% colnames(master[[i]]))){
      stop (c("one of the following headers missing: ", paste(headers,collapse = " ")))}
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

#read (and write) a processing sheet 
google.make.exp.sheet.illumina<-function(outDir,sheeturls,experiment_id){
  master<-list()
  headers<-c("Primer_set","Primer_F","Primer_R","Min_length","Max_length","ss_sample_id","experiment_id")
  for(i in 1:length(sheeturls)){
    master[[i]]<-google.read.master.url(sheeturls[i])
    if(length(headers)!=sum(headers %in% colnames(master[[i]]))){
      stop (c("one of the following headers missing: ", paste(headers,collapse = " ")))}
    master[[i]]<-master[[i]][,headers]
  }
  
  #make a processing sheet
  experimentsheet<-as.data.frame(data.table::rbindlist(master))
  experimentsheet<-experimentsheet[experimentsheet$experiment_id %in% experiment_id,]
  
  #write file
  write.table(experimentsheet,paste0(outDir,paste(experiment_id,collapse = "_"),"_experiment_sheet.txt"),sep = "\t",quote = F,row.names = F)
  message(paste("file saved as",paste0(outDir,paste(experiment_id,collapse = "_"),"_experiment_sheet.txt")))
  
  return(experimentsheet)
}

#workaround for getting filenames of fastas to blast
google.get.startingfastas<-function(outDir,sheeturls,experiment_id,usingobiuniq){
experimentsheet<-google.make.experiment.sheet(outDir,sheeturls,experiment_id) #in process, writes a sheet to file 

#get barcodes used  
barcodes.used<-unique(experimentsheet$barcode_id)
barcodes.used <- barcodes.used[!is.na(barcodes.used)]
barcodes.used<-gsub("BC","barcode",barcodes.used)

#size select, for each fragment, I checked and seqs appear to have primers plus one base (at each end)
experimentsheet$primer_combo<-paste0(experimentsheet$Primer_F,experimentsheet$Primer_R)
primer_combo<-unique(experimentsheet$primer_combo)

#barcodes in each primer combo
primer_combo.bcs<-list()
for(i in 1:length(primer_combo)){
  primer_combo.bcs[[i]]<-unique(experimentsheet[experimentsheet$primer_combo==primer_combo[i],"barcode_id"])
  primer_combo.bcs[[i]]<-gsub("BC","barcode",primer_combo.bcs[[i]])
  names(primer_combo.bcs[[i]])<-gsub(" ","",
                                     experimentsheet[experimentsheet$primer_combo==primer_combo[i],"Primer_set"][1])
}
  if(usingobiuniq){
    for(i in 1:length(primer_combo.bcs)){
      print(paste0(experiment_id,"_",names(primer_combo.bcs[[i]][1]),".uniq.filtlen.wlen.obi.fasta"))
    }}

  if(usingobiuniq==F){
    for(i in 1:length(primer_combo.bcs)){
      print(paste0(experiment_id,"_",names(primer_combo.bcs[[i]][1]),".filtlen.wlen.obi.fasta"))
  }}
  
}

#subset a mastersheet based on required values
# e.g
# ms_ss<-subset_mastersheet(master_sheet,
#                           list(experiment_id=c("2018_02")
#                           ,Primer_set=c("12SV5.2","12SV51"),
#                           Sample_Type=c("Field")))
subset_mastersheet<-function(master_sheet,...){
  
  a<-as.list(...)
  
  print(a)
  
  m2<-master_sheet
  
  for(i in 1:length(a)){
    coln<-grep(paste0(names(a)[i],"$"), colnames(m2))
    m2<-m2[m2[,coln] %in% a[[i]],]
  }
  return(m2)
}

#do xtabs to 1) find levels in factors and 2) see counts
master_xtabs<-function(master_sheet,columns){
  
  coln<-as.numeric(0)
  for(i in 1:length(columns)){
    coln[i]<-grep(paste0(columns[i],"$"), colnames(master_sheet))
  }
  
  xtabs(data = master_sheet[,coln],addNA = T)
}  




  

