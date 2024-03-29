#Run all code below this line
########################################CODE
setwd(outdir)

if(ms_option==1) packages.needed<-c("googlesheets4","httpuv","data.table")

if(ms_option==2) packages.needed<-c("data.table")

for(i in 1:length(packages.needed)){
  if(!packages.needed[i] %in% rownames(installed.packages())) {
    message(paste(packages.needed[i],"not found. Installing..."))
    install.packages(packages.needed[i])
  }
}

#source(paste0(bastoolsDir,"master_functions.R"))
#bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/"
#setwd(bastoolsDir)
#googlesheets4::gs4_auth(email = email)


######################################load functions
google.read.master.url<-function(sheeturl,out=NULL,ws="Master_Samplesheet"){
  url2<-stringr::str_split(sheeturl,"/d/")[[1]][2]
  url2<-stringr::str_split(url2,"/")[[1]][1]
  ss_info<-googlesheets4::gs4_get(ss = url2)
  if(ws == "ENA_sample_data") {
    ss_data<-googlesheets4::read_sheet(ss = url2,sheet = ws,col_types = "c") 
    colnames(ss_data)<-ss_data[2,]
    ss_data<-ss_data[3:length(ss_data$sample_alias),]
    ss_data<-as.data.frame(ss_data[!is.na(ss_data$sample_alias),])
  } else{
    if(ws == "ENA_library_data" | ws == "ENA_library_data_Single_End") { 
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

master_xtabs<-function(master_sheet,columns){
  
  coln<-as.numeric(0)
  for(i in 1:length(columns)){
    coln[i]<-grep(paste0(columns[i],"$"), colnames(master_sheet))
  }
  
  xtabs(data = master_sheet[,coln],addNA = T)
}  


########################################SAMPLE DATA
#first read ENA_sample_data sheet
if(ms_option==1) sample_sheet<-google.read.master.url(sheeturl,ws = "ENA_sample_data")
if(ms_option==2) sample_sheet<-data.table::fread(sample_sheet_file,skip = 2,data.table = F)  

#check on number of samples in each category
if(!is.null(sample.subsetlist)){
  message("Showing xtabs for entire datasheet, including Sample_Type by default")
  master_xtabs(master_sheet = sample_sheet,columns=c("sample_type",names(sample.subsetlist)))
  
  #subset sample sheet for this study
  message("Subsetting datasheet")
  sub_sample_sheet<-subset_mastersheet(sample_sheet, sample.subsetlist)
  
  #check again to see the subset made sense
  print(master_xtabs(sub_sample_sheet,columns=c("sample_type",names(sample.subsetlist))))
} else sub_sample_sheet<-sample_sheet

#write 
write.table(data.frame(V1=c("#checklist_accession","#unique_name_prefix"),V2=c("ERC000011","")),
            file = sample_sheet_name,quote = F,row.names = F,sep = "\t",col.names = F)
write.table(sub_sample_sheet,file = sample_sheet_name,quote = F,row.names = F,sep = "\t",append = T)
message("warning about appending column names to file is ok")

########################################LIBRARY DATA
#first read mastersheet to get libraries wanted
if(ms_option==1) master_sheet<-google.read.master.url(sheeturl)
if(ms_option==2) master_sheet<-data.table::fread(master_sheet_file,data.table = F)  

if(!"library_name" %in% colnames(master_sheet)) stop("Must have column library_name in Master_Samplesheet")

message("Count libraries before subsetting=",length(master_sheet$library_name))
#subset 
if(!is.null(library.subsetlist)){
  sub_master_sheet<-subset_mastersheet(master_sheet, library.subsetlist)
} else sub_master_sheet<-master_sheet

message("Count libraries after subsetting master sheet=",length(sub_master_sheet$library_name))

#read library sheet
if(ms_option==1) library_sheet<-google.read.master.url(sheeturl,ws = "ENA_library_data")
if(ms_option==2) library_sheet<-data.table::fread(library_sheet_file,data.table = F)  

#subset by sample sheet
sub_library_sheet<-library_sheet[library_sheet$library_name %in% sub_master_sheet$library_name,]
message("Count libraries after subsetting from library sheet=",length(sub_library_sheet$library_name))

########################################FILE LIST 
#get filenames and check they exist
setwd(fileDir)
if(minion) files<-sub_library_sheet$forward_file_name else files<-c(sub_library_sheet$forward_file_name,sub_library_sheet$reverse_file_name)
message(paste(sum(files %in% list.files()), "of", length(files), "files found, of which",
              length(files[duplicated(files)]), "are duplicates", files[duplicated(files)]))

write.table(files,file = filelist,quote = F,row.names = F,sep = "\t",col.names = F)

########################################MD5SUMS
checksum<-data.frame(md5=rep(0,length(files)))

for(i in 1:length(files)){
  checksum[i,1]<-system2("md5sum", args=files[i],stdout = T)
}
checksum2<-as.data.frame(do.call(rbind,stringr::str_split(checksum$md5," ")))
checksum2$V2<-NULL

#add md5 to library data
if(minion) {sub_library_sheet$forward_file_md5<-checksum2$V1[1:length(sub_library_sheet$forward_file_name)]
} else {
  sub_library_sheet$forward_file_md5<-checksum2$V1[1:length(sub_library_sheet$forward_file_name)]
  sub_library_sheet$reverse_file_md5<-checksum2$V1[length(sub_library_sheet$reverse_file_name)+1:
                                                   length(sub_library_sheet$reverse_file_name)]
}

#change NAs to empty
sub_library_sheet[is.na(sub_library_sheet)]<-""

#remove rev filename if minion
if(minion) {
  colnames(sub_library_sheet)<-gsub("forward_file_name","file_name",colnames(sub_library_sheet))
  colnames(sub_library_sheet)<-gsub("forward_file_md5","file_md5",colnames(sub_library_sheet))
  colnames(sub_library_sheet)<-gsub("forward_file_unencrypted_md5","file_unencrypted_md5",colnames(sub_library_sheet))
  sub_library_sheet<-sub_library_sheet[,-grep("reverse_.*",colnames(sub_library_sheet))]
}
  
 # sample_alias	instrument_model	library_name	library_source	library_selection	library_strategy	design_description	library_construction_protocol	file_name	file_md5


#write library sheet
write.table(x = sub_library_sheet,file = paste0(outdir,library_sheet_name),quote = F,row.names = F,sep = "\t",col.names = T)

#write password file
write.table(passwd,file = "ENA_password.txt",quote = F,row.names = F,sep = "\t",col.names = F)

message("warning about appending column names to file is ok")

message("Next do:

source /opt/env.sh

Then:

cd ",fileDir,"

Then: 

nohup emg_ena_upload.sh -l ", filelist," -u Webin-XXXX -p ENA_password.txt 2>upload.log 1>&2 & 

change name of log file if you wish
replace Webin-XXXX with your username

to complete metadata entry, use the output files that were created:","
",
paste0(outdir,sample_sheet_name),"
",
paste0(outdir,library_sheet_name)," 

and enter them here https://wwwdev.ebi.ac.uk/ena/submit/sra/#home for testing
or here for real https://www.ebi.ac.uk/ena/submit/sra/

if you havent done so already you will need to create a project first
read my notes for more
https://docs.google.com/document/d/1DOHRvA9cjwfrsLDkeZ4SRfIP7GIo9d0dxl5vjTBqfFY/edit#")
