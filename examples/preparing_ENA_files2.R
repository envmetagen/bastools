#making sample sheet, library sheet (adding md5sums), file list and password file


#first make a copy this file, do not edit this directly

#CONFIG
ms_option=1 #choose 1 or 2, see below

#full path. Must already exist. If using Option 2 below, this folder must contain the 3 tsv files required.
#Final ENA upload files will also be added here
outdir<-"/home/tutorial/test/"
#directory with files for upload
fileDir<-"/mnt/Disk2/MISEQ_RUNS/2017_12_AZORES-IRANVERTS-NZFROG-OZ-GUELTA/PP121217-56327499/BEFORE_PAIRED_END/FIRST_ATTEMPT/"


####Option1
#use google. This is ok, but can be problematic because 
#it will try to open browser if you havent stored credentials before and this is too slow if accessing emg machines . 
#If using it set the following options, otherwise leave as they are
#datasheet
sheeturl<-"https://docs.google.com/spreadsheets/d/1ru01umKc4tCvNxUlxUJPcUyz9jon1L-QJUm67etL3aE/edit?usp=sharing"
#email  
email="basegeter@gmail.com"

####Option 2
#download the tsv files from the google doc and put them in 'outdir'. 
#You will need Master_Samplesheet, ENA_sample_data, ENA_library_data
master_sheet_file<-"EMG_Project_EtOH_biases_datasheet_v2 - Master_Samplesheet.tsv" 
sample_sheet_file<-"EMG_Project_EtOH_biases_datasheet_v2 - ENA_sample_data.tsv"
library_sheet_file<-"EMG_Project_EtOH_biases_datasheet_v2 - ENA_library_data.tsv"

#options for subsetting sample sheet. Used if you only want to include a subset of the samples
#a named list e.g. where Site="Atar" #sample.subsetlist<-list(Site=c("Atar"))
#if none use sample.subsetlist<-NULL
sample.subsetlist<-NULL

#options for subsetting library sheet. Used if you only want to include a subset of the libraries
#Note that the column being specified is one on the master_sheet. e.g. library.subsetlist<-list(study="NZ1")
#where the master_sheet has a column called study, and we only want to include libraries where study=NZ1
#if none use library.subsetlist<-NULL
library.subsetlist<-NULL



#data.table, headers thing

#name to give to final sample sheet filename
sample_sheet_name<-"filipa_test_ENA_SAMPLE_SHEET.txt"
#name to give to final library sheet filename
library_sheet_name<-"filipa_test_ENA_LIBRARY_SHEET.txt"
#name to give to list of files for upload
filelist<-"filipa_test_ENA_FILES.txt"
#ENA password
passwd<-"xxx4452"

#Run all code below this line
########################################CODE
setwd(outdir)

if(ms_option=1) packages.needed<-c("googlesheets4","httpuv","data.table")

if(ms_option=2) packages.needed<-c("data.table")

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
  master_xtabs(sub_sample_sheet,columns=c("sample_type",names(sample.subsetlist)))
} else sub_sample_sheet<-sample_sheet

#write 
write.table(data.frame(V1=c("#checklist_accession","#unique_name_prefix"),V2=c("ERC000011","")),
            file = sample_sheet_name,quote = F,row.names = F,sep = "\t",col.names = F)
write.table(sub_sample_sheet,file = sample_sheet_name,quote = F,row.names = F,sep = "\t",append = T)
message("warning appending column names to file is ok")

########################################LIBRARY DATA
#first read mastersheet to get libraries wanted
if(ms_option==1) master_sheet<-google.overlord(url = sheeturl)
if(ms_option==2) master_sheet<-data.table::fread(master_sheet_file,data.table = F)  

if(!"library_name" %in% colnames(master_sheet)) stop("Must have column library_name in Master_Samplesheet")

message("Count libraries before subsetting=",length(master_sheet$library_name))
#subset 
if(!is.null(sample.subsetlist)){
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
files<-c(sub_library_sheet$forward_file_name,sub_library_sheet$reverse_file_name)
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
sub_library_sheet$forward_file_md5<-checksum2$V1[1:length(sub_library_sheet$forward_file_name)]
sub_library_sheet$reverse_file_md5<-checksum2$V1[length(sub_library_sheet$reverse_file_name)+1:
                                                length(sub_library_sheet$reverse_file_name)]

#change NAs to empty
sub_library_sheet[is.na(sub_library_sheet)]<-""

#write library sheet
write.table(x = sub_library_sheet,file = library_sheet_name,quote = F,row.names = F,sep = "\t",col.names = T)

#write password file
write.table(passwd,file = "ENA_password.txt",quote = F,row.names = F,sep = "\t",col.names = F)

#run from terminal
# emg_ena_upload.sh -l filelist -u Webin-51203 -p ENA_password.txt
        
#go to https://wwwdev.ebi.ac.uk/ena/submit/sra/#home #FOR TESTING
#or here for real https://www.ebi.ac.uk/ena/submit/sra/

####need to see if can upload md5s directly to google sheet
