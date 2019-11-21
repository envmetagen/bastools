#uploading ENA data via Aspera

#CONFIG
bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/"

#filturb datasheet
sheeturl<-"https://docs.google.com/spreadsheets/d/1FUSaeVaYzms2EOGUoCAB4jaRKzguD3AKTsC8lYwaKP4/edit?ts=5dae01be#gid=1531090624"

#options for subsetting sample sheet 
sample.subsetlist<-list(study=c("North"))

#options for subsetting library sheet 
library.subsetlist<-list(experiment_id=c("2018_07"), Primer_set=c("12SV51"), study=c("North"))

#email
email="basegeter@gmail.com"

#directory with files for upload
fileDir<-"/mnt/Disk2/MISEQ_RUNS/2018_07_FILTURB_IRANVERTS/FASTQ_GENERATION_3_18/"
##note that the ENA upload sheets will be saved here

#final sample sheet filename
sample_sheet<-"FILTURB_NORTH_ENA_SAMPLE_SHEET.txt"

#name to give list of files
filelist<-"FILTURB_NORTH_ENA_FILES.txt"

#CODE

setwd(bastoolsDir)
googlesheets4::sheets_auth(email = email)
source("googlesheet.foos.R")
setwd(fileDir)

########################################SAMPLE DATA
#first read ENA_sample_data sheet
sample_sheet<-google.read.master.url(sheeturl,ws = "ENA_sample_data")

#check on number of samples in each category
message("Showing xtabs for entire datasheet, including Sample_Type by default")
master_xtabs(sample_sheet,columns=c("sample_type",names(sample.subsetlist)))

#subset sample sheet for this study
message("Subsetting datasheet")
sub_sample_sheet<-subset_mastersheet(sample_sheet, sample.subsetlist)

#check again to see the subset made sense
master_xtabs(sub_sample_sheet,columns=c("sample_type",names(sample.subsetlist)))

#write - Needed here?
write.table(data.frame(V1=c("#checklist_accession","#unique_name_prefix"),V2=c("ERC000011","")),
            file = sample_sheet,quote = F,row.names = F,sep = "\t",col.names = F)
write.table(sub_sample_sheet,file = sample_sheet,quote = F,row.names = F,sep = "\t",append = T)


########################################LIBRARY DATA
#first read mastersheet to get libraries wanted
master_sheet<-google.read.master.url(sheeturl)

#subset 
sub_master_sheet<-subset_mastersheet(master_sheet, library.subsetlist)

#read library sheet
library_sheet<-google.read.master.url(sheeturl,ws = "ENA_library_data")

#subset by sample sheet
library_sheet<-library_sheet[library_sheet$library_name %in% sub_master_sheet$library_name,]

########################################UPLOAD SEQUENCES
#get filenames and check they exist
files<-c(library_sheet$forward_file_name,library_sheet$reverse_file_name)
message(paste(sum(files %in% list.files()), "of", length(files), "files found"))
write.table(files,file = filelist,quote = F,row.names = F,sep = "\t",col.names = F)

#cp files to temp folder

#run Aspera
            
