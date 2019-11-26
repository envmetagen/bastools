#uploading ENA data via Aspera

#CONFIG
#bastoolsDir<-"/Users/basti/Documents/WORK/CIBIO/STATS_AND_CODE/git-bastools/bastools/"
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

#outdir ENA upload sheets will be saved here
outdir<-"/mnt/Disk2/MISEQ_RUNS/2018_07_FILTURB_IRANVERTS/FASTQ_GENERATION_3_18/"

#final sample sheet filename
sample_sheet_name<-"FILTURB_NORTH_ENA_SAMPLE_SHEET.txt"

#final library sheet filename
library_sheet_name<-"FILTURB_NORTH_ENA_LIBRARY_SHEET.txt"

#name to give list of files
filelist<-"FILTURB_NORTH_ENA_FILES.txt"

#make md5 files?
makemd5s=T

#CODE

setwd(bastoolsDir)
googlesheets4::sheets_auth(email = email)
source("googlesheet.foos.R")
setwd(outdir)

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

#write 
setwd(outdir)
write.table(data.frame(V1=c("#checklist_accession","#unique_name_prefix"),V2=c("ERC000011","")),
            file = sample_sheet_name,quote = F,row.names = F,sep = "\t",col.names = F)
write.table(sub_sample_sheet,file = sample_sheet_name,quote = F,row.names = F,sep = "\t",append = T)


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
setwd(fileDir)
files<-c(library_sheet$forward_file_name,library_sheet$reverse_file_name)
message(paste(sum(files %in% list.files()), "of", length(files), "files found, of which",
              length(files[duplicated(files)]), "are duplicates", files[duplicated(files)]))

write.table(files,file = filelist,quote = F,row.names = F,sep = "\t",col.names = F)


#creating md5s
if(makemd5s){
  
  checksum<-data.frame(md5=rep(0,length(files)))
  
for(i in 1:length(files)){
  checksum[i,1]<-system2("md5sum", args=files[i],stdout = T)
  #a<-read.table(paste0(files[i],".md5"))
  #a$V2=NULL
  #write.table(a,file = paste0(files[i],".md5"),quote = F,row.names = F,sep = "\t",col.names = F)
}
  checksum2<-as.data.frame(do.call(rbind,stringr::str_split(checksum$md5," ")))
  checksum2$V2<-NULL

#add md5 to library data
library_sheet$forward_file_md5<-checksum2$V1[1:length(library_sheet$forward_file_name)]
library_sheet$reverse_file_md5<-checksum2$V1[length(library_sheet$reverse_file_name)+1:
                                                length(library_sheet$reverse_file_name)]
}

#write 
write.table(x = library_sheet,file = library_sheet_name,quote = F,row.names = F,sep = "\t",col.names = T)


write.table("tigress33",file = "ENA_password.txt",quote = F,row.names = F,sep = "\t",col.names = F)

#run from terminal
#emg_ena_upload.sh -l filelist -u Webin-51203 -p ENA_password.txt
        
#go to https://wwwdev.ebi.ac.uk/ena/submit/sra/#home #FOR TESTING
#or here for real https://www.ebi.ac.uk/ena/submit/sra/

####need to see if can upload md5s directly to google sheet
