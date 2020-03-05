#making sample sheet, library sheet (adding md5sums), file list and password file

#CONFIG
#bastoolsDir<-"/Users/basti/Documents/WORK/CIBIO/STATS_AND_CODE/git-bastools/bastools/"
bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/"
#datasheet
sheeturl<-"https://docs.google.com/spreadsheets/d/1LUN1siU0D7ZAkj7aq4p98rAoZaQJ8Cg44MzGuGN91cU/edit#gid=0"
#options for subsetting sample sheet - if none use sample.subsetlist<-NULL
#sample.subsetlist<-list(Site=c("Atar"))
sample.subsetlist<-list(notes="Used for A case study of introduced rats consuming endemic New Zealand frogs")
#options for subsetting library sheet (uses master_sheet)
library.subsetlist<-list(study="NZ1")
#email
email="basegeter@gmail.com"
#directory with files for upload
fileDir<-"/mnt/Disk2/MISEQ_RUNS/2017_12_AZORES-IRANVERTS-NZFROG-OZ-GUELTA/PP121217-56327499/BEFORE_PAIRED_END/FIRST_ATTEMPT/"
#outdir ENA upload sheets will be saved here
outdir<-"/mnt/Disk2/MISEQ_RUNS/2017_12_AZORES-IRANVERTS-NZFROG-OZ-GUELTA/"
#final sample sheet filename
sample_sheet_name<-"NZFROG_CASE_STUDY_ENA_SAMPLE_SHEET.txt"
#final library sheet filename
library_sheet_name<-"NZFROG_CASE_STUDY_ENA_LIBRARY_SHEET.txt"
#name to give list of files
filelist<-"NZFROG_CASE_STUDY_ENA_FILES.txt"
#ENA password
passwd<-"tigress33"


########################################CODE

source(paste0(bastoolsDir,"master_functions.R"))

setwd(bastoolsDir)
googlesheets4::sheets_auth(email = email)
setwd(outdir)

########################################SAMPLE DATA
#first read ENA_sample_data sheet
sample_sheet<-google.read.master.url(sheeturl,ws = "ENA_sample_data")

#check on number of samples in each category
message("Showing xtabs for entire datasheet, including Sample_Type by default")
if(!is.null(sample.subsetlist)){
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
#warning appending column names to file is ok

########################################LIBRARY DATA
#first read mastersheet to get libraries wanted
master_sheet<-google.overlord(url = sheeturl)

if(!"library_name" %in% colnames(master_sheet)) stop("Must have column library_name in Master_Samplesheet")

message("Count libraries before subsetting=",length(master_sheet$library_name))
#subset 
sub_master_sheet<-subset_mastersheet(master_sheet, library.subsetlist)
message("Count libraries after subsetting master sheet=",length(sub_master_sheet$library_name))

#read library sheet
library_sheet<-google.read.master.url(sheeturl,ws = "ENA_library_data")

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
