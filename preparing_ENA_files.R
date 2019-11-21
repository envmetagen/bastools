#uploading ENA data via Aspera

#CONFIG
bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/"

#filturb datasheet
sheeturl<-"https://docs.google.com/spreadsheets/d/1FUSaeVaYzms2EOGUoCAB4jaRKzguD3AKTsC8lYwaKP4/edit?ts=5dae01be#gid=1531090624"

#options for subsetting sample sheet 
subsetlist<-list(study=c("North"))

#email
email="basegeter@gmail.com"


#CODE

setwd(bastoolsDir)
googlesheets4::sheets_auth(email = email)
source("/home/bastian.egeter/git_bastools/bastools/googlesheet.foos.R")

#first read ENA_sample_data sheet
master_sheet<-google.read.master.url(sheeturl,ws = "ENA_sample_data")

#check on number of samples in each category
message("Showing xtabs for entire datasheet, including Sample_Type by default")
master_xtabs(master_sheet,columns=c("sample_type",names(subsetlist)))

#subset sample sheet for this study
message("Subsetting datasheet")
sub_sample_sheet<-subset_mastersheet(master_sheet, subsetlist)

#check again to see the subset made sense
master_xtabs(sub_sample_sheet,columns=c("sample_type",names(subsetlist)))

#write
write.table(sub_sample_sheet,file = paste0(gsub(" ","_",ss_info$sheet_title),".txt"),quote = F,row.names = F,sep = "\t")