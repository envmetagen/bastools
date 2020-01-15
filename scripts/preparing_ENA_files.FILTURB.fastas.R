#making sample sheet, library sheet (adding md5sums), file list and password file

#CONFIG
#bastoolsDir<-"/Users/basti/Documents/WORK/CIBIO/STATS_AND_CODE/git-bastools/bastools/"
bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/"
#datasheet
sheeturl<-"https://docs.google.com/spreadsheets/d/1FUSaeVaYzms2EOGUoCAB4jaRKzguD3AKTsC8lYwaKP4/edit?ts=5de15df2#gid=1531090624"
#options for subsetting sample sheet - if none use sample.subsetlist<-NULL
sample.subsetlist<-list(experiment_id=c("SANG_1"))
#sample.subsetlist<-NULL
#options for subsetting library sheet 
library.subsetlist<-list(experiment_id=c("SANG_1"))
#email
email="basegeter@gmail.com"

#outdir ENA upload sheets will be saved here
outdir<-"/home/bastian.egeter/DATABASES/"

#name to give list of files
filelist<-"GUELTA1_ENA_FILES.txt"
#ENA password
passwd<-"your_password"


########################################CODE

setwd(bastoolsDir)
googlesheets4::sheets_auth(email = email)
source("googlesheet.foos.R")
setwd(outdir)

########################################SAMPLE DATA
#first read library_data sheet
library_sheet<-google.read.master.url(sheeturl,ws = "Library_data")

sub.library_sheet<-subset_mastersheet(library_sheet,library.subsetlist)

