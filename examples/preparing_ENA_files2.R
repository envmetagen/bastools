#making sample sheet, library sheet (adding md5sums), file list and password file

#first make a copy this file, do not edit this directly
#change settings as desired
#run in R, e.g. from terminal Rscript this_file.R

#CONFIG
ms_option=1 #choose 1 or 2, see below

#full path. Must already exist. If using Option 2 below, this folder must contain the 3 tsv files required.
#Final ENA upload files will also be added here
outdir<-"/home/tutorial/test/"
#directory with files for upload
fileDir<-"/mnt/Disk2/MISEQ_RUNS/2017_12_AZORES-IRANVERTS-NZFROG-OZ-GUELTA/PP121217-56327499/BEFORE_PAIRED_END/FIRST_ATTEMPT/"

####ms_option=1
#use google. This is ok, but can be problematic because 
#it will try to open browser if you havent stored credentials before and this is too slow if accessing emg machines . 
#If using it set the following options, otherwise leave as they are
#datasheet
sheeturl<-"https://docs.google.com/spreadsheets/d/1ru01umKc4tCvNxUlxUJPcUyz9jon1L-QJUm67etL3aE/edit?usp=sharing"
#email  
email="basegeter@gmail.com"

####ms_option=2
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

#name to give to final sample sheet filename
sample_sheet_name<-"filipa_test_ENA_SAMPLE_SHEET.txt"
#name to give to final library sheet filename
library_sheet_name<-"filipa_test_ENA_LIBRARY_SHEET.txt"
#name to give to list of files for upload
filelist<-"filipa_test_ENA_FILES.txt"
#ENA password
passwd<-"xxx4452"

source("/home/bastian.egeter/git_bastools/bastools/scripts/ENA_prep_script.R")