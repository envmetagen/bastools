#making sample sheet, library sheet (adding md5sums), file list and password file

##############INSTRUCTIONS
#works on emg1 and emg2
#first make a copy this file, do not edit this directly
#change settings as desired and save this file
#then, from terminal 
#"Rscript this_file.R"

#to upload fastq.gz files to ena:

#source /opt/env.sh
#cd to the filelist directoy (here /mnt/Disk2/FilipaMSMartins/7.run_180201/ena/upload/)
# nohup emg_ena_upload.sh -l filelist -u Webin-XXXX -p ENA_password.txt 2>upload.log 1>&2 & 
#(change name of log file if you wish)
#replacing 'filelist' with the filelist output below (here "filipa_test_ENA_FILES.txt") 
#replace Webin-XXXX with your username

#to complete metadata entry, use the output files (sample_sheet_name and library_sheet_name) that were created and enter them here:
#go to https://wwwdev.ebi.ac.uk/ena/submit/sra/#home #FOR TESTING
#or here for real https://www.ebi.ac.uk/ena/submit/sra/

#if you havent done so already you will need to create a project first
#read my notes for more
#https://docs.google.com/document/d/1DOHRvA9cjwfrsLDkeZ4SRfIP7GIo9d0dxl5vjTBqfFY/edit#

##################################

#CONFIG
ms_option=1 #choose 1 or 2, see below

#full path, including trailing '/'. Must already exist. If using Option 2 below, this folder must contain the 3 tsv files required.
#Final ENA upload files (sample_sheet_name and library_sheet_name) will also be added here
#the filelist and password file will be written to fileDir
outdir<-"/mnt/Disk2/FilipaMSMartins/7.run_180201/ena"
#directory with files for upload
fileDir<-"/mnt/Disk2/FilipaMSMartins/7.run_180201/ena/upload/"

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