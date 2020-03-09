#Minion Pipeline 2

#assumes fastq.gz files as starting, named barcode01.fastq.gz, barcode02.fastq.gz ...

bastoolsDir<-"/home/tutorial/TOOLS/bastools/" #change to your bastools directory
outDir<-"/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/MUSSELS/March_2020/2019_September_001_Run8/" #must exist
filesDir<-"/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/MUSSELS/March_2020/2019_September_001_Run8/fastq_files/"
# get sheet first separately (better than including in pipe)
# source(paste0(bastoolsDir,"master_functions.R"))
# mastersheet<-google.overlord("https://docs.google.com/spreadsheets/d/1k1mAGogWq9rXcwBKDyxG9oZ0OrWreBRRcSEUw0RGwyk/edit?ts=5d776492#gid=1377121809")
# write.table(mastersheet,paste0(outDir,"mastersheet.txt"),append = F,quote = F,row.names = F,sep = "\t")
mastersheet<-paste0(outDir,"mastersheet.txt")
#options for subsetting master sheet. This functions to sleect the samples you want to analyse.
#Each item in list is a column heading in master sheet and each character within the item should be what you want to include 
#(sample_type should always be lower case, even if it is not so on google)
subsetlist<-list(experiment_id="2019_September_001",primer_set=c("VENE","UNIO"))
##set taxidlimit and taxidname if desired. 
taxidlimit=33208 #set to NULL if not required
taxidname="metazoa"

                 