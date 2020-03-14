#Minion Pipeline 2

#assumes fastq.gz files as starting, named barcode01.fastq.gz, barcode02.fastq.gz ...

#usually, one run, one primer

stepstotake=c("step1","step2","step3","step4","step5")

#all dirs must have trailing /
bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools//" #bastools directory
outDir<-"/home/bastian.egeter/Minion_data/MINION.PIPE2/2019_September_001/" #must exist
filesDir<-"/home/bastian.egeter/Minion_data/MINION.PIPE2/2019_September_001/"
blast_exec<-"/home/bastian.egeter/Tools/ncbi-blast-2.9.0+/bin/blastn"
# get sheet first separately (better than including in pipe)
# source(paste0(bastoolsDir,"master_functions.R"))
# mastersheet<-google.overlord("https://docs.google.com/spreadsheets/d/1k1mAGogWq9rXcwBKDyxG9oZ0OrWreBRRcSEUw0RGwyk/edit?ts=5d776492#gid=1377121809",email="basegeter@gmail.com")
# write.table(mastersheet,paste0(outDir,"mastersheet.txt"),append = F,quote = F,row.names = F,sep = "\t")
mastersheet<-paste0(outDir,"mastersheet.txt")
#options for subsetting master sheet. This functions to sleect the samples you want to analyse.
#Each item in list is a column heading in master sheet and each character within the item should be what you want to include 
#(sample_type should always be lower case, even if it is not so on google)
subsetlist<-list(experiment_id="2019_September_001",primer_set=c("VENE"))

#cutadapt error rate
ca.error.rate=0.3

##set taxidlimit and taxidname if desired. 
taxidlimit=NULL #set to NULL if not required
taxidname=NULL
refdb="/home/bastian.egeter/DATABASES/16S_Bivalves/16S.s2.Bivalvia"

source(paste0(bastoolsDir,"scripts/bas.minion.pipeline2.R"))      