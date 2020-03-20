#Minion Pipeline 2

#assumes fastq.gz files as starting, named barcode01.fastq.gz, barcode02.fastq.gz ...

#usually, one run, one primer

#step1 - cutadapt
#step2 - print length data
#step3 - size select
#step4 - cat files
#step5 - blast
#step6 - filter blast results
#step7 - bin blast results
#step8 - make otutab
#step9 - merge with otutab
#step10 - make contributor files
#step11 - make krona plots

stepstotake=c("step4","step5","step6","step7","step8","step9","step10","step11")

#all dirs must have trailing /
bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/" #bastools directory
filesDir<-"/home/bastian.egeter/Minion_data/MINION.PIPE2/2019_August_002/"
blast_exec<-"/home/bastian.egeter/Tools/ncbi-blast-2.9.0+/bin/blastn"
ncbiTaxDir="/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/taxdump/Jan2020/"
KronaPath<-"/home/bastian.egeter/Tools/Krona.install/bin/ktImportText"

preffered_rlibsfolder="~/MBC_pipelines/MBC/Rlibs" #had to add this for emg1

# get sheet first separately (better than including in pipe)
# source(paste0(bastoolsDir,"master_functions.R"))
# mastersheet<-google.overlord("https://docs.google.com/spreadsheets/d/1k1mAGogWq9rXcwBKDyxG9oZ0OrWreBRRcSEUw0RGwyk/edit?ts=5d776492#gid=1377121809",email="basegeter@gmail.com")
# write.table(mastersheet,paste0(outDir,"mastersheet.txt"),append = F,quote = F,row.names = F,sep = "\t")
mastersheet<-"mastersheet.txt"
#options for subsetting master sheet. This functions to select the samples you want to analyse.
#Each item in list is a column heading in master sheet and each character within the item should be what you want to include 
#(sample_type should always be lower case, even if it is not so on google)
subsetlist<-list(experiment_id="2019_August_002",primer_set=c("UNIO"))

#cutadapt error rate
ca.error.rate=0.3

#add a suffix to the final fasta, useful for running different pipelines from step4 onwards. put as NULL to ignore
catted_suffix<-"SC3"

##set BLAST taxidlimit and taxidname if desired. 
taxidlimit=NULL #set to NULL if not required
taxidname=NULL
refdb="/home/bastian.egeter/DATABASES/16S_Bivalves/16S.s2.Bivalvia.7species"
max_target_seqs = 50

#filtering/binning options #based on results, scenario3
#https://docs.google.com/presentation/d/1NnE9yQJzuRChQaj0GnSuZWq1QmZQO95SBb9oKluD1vI/edit#slide=id.g71516c64c3_0_190
min_qcovs=90
max_evalue = 10

spident=97
topS=30
gpident=90
topG=30
fpident=70
topF=30
abspident=60 #arbitrary
topAbs=30

source(paste0(bastoolsDir,"scripts/bas.minion.pipeline2.R"))
