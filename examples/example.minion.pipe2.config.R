#Minion Pipeline 2

#assumes fastq.gz files as starting, named barcode01.fastq.gz, barcode02.fastq.gz ...

#once settings are changed, to run, and record time properly, do the following in terminal:
#cmd='nohup Rscript example.minion.pipe2.config.R'
#datetime=`date "+%F %R"`
#/usr/bin/time -q  -a --format "%E\t%e\t%M\t$datetime\t$*\t%x\t%I\t%O\t%W" bash -c "$cmd" 2>example.minion.pipe2.log 1>&2 &

#usually, one run, one primer

#step1 - cutadapt
#step2 - print length data
#step3 - size select
#step4 - count reads up to blast step
#step5 - cat files
#step6 - blast
#step7 - plot blast hits
#step8 - filter blast results
#step9 - bin blast results
#step10 - make otutab
#step11 - merge with otutab
#step12 - apply taxon filter
#step13 - make contributor files
#step14 - make krona plots

stepstotake=c("step1","step2","step3","step4","step5","step6","step7","step8","step9","step10","step11","step12","step13","step14")

#all dirs must have trailing /
bastoolsDir<-"/home/tutorial/TOOLS/bastools/" #bastools directory
blast_exec<-"/home/tutorial/TOOLS/ncbi-blast-2.9.0+/bin/blastn"
ncbiTaxDir="/home/tutorial/TOOLS/DBS/ncbi_taxonomy/taxdump/"
KronaPath<-"/home/tutorial/TOOLS/Krona.install/bin/ktImportText"

filesDir<-"/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/MUSSELS/March_2020/test_polishing/"
#specifiy location of unpolished raw reads, must be separate folder
origfilesDir<-"/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/MUSSELS/March_2020/test_polishing/raw/"


preffered_rlibsfolder="~/MBC/Rlibs/" #had to add this for emg1

#used polishing script (then start at step 2)
polished=T


# get sheet first separately (better than including in pipe)
# source(paste0(bastoolsDir,"master_functions.R"))
# mastersheet<-google.overlord("https://docs.google.com/spreadsheets/d/1k1mAGogWq9rXcwBKDyxG9oZ0OrWreBRRcSEUw0RGwyk/edit?ts=5d776492#gid=1377121809",email="basegeter@gmail.com")
# write.table(mastersheet,paste0(filesDir,"mastersheet.txt"),append = F,quote = F,row.names = F,sep = "\t")
mastersheet<-"MUSSELS_mastersheet.txt"
#options for subsetting master sheet. This functions to select the samples you want to analyse.
#Each item in list is a column heading in master sheet and each character within the item should be what you want to include 
#(sample_type should always be lower case, even if it is not so on google)
subsetlist<-list(experiment_id="2019_August_002",primer_set=c("VENE"))

#cutadapt error rate
ca.error.rate=0.3

#add a suffix to the final fasta, useful for running different pipelines from step4 onwards. put as NULL to ignore
catted_suffix<-"SC3"

##BLAST settings 
taxidlimit=NULL #set to NULL if not required
taxidname=NULL #only required if taxidlimit is not NULL
refdb="/home/tutorial/TOOLS/mussels/databases/16S_Bivalves_7_species/16S.s2.Bivalvia.7species"
opts=c("-task","blastn","-outfmt", "6 qseqid pident qcovs saccver staxid ssciname","-num_threads", 64,
       "-max_target_seqs", 100, "-max_hsps",1,"-word_size", 11, "-perc_identity", 50,
       "-qcov_hsp_perc", 98, "-gapopen", 0, "-gapextend", 2, "-reward", 1, "-penalty", -1)

#filtering/binning options
spident=97
topS=30
gpident=90
topG=30
fpident=70
topF=30
abspident=60 #arbitrary
topAbs=30

disabledTaxaFiles=NULL #otherwise full path to file

#taxonfilter
taxonpc=0.001

source(paste0(bastoolsDir,"scripts/bas.minion.pipeline2.R"))
