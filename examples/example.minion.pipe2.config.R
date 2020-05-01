#Minion Pipeline 2
#assumes fastq.gz files as starting, named barcode01.fastq.gz, barcode02.fastq.gz ...
#or if polished pipeline a single file called "results.fasta.gz"
#usually, one run, one primer

#once settings are changed, to run, and record time properly, do the following in terminal:
#cmd='nohup Rscript example.minion.pipe2.config.R'
#datetime=`date "+%F %R"`
#/usr/bin/time -q  -a --format "%E\t%e\t%M\t$datetime\t$*\t%x\t%I\t%O\t%W" bash -c "$cmd" 2>example.minion.pipe2.log 1>&2 &

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
#step12 - make contributor files
#step13 - make krona plots

#stepstotake=c("step1","step2","step3","step4","step5","step6","step7","step8","step9","step10","step11","step12","step13")
stepstotake="all"

#common settings:

  #all dirs must have trailing /
  bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/" #bastools directory
  blast_exec<-"/home/bastian.egeter/Tools/ncbi-blast-2.9.0+/bin/blastn"
  ncbiTaxDir="/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/taxdump/Jan2020/"
  obitaxdb<-"/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/obitaxdump/Jan2020/obitaxdb_jan_2020"
  KronaPath<-"/home/bastian.egeter/Tools/Krona.install/bin/ktImportText"
  preffered_rlibsfolder="~/MBC/Rlibs/" #silly workaround for emg1

#specific settings:
  
  #specifiy location of unpolished raw reads, must be separate folder
  origfilesDir<-"/home/bastian.egeter/Minion_data/MINION.PIPE2/2019_August_002/raw/"
  
  #location of mastersheet, results.fasta.gz (if applicable), all outputs will be put here 
  filesDir<-"/home/bastian.egeter/Minion_data/MINION.PIPE2/2019_August_002/polished/"
  
  #used polishing script?
  polished=T
  
  # get sheet first separately (better than including in pipe)
  # source(paste0(bastoolsDir,"master_functions.R"))
  # mastersheet<-google.overlord("https://docs.google.com/spreadsheets/d/1k1mAGogWq9rXcwBKDyxG9oZ0OrWreBRRcSEUw0RGwyk/edit?ts=5d776492#gid=1377121809",email="basegeter@gmail.com")
  # write.table(mastersheet,paste0(filesDir,"mastersheet.txt"),append = F,quote = F,row.names = F,sep = "\t")
  mastersheet<-"MUSSELS_mastersheet.txt"
  #options for subsetting master sheet. This functions to select the samples you want to analyse.
  #Each item in list is a column heading in master sheet and each character within the item should be what you want to include 
  #(sample_type should always be lower case, even if it is not so on google)
  subsetlist<-list(experiment_id="2019_August_002",primer_set=c("UNIO"))
  
  #cutadapt error rate
  ca.error.rate=0.2
  
  #add a suffix to the final fasta, useful for running different pipelines from step4 onwards. put as NULL to ignore
  catted_suffix<-"SC1"
  
  ##BLAST settings 
  taxidlimit=6544 #set to NULL if not required
  refdb="/mnt/Disk2/nt_April2020/nt"
  opts=c("-task","blastn","-outfmt", "6 qseqid pident qcovs saccver staxid ssciname","-num_threads", 64,
         "-max_target_seqs", 50, "-max_hsps",1,"-word_size", 11, "-perc_identity", 50,
         "-qcov_hsp_perc", 98, "-gapopen", 0, "-gapextend", 2, "-reward", 1, "-penalty", -1)
  
  #filtering/binning options
  spident=99
  topS=2
  gpident=97
  topG=4
  fpident=95
  topF=5
  abspident=93
  topAbs=7
  
  disabledTaxaFiles="/home/bastian.egeter/Minion_data/MINION.PIPE2/2019_August_002/disabled.taxa.txt" #otherwise NULL

source(paste0(bastoolsDir,"scripts/bas.minion.pipeline2.R"))
