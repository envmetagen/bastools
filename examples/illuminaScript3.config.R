#############################################################################
#test.config

#first create processing sheet using create.mastersheet.R

#put all final fastas (usually one per primer) in a folder called "final_fastas", with no other fastas.
#put all final otutabs (usually one per primer) in a folder called "final_otutabs", with no other tabs.
#the parent directory should be empty also, as the outputs will be put here
#change setting below as necessary

#once settings are changed, from terminal run
# nohup Rscript test.config.R  2>run.processing.log 1>&2 &

  
###################
#SETTINGS

#set starting fastas, taxidlimits, and taxidnames if desired. 
startingfastas<-"12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.fasta"
####To blast entire nt just use
#startingfastas<-c("file1.fasta","file2.fasta")
#stepstotake<-c("step7","step8","step8a","step9","step10","step11","step12","step13","step14")
stepstotake<-c("step9","step10","step11","step12","step13","step14")

#all directories must have trailing "/"

bastoolsDir<-"/home/tutorial/TOOLS/bastools/"
#full path to rootir
outDir<-"/home/tutorial/test/"
#full path to ncbi taxonomy directory
ncbiTaxDir<-"/home/tutorial/TOOLS/metabinkit.install/db"
#full path to obitools taxonomy
#obitaxdb<-"/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/obitaxdump/October-2019/obitaxdb_oct_2019"
#full path to blast exec (or just # it)
blast_exec<-"/home/tutorial/TOOLS/metabinkit.install/bin/blastn"
#full path to blast db
refdb = "/mnt/Disk1/Tools/BLAST+/DBs/nt_v5/nt"
#sheet urls (here using "processing sheet")
#sheeturls<-"https://docs.google.com/spreadsheets/d/1Eu8EDcGgeGh3yv-q3ATq4uMrRsA3iL9Gqy4HfYhRCoU/edit#gid=0"
#experiment_id<-"2018_04"
taxidlimit=NULL # taxidlimit=c("2210","9606") in same order as starting fastas  
opts=c("-task","blastn","-outfmt", "6 qseqid pident qcovs saccver staxid ssciname sseq","-num_threads", 64,
       "-max_target_seqs", 50, "-max_hsps",1,"-word_size", 11, "-perc_identity", 50,
       "-qcov_hsp_perc", 98, "-gapopen", 0, "-gapextend", 2, "-reward", 1, "-penalty", -1)

#minimum query coverage for filtering blast results
min_qcovs=70
#maximum evalue for filtering blast results
max_evalue=0.001

use.metabin=T #use metabin or older script for processing?
KronaPath="/home/tutorial/TOOLS/Krona.install/bin/ktImportText" #full path
known_flags="/home/tutorial/REPTILE/auto_and_manual_flags_reptile.txt" #full path
#get from google?
use_flagged_accessions_mbk=T
TaxlevelTestall<-c("K","P","C","O","F") #TaxlevelTestall<-c("K","P","C","O","F","G","S")

#percentage identity threshold to consider (1=1%)
#top=1 #this will overwrite all other tops! hash it out to not use it!
topS=100
topG=100
topF=100
topAF=2

#percentage identity for species level
spident=98
#percentage identity for genus level
gpident=95
#percentage identity for family level
fpident=92
#percentage identity for higher-than-family level
abspident=80
#taxon filter percent (0.1=0.1%)
filterpc=0.1

## allow taxa_disabling
SpeciesBL="species_blacklist_test.txt"
GenusBL="genus_blacklist_test.txt"
FamilyBL=NULL


source("/home/bastian.egeter/git_bastools/bastools/scripts/illumina_post_MBC_script3.R")