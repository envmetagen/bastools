#############################################################################
#test.config

#first create processing sheet using create.mastersheet.R

#the parent directory should be empty also, as the outputs will be put here
#change setting below as necessary

#once settings are changed, from terminal run
# nohup Rscript test.config.R  2>run.processing.log 1>&2 &


###################
#SETTINGS

#set starting fastas, taxidlimits, and taxidnames if desired. 
startingfastas<-c(
  "2018_02_12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.fasta"
  ,"2018_04_12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.fasta"
  ,"2018_07_12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.fasta"
  )
####To blast entire nt just use
#startingfastas<-c("file1.fasta","file2.fasta")
#stepstotake<-c("step0","step7","step8","step9","step10","step11","step12","step13","step14")


stepstotake<-c("step9","step10","step11","step12","step13","step14")

#all directories must have trailing "/"

bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/"
#full path to rootir
outDir<-"/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/FILTURB_12S_BLAST/"
#full path to ncbi taxonomy directory
ncbiTaxDir<-"/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/taxdump/October-2019/"
#full path to obitools taxonomy
obitaxdb<-"/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/obitaxdump/October-2019/obitaxdb_oct_2019"
#full path to blast exec (or just # it)
blast_exec<-"/home/bastian.egeter/Tools/ncbi-blast-2.9.0+/bin/blastn"
#full path to blast db
refdb = "/home/bastian.egeter/DATABASES/12S_Portugal_amphibians_DB_for_north_paper"
#sheet urls (here using "master sheet")
sheeturls<-c("https://docs.google.com/spreadsheets/d/1FUSaeVaYzms2EOGUoCAB4jaRKzguD3AKTsC8lYwaKP4/edit?ts=5dae01be#gid=0")
experiment_id<-c("2018_02","2018_04","2018_07")
#minimum query coverage for filtering blast results
min_qcovs=70
#maximum evalue for filtering blast results
max_evalue=0.001
#percentage identity threshold to consider (1=1%)
top=1
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


source("/home/bastian.egeter/git_bastools/bastools/scripts/illumina_post_MBC_script2.R")



