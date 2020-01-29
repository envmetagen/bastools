#############################################################################

#once settings are changed, from terminal run
# nohup Rscript example.rebin.illumina.config.R  2>run.processing.log 1>&2 &

#THIS IS FOR IRAN ILLUMINA RUNS

###################
#SETTINGS

#step9 - rebin
#step10 - merge with MBC otutab
#step11 - apply minor taxon filter
#step12 - cut out project taxatables
#step13 - make new contributor files for later inspection
#step14 - make krona plots

stepstotake<-c("step9","step10","step11","step12","step13","step14")

#full paths to filtered blastfiles, can be multiple
files<-c("/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/2020_01/none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.blast.filt.txt" 
    )
#exp ids in same order as files
experiment_id<-c("2020_01")

#full path to disabled taxa files, can be multiple (must have "contributors",taxids", "disable_species","disable_genus","disable_family"
#as colnames, can have other cols, tab-sep)
disabledTaxaFiles<-c("/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/POST_TAXONOMY_CHECK/IRAN_REBIN/Disabled_taxa_All_Table.txt")

#dont need to change, just an output of final disabled taxa table for the record (possibly obsolete now)
disabledTaxaOut<-"final.disabled.taxa.after.binfunc.txt"

#all directories must have trailing "/"

#full path to bastools directory
bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/" 

#full path to out dir: this must exist and should be more or less empty
outDir<-"/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/2020_01/rebins/" 

#full path to ncbi taxonomy directory
ncbiTaxDir<-"/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/taxdump/Jan2020/"

#full path to obitools taxonomy - this is the stem filename of the obitaxonomy files
obitaxdb<-"/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/obitaxdump/Jan2020/obitaxdb_jan_2020" 

#full path to Krona exec
KronaPath<-"/home/tutorial/TOOLS/Krona.install/bin/ktImportText"

project_name<-"IRAN-" #INCLUDE THE DASH HERE

#Use force=T to ignore any contributor entries where no levels were disabled when consistency checking. 
#Use force=F to do more thorough consistency checks
force=F

#THE FOLLOWING SETTINGS MUST BE IDENTICAL TO THOSE USED FOR ORIGINAL BINNING 

#percentage identity for species level
spident=98
#percentage identity for genus level
gpident=95
#percentage identity for family level
fpident=92
#percentage identity for higher-than-family level
abspident=80
#to filter taxatable by %taxon (0.1=0.1%) 
filterpc=0.1

source("/home/tutorial/TOOLS/bastools/scripts/rebin.pipe.illumina.R")
