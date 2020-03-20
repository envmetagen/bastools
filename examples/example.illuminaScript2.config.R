#############################################################################

#change setting below as necessary

#once settings are changed, from terminal run
# nohup Rscript example.config.R  2>run.processing.log 1>&2 &

###################
#SETTINGS

#step7 - blast
#step8 - filter blast results
#step9 - bin
#step10 - merge with MBC otutab
#step11 - apply minor taxon filter
#step12 - split taxatables by project
#step13 - make contributor files for later inspection
#step14 - make krona plots

#delete any steps you do not want to take
stepstotake<-c("step7","step8","step9","step10","step11","step12","step13","step14")

#set starting fastas, as well as taxidlimits and taxidnames if desired. 
#can be multiple. Taxidlimits and names must be in same order as fastas

#example 1
# startingfastas<-data.frame(
#   infasta=c("12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.fasta",
#             "16S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.fasta"),
#                            taxidlimit=c(7711,40674),
#                            taxidname=c("chordata","mammalia"))

####To blast entire database just use
#startingfastas<-c("file1.fasta","file2.fasta")

startingfastas<-data.frame(infasta=
  "none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.fasta",
  taxidlimit=c(7711),
  taxidname=c("chordata"))

#all directories must have trailing "/"

#full path to bastools directory
bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/" 
#full path to out dir: this must exist and already contain the final fasta and corresponding otu table 
outDir<-"/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/HSJUN19BAS/NT/" 
#full path to ncbi taxonomy directory
ncbiTaxDir<-"/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/taxdump/Jan2020/"
#full path to blast exec
blast_exec<-"/home/bastian.egeter/Tools/ncbi-blast-2.9.0+/bin/blastn"
#full path to Krona exec
KronaPath<-"/home/bastian.egeter/Tools/Krona.install/bin/ktImportText"
#full path to blast db - this is the stem filename of the db files
refdb = "/mnt/Disk1/Tools/BLAST+/DBs/nt_v5/nt_v5"
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



