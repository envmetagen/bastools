#transform the output fastas from msi (fastas in multiple folders) into usual for illscript3,
#i.e. the same as the output from mbc
#output fastas and otutabs will be split by primer
#reads with no_adapter are removed

root_folder<-"/media/sf_Documents/WORK/CIBIO/temp/MINION_TEST/2020_July_NM_B/"
master_sheet<-"/media/sf_Documents/WORK/CIBIO/temp/MINION_TEST/2020_July_NM_B/master_sheet.txt" 
#master sheet for single run only
fasta.out<-"/media/sf_Documents/WORK/CIBIO/temp/MINION_TEST/2020_July_NM_B/2020_July_NM_B.fasta"
otutab.out<-"/media/sf_Documents/WORK/CIBIO/temp/MINION_TEST/2020_July_NM_B/2020_July_NM_B.otutab.tsv"
#if multiple primers, the files will be renamed accordingly

#run
source("/home/bastian.egeter/git_bastools/bastools/scripts/msi2usual.R")
#source("/home/tutorial/TOOLS/bastools/scripts/msi2usual.R")
