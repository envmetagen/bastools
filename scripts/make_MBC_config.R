#make MBC processing sheet
source("/home/bastian.egeter/git_bastools/bastools/master_functions.R")

url<-("https://docs.google.com/spreadsheets/d/1Kw8ONKtlX1hGr0LM-FvkowXYElEJ_taFNu2ZSXqp_n8/edit#gid=0")
tokenDir<-"/home/bastian.egeter/git_bastools/bastools/"
email="basegeter@gmail.com"
subsetlist<-list(experiment_id=c("2020_01"), Primer_set=c("12SV5.1"))
out<-"IRAN_experiment_sheet.tsv"

#first setwd to target directory for sheet
##then:
google.make.MBC.sheet(tokenDir,email,url,out,subsetlist)

#command for MBC
#source path_where_mbc_was_installed/env.sh
#mbc conf=12SV1.conf step6
  