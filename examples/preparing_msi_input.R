#preparing msi
url="https://docs.google.com/spreadsheets/d/14QeZ0KyEHpqCZ3dsj4a6iojxnuLbasR7Z6y68FK9vn8/edit#gid=1377121809"
METADATAFILE="/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/MINION/2020_July_NM/master_sheet.txt"
subsetList=list(experiment_id="2020_July_NM")

##############################################################

####script
source("/home/bastian.egeter/git_bastools/bastools/master_functions.R")

a<-google.overlord(url,for.post.illscript2 = F,for.msi = T,subsetList = subsetList)

write.table(a,METADATAFILE,quote = F,sep = "\t",row.names = F)
