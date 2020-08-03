#preparing msi
url="https://docs.google.com/spreadsheets/d/14QeZ0KyEHpqCZ3dsj4a6iojxnuLbasR7Z6y68FK9vn8/edit#gid=1377121809"
METADATAFILE="/home/bastian.egeter/msi_output/2020_July_NM/master_sheet.txt"

##############################################################

####script
a<-google.overlord(url,for.post.illscript2 = F,for.msi = T,subsetList = subsetList)

write.table(a,METADATAFILE,quote = F,sep = "\t",row.names = F)