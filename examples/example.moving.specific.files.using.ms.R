#source("/home/bastian.egeter/git_bastools/bastools/master_functions.R")

source("/home/tutorial/TOOLS/bastools/master_functions.R")

#set to original file directory
setwd("/mnt/Disk2/MISEQ_RUNS/2018_02_MICROSATS-SOILPHOS-FILTURB-NZFROG-IRANVERTS-ICVERTS-GUELTA/original_miseq_files/")

a<-google.overlord("https://docs.google.com/spreadsheets/d/1E0LjEqJx4VJQ2pZEPf5EKrhHdkGeo99ZviKjaXxZCPs/edit#gid=0",
                subsetList = list(experiment_id="2018_02"))

files<-unlist(stringr::str_split(a$sequencing_filenames," "))

for(i in 1:length(files)){
  file.copy(from = files[i],to = "./temp/") #must exist
}

write.table(a,"./temp.copy/SOILPHOS_2018_02_PHOD_master.txt",append = F,quote = F,row.names = F,sep = "\t")

file.copy(from = "test.binthresh.out.tsv",to = "./temp/")
file.remove("./temp/test.binthresh.out.tsv")
