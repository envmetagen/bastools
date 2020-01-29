
packages.needed<-c("dplyr","googlesheets4","processx","httpuv","treemap","igraph","XML","data.tree","DT","ggplot2","vegan"
                   ,"stringr")

for(i in 1:length(packages.needed)){
  if(!packages.needed[i] %in% rownames(installed.packages())) {
    message(paste(packages.needed[i],"not found. Installing..."))
    install.packages(packages.needed[i])
  }
}

#check robitools
robitools<-c("ROBIBarcodes","ROBITools","ROBITaxonomy","ROBIUtils")
for(i in 1:length(robitools)){
  if(!robitools[i] %in% rownames(installed.packages())) message(paste(robitools[i],"not found. Follow next steps"))
}

#bas.install.robitools
path_to_parent<-"/home/tutorial/TOOLS/bastools/ROBITOOLS"
file.remove(list.files(path = path_to_parent,recursive = T,full.names = T,pattern = "\\.o$"))

folders<-list.files(path = path_to_parent,full.names = T)
if(length(grep(".zip",folders))>0) folders<-folders[-grep(".zip",folders)]
folders<-folders[c(grep("ROBITaxonomy",folders),grep("ROBIUtils",folders),grep("ROBIBarcodes",folders),grep("ROBITools",folders))]
folders[4]<-paste0(grep("ROBITools",folders,value = T),"/ROBITools")

for(i in 1:length(folders)){
  system2(command = "R",args = c("CMD", "INSTALL","--no-multiarch", "--with-keep.source",folders[i]), wait = T)
}

#use ctrl+shift+F10 to restart R

#get taxonomy dbs
source("/home/tutorial/TOOLS/bastools/master_functions.R")
getNCBItaxonomy(path = "/home/tutorial/DATABASES/nt_taxdump")
NCBI2obitaxonomy(path = "/home/tutorial/DATABASES/nt_taxdump",out = "/home/tutorial/DATABASES/obitaxdb_Jan2020")
