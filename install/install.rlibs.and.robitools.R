
bastoolsDir<-"/home/tutorial/TOOLS/bastools/" #change to your bastools dir. must have trailing "/"

update.taxonomies=F #change to T if needed and amend directories below to desired locations

if(update.taxonomies==T) {
 source(paste0(bastoolsDir,"master_functions.R"))
   getNCBItaxonomy(path = "/home/tutorial/DATABASES/nt_taxdump")
   NCBI2obitaxonomy(path = "/home/tutorial/DATABASES/nt_taxdump",out = "/home/tutorial/DATABASES/obitaxdb_Jan2020")
}

packages.needed<-c("dplyr","googlesheets4","processx","httpuv","treemap","igraph","XML","data.tree","DT","ggplot2","vegan"
                   ,"stringr","phylotools","htmlwidgets","tidyverse","cowplot","ggpubr","insect","plotly","seqRFLP","rlang","fs","purrr")

for(i in 1:length(packages.needed)){
  if(!packages.needed[i] %in% rownames(installed.packages())) {
    message(paste(packages.needed[i],"not found. Installing..."))
    install.packages(packages.needed[i])
  }
}

#check robitools
robitools<-c("ROBIBarcodes","ROBITools","ROBITaxonomy","ROBIUtils")
robi.installed<-list()
for(i in 1:length(robitools)){
  if(!robitools[i] %in% rownames(installed.packages())) {robi.installed[[i]]<-F} else {robi.installed[[i]]<-T}  
}

if(sum(do.call(rbind,robi.installed))<4){
#bas.install.robitools
file.remove(list.files(path = paste0(bastoolsDir,"ROBITOOLS"),recursive = T,full.names = T,pattern = "\\.o$"))

folders<-list.files(path = paste0(bastoolsDir,"ROBITOOLS"),full.names = T)
if(length(grep(".zip",folders))>0) folders<-folders[-grep(".zip",folders)]
folders<-folders[c(grep("ROBITaxonomy",folders),grep("ROBIUtils",folders),grep("ROBIBarcodes",folders),grep("ROBITools",folders))]
folders[4]<-paste0(grep("ROBITools",folders,value = T),"/ROBITools")

for(i in 1:length(folders)){
  system2(command = "R",args = c("CMD", "INSTALL","--no-multiarch", "--with-keep.source",folders[i]), wait = T)
}
}

#use ctrl+shift+F10 to restart R


