
packages.needed<-c("dplyr","googlesheets4","processx","httpuv","treemap","igraph","XML","data.tree","DT","ggplot2")

for(i in 1:length(packages.needed)){
  if(!packages.needed[i] %in% rownames(installed.packages())) {
    message(paste(packages.needed[i],"not found. Installing..."))
    install.packages(packages.needed[i],dependencies = T)
  }
}

#check robitools
robitools<-c("ROBIBarcodes","ROBITools","ROBITaxonomy","ROBIUtils")
for(i in 1:length(robitools)){
  if(!robitools[i] %in% rownames(installed.packages())) message(paste(robitools[i],"not found. Follow next steps"))
}

#bas.install.robitools
path_to_parent<-"/home/bastian.egeter/ROBITOOLS"
file.remove(list.files(path = path_to_parent,recursive = T,full.names = T,pattern = "\\.o$"))

folders<-list.files(path = path_to_parent,full.names = T)
folders<-folders[-grep(".zip",folders)]
folders<-folders[c(grep("ROBITaxonomy",folders),grep("ROBIUtils",folders),grep("ROBIBarcodes",folders),grep("ROBITools",folders))]
folders[4]<-paste0(grep("ROBITools",folders,value = T),"/ROBITools")

for(i in 1:length(folders)){
  system2(command = "R",args = c("CMD", "INSTALL","--no-multiarch", "--with-keep.source",folders[i]), wait = T)
}

#use ctrl+shift+F10 to restart R