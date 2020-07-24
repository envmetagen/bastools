source("/home/tutorial/TOOLS/bastools/master_functions.R")
#convert msi output to taxatab

#msi<-data.table::fread("/home/tutorial/TOOLS/iranvert/results/run3/results.tsv.gz",data.table = F,fill = T)

msi<-data.table::fread("/home/tutorial/TOOLS/iranvert/results/run3/binres.tsv.gz",data.table = F,fill = T)
msi$taxon<-paste(msi$K,msi$P,msi$C,msi$O,msi$F,msi$G,msi$S,sep = ";")


#msi$taxon<-paste(msi$kingdom,msi$phylum,msi$class,msi$order,msi$family,msi$genus,msi$species,sep = ";")

msi$taxon<-gsub(";;;;;;","no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",msi$taxon)
msi$taxon<-gsub(";;;",";NA;NA;",msi$taxon)

msi<-msi[msi$adapter!="no_adapter",]

msi.split<-split(msi,f=msi$adapter)
message(unique(paste(msi$adapter," ")))

taxatab.list<-list()
for(i in 1:length(msi.split)){
taxatab.list[[i]]<-as.data.frame(tidyr::pivot_wider(msi.split[[i]][,c("sample","taxon","nreads")],names_from = sample,values_from=nreads,
                                       values_fill = 0,values_fn=sum))
}


a<-taxatab.list[[1]]
#a<-rm.0readtaxSam(a)    
master_sheet<-data.frame(ss_sample_id=colnames(a[,-1]))
hm<-taxatab.heatmap(taxatab = a,master_sheet = master_sheet,inc.values = F,taxafontsize=8,colfontsize=8)
#saving image
jpeg(filename="/home/tutorial/TOOLS/iranvert/results/taxa_tables/iran.hm.jpg",
     unit="in",
     width=27,
     height=13,
     pointsize=12,
     res=200)
ComplexHeatmap::draw(hm)
dev.off()
