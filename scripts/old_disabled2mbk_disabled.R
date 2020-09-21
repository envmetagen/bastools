#convert old disabled file to new format for mbk

orig<-data.table::fread("temp/Disabled_taxa_All_Table_March.txt",data.table = F)
ncbiTaxDir<-"/home/tutorial/TOOLS/metabinkit.install/db/"
disabled.sp.out<-"/home/tutorial/temp/disabled.sp.iran.txt"
disabled.g.out<-"/home/tutorial/temp/disabled.g.iran.txt"
disabled.f.out<-"/home/tutorial/temp/disabled.f.iran.txt"


##SPECIES
origS<-orig[,c("taxids","disable_species")]
#remove NA and FALSE
origS<-origS[!is.na(origS$disable_species),]
origS<-origS[!origS$disable_species==F,]
#remove dups
origS<-origS[!duplicated(origS$taxids),]
#done, save file
write.table(origS$taxids,file = disabled.sp.out,quote = F,row.names = F,col.names = F)

##GENUS
origG<-orig[,c("taxids","disable_genus")]
#remove NA and FALSE
origG<-origG[!is.na(origG$disable_genus),]
origG<-origG[!origG$disable_genus==F,]
#remove dups
origG<-origG[!duplicated(origG$taxids),]
#add taxid lineage
origG<-add.lineage.df(origG,ncbiTaxDir = ncbiTaxDir,as.taxids = T)
origG<-origG[,c("G"),drop=F]
#remove dups
origG<-origG[!duplicated(origG$G),,drop=F]
#done, save file
write.table(origG$G,file = disabled.g.out,quote = F,row.names = F,col.names = F)

##FAMILY
origF<-orig[,c("taxids","disable_family")]
#remove NA and FALSE
origF<-origF[!is.na(origF$disable_family),]
origF<-origF[!origF$disable_family==F,]
#remove dups
origF<-origF[!duplicated(origF$taxids),]
#add taxid lineage
origF<-add.lineage.df(origF,ncbiTaxDir = ncbiTaxDir,as.taxids = T)
origF<-origF[,c("F"),drop=F]
#remove dups
origF<-origF[!duplicated(origF$F),,drop=F]
#done, save file
write.table(origF$F,file = disabled.f.out,quote = F,row.names = F,col.names = F)

