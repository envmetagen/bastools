pathTobowtieScript ="/home/tutorial/SCRIPTS/emg_setup/scripts/emg_fasta_galign"
fasta="example_modecopcroutput.fasta"
out="bowtie2.test.txt"

#bowtie2-build --large-index -f --threads 2 COI.fullgb.ecopcrResults.txt.clean.fasta COI.fullgb.ecopcrResults.txt.clean.fasta.db2
#bowtie2 -a -f --large-index --end-to-end  --sam-no-qname-trunc -p 2 --no-unal -x COI.fullgb.ecopcrResults.txt.clean.fasta.db2 
#-U COI.fullgb.ecopcrResults.txt.clean.fasta --min-score L,0,0.00001 -S COI.fullgb.ecopcrResults.txt.clean.sam

#cat COI.fullgb.ecopcrResults.txt.clean.sam | grep "XM:" | cut -f 1,3-|sed -E "s/([^\t]+)\t([^\t]+)\t.*XM:i:([0-9]+)\t.*XO:i:([0-9]+).*/\1\t\2\t\3\t\4/" | awk 'BEGIN{FS="\t";OFS="\t";} $3<4 { print;} '> COI.fullgb.ecopcrResults.txt.clean.sam2

bowtie.res<-data.table::fread("COI.fullgb.ecopcrResults.txt.clean.sam2",sep = "\t",data.table = F)

ecopcroutput<-"COI.fullgb.ecopcrResults.txt.clean"

add.res.bowtie<-function(ecopcroutput){
  
  ecopcr<-data.table::fread(ecopcroutput,sep = "\t",data.table = F)
  
  ecopcr$seq.name<-paste0(ecopcr$AC, " taxid=", ecopcr$taxid,"; organism=", gsub(" ","_",ecopcr$species_name))
  ecopcr$seq.text<-ecopcr$sequence #only using insert here, not fullseq
  
  phylotools::dat2fasta(ecopcr[,c("seq.name","seq.text")],outfile = paste0(ecopcroutput,".fasta"))
  
  system2(paste0(bastoolsDir,"scripts/emg_fasta_galign_ALL_PERFECT"),args=c("-B","-c","-q",paste0(ecopcroutput,".fasta")
                                                                ,"-r",paste0(ecopcroutput,".fasta"),"-o",paste0(ecopcroutput,".bowtie.res")),wait = T)
  
  bowtie.res<-data.table::fread(paste0(ecopcroutput,".bowtie.res"),data.table = F)
  
  colnames(bowtie.res)<-c("query","subject","mismatches","gaps")
  
  #lca
  bowtie.res$taxids<-do.call(rbind,stringr::str_split(string = bowtie.res$query,pattern = " taxid="))[,2]
  bowtie.res$taxids<-do.call(rbind,stringr::str_split(string = bowtie.res$taxids,pattern = ";"))[,1]
  bowtie.res$AC<-do.call(rbind,stringr::str_split(string = bowtie.res$query,pattern = " "))[,1]
  
  #this is mainly to update any old taxids
  bowtie.res<-add.lineage.df(df = bowtie.res,ncbiTaxDir = ncbiTaxDir)
  
  #lca
  lcasp = aggregate(bowtie.res$taxids, by=list(bowtie.res$query),function(x) ROBITaxonomy::lowest.common.ancestor(obitaxoR,x))
  
  #get lca names
  colnames(lcasp)<-gsub("x","taxids",colnames(lcasp))
  
    if(sum(is.na(lcasp$taxids))>0){
      message("************
                ERROR: Some taxids were not recognized by ROBITaxonomy::lowest.common.ancestor, probably need to update obitaxdb using NCBI2obitaxonomy
                *************")
    }
  
  #get new lineage
  lcasp<-add.lineage.df(df = lcasp,ncbiTaxDir = ncbiTaxDir)
  
  lcasp$path<-paste(lcasp$K,lcasp$P,lcasp$C,lcasp$O,lcasp$F,lcasp$G,lcasp$S,sep = ";")
  
  lcasp$fam.res<-TRUE
  lcasp$gen.res<-TRUE
  lcasp$sp.res<-TRUE
  
  lcasp$fam.res<- !lcasp$path %in% lcasp$path[grep(";unknown;unknown;unknown$",lcasp$path)]
  lcasp$gen.res<- !lcasp$path %in% lcasp$path[grep(";unknown;unknown$",lcasp$path)]
  lcasp$sp.res<- !lcasp$path %in% lcasp$path[grep(";unknown$",lcasp$path)]
  
  lcasp$res<-"htf" #higher than family
  lcasp$res[lcasp$fam.res==T]<-"fam"
  lcasp$res[lcasp$gen.res==T]<-"gen"
  lcasp$res[lcasp$sp.res==T]<-"sp"
}

 