#workaround for getting filenames of fastas to blast
google.get.startingfastas<-function(outDir,sheeturls,experiment_id,usingobiuniq){
  experimentsheet<-google.make.experiment.sheet(outDir,sheeturls,experiment_id) #in process, writes a sheet to file 
  
  #get barcodes used  
  barcodes.used<-unique(experimentsheet$barcode_id)
  barcodes.used <- barcodes.used[!is.na(barcodes.used)]
  barcodes.used<-gsub("BC","barcode",barcodes.used)
  
  #size select, for each fragment, I checked and seqs appear to have primers plus one base (at each end)
  experimentsheet$primer_combo<-paste0(experimentsheet$Primer_F,experimentsheet$Primer_R)
  primer_combo<-unique(experimentsheet$primer_combo)
  
  #barcodes in each primer combo
  primer_combo.bcs<-list()
  for(i in 1:length(primer_combo)){
    primer_combo.bcs[[i]]<-unique(experimentsheet[experimentsheet$primer_combo==primer_combo[i],"barcode_id"])
    primer_combo.bcs[[i]]<-gsub("BC","barcode",primer_combo.bcs[[i]])
    names(primer_combo.bcs[[i]])<-gsub(" ","",
                                       experimentsheet[experimentsheet$primer_combo==primer_combo[i],"Primer_set"][1])
  }
  if(usingobiuniq){
    for(i in 1:length(primer_combo.bcs)){
      print(paste0(experiment_id,"_",names(primer_combo.bcs[[i]][1]),".uniq.filtlen.wlen.obi.fasta"))
    }}
  
  if(usingobiuniq==F){
    for(i in 1:length(primer_combo.bcs)){
      print(paste0(experiment_id,"_",names(primer_combo.bcs[[i]][1]),".filtlen.wlen.obi.fasta"))
    }}
  
}