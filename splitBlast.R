splitBlast<-function(query,out,refdb){
message("Reminders: must not be any files in the folder with '*.0[0-9].fasta' in filename")
#split first
f<-process$new(command = "pyfasta", args=c("split","-n","10",query), echo_cmd = T)
f$wait()
#blast split files
a<-list.files(pattern = "*.0[0-9].fasta")
h<-list()
for(i in 1:length(a)){
  h[[i]]<-process$new(command = "blastn", args=c("-query", a[i], "-task", "blastn","-db",refdb,"-outfmt",
                                                 "7 qseqid qlen qstart qend slen sstart send length pident qcovs sstrand",
                                                 "-num_threads", "16", "-subject_besthit", "-max_hsps", "1","-max_target_seqs", "1"),
                      echo_cmd = T,stdout=paste0(gsub(x = a[i],pattern = "fasta",replacement = "blastTemp.txt")))
}
for(i in 1:length(h)){
  h[[i]]$wait()
}

#concatenate blast results
b<-list.files(pattern = "*blastTemp.txt")
g<-process$new(command = "cat", args=c(b), echo_cmd=T,stdout=out)
g$wait()

#remove intermediate files
d<-list.files(pattern = "*flat")
e<-list.files(pattern = "*gdx")
unlink(c(a,b,d,e))
}