blast.min.bas<-function(infastas,refdb,blast_exec="blastn",wait=F){
  for(i in 1:length(infastas)){
  system2(command = blast_exec, args=c("-query", infastas[i],
          "-task", "megablast","-db",refdb,"-outfmt", '"6 qseqid evalue staxid pident qcovs"', 
          "-num_threads", "16", "-max_target_seqs", "100", "-max_hsps","1", "-out",
          paste0(gsub(x = infastas[i],pattern = ".fasta",replacement = ".blast.txt"))),stdout = NULL,wait = wait)
  }
}

