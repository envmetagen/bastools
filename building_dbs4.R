map2targets<-function(queries.to.map,refs,out){
  #use blast to align seqs to ref
  message("mapping sequences to reference")
  system2(command = "makeblastdb", args=c("-in", refs, "-dbtype", "nucl", "-parse_seqids","-out","refdb"),wait=T)

  system2(command = "blastn", args=c("-query", queries.to.map, "-task", "megablast","-db","refdb",
                "-outfmt",'"7 qseqid qlen qstart qend slen sstart send length pident qcovs sstrand"',
                "-num_threads", "16","-max_target_seqs", "1"),stdout=out,wait = T)
}

