mapTrim<-function(query,buffer,blast.results.file,qc=0.7,out){
  message("reading query file")
  #read in query file and count
  n<-phylotools::read.fasta(query)
  n$qseqid = sub(" .*", "", x = n$seq.name)
  n$seq.text<-as.character(n$seq.text)
  #count
  count_queries=length(n$qseqid)
  
  message("reading mapping results")
  #turn result into table, ignoring comment lines (so ignoring defintions). The table length excludes queries with
  #non-hits
  j<-read.table(file = blast.results.file)
  colnames(j)<-c("qseqid", "qlen", "qstart", "qend",
                 "slen", "sstart", "send", "length", "pident", "qcovs","sstrand")
  count_hits<-length(j$qseqid)
  
  #remove hits less than specified subject cover
  j2<-j[j$length>j$slen*qc,]
  count_qc<-length(j2$qseqid)
  
  #merge query and blast results
  k<-merge(x = j2,y = n,by = "qseqid",all.y = F)
  
  #find query position that matches first subject position (i.e. first base of primer-binding site)
  ##the strand problem...
  kplus<-k[k$sstrand=="plus",]
  kminus<-k[k$sstrand=="minus",]
  kplus$q_start_base<-kplus$qstart-kplus$sstart
  kminus$q_start_base<-kminus$qstart-kminus$send
  kboth<-rbind(kplus,kminus)
  
  #subtract the buffer from this (to only accept queries that extend left of primer)
  kboth$q_left_buff<-kboth$q_start_base-buffer
  #find final query base required
  kboth$q_right_buff<-kboth$q_left_buff+kboth$slen+buffer*2
  
  #remove queries with negative starting base
  k2<-kboth[kboth$q_left_buff>-1,]
  #count
  count_left_buffer<-length(k2$qseqid)
  
  #extract the sequence within thresholds
  k2$ex.seq<-substr(k2$seq.text,start = k2$q_left_buff,stop = k2$q_right_buff)
  #remove queries shorter than expected
  k2$newlen<-nchar(k2$ex.seq)
  #caLculate desired length, with slight leeway (to account for 1-3bp difference in calculation)
  k2$desiredlen<-k2$q_right_buff-k2$q_left_buff-3
  k3<-k2[k2$newlen>k2$desiredlen,]
  #count
  count_right_buffer<-length(k3$qseqid)
  
  message("outputting as fasta")
  k3_export<-k3[,c("seq.name","ex.seq")]
  colnames(k3_export)<-c("seq.name","seq.text")
  count_final_db<-length(k3_export$seq.name)
  phylotools::dat2fasta(k3_export,outfile = out)
  
  ###########################
  message(c("Done. ", "From ", count_queries," sequences, ", count_hits, " mapped to a reference, ",count_qc,
            " of which had > ",qc*100,"% coverage, ",
            count_left_buffer," of which passed left buffer and, of those, ",count_right_buffer," passed right buffer."))
}
