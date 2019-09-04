build.db3<-function(query,refs,blast.results.file){
  #use blast to align seqs to ref
  message("mapping sequences to reference")
  g<-processx::run(command = "makeblastdb", args=c("-in", refs, "-dbtype", "nucl", "-parse_seqids","-out","refdb"),
                   echo=F,echo_cmd = T)


  h<-process$new(command = "blastn", args=c("-query", query, "-task", "blastn","-db","refdb",
                "-outfmt",
                "7 qseqid qlen qstart qend slen sstart send length pident qcovs sstrand",
                "-num_threads", "16", "-subject_besthit", "-max_hsps", "1","-max_target_seqs", "1"),
                 echo_cmd = T,stdout=blast.results.file)
  h$wait()
  #just to report no. of hits
  j<-read.table(file = blast.results.file)
  colnames(j)<-c("qseqid", "qlen", "qstart", "qend",
                 "slen", "sstart", "send", "length", "pident", "qcovs","sstrand")
  count_hits<-length(j$qseqid)
  message(c(count_hits, "hits"))
  d<-"SUCCESS"
}

build.db4<-function(query,buffer,blast.results.file,out.final.db,obitaxo){
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

  #remove hits less than 70% subject cover (#easier yo use megablast algorithm?)
  j2<-j[j$length>j$slen*0.7,]
  count_70<-length(j2$qseqid)

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

  #select one seq per species
  k3$seq.name<-as.character(k3$seq.name)
  a<-stringr::str_split(string = k3$seq.name,pattern = "; ")
  b<-sapply(a, function(x){x[3]})
  b<-gsub(pattern = "species=",replacement = "",x = b)
  k3$species<-b
  k4<-k3[!duplicated(k3$species),]
  #count
  count_final_db<-length(k4$qseqid)

  message("outputting as fasta")
  #combine nominal accession, species tag, taxid tag and defintion into new defintion and output fasta
  #a<-stringr::str_split(string = k4$seq.name,pattern = "; ")
  #b<-sapply(a, function(x){x[5]})
  #b<-gsub(pattern = "merged_taxid=\\{",replacement = "",x = b)
  #b<-gsub(pattern = ": *.}",replacement = "",x = b)
  #k4$taxid<-b
  #k4$newdef<-paste0(k4$qseqid, " species=",k4$species,"; taxid=",k4$taxid,";")
  k4_export<-k4[,c("seq.name","ex.seq")]
  invisible(seqRFLP::dataframe2fas(k4_export,file = gsub(x=out.final.db,pattern = ".fasta",replacement = ".1.fasta")))
  #####why is obiconvert to ecopcr failing at next step?####
  #no idea but if I run obiaddtaxids first then it seems to work
  obiaddtaxids.Bas(infile = gsub(x=out.final.db,pattern = ".fasta",replacement = ".1.fasta"),
                   taxo = obitaxo,k = "species",out = out.final.db)

  ###########################
  message(c("Done. ", "From ", count_queries," sequences, ", count_hits, " mapped to a reference, ",count_70,
            " of which had > 70% coverage, ",
            count_left_buffer," of which passed left buffer and ",count_right_buffer," passing right buffer, corresponding to ",
            count_final_db," unique species."))
}

build.db4.uniqFam<-function(query,buffer,blast.results.file,out.final.db,obitaxo){
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
                 "slen", "sstart", "send", "length", "pident", "qcovs","qseq","sseq","sstrand")
  count_hits<-length(j$qseqid)

  #remove hits less than 70% subject cover
  j2<-j[j$length>j$slen*0.7,]
  count_70<-length(j2$qseqid)

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

  #select one uniq seq per family
  k3$seq.name<-as.character(k3$seq.name)
  k3$taxid<-stringr::str_match(k3$seq.name, "merged_taxid=\\{(.*?):")[,2]
  d<-bastools::get.classic.taxonomy.Bas(k3,obitaxdb = obitaxoR)
  k3$family<-d$family_name_ok
  k4<-k3[!duplicated(k3[,c("family","ex.seq")]),]

  #make new name
  k4$new.seq.name<-paste0(k4$qseqid," ",k4$family,"_sp.;")
  count_final_db<-length(k4$qseqid)

  message("outputting as fasta")
  k4_export<-k4[,c("new.seq.name","ex.seq")]
  colnames(k4_export)<-c("seq.name","seq.text")
  phylotools::dat2fasta(k4_export,outfile = gsub(x=out.final.db,pattern = ".fasta",replacement = ".1.fasta"))

  #####why is obiconvert to ecopcr failing at next step?####
  #no idea but if I run obiaddtaxids first then it seems to work
  #obiaddtaxids.nodump.Bas(infile = gsub(x=out.final.db,pattern = ".fasta",replacement = ".1.fasta"),
     #              taxo = obitaxo,k = "family",out = out.final.db)

  ###########################
  message(c("Done. ", "From ", count_queries," sequences, ", count_hits, " mapped to a reference, ",count_70,
            " of which had > 70% coverage, ",
            count_left_buffer," of which passed left buffer and ",count_right_buffer," passing right buffer, corresponding to ",
            count_final_db," unique species."))
}
