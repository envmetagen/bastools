out="18S.nseqs.with.insert.txt"

#check that sequences which mapped overlap target region  
##read in blast results
a<-data.table::fread(gsub(".fasta",".checked.NoPri.blast",catted_DLS),data.table = F)
colnames(a)<-c("qseqid","sseqid", "stitle", "qlen", "qstart", "qend", "slen", "sstart", "send", "length", "pident", "qcovs", "sstrand")
###note that there are multiple hsps using blast settings above

##correct minus strands, change sstart to send and send to sstart
kplus<-a[a$sstrand=="plus",]
kminus<-a[a$sstrand=="minus",]
kminus$sstart.tmp<-kminus$send
kminus$send<-kminus$sstart
kminus$sstart<-kminus$sstart.tmp
kminus$sstart.tmp=NULL
kboth<-rbind(kplus,kminus)

#extract target start/stop to new cols
kboth$target.start<-as.numeric(gsub(";","",gsub("target.start=","",stringr::str_extract(kboth$stitle,"target.start=(.*?);"))))
kboth$target.end<-as.numeric(gsub(";","",gsub("target.end=","",stringr::str_extract(kboth$stitle,"target.end=(.*?);"))))

#get insert start and end 
kboth$insert.start<-kboth$target.start+buffer+nchar(Pf)
kboth$insert.end<-kboth$target.end-buffer-nchar(Pr)

#check that alignments cover whole insert
kboth$aln.left.of.insert.start<-kboth$sstart<kboth$insert.start+1
kboth$aln.right.of.insert.end<-kboth$send>kboth$insert.end-1
kboth.contains.insert<-kboth[kboth$aln.left.of.insert.start==T & kboth$aln.right.of.insert.end==T,]

#join kboth.contains.insert and gsub(".fasta",".checked.wPos.fasta",catted_DLS) into single fasta
##extract those seqs first
write.table(kboth.contains.insert$qseqid,quote = F,row.names = F,append = F,file = "ref.ids.tmp")
system2(command = "seqkit", args=c("grep","-f","ref.ids.tmp",gsub(".fasta",".checked.fasta",catted_DLS))
        ,stdout = gsub(".fasta",".checked.mapcaught.insert.fasta",catted_DLS),wait = T)
unlink("ref.ids.tmp")

#join 
system2(command = "cat", args=c(gsub(".fasta",".checked.wPos.fasta",catted_DLS),gsub(".fasta",".checked.mapcaught.insert.fasta",catted_DLS))
        ,stdout = gsub(".fasta",".checked.wPos.plus.mapcaught.insert.fasta",catted_DLS),wait = T)

##ntaxa 
taxa.afterchecks<-count.taxa.in.fams.in.fasta(fasta = gsub(".fasta",".checked.wPos.plus.mapcaught.insert.fasta",catted_DLS),ncbiTaxDir)

nseqs.afterchecks<-count.nseqs.in.fams.in.fasta(fasta = gsub(".fasta",".checked.wPos.plus.mapcaught.insert.fasta",catted_DLS),ncbiTaxDir)

write.table(taxa.afterchecks,file = out,append = F,quote = F,sep = "\t",row.names = F)
write.table(nseqs.afterchecks,file = out,append = F,quote = F,sep = "\t",row.names = F)
