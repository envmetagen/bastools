#options for mapping seqeunces to refs

#blast 
#-map refs against seqs using min_qcovs 100%
#-would require to report all hits (as many as are in seqs). But might not be many hits because requires 100% query cover
#-refs should be buffer-primer-insert-primer-buffer
#-would it be worth excluding refs from seqdb?
#pident threshold!?

make.blastdb.bas(infasta = "COI.fullgb.checked.headGenious.fasta",ncbiTaxDir = ncbiTaxDir)

system2(command = "blastn", args=c("-query", "COI.fullgb.refs.headGenious.fasta", "-task", "blastn","-db", "COI.fullgb.checked.headGenious",
                                   "-outfmt",'"6 qseqid sseqid qlen qstart qend slen sstart send length pident qcovs sstrand"',
                                   "-num_threads", "16","-max_target_seqs", "10000","-max_hsps","1","-qcov_hsp_perc","100"),stdout="test1.txt",wait = T)

#blast 
#-map seqs against refs
#-cannot do 100% cover, so not possible...
#-refs should be buffer-primer-insert-primer-buffer
#-would it be worth excluding refs from seqdb?
#pident threshold!?

#bowtie
#map seqs against refs
#-local, same issue as blast regardin 100% cover?
#-accepted alignment = ??
#-would only need best match, if excluding refs from db

#bowtie
#map refs against seqs
#-accepted alignment = ??
#-end to end
#-but maybe have more control over what gets through, based on hits? hard to know

system2(command = "bowtie2-build",args=c("--large-index", "-f", "--threads", "16"
                                         ,"COI.fullgb.checked.over1500.fasta", "COI.fullgb.checked.over1500.fasta.db2"),wait = T)

system2(command = "bowtie2", args = c("-f", "--large-index", "--local","--sam-no-qname-trunc", "-p", "16", "-x"
                                      , "COI.fullgb.checked.over1500.fasta.db2", "-U",  "COI.fullgb.checked.headGenious.fasta"
                                      ,"--very-sensitive-local", "--min-score", "L,1,0", "-S", "COI.fullgb.checked.headGenious.fasta.sam"),wait=T)

a<-data.table::fread(cmd =paste("grep","XM:","COI.fullgb.checked.headGenious.fasta.sam"),data.table = F,sep = "\t")
a<-a[,c(1,3,12,15,16,18,10)] #likely need to change if chaging bowtie2 command
colnames(a)<-c("query","subject","score","mismatches","gaps","editD","sequence")

a$score<-gsub("AS:i:","",a$score)
a$mismatches<-gsub("XM:i:","",a$mismatches)
a$gaps<-gsub("XO:i:","",a$gaps)
a$editD<-gsub("NM:i:","",a$editD)




#--min-score L,0,0.00001  "--no-unal"

#need to come up with test files and known outcomes!



#find longest seqs first to get whole gene...
full.gene<-phylotools::read.fasta("COI.fullgb.checked.over1500.fasta")

make.blastdb.bas(infasta = "COI.fullgb.checked.over1500.fasta",ncbiTaxDir = ncbiTaxDir,dbversion = 4)

system2(command = "blastn", args=c("-query", "COI.fullgb.refs.fasta", "-task", "megablast","-db", "COI.fullgb.checked.over1500",
                                   "-outfmt",'"6 qseqid sseqid qlen qstart qend slen sstart send length pident qcovs sstrand"',
                                   "-num_threads", "16","-max_target_seqs", "1","-subject_besthit","-max_hsps","1","-qcov_hsp_perc","100"),stdout="test1.txt",wait = T)

a<-data.table::fread("test1.txt",data.table = F)

time1<-Sys.time()
system2(command = "blastn", args=c("-query", "COI.fullgb.checked.fasta", "-task", "blastn","-db", "COI.fullgb.checked.over1500",
                                   "-outfmt",'"6 qseqid sseqid qlen qstart qend slen sstart send length pident qcovs sstrand"',
                                   "-num_threads", "16","-max_target_seqs", "1","-max_hsps","1","-subject_besthit"),stdout="test2.txt",wait = T)
time2<-Sys.time()

b<-data.table::fread("test2.txt",data.table = F)
unique(b$V1)
