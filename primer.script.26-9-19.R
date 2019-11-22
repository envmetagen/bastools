###########################################################################################
#remove dashes
outfasta = gsub(".fasta",".fixed.fasta",starting_file)

remove.dashes.fasta(infasta = starting_file,outfasta = outfasta)
count.download<-bascount.fasta(outfasta)

stepcounts<-as.data.frame(count.download)

#to save time later, remove all seqs less than minimum required length
f<-process$new(command = "obigrep", args = c("-l",sum(min_length,nchar(Pf),nchar(Pr),buffer,buffer),
                                             outfasta), echo_cmd = T,
               stdout=gsub(".fixed.fasta","formatted.minL.fasta",outfasta))
f$wait()
stepcounts$count.after.rm.minL<-bascount.fasta(gsub(".fixed.fasta","formatted.minL.fasta",outfasta))

#add lineage
add.lineage.fasta.BAS(infasta=gsub(".fixed.fasta","formatted.minL.fasta",outfasta),
                      ncbiTaxDir=ncbiTaxDir,out="formatted.minL.lineage.fasta")
bascount.fasta("formatted.minL.lineage.fasta")

#remove family unknown
fastatemp<-phylotools::read.fasta("formatted.minL.lineage.fasta")
if(length(grep("family=unknown",fastatemp$seq.name))>0){
  fastatemp<-fastatemp[-grep("family=unknown",fastatemp$seq.name),]}
phylotools::dat2fasta(dat = fastatemp,outfile = "formatted.minL.lineage.goodfam.fasta")


############actually should derep within taxid, as I was doing before

stepcounts$count.after.rmNoFams<-bascount.fasta("formatted.minL.lineage.goodfam.fasta")
stepcounts$count.familes.after.rmNoFams<-length(unique(stringr::str_match(fastatemp$seq.name, "family=(.*?);")[,2]))

#ensure uniq ids
system2(command = "obiannotate", args=c("--uniq-id","formatted.minL.lineage.goodfam.fasta"), 
        stdout="formatted.minL.lineage.goodfam.uid.fasta",wait=T)
bascount.fasta("formatted.minL.lineage.goodfam.uid.fasta")


######################should add this chunk to obiconvert.BAS (with option for adding to existing)
#remove "|"s to stop ecopcr errors
f<-process$new(command = "sed", args = c("s/|/_/g",
                                         "formatted.minL.lineage.goodfam.uid.fasta"), echo_cmd = T,
               stdout="formatted.minL.lineage.goodfam.uid2.fasta")
f$wait()
bascount.fasta("formatted.minL.lineage.goodfam.uid2.fasta")

#convert to ecopcrdb
unlink(x = "formatted.minL.lineage.goodfam.uid2.ecopcrdb*")
####################################

obiconvert.Bas(infile = "formatted.minL.lineage.goodfam.uid2.fasta",in_type = "fasta",
               out = "formatted.minL.lineage.goodfam.uid2.ecopcrdb",taxo = obitaxo,
               out_type = "--ecopcrdb-output" )

#run ecopcr with buffer option
ecoPCR.Bas(Pf,Pr,ecopcrdb = "formatted.minL.lineage.goodfam.uid2.ecopcrdb",max_error = max_error_buildrefs,
           min_length,max_length,out = "all.ecopcr.hits.txt",  buffer = buffer)

#build refs fasta from results
ecopcr2refs2(ecopcrfile="all.ecopcr.hits.txt",outfile = "mapping.reference.fasta",bufferecopcr = buffer)
stepcounts$mapping.reference<-bascount.fasta("mapping.reference.fasta")

###########################################################################################
#MAP ALL SEQUENCES AGAINST TARGET REFERENCE DATABASE
map2targets(queries.to.map = "formatted.minL.lineage.goodfam.uid2.fasta",
            refs = "mapping.reference.fasta",out = "formatted.minL.lineage.goodfam.uid2.mapped.txt")

#trim
mapTrim2.simple(query = "formatted.minL.lineage.goodfam.uid2.fasta",
         blast.results.file = "formatted.minL.lineage.goodfam.uid2.mapped.txt",
         out = "formatted.minL.lineage.goodfam.uid2.mapped.fasta",qc=0.97)

stepcounts$count.aftermap<-bascount.fasta("formatted.minL.lineage.goodfam.uid2.mapped.fasta")
#count.families
tempfasta<-phylotools::read.fasta("formatted.minL.lineage.goodfam.uid2.mapped.fasta")
stepcounts$count.familes.after.map<-length(unique(stringr::str_match(tempfasta$seq.name, "family=(.*?);")[,2]))

#create final ecopcrdb
unlink(x = "final.ecopcrdb*")
obiconvert.Bas(infile = "formatted.minL.lineage.goodfam.uid2.mapped.fasta",
               in_type = "fasta",out = "final.ecopcrdb",taxo = obitaxo,
               out_type = "--ecopcrdb-output")
###########################################################################################
#run final ecopcr without buffer option
ecoPCR.Bas(Pf,Pr,ecopcrdb = "final.ecopcrdb",max_error = max_error_ecopcr,
           min_length,max_length,out = "final.ecopcr.hits.txt")
###########################################################################################
#convert final ecopcrdb
system2(command = "obitab", args=c("-o","final.ecopcrdb"), stdout="final.ecopcrdb.tab", wait = T)
###########################################################################################
#DO STATS
make.primer.bias.tables(originaldbtab = "final.ecopcrdb.tab",ecopcrfile = "final.ecopcr.hits.txt",
                        out_bias_file = out_bias_file,
                        out_mod_ecopcrout_file = out_mod_ecopcrout_file,Pf = Pf, Pr = Pr,
                        obitaxoR = obitaxoR,min_length = min_length,max_length = max_length)

biastemp<-data.table::fread(out_bias_file,sep = "\t")

stepcounts$count.amped<-sum(biastemp$nseqs.amped,na.rm = T)
stepcounts$families.amped<-sum(biastemp$amplified)
stepcounts$count.uniq.brcds.amped<-sum(biastemp$n.uniq.brcds.amped,na.rm = T)

write.table(stepcounts,stepcountfile,quote = F,sep = "\t",row.names = F)

