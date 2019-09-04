#PIPELINE/CONFIG/SETUP
library(PrimerMiner)
library(taxize)
library(processx)
library(bastools)
library(bold)
library(dplyr)
setwd("/media/sf_Documents/WORK/CIBIO/STATS_AND_CODE/bastools_dir/bastools/R/")
file.sources<-list.files()
sapply(file.sources,source,.GlobalEnv)

#start in fresh, empty folder
setwd("/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/Primer_tests/bastools_results/Filipa_1st_test/COI_ALL_GROUPS/testfamsuniq/")

groups=c("Clitellata","Insecta","Malacostraca","Mollusca","Nematoda","Ostracoda","Platyhelminthes",
        "Trombidiformes")

#weird, Leiodidae causing bold_seq to crash. no idea why!

target_gene=c("COI","COX1")
obitaxo<-system.file("extdata", "taxa.csv", package = "bastools")
obitaxo<-gsub(pattern = "taxa.csv",replacement = "obitax_26-4-19",x = obitaxo)
Pf="GGNTGAACHGTHTAYCCHCC"
Pr="TCDGGRTGNCCRAARAAYCA"
max_error_finding_references=4
max_error_final_ecopcr=6
min_length=305 #insert
max_length=323 #insert
buffer=10
BOLDTF = T
genbankTF = F
mitoTF = F

#########CHECK THAT BOLD DOWNLOAD IS NOT GETTING OVERWRITTEN FOR MULTIPLE GROUPS

############################
#CREATE DATABASE
#download sequences
build.db1(groups = groups2,target_gene = target_gene,BOLD = BOLDTF,genbank = genbankTF,mito = mitoTF)
#had to restart a couple of times due to connection breaking, but completed ok
count_build1<-bascount.fastas.recur() #run in directory with all seqs - probably overcounts gb files,
#only perform this on fastas
#format species names
build.db2(out = "build2.fasta")
count_build2<-bascount.fasta("build2.fasta")
#to save time later, remove all seqs less than minimum required length
l<-sum(min_length,nchar(Pf),nchar(Pr),buffer,buffer)
f<-process$new(command = "obigrep", args = c("-l",l,"build2.fasta"), echo_cmd = T,stdout="build2A.fasta")
f$wait()
count_build2A<-bascount.fasta("build2A.fasta")
#add taxids, only sequences with a species level taxid are kept (genus level goes to dump_genus.fasta)

###############this has major error, sequences such as genus_sp. or family_sp. pass into the next file! need to fix
obiaddtaxids.Bas(infile = "build2A.fasta",taxo = obitaxo,k = "species",out = "build2B.fasta")
count_build2B<-bascount.fasta("build2B.fasta")
# note, could add taxids and uniq-id to build.db2
f<-process$new(command = "obiannotate", args=c("--uniq-id","build2B.fasta"), echo_cmd = T,stdout="build2C.fasta")
f$wait()
count_build2C<-bascount.fasta("build2C.fasta")
#######################
#DEREPLICATE AT THIS POINT?
f<-process$new(command = "obiuniq", args=c("-m","taxid","-i","build2C.fasta"), echo_cmd = T,stdout="build2C.derep.fasta")
f$wait()
count_build2C.derep<-bascount.fasta("build2C.derep.fasta")


############################
#CREATE TARGET REFERENCE SEQUENCES
####update to match derepped file

#make ecopcrdb
obiconvert.Bas(infile = "build2C.derep.fasta",in_type = "fasta",out = "build2C.derep.ecopcr",taxo = obitaxo,
               out_type = "--ecopcrdb-output" )
#run ecopcr
build2C.ecopcr.Results<-ecoPCR.Bas(Pf = Pf,Pr = Pr,ecopcrdb = "build2C.derep.ecopcr",
                                   max_error = max_error_finding_references,
                                   min_length = min_length,max_length = max_length)
#build reference sequence set (select one amplicon, including pbs, per genus)
build.refs(input.ecopcr.results = build2C.ecopcr.Results,output = "build2D.fasta")
count_build2D<-bascount.fasta("build2D.fasta")
#############################
#map seqs to refs to make final db
#split first
f<-process$new(command = "pyfasta", args=c("split","-n","10","build2C.derep.fasta"), echo_cmd = T)
f$wait()

#make refdb
g<-processx::run(command = "makeblastdb", args=c("-in", "build.2D.fasta", "-dbtype", "nucl",
                "-parse_seqids","-out","refdb"), echo=F,echo_cmd = T)
g$wait()

#blast split files
a<-list.files(pattern = "*derep.0")
h<-list()
for(i in 1:length(a)){
h[[i]]<-process$new(command = "blastn", args=c("-query", a[i], "-task", "megablast","-db","refdb",
                                          "-outfmt",
                                          "7 qseqid qlen qstart qend slen sstart send length pident qcovs sstrand",
                                          "-num_threads", "16", "-subject_besthit", "-max_hsps", "1","-max_target_seqs", "1"),
               echo_cmd = T,stdout=paste0(gsub(x = a[i],pattern = "fasta",replacement = "blast.txt")))
}
for(i in 1:length(h)){
  h[[i]]$wait()
}

#concatenate blast results
g<-process$new(command = "cat", args=c("*blast.txt"), echo_cmd=T,stdout="x.txt")
g$wait()

build.db4.uniqFam(query = "../build2C.derep.fasta",buffer = 10,blast.results.file = "../x.txt",out.final.db = "build4.fasta")

###derep by family
#DEREPLICATE AT THIS POINT?
f<-process$new(command = "obiuniq", args=c("-m","family","-i","build4.1.fasta"), echo_cmd = T,stdout="build4.1.derep.fasta")
f$wait()
count_build4.derep<-bascount.fasta("build4.1.derep.fasta")

#count_build4<-bascount.fasta("build4.1.fasta")
###########################
##create final ecopcrdb
obiconvert.Bas(infile = "build4.1.derep.fasta",in_type = "fasta",out = "build4.1.derep.ecopcr",taxo = obitaxo,
               out_type = "--ecopcrdb-output")
############################
message(c("downloaded sequences=",count_build1,", formatted=",count_build2,", passing required length filter=",
count_build2A,", with species level taxonomy=",count_build2B,
", mapped to ref and passed buffers and were unique species=",count_build4))
#Move to final.tables.R script
