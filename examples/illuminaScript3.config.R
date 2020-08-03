#############################################################################
#test.config

#first copy this file, do not edit directly

#This config takes as input a fasta (or multiple fastas) and an otu table (or multiple) as produced by MBC 
#All files must be unzipped
#It can process data coming from multiple projects, even if the same primers were used.
#Requirement, the first part of the sample names in the OTU table, the part up to the first dash, 
#must be the project name (i.e. NZFROG-XXX-XX-XX-XX, FILTURB-XX-XX-XX)
#The rest of the sample name is not important
#The otu table and fasta file names must correspond as in the following example: my_sequences.fasta , my_sequences.otutab.tsv
#The otu tables and fasta files must be put in the outDir (as specified below)

#once settings are changed, save the file and from terminal run one line at a time:
#source /home/bastian.egeter/metabinkit.install/metabinkit_env.sh #(or your path to metabinkit)
# nohup Rscript test.config.R  2>run.processing.log 1>&2 & #change log file name if you wish

###################
#SETTINGS
#change settings below as necessary (this directory works on emg1 and emg2)
bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/" #change to your bastools directory, must have trailing "/"
#full path to outDir (containing input fasta and otutabs)
outDir<-"/mnt/Disk1/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/test/" # must have trailing "/"
#full path to ncbi taxonomy directory
ncbiTaxDir<-"/home/bastian.egeter/metabinkit.install/db" #(or your path to metabinkit/db)
#full path to blast exec
blast_exec<-"/home/bastian.egeter/metabinkit.install/bin/blastn" #(or your path to metabinkit/bin/blastn)
#full path to blast database
refdb = "/mnt/Disk1/Tools/BLAST+/DBs/nt_v5/nt"
#full path to kronatools
KronaPath="/home/bastian.egeter/Tools/Krona.install/bin/ktImportText" #(or your path to krona/bin/ktImportText)
#Options for blasting. I recommend keeping the qcov_hsp_perc, gapopen, gapextend, reward, penalty as they are (forces full length alignments)
opts=c("-task","megablast","-outfmt", "6 qseqid pident qcovs saccver staxid ssciname sseq","-num_threads", 64,
       "-max_target_seqs", 50, "-max_hsps",1,"-word_size", 20, "-perc_identity", 70,
       "-qcov_hsp_perc", 98, "-gapopen", 0, "-gapextend", 2, "-reward", 1, "-penalty", -1,"-dust","no")

#starting fastas
startingfastas<-c("12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.fasta",
                  "16S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.fasta")

#set a taxonomic limit to search blast database, one per fasta (here 7742 is vertebrates)
taxidlimit<-c("7742","7742") #To blast entire nt just use taxidlimit=NULL

#which steps do you want to take, see below, remove any you do not want to include
stepstotake<-c("step7","step8","step8a","step9","step10","step11","step12","step13","step14","step15")
#e.g. partial run: stepstotake<-c("step9","step10","step11","step12","step13","step14")

#step7: BLAST
  #output : *blast.txt
#step8: format BLAST results, add taxonomy
  #output : *blast.filt.txt 
#step8a: produce barcode gap analysis report (html) and look for database taxonomic errors
  #there are a few files produced here, not fully detailed here
    #*barcode.gap.report.html
    #*flagged.tsv : the accession numbers of the flagged entries
    #*flagged.details.tsv : details of the flagged entries, for manual checking
    #*non-flagged.potential.errors.tsv : sequences that had >=99% a match to a different taxon
#step9: binning BLAST results (settings below)
  #output : *bins.txt
#step10: merge with otutab
  #output : *taxatable.txt
#step11: apply a % filter by rows 
  #(i.e. remove reads of taxonA from sample1 if number of reads in sample1 are less than the sum of reads for taxonA times 'filterpc')
  #filterpc (0.1=0.1%)
  filterpc=0.1 #to not use this set to 0
  #output : *taxatable.tf.txt
#step12: splice taxa tables by project and produce OTU.per.bin file (number of OTUs in each taxonomic bin)
  #output : *taxatable.tf.spliced.txt
  #output : *otus.per.bin.tsv
#step13: make contributor files
  #Not fully detailed here, this is a file that shows the contributors for each taxonomic bin (to help choose taxa for disabling)
  #output : *taxatable.tf.spliced.contr.txt
#step14: make krona plots of the taxatables
  #output : *taxatable.tf.spliced.krona.html
#step15: make standard heatmaps
  #output : *taxatable.tf.spliced.heatmap.jpg

##BINNING SETTINGS (step9) see metabinkit for details (https://github.com/envmetagen/metabinkit)
use.metabin=T #use metabin or older script for processing? always use.metabin=T
#percentage identity for species level
spident=98
#percentage identity for genus level
gpident=95
#percentage identity for family level
fpident=92
#percentage identity for higher-than-family level
abspident=80
#TOP threshold for each level (1=1%)
topS=2
topG=2
topF=2
topAF=2
#disable specific taxa?
SpeciesBL="species_blacklist_test.txt"
GenusBL="genus_blacklist_test.txt"
FamilyBL=NULL
#exclude specific accessions, if none use known_flags=NULL 
known_flags="/home/bastian.egeter/git_bastools/bastools/known_flags_downloaded_22-7-20.txt" 
#use the flags that were newly identified by step 8a?
use_flagged_accessions_mbk=T 

#Settings for step 8a, only change these settings if you want to run parts of this script alone (havent written user notes on this yet)
steps<-c("selfblast", "find_db_errors", "calc_bar_gap","thresher") #options: "selfblast", "find_db_errors", "calc_bar_gap","thresher"
threshersteps<-c("threshblast","threshbin","threshplots") #"threshblast","threshbin","threshplots"
TaxlevelTestall<-c("K","P","C","O","F") #TaxlevelTestall<-c("K","P","C","O","F","G","S") 
#thresher settings
plot.at.level<-"O" 
limit.plot.to.taxon<-NULL #c("Mammalia","C") #the taxon name and level can be NULL
#binning settings to loop through
tops<-c(0,1,100)
#order=S,G,F,AF
pidents.list<-list(one=c(99,97,95,90),two=c(98,94,92,88),three=c(93,85,75,60)) #can be more than three

source("/home/bastian.egeter/git_bastools/bastools/scripts/illumina_post_MBC_script3.R")