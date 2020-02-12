#config for new pipeline
######################################
#CONFIG COMMON
ncbiTaxDir="/home/tutorial/TOOLS/DBS/ncbi_taxonomy/taxdump/"
obitaxo<-"/home/tutorial/TOOLS/DBS/ncbi_taxonomy/Jan2020/obitaxdb_jan_2020"
makeblastdb_exec = "/home/tutorial/TOOLS/ncbi-blast-2.9.0+/bin/makeblastdb"
blast_exec = "/home/tutorial/TOOLS/ncbi-blast-2.9.0+/bin/blastn"
bastoolsDir = "/home/tutorial/TOOLS/bastools/"

#####################################
#CONFIG 18S Primer
#stepstotake<-c("step3","step4","step5","step6","step7","step8","step9","step10","step11")
stepstotake<-c("step8a","step9","step10","step11") #testing with upper limit of the miseq run used in the final ecopcr (step 6) only
outDir<-"/media/sf_Documents/WORK/CIBIO/temp/primer_testing/"

Pf="GGNTGAACHGTHTAYCCHCC"
Pr="TCDGGRTGNCCRAARAAYCA"
max_error_buildrefs=2
max_error_ecopcr=4
min_length=50 #insert
max_length=500-nchar(Pf)-nchar(Pr) #testing with upper limit of the miseq run used in the final ecopcr 
long_length=max_length+400 #used to count how many seqs are being missed by first ecopcr due to length rather than mismatches
Ta=50
buffer=5
#top=1 #100-top = the percentage identity to consider when calculating resolution (bas2)
out_bias_file = "example_bias.txt"  
out_mod_ecopcrout_file = "example_modecopcroutput.txt"
catted_DLS<-"COI.fullgb.fasta" #this is the starting fasta
script<-paste0(bastoolsDir,"scripts/primer_script_FULLGB_Feb_2020.R")

#do not modify fasta, skip coverage steps
simplepipe<-T

source(script)