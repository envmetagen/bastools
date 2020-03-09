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
#stepstotake<-c("step1","step2","step3","step4","step5","step6","step7","step8")
stepstotake<-c("step1","step2","step3","step4","step5","step6","step7") 
outDir<-"/media/sf_Documents/WORK/CIBIO/temp/primer_testing/"

Pf="TTTGTCTGSTTAATTSCG"
Pr="CACAGACCTGTTATTGC"
max_error_buildrefs=3
max_error_ecopcr=5
min_length=50 #insert
max_length=300-nchar(Pf)-nchar(Pr) #testing with upper limit of the miseq run used in the final ecopcr 
Ta=50
buffer=5

out_bias_file = "18S.s2.metazoa.bias.txt"  
out_mod_ecopcrout_file = "18S.s2.metazoa.modecopcroutput.txt"
catted_DLS<-"18S.s2.metazoa.fasta" #this is the starting fasta...usually made from nunos fasta
cattedDLS.s2<-"18S.s2.tsv" #this is nunos table
taxonlimit<-"Metazoa" #this is an inital limited of the table, helps a lot to reduce size, set to NULL to ignore
script<-paste0(bastoolsDir,"scripts/primer_script_FULLGB_Feb_2020.R")

source(script)