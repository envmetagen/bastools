#config for new pipeline
######################################
#CONFIG COMMON
ncbiTaxDir="/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/taxdump/October-2019/"
obitaxo<-"/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/obitaxdump/obitaxdb"
obitaxoR<-ROBITaxonomy::read.taxonomy(dbname = obitaxo)

#####################################
#CONFIG 16S
stepstotake<-c("step2a","step3","step4","step5","step6")
outDir<-"/home/bastian.egeter/FILIPA_NEW/bas_nuccore_testing/16S/"
gene="16S"
groups=c("Clitellata","Insecta","Malacostraca","Mollusca","Platyhelminthes","Trombidiformes")
group.taxids=c("42113","50557","6681","6447","6157","83136")
group.ranks=c("class","class","class","phylum","phylum","order")
Pf="RGACGAGAAGACCCTATARA"
Pr="ACGCTGTTATCCCTAARGTA"
max_error_buildrefs=4
max_error_ecopcr=7
min_length=125 #insert
max_length=173 #insert
Ta=52
buffer=10
out_bias_file = "16S.test.bias.txt"
out_mod_ecopcrout_file = "16S.test.mod_ecopcrout.txt"
stepcountfile<-"16S.stepcounts.txt"


source("/home/bastian.egeter/git_bastools/bastools/primer_script_oct_2019.R")