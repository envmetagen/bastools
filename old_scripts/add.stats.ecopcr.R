
#make fasta from ecopcr inserts 

mod_ecopcrout = "COI.modecopcroutput.LIMIT50-500.FULLGB.txt"
makeblastdb_exec = "/home/bastian.egeter/Tools/ncbi-blast-2.9.0+/bin/makeblastdb"
ncbiTaxDir = "/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/taxdump/Jan2020/"
blast_exec = "/home/bastian.egeter/Tools/ncbi-blast-2.9.0+/bin/blastn"
obitaxdb = "/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/obitaxdump/Jan2020/obitaxdb_jan_2020"

source("/home/bastian.egeter/git_bastools/bastools/bin.blast.R")

add.res.bas2("COI.modecopcroutput.LIMIT50-500.FULLGB.txt",makeblastdb_exec,ncbiTaxDir,blast_exec,obitaxdb)

#add stats to ecopcroutput


if(add.res.ecopcroutput(add.res.bas2(mod_ecopcrout)))