

bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/" 

###############

source(paste0(bastoolsDir,"building_dbs.R"))
source(paste0(bastoolsDir,"manip_fasta.R"))
setwd(paste0(bastoolsDir,"tests"))

#I am parsing genbank files to extract the gene of interest
#this is done by first splitting the genbank file by record (no issues there)
#then, matching a query string within the primary field "rRNA" and subfield "product"

# the issue I have is that I like to search for a couple of different strings
# e.g. "18S*" (the default) OR 

# view example files
readLines("Xenidae_41039_18S.gb") #works, see below
readLines("Xyronotidae_62795_18S.gb") #partially fails, see below

# currently, for 18S gene, I look for the string "18S"

# this works (emg1)
bascount.gb("Xenidae_41039_18S.gb")
extract.gene.gb(gbfile = "Xenidae_41039_18S.gb",gene="18S")
bascount.fasta("Xenidae_41039_18S.extract.fasta")

# here one out of the two records fails, because "/product="small subunit ribosomal RNA"" does not contain "18S"    
bascount.gb("Xyronotidae_62795_18S.gb")
extract.gene.gb(gbfile = "Xyronotidae_62795_18S.gb",gene="18S")
bascount.fasta("Xyronotidae_62795_18S.extract.fasta")

###NOTES
# the function calls a script "parse-genbank-18S.py" (which Ive pasted below). 

# I have also pasted the primary function here (you can see I had a couple of workarounds, but
#these required having multiple py scripts)

# On a separate issue, regardless of this I have a separate script for each gene, which seems wasteful
# would be nice if I could pass these from R


#parse-genbank.py modified from https://github.com/adina/scripts-for-ngs/blob/master/parse-genbank.py

import sys
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


genome=SeqIO.read(sys.argv[1], 'genbank')


l = []
n = 0
for record in list(SeqIO.parse(sys.argv[1], 'genbank')):
  name = record.name
for feature in record.features:
  if 'source' in feature.type:
  taxid=''.join(feature.qualifiers['db_xref'])
taxid=re.sub(r'.*taxon:','',taxid)
org=''.join(feature.qualifiers['organism'])
org=re.sub(r'.*organism=','',org)
for feat in genome.features:
  if feat.type == "rRNA":
  if '18S' in feat.qualifiers['product'][0]:      
  
  
  # here is where I need an OR statement: "18S" OR "small subunit ribosomal RNA" 
      # ideally these should be a character vector in R and fed to this script. 
      #There can be up to 6 variations to search for
  # furthermore, I would like to have an OR statement for feat.qualifiers: ['product'] OR ['note']
  
  
  start = feat.location.start.position
end = feat.location.end.position
pos = [start, end]
l.append(pos)
print '>' + name + ' organism=' + org + '; taxid=' + taxid + ";"

print feat.extract(genome.seq)
n = n + 1


#primary function
extract.gene.gb<-function(gbfile,gene){
  #split gb file by record
  system2(command = "cat", args=c(gbfile, "|", "sh",paste0(bastoolsDir,"split_gb.sh")),
          wait=T) 
  #remove last file cause its always empty
  files<-list.files(pattern = "^outTemp.*")
  a<-suppressWarnings(max(as.numeric(do.call(rbind,stringr::str_split(files,"outTemp"))[,2])))
  unlink(paste0("outTemp",a))
  
  #extract gene for each file
  for(i in 1:length(list.files(pattern = "^outTemp.*"))){
    a<-list.files(pattern = "^outTemp.*")[i]
    if(gene=="18S") script<-paste0(bastoolsDir,"parse-genbank-18S.py")
    if(gene=="16S") script<-paste0(bastoolsDir,"parse-genbank-16S.py")
    if(gene=="COI") script<-paste0(bastoolsDir,"parse-genbank-COI.py")
    system2("python",args = c(script, a),wait = T,
            stdout = gsub("outTemp","extract.outTemp",a),stderr = F)
    
    # #some files have /note instead of /product
    # count<-system2("wc",args = c("-l",gsub("outTemp","extract.outTemp",a)),wait = T,stdout = T)
    # if(as.numeric(do.call(rbind,stringr::str_split(count," "))[,1])==0){
    #   system2("python",args = c(gsub(".py","_note.py",script), a),wait = T,
    #           stdout = gsub("outTemp","extract.outTemp",a),stderr = F)
    # }
    # #some files have "small subunit ribosomal RNA" in /product (no "16S")
    # count<-system2("wc",args = c("-l",gsub("outTemp","extract.outTemp",a)),wait = T,stdout = T)
    # if(as.numeric(do.call(rbind,stringr::str_split(count," "))[,1])==0){
    #   system2("python",args = c(gsub(".py","_ssrrna.py",script), a),wait = T,
    #           stdout = gsub("outTemp","extract.outTemp",a),stderr = F)
    # }
    
    #some files fail for other reasons (usually not having field rRNA, dont mind missing these)
    count<-system2("wc",args = c("-l",gsub("outTemp","extract.outTemp",a)),wait = T,stdout = T)
    if(as.numeric(do.call(rbind,stringr::str_split(count," "))[,1])==0){
      message(paste("One record from",gbfile,"failed and was excluded"))
      unlink(gsub("outTemp","extract.outTemp",a))
    }
  }
  
  #cat files
  if(length(list.files(pattern = "^extract.outTemp.*"))!=0){
    
    system2("cat",args=c(list.files(pattern = "^extract.outTemp.*")),
            stdout = gsub(".gb",".extract.fasta",gbfile),wait = T)
    
    #remove extraneuos files
    unlink(list.files(pattern = "^extract.outTemp.*"))
    unlink(list.files(pattern = "^outTemp.*"))
  }
}



