#taxatable filtering and report config

#usually will be for one primer/run, but if all settings can be applied globally then multiple tables can be provided

#############################################################################

#change setting below as necessary

#before running for the very first time run this, then hash them out again:
# email="your.name@email.com"
# bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/" #change for your path to bastoolDir
# setwd(bastoolsDir) 
# googlesheets4::sheets_auth(email = email)

#colnames do not match current for hiseq so do
#test<-data.table::fread(taxatabs,data.table = F)
#colnames(test)<-c("taxon",paste(colnames(test[,-1]),"-HSBAS",sep = ""))
#write.table(test,file=taxatabs,append = F,row.names = F,quote = F,sep = "\t")

bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/" #change to your bastools directory

#can be multiple
taxatabs<-c("/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/POST_TAXONOMY_CHECK/CRAY_REBIN/CRAY-HSJUN19BAS_COI.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.rebins.taxatable.tf.spliced.txt")

#datasheet url
ss_url<-"https://docs.google.com/spreadsheets/d/1KZLoXHTgtkD0btSWjyAmFiGJ_cPcYITyfFSlzehisRI/edit#gid=1531090624"

#options for subsetting master sheet. This functions to sleect the samples you want to analyse.
#Each item in list is a column heading in master sheet and 
#each character within the item should be what you want to include (sample_type should always be lower case,
#even if it is not so on google)
subsetlist<-list(experiment_id="HSJUN19BAS",primer_set="LERAY-XT",MPLX="N",sample_type=c("Extraction_Negative","GIT_contents","PCR_negative"),
                 Replicate_Name=c("ZYMO"))

#Detection below taxonpc % of taxon read count will be removed (0.1=0.1%)
taxonpc = 0.1
#Detection below samplepc % of sample read count will be removed (0.1=0.1%)
samplepc=0.1

#the absolute value for removing detections
filter_dxn = 100

#sample_type used to descirbe your real samples (not negatives)
real = "GIT_contents"

#this should be a dataframe with the taxon to group in column 1 and taxon to group to in column 2. 
#Can be made separately, this is just example. Ensure the order is correct
taxa.to.group<-data.frame(taxa.to.group=c("Eukaryota;Chordata;Amphibia;Anura;Alytidae;Discoglossus;NA"
                                           ,"Eukaryota;Chordata;Amphibia;Anura;Bufonidae;Epidalea;NA"),
                           group.to=c("Eukaryota;Chordata;Amphibia;Anura;Alytidae;Discoglossus;Discoglossus galganoi"
                                      ,"Eukaryota;Chordata;Amphibia;Anura;Bufonidae;Epidalea;Epidalea calamita"))

xLevel<-"family"

#variables to use as x axis for barplots
stackplot.vars<-c("Exact Site")

#grouping for removing detections in <1 rep
rep.rm<-"Sample_Name"

#sum reps by
sumrepsby<-"biomaterial"
