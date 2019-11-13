#merge taxatables

source("/home/bastian.egeter/git_bastools/bastools/merge_MBC_otutab_with_bin_blast.R")
source("/home/bastian.egeter/git_bastools/bastools/googlesheet.foos.R")
source("/home/bastian.egeter/git_bastools/bastools/taxatab.filter.R")
source("/home/bastian.egeter/git_bastools/bastools/plotting.R")

library(ggplot2)

# #taxatabs<-c("/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/2018_07/FILTURB-12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.taxatable.tf.spliced.txt"
#             ,"/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/2018_04/FILTURB-12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.taxatable.tf.spliced.txt"
#             ,"/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/2018_02/FILTURB-12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.taxatable.tf.spliced.txt")


#testng blast against local 12S amphibian database
taxatabs<-c("/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/FILTURB_12S_BLAST/FILTURB-2018_02_12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.taxatable.tf.spliced.txt"
            ,"/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/FILTURB_12S_BLAST/FILTURB-2018_04_12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.taxatable.tf.spliced.txt"
            ,"/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/FILTURB_12S_BLAST/FILTURB-2018_07_12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.taxatable.tf.spliced.txt"
            )


#filturb datasheet
ss_url<-"https://docs.google.com/spreadsheets/d/1FUSaeVaYzms2EOGUoCAB4jaRKzguD3AKTsC8lYwaKP4/edit?ts=5dae01be#gid=1531090624"

##########################################
#organise mastersheet
master_sheet<-google.read.master.url(ss_url)

#check on number of samples in each category
master_xtabs(master_sheet,columns=c("experiment_id","Sample_Type","Primer_set","study"))

#subset master sheet for this study
ms_ss<-subset_mastersheet(master_sheet, experiment_id=c("2018_07"),Primer_set=c("12SV51"),study=c("Both","North"))

#check again to see the subset made sense
master_xtabs(ms_ss,columns=c("experiment_id","Sample_Type","Primer_set","study"))

##########################################
#import taxtabs
all.taxatabs<-bas.merge.taxatables(taxatabs)

#make new taxa table based on subsetted master sheet
all.taxatabs.ss<-cbind(taxon=all.taxatabs$taxon,all.taxatabs[colnames(all.taxatabs) %in% ms_ss$ss_sample_id])

#remove samples/taxa with 0 reads
all.taxatabs.ss<-rm.0readtaxSam(all.taxatabs.ss)

#check negatives
negatives<-negs.stats(taxatab=all.taxatabs.ss,ms_ss = ms_ss,real=c("Field","Tissue"),ex_hominidae=F)

##########################################
#group taxa
all.taxatabs.ss<-bas.group.taxa(taxatab = all.taxatabs.ss,taxon="Eukaryota;Chordata;Amphibia;Caudata;Salamandridae;Salamandra;NA",
                                       jointo="Eukaryota;Chordata;Amphibia;Caudata;Salamandridae;Salamandra;Salamandra salamandra")

all.taxatabs.ss<-bas.group.taxa(taxatab = all.taxatabs.ss,taxon="Eukaryota;Chordata;Amphibia;Anura;Bufonidae;Epidalea;NA",
                                       jointo="Eukaryota;Chordata;Amphibia;Anura;Bufonidae;Epidalea;Epidalea calamita")

all.taxatabs.ss<-bas.group.taxa(taxatab = all.taxatabs.ss,taxon="Eukaryota;Chordata;Amphibia;Anura;Alytidae;Discoglossus;NA",
                                       jointo="Eukaryota;Chordata;Amphibia;Anura;Alytidae;Discoglossus;Discoglossus galganoi")

#Filter samples before removing the NAs and no_hits
all.taxatabs.ss<-taxon.filter.solo.df(taxatab = all.taxatabs.ss,taxonpc = 0.1)
all.taxatabs.ss<-sample.filter.solo(taxatab = all.taxatabs.ss,samplepc=0.1)

#remove the NAs and no_hits
all.taxatabs.ss<-all.taxatabs.ss[all.taxatabs.ss$taxon!="NA;NA;NA;NA;NA;NA;NA",]
all.taxatabs.ss<-all.taxatabs.ss[all.taxatabs.ss$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]

##########################################

#check negatives
negatives<-negs.stats(taxatab=all.taxatabs.ss,ms_ss = ms_ss,real=c("Field","Tissue"),ex_hominidae=F)

#at this point need a few key stats to help decisions
count.dxns.by.taxon(all.taxatabs.ss)
range.dxns.by.taxon(all.taxatabs.ss)
taxatab.stackplot(all.taxatabs.ss)

#check with removing detections (here using 100 becuase it sounds high and doesnt remove any slamander detections)
all.taxatabs.ss.joined50<-filter.dxns(all.taxatabs.ss,filter_dxn = 100)

#check negatives
negatives<-negs.stats(taxatab=all.taxatabs.ss.joined50,ms_ss = ms_ss,real=c("Field","Tissue"),ex_hominidae=F)

count.dxns.by.taxon(all.taxatabs.ss.joined50)
range.dxns.by.taxon(all.taxatabs.ss.joined50)
taxatab.stackplot(all.taxatabs.ss.joined50)

taxtab.pca.plot(all.taxatabs.ss.joined50,master_sheet,"field_method")
#need to fix colours

  

#write
write.taxatab(all.taxatabs.ss.joined50,"all.taxatabs.ss.joined50.test.txt")

#look at read distribution for each taxon

# #remove easy contaminants - NOT NECESSARY FOR AMPHIBIAN BLAST RESULTS
# contaminants<-(c("Suidae","Primates","Felidae","Artiodactyla","Microhylidae","Oryctolagus","Alouatta","Passeriformes","Tarentola boettgeri","Phasianidae","Ardeidae"))
# all.taxatabs.filt.rmCont1<-all.taxatabs.filt
# for(i in 1:length(contaminants)){
#   if(length(grep(contaminants[i],all.taxatabs.filt.rmCont1$taxon))!=0){
#     all.taxatabs.filt.rmCont1<-all.taxatabs.filt.rmCont1[-grep(contaminants[i],all.taxatabs.filt.rmCont1$taxon),]
#   }
# }


  