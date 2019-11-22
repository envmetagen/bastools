

count.MBC<-function(MBCtsvDir,ms_ss,otutab,illumina_script_taxatab,illumina_script_taxatab_tf){
##########################################
  message("Currently set up for one run, one primer, modify later if needed")
#import counts step1
step1<-data.table::fread(paste0(MBCtsvDir,"step1_stats.tsv"),data.table = F)
step1$ss_sample_id<-gsub(".none.*$","",step1$Stats)
demuliplexed_files<-sum(step1[step1$ss_sample_id %in% ms_ss$ss_sample_id,"Number of reads"])/2

#import counts step2
step2<-data.table::fread(paste0(MBCtsvDir,"step2_stats.tsv"),data.table = F)
step2$ss_sample_id<-gsub(".none.*$","",step2$Stats)
after_paired_end<-sum(step2[step2$ss_sample_id %in% ms_ss$ss_sample_id,"Number of reads"])

#import counts step3
step3<-data.table::fread(paste0(MBCtsvDir,"step3_stats.tsv"),data.table = F)
step3$ss_sample_id<-gsub(".none.*$","",step3$Stats)
after_cutadapt<-sum(step3[step3$ss_sample_id %in% ms_ss$ss_sample_id,"nseqs"])

# #import counts step4 - THIS IS UNIQ, SO I GUESS IT DOESNT AFFECT READ COUNT?
# step4<-data.table::fread(paste0(MBCtsvDir,"step4_stats.tsv"),data.table = F)
# step4$ss_sample_id<-gsub(".none.*$","",step4$Stats)
# OTUs_after_cutadapt<-sum(step4[step4$ss_sample_id %in% ms_ss$ss_sample_id,"nseqs"])
# 
# sum(step4[,"nseqs"])

#import counts step5
#command used to make "step5_stats_BAS.tsv"

# find . -name \*none.flash2_merged.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.fasta.gz -print0 
# | xargs -0 zgrep ">" | sed 's/\.\///;s/.*\///;s/.none.*size=/ /' > step5_stats_BAS.tsv         
# 
# step5<-data.table::fread(paste0(MBCtsvDir,"step5_stats_BAS.tsv"),data.table = F)
# step5_A<-step5[step5$V1 %in% ms_ss$ss_sample_id,]
# sum(step5_A$V2)

#huh, well after all that, this is the same read count as in otutab, so can skip it


#import otu tab
otutab_A<-data.table::fread(otutab,data.table = F)
#subset 
otutab_B<-cbind(OTU=otutab_A$`#OTU ID`,otutab_A[,colnames(otutab_A) %in% ms_ss$ss_sample_id])
#remove 0-read OTUs and samples
otutab_C<-rm.0readOTUSam(taxatab = otutab_B)
#OTUs.in.otutab<-length(otutab_C$OTU) #NOT REALLY NECESSARY
after_size_select<-sum(otutab_C[,-1])

#import first taxa tab
taxatab_A<-data.table::fread(illumina_script_taxatab,data.table = F)
#subset
taxatab_B<-cbind(taxon=taxatab_A$taxon,taxatab_A[,colnames(taxatab_A) %in% ms_ss$ss_sample_id])
#remove 0-read OTUs and samples
taxatab_C<-rm.0readtaxSam(taxatab = taxatab_B)
##this is the same as otutab

after_blast<-sum(taxatab_C[taxatab_C$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",-1])
after_blast_filt<-after_blast-sum(taxatab_C[taxatab_C$taxon=="NA;NA;NA;NA;NA;NA;NA",-1])

#import filtered taxa tab
taxatab_A<-data.table::fread(illumina_script_taxatab_tf,data.table = F)
#subset 
taxatab_B<-cbind(taxon=taxatab_A$taxon,taxatab_A[,colnames(taxatab_A) %in% ms_ss$ss_sample_id])
#remove 0-read OTUs and samples
taxatab_C<-rm.0readtaxSam(taxatab = taxatab_B)
#remove no hits and NA
taxatab_D<-taxatab_C[taxatab_C$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]
taxatab_D<-taxatab_D[taxatab_D$taxon!="NA;NA;NA;NA;NA;NA;NA",]
after.taxon.filter<-sum(taxatab_D[,-1])


out<-data.frame("Demuliplexed files"=demuliplexed_files,
                "After paired end alignment"=after_paired_end,
                "After primer trimming"=after_cutadapt,
                "After size selection"=after_size_select,
                "After blast"=after_blast,
                "After blast filters"=after_blast_filt,
                "After initial taxon filter"=after.taxon.filter)

colnames(out)<-gsub("\\."," ",colnames(out))

return(out)
}
