
#DO CALCULATIONS AND STATS
make.primer.bias.tables<-function(originaldbtab,ecopcrfile,
                                  out_bias_file,out_mod_ecopcrout_file,Pf,Pr, obitaxoR,min_length,
                                  max_length){
  
  ##########################################
  #GENERAL CLEANING 
  
  #read results
  ecopcroutput<-data.table::fread(ecopcrfile,sep = "\t")
  #remove hits outside desired lengths
  ecopcroutput<-ecopcroutput[!ecopcroutput$amplicon_length<min_length,]
  ecopcroutput<-ecopcroutput[!ecopcroutput$amplicon_length>max_length,]
  #remove duplicates (i.e. pick one entry per AC, based on lowest mismatches)
  ecopcroutput$total_mismatches<-as.numeric(ecopcroutput$forward_mismatch)+as.numeric(ecopcroutput$reverse_mismatch)
  ecopcroutput <- ecopcroutput[order(ecopcroutput$AC,ecopcroutput$total_mismatches),]
  ecopcroutput<-ecopcroutput[!duplicated(ecopcroutput$AC),]
  #remove weird primer mismatches (only a few usually)
  ecopcroutput<-ecopcroutput[!nchar(ecopcroutput$forward_match,allowNA = T)<nchar(Pf),]
  ecopcroutput<-ecopcroutput[!nchar(ecopcroutput$reverse_match,allowNA = T)<nchar(Pr),]
  ##########################################
  #READ ORIGINAL DB
  originaldb<-as.data.frame(data.table::fread(originaldbtab,header = TRUE,sep = "\t"))
  colnames(originaldb)<-gsub("taxid","taxids",colnames(originaldb))
  originaldb<-add.lineage.df(originaldb,ncbiTaxDir)
  ##########################################
  #LIST FAMILIES IN ODB
  all_primer_bias<-data.frame(row.names = 1:length(unique(originaldb$F)))
  all_primer_bias$in.odb<-unique(originaldb$F)[order(as.character(unique(originaldb$F)))]
  ##########################################
  #COUNT NO. SEQS IN ORIGINALDB 
  originaldb$n<-1
  nseqs.odb<-aggregate(originaldb$n,by=list(originaldb$F),FUN = sum)
  colnames(nseqs.odb)<-c("in.odb","nseqs.odb")
  all_primer_bias<-merge(all_primer_bias,nseqs.odb,by = "in.odb",all.x = T)
  ##########################################
  #COUNT No. Taxa IN ORIGINALDB 
  originaldb$path<-paste0(originaldb$F,";",originaldb$G,originaldb$S)
  originaldb.taxa<-originaldb[!duplicated(originaldb$path),c("path","n")]
  nseqs.odb<-aggregate(originaldb.taxa$n,by=list(originaldb.taxa$path),FUN = sum)
  colnames(nseqs.odb)<-c("path","nseqs.odb.taxa")
  nseqs.odb$F<-do.call(rbind,strsplit(nseqs.odb$path,";"))[,1]
  nseqs.odb<-aggregate(nseqs.odb$n,by=list(nseqs.odb$F),FUN = sum)
  colnames(nseqs.odb)<-c("in.odb","ntaxa.odb")
  all_primer_bias<-merge(all_primer_bias,nseqs.odb,by = "in.odb",all.x = T)
  ##########################################
  #ADD LOGICAL IF FAMILY AMPED
  all_primer_bias$amplified<-all_primer_bias$in.odb %in% unique(ecopcroutput$family_name)
  ##########################################
  #COUNT NO. SEQS THAT AMPED 
  ecopcroutput$n<-1
  nseqs.amped<-aggregate(ecopcroutput$n,by=list(ecopcroutput$family_name),FUN = sum)
  colnames(nseqs.amped)<-c("in.odb","nseqs.amped")
  all_primer_bias<-merge(all_primer_bias,nseqs.amped,by = "in.odb",all.x = T)
  ##########################################
  #COUNT No. Taxa THAT AMPED
  ecopcroutput$path<-paste0(ecopcroutput$family_name,";",ecopcroutput$genus_name,ecopcroutput$species_name)
  ecopcroutput.taxa<-ecopcroutput[!duplicated(ecopcroutput$path),c("path","n")]
  nseqs.amped<-aggregate(ecopcroutput.taxa$n,by=list(ecopcroutput.taxa$path),FUN = sum)
  colnames(nseqs.amped)<-c("path","nseqs.amped.taxa")
  nseqs.amped$F<-do.call(rbind,strsplit(nseqs.amped$path,";"))[,1]
  nseqs.amped<-aggregate(nseqs.amped$n,by=list(nseqs.amped$F),FUN = sum)
  colnames(nseqs.amped)<-c("in.odb","ntaxa.amped")
  all_primer_bias<-merge(all_primer_bias,nseqs.amped,by = "in.odb",all.x = T)
  ##########################################
  #percentage seqs amped
  all_primer_bias$pc_seqs_amped<-round(all_primer_bias$nseqs.amped/all_primer_bias$nseqs.odb*100,digits = 2)
  #percentage taxa amped
  all_primer_bias$pc_taxa_amped<-round(all_primer_bias$ntaxa.amped/all_primer_bias$ntaxa.odb*100,digits = 2)
  ##########################################
  
  
  ##########################################
  #for rest of stats keep only unique barcodes for each family
  message("dereplicating by family")
  ecopcroutput$fullseq<-paste0(ecopcroutput$forward_match,ecopcroutput$sequence,ecopcroutput$reverse_match)
  ecopcroutput<-ecopcroutput[!duplicated(ecopcroutput[,c("family_name","fullseq")]),]
  message("Total_unique_barcodes=",length(ecopcroutput$AC))
  message("Total_families_odb=",length(all_primer_bias$in.odb))
  message("Total_families_amped=",sum(all_primer_bias$amplified))
  ##########################################
  ##########################################
  ##########################################
  ##########################################
  #COUNT NO. UNIQUE BARCODES THAT AMPED 
  nseqs.amped<-aggregate(ecopcroutput$n,by=list(ecopcroutput$family_name),FUN = sum)
  colnames(nseqs.amped)<-c("in.odb","n.uniq.brcds.amped")
  all_primer_bias<-merge(all_primer_bias,nseqs.amped,by = "in.odb",all.x = T)
  ##########################################
  #ADD LINEAGES TO ALL_PRIMER_BIAS AND ECOPCROUTPUT
  y=originaldb[,c("K","P","C","O","F")]
  y=y[!duplicated(y),]
  all_primer_bias<-merge(x = all_primer_bias,y = y,by.x = "in.odb",by.y = "F",all.y = F)
  ecopcroutput<-merge(x = ecopcroutput,y = y,by.x = "family_name",by.y = "F",all.y = F)
  
  ##########################################################################################
  #add 3 prime mms to ecopcroutput
  ecopcroutput<-add.3pmms(ecopcroutput,Pf,Pr) 
  #add tm
  ecopcroutput<-add.tm.ecopcroutput(ecopcroutput)
  #diff tm
  ecopcroutput$diff_tm<-ecopcroutput$fTms-ecopcroutput$rTms
  #add gc content
  ecopcroutput<-add.gc.ecopcroutput(ecopcroutput)
  #gc clamp present?###########
  #add diff ta tm
  ecopcroutput$tm.ta_fw<-ecopcroutput$fTms-Ta
  ecopcroutput$tm.ta_rv<-ecopcroutput$rTms-Ta  
  #Tm of last 6 bp of fw primer divided by overall Tm (%)
  ecopcroutput$tm_fw_3p6_perc<-ecopcroutput$fTms3prime6/ecopcroutput$fTms*100
  ecopcroutput$tm_rv_3p6_perc<-ecopcroutput$rTms3prime6/ecopcroutput$rTms*100
  #add taxonomic resolution
  ecopcroutput<-add.res.ecopcroutput(ecopcroutput)
  ###########################################################################################
  #mean no. primer mismatches fw
  mean_mms_fw<-calc.stat.ecopcroutput(ecopcroutput,variable="forward_mismatch","mean")
  #mean no. primer mismatches rv
  mean_mms_rv<-calc.stat.ecopcroutput(ecopcroutput,variable="reverse_mismatch","mean")
  #mean no. primer mismatches total
  mean_mms_total<-calc.stat.ecopcroutput(ecopcroutput,variable="total_mismatches","mean")
  #Mean no. 3 prime mismatches (last half and last 6 bases) 
  mean.3pmms<-calc.3pmms.fam(ecopcroutput,Pf,Pr) ######LIST
  mean_fmms3Phalf<-mean.3pmms[[1]]
  mean_rmms3Phalf<-mean.3pmms[[2]]
  mean_fmms3P6<-mean.3pmms[[3]]
  mean_rmms3P6<-mean.3pmms[[4]]
  #mean ftm 
  mean_ftm<-calc.stat.ecopcroutput(ecopcroutput,variable="fTms","mean")
  #mean rtm 
  mean_rtm<-calc.stat.ecopcroutput(ecopcroutput,variable="rTms","mean")
  #mean fTm 3' half
  mean_ftm3Phalf<-calc.stat.ecopcroutput(ecopcroutput,variable="fTms3primehalf","mean")
  #mean rTm 3' half
  mean_rtm3Phalf<-calc.stat.ecopcroutput(ecopcroutput,variable="rTms3primehalf","mean")
  #mean fTm for 3' half - 6bp  
  mean_ftm3P6<-calc.stat.ecopcroutput(ecopcroutput,variable="fTms3prime6","mean")  
  #mean rTm for 3' half - 6bp  
  mean_rtm3P6<-calc.stat.ecopcroutput(ecopcroutput,variable="rTms3prime6","mean")  
  #% of unqiue seqs that get to fam or better
  mean_fam.res<-calc.fam.res(ecopcroutput)
  #mean amplicon length
  mean_amplicon.len<-aggregate(x =  nchar(ecopcroutput$sequence), by = list(ecopcroutput$family_name),FUN = mean)
  colnames(mean_amplicon.len)<-gsub("x","mean_amplicon.len",colnames(mean_amplicon.len))
  #mean gc_fw
  mean_fgc<-calc.stat.ecopcroutput(ecopcroutput,variable="fgc","mean")
  #mean gc_rv
  mean_rgc<-calc.stat.ecopcroutput(ecopcroutput,variable="rgc","mean")
  #mean tm.ta.fw
  mean_tm.ta.fw<-calc.stat.ecopcroutput(ecopcroutput,variable="tm.ta_fw","mean")
  #mean tm.ta.rv
  mean_tm.ta.rv<-calc.stat.ecopcroutput(ecopcroutput,variable="tm.ta_rv","mean")
  #mean diff tm
  mean_diff_tm<-calc.stat.ecopcroutput(ecopcroutput,variable="diff_tm","mean")
  #mean tm_fw_3p6_perc
  mean_tm_fw_3p6_perc<-calc.stat.ecopcroutput(ecopcroutput,variable="tm_fw_3p6_perc","mean")
  #mean tm_rv_3p6_perc
  mean_tm_rv_3p6_perc<-calc.stat.ecopcroutput(ecopcroutput,variable="tm_rv_3p6_perc","mean")
  ##########################################################################################
  #VARIANCES
  #var no. primer mismatches fw
  var_mms_fw<-calc.stat.ecopcroutput(ecopcroutput,variable="forward_mismatch","var")
  #var no. primer mismatches rv
  var_mms_rv<-calc.stat.ecopcroutput(ecopcroutput,variable="reverse_mismatch","var")
  #var no. primer mismatches total
  var_mms_total<-calc.stat.ecopcroutput(ecopcroutput,variable="total_mismatches","var")
  #var amplicon length
  var_amplicon.len<-aggregate(x =  nchar(ecopcroutput$sequence), by = list(ecopcroutput$family_name),FUN = var)
  colnames(var_amplicon.len)<-gsub("x","var_amplicon.len",colnames(var_amplicon.len))
  #var gc_fw
  var_fgc<-calc.stat.ecopcroutput(ecopcroutput,variable="fgc","var")
  #var gc_rv
  var_rgc<-calc.stat.ecopcroutput(ecopcroutput,variable="rgc","var")
  #var tm.ta_fw
  var_tm.ta.fw<-calc.stat.ecopcroutput(ecopcroutput,variable="tm.ta_fw","var")
  #var tm.ta_rv
  var_tm.ta.rv<-calc.stat.ecopcroutput(ecopcroutput,variable="tm.ta_rv","var")
  #diff tm
  var_diff_tm<-calc.stat.ecopcroutput(ecopcroutput,variable="diff_tm","var")
  #var tm_fw_3p6_perc
  var_tm_fw_3p6_perc<-calc.stat.ecopcroutput(ecopcroutput,variable="tm_fw_3p6_perc","var")
  #var tm_rv_3p6_perc
  var_tm_rv_3p6_perc<-calc.stat.ecopcroutput(ecopcroutput,variable="tm_rv_3p6_perc","var")
  #mms_fw_3p6
  var_mms_fw_3p6<-calc.stat.ecopcroutput(ecopcroutput,variable="f_mismatches_3prime6","var")
  #mms_rv_3p6
  var_mms_rv_3p6<-calc.stat.ecopcroutput(ecopcroutput,variable="r_mismatches_3prime6","var")
  
  ##########################################################################################
  #compile
  h<-list()
  for(i in 1:length(ls(pattern = "mean\\_.*"))){
    h[[i]]<-get(ls(pattern = "mean\\_.*")[i])
  }
  all_means<-do.call(cbind,h)
  h<-list()
  for(i in 1:length(ls(pattern = "var\\_.*"))){
    h[[i]]<-get(ls(pattern = "var\\_.*")[i])
  }
  all_vars<-do.call(cbind,h)
  
  all_mean_and_var<-cbind(all_means,all_vars)
  all_mean_and_var<-all_mean_and_var[,-grep("family",colnames(all_mean_and_var))]
  all_mean_and_var$Group.1.1=NULL
  
  #compile all
  all_primer_bias<-merge(all_primer_bias,all_mean_and_var,all.x = T,by.x = "in.odb", by.y = "Group.1")
  
  #remove extraneous columns
  ecopcroutput$path=NULL
  ecopcroutput$n=NULL
  
  #write primer bias file
  write.table(x=all_primer_bias,file = out_bias_file,quote = F,sep = "\t",row.names = F)
  #write final, modified ecopcroutput file
  write.table(x=ecopcroutput, file = out_mod_ecopcrout_file,quote = F,sep = "\t",row.names = F)
}
