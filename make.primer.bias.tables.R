
#DO CALCULATIONS AND STATS
# sp_amplifications_per_family<-{
#   
#   #Amped or not
#   in_original_db<-as.character(unique(originaldb$F))
#   Families_amp_success<-as.data.frame(in_original_db)
#   Families_amp_success$amplified<-Families_amp_success$in_original_db %in% unique(ecopcroutput$family_name)
#   
#   #what % of species within each family were amplified?
#   length(unique(ecopcroutput$species_name))
#   length(ecopcroutput$species_name)
#   
#   amped_species_count<-aggregate(ecopcroutput$species_name[!duplicated(ecopcroutput$species_name)],
#                                  by = list(ecopcroutput$family_name[!duplicated(ecopcroutput$species_name)]),FUN=length)
#   orig_species_count<-aggregate(originaldb$S[!duplicated(originaldb$S)],
#                                 by = list(originaldb$F[!duplicated(originaldb$S)]),FUN=length)
#   
#   
#   
#   ##No. of species within each family in original
#   unique_species_ori<-originaldb[!duplicated(originaldb$species_name),]
#   pc_sp_ori<-as.data.frame(table(unique_species_ori[,c("family_name","species_name")]))
#   pc_sp_ori<-pc_sp_ori[pc_sp_ori$Freq>0,]
#   pc_sp_ori<-as.data.frame(table(pc_sp_ori$family_name,pc_sp_ori$Freq))
#   pc_sp_ori<-pc_sp_ori[pc_sp_ori$Freq>0,]
#   
#   #combine
#   pc_sp_ori_eco<-merge(pc_sp_ori[,c(1,3)],pc_sp_eco[c(1,3)],by = "Var1",all.x = T,all.y = F)
#   pc_sp_ori_eco$prop_amped<-pc_sp_ori_eco$Freq.y/pc_sp_ori_eco$Freq.x
#   colnames(pc_sp_ori_eco)<-c("family","no. sp in ecopcrdb","no. amplified","prop amped")

make.primer.bias.tables<-function(originaldb,ecopcroutput,out_bias_file,out_mod_ecopcrout_file,Pf,Pr, obitaxoR){
  
  in_original_db<-as.character(unique(originaldb$F))
  amplified<-unique(ecopcroutput$family_name)
  not_amplified<-in_original_db[!in_original_db %in% amplified]
  
  #combine results
  Families_amp_success<-data.frame(row.names = 1:length(in_original_db))
  #Families_amp_success$all.known.families<-all_families
  Families_amp_success$in.ecopcrdb<-in_original_db
  Families_amp_success$amplified<-in_original_db %in% amplified
  
  ################SHOULD ENSURE THIS WORKS FOR SPECIES COUNTS WITHIN FAMILES ALSO
  
  ###########################################################################################
  
  #PRIMER MISMATCHES BY FAMILY
  
  forward_mismatches<-aggregate(ecopcroutput[, "forward_mismatch"],list(ecopcroutput$family_name), mean)
  colnames(forward_mismatches)<-c("family","mean no. mismatches")
  ##reverse primer
  reverse_mismatches<-aggregate(ecopcroutput[, "reverse_mismatch"],list(ecopcroutput$family_name), mean)
  colnames(reverse_mismatches)<-c("family","mean no. mismatches")
  ##total mismatches
  ecopcroutput$total_mismatches<-ecopcroutput$forward_mismatch+ ecopcroutput$reverse_mismatch
  total_mismatches<-aggregate(ecopcroutput[, "total_mismatches"], list(ecopcroutput$family_name), mean)
  colnames(total_mismatches)<-c("family","mean no. mismatches")
  #combine
  merged_mismatches<-merge(forward_mismatches,reverse_mismatches,by="family")
  merged_mismatches<-merge(merged_mismatches,total_mismatches,by="family")
  colnames(merged_mismatches)<-c("family","mean_f_mms","mean_r_mms","mean_total_mms")
  
  ##########################################################################################
  
  #add 3' mismatches to ecopcroutput
  
  f_mismatch_table<-mismatch.table(ecopcroutput,Pf,"f")
  f_mismatches_3prime<-as.data.frame(rowSums(f_mismatch_table[,as.integer(nchar(Pf)/2):nchar(Pf)]))
  colnames(f_mismatches_3prime)<-"f_mismatches_3prime"
  r_mismatch_table<-mismatch.table(ecopcroutput, Pr,"r")
  r_mismatches_3prime<-as.data.frame(rowSums(r_mismatch_table[,as.integer(nchar(Pr)/2):nchar(Pr)]))
  colnames(r_mismatches_3prime)<-"r_mismatches_3prime"
  ecopcroutput<-cbind(ecopcroutput,f_mismatches_3prime,r_mismatches_3prime)
  
  #add 3' mismatches to ecopcroutput - 6bp
  f_mismatches_3prime6<-as.data.frame(rowSums(f_mismatch_table[,as.integer(nchar(Pf)-5):nchar(Pf)]))
  colnames(f_mismatches_3prime6)<-"f_mismatches_3prime6"
  r_mismatches_3prime6<-as.data.frame(rowSums(r_mismatch_table[,as.integer(nchar(Pr)-5):nchar(Pr)]))
  colnames(r_mismatches_3prime6)<-"r_mismatches_3prime6"
  ecopcroutput<-cbind(ecopcroutput,f_mismatches_3prime6,r_mismatches_3prime6)
  
  #mean primer mismatch tables for each base by family
  fam_f_mismatch_table<-family.mean.mismatch(ecopcroutput,Pf,"f")
  fam_r_mismatch_table<-family.mean.mismatch(ecopcroutput,Pr,"r")
  
  #mean 3' mismatches (last half of bases)
  fam_f_mismatches_3prime<-as.data.frame(rowMeans(fam_f_mismatch_table[,as.integer(nchar(Pf)/2+1):nchar(Pf)]))
  fam_f_mismatches_3prime$family<-rownames(fam_f_mismatches_3prime)
  colnames(fam_f_mismatches_3prime)<-c("mean_3prime_mms_f","family")
  fam_r_mismatches_3prime<-as.data.frame(rowMeans(fam_r_mismatch_table[,as.integer(nchar(Pr)/2+1):nchar(Pr)]))
  fam_r_mismatches_3prime$family<-rownames(fam_r_mismatches_3prime)
  colnames(fam_r_mismatches_3prime)<-c("mean_3prime_mms_r","family")
  
  #mean 3' mismatches (last 6 bases)
  fam_f_mismatches_3prime6<-as.data.frame(rowMeans(fam_f_mismatch_table[,as.integer(nchar(Pf)-5):nchar(Pf)]))
  fam_f_mismatches_3prime6$family<-rownames(fam_f_mismatches_3prime6)
  colnames(fam_f_mismatches_3prime6)<-c("mean_3prime6bp_mms_f","family")
  fam_r_mismatches_3prime6<-as.data.frame(rowMeans(fam_r_mismatch_table[,as.integer(nchar(Pr)-5):nchar(Pr)]))
  fam_r_mismatches_3prime6$family<-rownames(fam_r_mismatches_3prime6)
  colnames(fam_r_mismatches_3prime6)<-c("mean_3prime6bp_mms_r","family")
  
  ########################################################################################
  
  #Tm for each family
  
  ecopcroutput$fTms<-Tm.calc(ecopcroutput$forward_match)
  ecopcroutput$rTms<-Tm.calc(ecopcroutput$reverse_match)
  meanftmFam<-aggregate(x =  ecopcroutput$fTms,by = list(ecopcroutput$family_name), FUN = mean)
  colnames(meanftmFam)<-c("family","mean_fTm")
  meanrtmFam<-aggregate(x =  ecopcroutput$rTms,by = list(ecopcroutput$family_name), FUN = mean)
  colnames(meanrtmFam)<-c("family","mean_rTm")
  
  #mean Tm for 3' half
  ecopcroutput$fTms3primehalf<-Tm.calc(substr(x = ecopcroutput$forward_match,
                                              start = as.integer(nchar(as.character(ecopcroutput$forward_match))/2+1),
                                              stop = nchar(as.character(ecopcroutput$forward_match))))
  ecopcroutput$rTms3primehalf<-Tm.calc(substr(x = ecopcroutput$reverse_match,
                                              start = as.integer(nchar(as.character(ecopcroutput$reverse_match))/2+1),
                                              stop = nchar(as.character(ecopcroutput$reverse_match))))
  meanftm3PFam<-aggregate(x =  ecopcroutput$fTms3primehalf, by = list(ecopcroutput$family_name), FUN = mean)
  colnames(meanftm3PFam)<-c("family","mean_fTm_3P")
  meanrtm3PFam<-aggregate(x =  ecopcroutput$rTms3primehalf, by = list(ecopcroutput$family_name),FUN = mean)
  colnames(meanrtm3PFam)<-c("family","mean_rTm_3P")
  
  #mean Tm for 3' half - 6bp
  ecopcroutput$fTms3prime6<-Tm.calc(substr(x = ecopcroutput$forward_match,
                                           start = as.integer(nchar(as.character(ecopcroutput$forward_match))-5),
                                           stop = nchar(as.character(ecopcroutput$forward_match))))
  ecopcroutput$rTms3prime6<-Tm.calc(substr(x = ecopcroutput$reverse_match,
                                           start = as.integer(nchar(as.character(ecopcroutput$reverse_match))-5),
                                           stop = nchar(as.character(ecopcroutput$reverse_match))))
  meanftm3P6Fam<-aggregate(x =  ecopcroutput$fTms3prime6, by = list(ecopcroutput$family_name),FUN = mean)
  colnames(meanftm3P6Fam)<-c("family","mean_fTm_3P6")
  meanrtm3P6Fam<-aggregate(x =  ecopcroutput$rTms3prime6, by = list(ecopcroutput$family_name),FUN = mean)
  colnames(meanrtm3P6Fam)<-c("family","mean_rTm_3P6")
  
  ##########################################################################################
  
  #taxonomic resolution
  
  #add taxonomic resolution to ecopcroutput
  ecopcroutput.res<-add.res.Bas(ecopcroutput,obitaxdb=obitaxoR)
  
  #remove extra columns
  ecopcroutput.res$genus=NULL
  ecopcroutput.res$genus_name=NULL
  ecopcroutput.res$species_name=NULL
  ecopcroutput.res$species=NULL
  ecopcroutput.res$forward_tm=NULL
  ecopcroutput.res$reverse_tm=NULL
  
  #calculate percentage species that have tax res to family or better
  famsplit<-split(ecopcroutput.res,f = ecopcroutput.res$family_name)
  b<-lapply(X = famsplit,FUN = res.fam.or.better)
  pc.res<-as.data.frame(t(as.data.frame(b)))
  pc.res$family<-names(famsplit)
  pc.res$pc_res_to_family_or_better<-pc.res$V1
  pc.res$V1=NULL
  #############################################
  #compile all
  all_primer_bias<-merge(Families_amp_success,merged_mismatches,all.x = T,by.x = "in.ecopcrdb", by.y = "family")
  all_primer_bias<-merge(all_primer_bias,fam_f_mismatches_3prime,all.x = T,by.x = "in.ecopcrdb", by.y = "family")
  all_primer_bias<-merge(all_primer_bias,fam_r_mismatches_3prime,all.x = T,by.x = "in.ecopcrdb", by.y = "family")
  all_primer_bias<-merge(all_primer_bias,fam_f_mismatches_3prime6,all.x = T,by.x = "in.ecopcrdb", by.y = "family")
  all_primer_bias<-merge(all_primer_bias,fam_r_mismatches_3prime6,all.x = T,by.x = "in.ecopcrdb", by.y = "family")
  all_primer_bias<-merge(all_primer_bias,meanftmFam,all.x = T,by.x = "in.ecopcrdb", by.y = "family")
  all_primer_bias<-merge(all_primer_bias,meanrtmFam,all.x = T,by.x = "in.ecopcrdb", by.y = "family")
  all_primer_bias<-merge(all_primer_bias,meanftm3PFam,all.x = T,by.x = "in.ecopcrdb", by.y = "family")
  all_primer_bias<-merge(all_primer_bias,meanrtm3PFam,all.x = T,by.x = "in.ecopcrdb", by.y = "family")
  all_primer_bias<-merge(all_primer_bias,meanftm3P6Fam,all.x = T,by.x = "in.ecopcrdb", by.y = "family")
  all_primer_bias<-merge(all_primer_bias,meanrtm3P6Fam,all.x = T,by.x = "in.ecopcrdb", by.y = "family")
  all_primer_bias<-merge(all_primer_bias,pc.res,all.x = T,by.x = "in.ecopcrdb", by.y = "family")
  
  #write primer bias file
  write.table(x=all_primer_bias,file = out_bias_file,quote = F,sep = "\t",row.names = F)
  #write final, modified ecopcroutput file
  write.table(x=ecopcroutput.res, file = out_mod_ecopcrout_file,quote = F,sep = "\t",row.names = F)
  
}
