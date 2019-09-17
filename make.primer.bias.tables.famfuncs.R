
#Mean No. 3 prime PRIMER MISMATCHES (last half and last 6 bases) BY FAMILY
calc.3pmms.fam<-function(ecopcroutput,Pf,Pr){
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

list(fam_f_mismatches_3prime,fam_r_mismatches_3prime,fam_f_mismatches_3prime6,fam_r_mismatches_3prime6)
}

calc.fam.res<-function(ecopcroutput){
#calculate percentage species that have tax res to family or better
famsplit<-split(ecopcroutput,f = ecopcroutput$family_name)
b<-lapply(X = famsplit,FUN = res.fam.or.better)
pc.res<-as.data.frame(t(as.data.frame(b)))
pc.res$family<-names(famsplit)
colnames(pc.res)<-gsub("V1","pc_res_fam_or_better",colnames(pc.res))
pc.res
}
variable="r_mismatches_3prime6"
appliedstat=var
calc.stat.ecopcroutput<-function(ecopcroutput,variable,appliedstat="mean"){
  if(appliedstat=="mean"){
  a<-as.data.frame(aggregate(x =  ecopcroutput[,..variable], by = list(ecopcroutput$family_name),FUN = mean))
  colnames(a)<-c("family_name",paste0("mean_",colnames(a)[2]))}
  if(appliedstat=="var"){
    a<-as.data.frame(aggregate(x =  ecopcroutput[,..variable], by = list(ecopcroutput$family_name),FUN = var))
    colnames(a)<-c("family_name",paste0("var_",colnames(a)[2]))}
  a
}

