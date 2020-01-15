build.refs<-function(input.ecopcr.results,output){
#select one hit per genus
a<-input.ecopcr.results[!duplicated(input.ecopcr.results$genus),]
#concatenate the primer forward binding site, sequence, reverse binding site
a$reverse_matchRC<-insect::rc(z = a$reverse_match)
a$newseq<-paste0(a$forward_match, a$sequence,a$reverse_matchRC)
a$definition<-paste0(a$AC," genus=",a$genus_name,"; taxid=",a$genus,";")
export<-a[,c("definition","newseq")]
invisible(seqRFLP::dataframe2fas(export,file = output))
}
