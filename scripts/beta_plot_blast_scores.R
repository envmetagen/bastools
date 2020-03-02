
#first show all paths of interest from taxatab, without removing "may not be improved taxa"

filtered.taxatab<-"/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/MSDEC18BAS/ICVERTS-12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.taxatable.tf.spliced.txt"
filtered_blastfile<-"/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/MSDEC18BAS/12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.blast.filt.txt"
binfile<-"/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/MSDEC18BAS/12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.bins.txt"

a<-check.low.res.df(filtered.taxatab,filtered_blastfile, binfile,rm.excess=F,out = F)

#find top hit for each pathofinterest
b<-aggregate(a$high.pident,by=list(a$pathofinterest),FUN=max)
colnames(b)<-c("pathofinterest","top_pident")

#add class and order cols
class<-as.data.frame(do.call(rbind,stringr::str_split(b$pathofinterest,";")))
class$path<-paste(class$V1,class$V2,class$V3,sep = ";")
b$class_path<-class$path

order<-as.data.frame(do.call(rbind,stringr::str_split(b$pathofinterest,";")))
order$path<-paste(order$V1,order$V2,order$V3,order$V4,sep = ";")
b$order_path<-order$path

#plot class
plot.class<-ggplot(b,aes(x=b$class_path,y=b$top_pident)) + geom_point() + theme(axis.text.x=element_text(size=8,angle=45, hjust=1))
#plot order
plot.order<-ggplot(b,aes(x=b$order_path,y=b$top_pident)) + geom_point() + theme(axis.text.x=element_text(size=8,angle=45, hjust=1))

#save files
ggsave(filename = "/home/bastian.egeter/top.blast.hits.class.pdf",plot = plot.class,device = "pdf",width = 15,height = 10)
ggsave(filename = "/home/bastian.egeter/top.blast.hits.order.pdf",plot = plot.order,device = "pdf",width = 15,height = 10)

write.table(b,file = "/home/bastian.egeter/top.blast.hits.txt",append = F,quote = F,row.names = F,sep = "\t")     

taxatab<-data.table::fread("/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/MSDEC18BAS/ICVERTS-12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.taxatable.tf.spliced.txt"
                           ,data.table = F)

stats.by.rank(taxatab)
                