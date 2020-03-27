

newbins<-data.table::fread("12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.bins.MSDEC18BAS_newbins.txt"
                           ,sep="\t",data.table = F)

oldbins<-data.table::fread("../MSDEC18BAS/12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.bins.txt"
                           ,sep="\t",data.table = F)

compd1d2 <- dataCompareR::rCompare(newbins, oldbins)
summary(compd1d2)
dataCompareR::generateMismatchData(compd1d2,newbins,oldbins)
