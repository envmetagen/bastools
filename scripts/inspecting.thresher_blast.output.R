#post thresher blast

source("/home/tutorial/TOOLS/bastools/master_functions.R")
setwd("/home/tutorial/temp/")


bin.thresh(blast.thresh.input = "/home/tutorial/temp/16S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.blast.filt.tempBLASTDBALL.tsv"
           ,tops=c(0,10),known_flags = "/home/tutorial/TOOLS/bastools/known_flags_downloaded_22-7-20.txt",final.table.out = "test.binthresh.out.tsv"
             ,pidents.list = list(one=c(99,98,97,95)))


plot.thresh(thresher.final.table = "test.binthresh.out.tsv",limit.plot.to.taxon = c("Mammalia","C"),plot.at.level = "S")
