#post thresher blast

source("/home/tutorial/TOOLS/bastools/master_functions.R")
setwd("/home/tutorial/temp/")


bin.thresh(
blast.thresh.input = "/home/tutorial/temp/16S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.blast.filt.tempBLASTDBALL.tsv"
           ,tops=c(0,1,100),known_flags = "/home/tutorial/TOOLS/bastools/known_flags_downloaded_22-7-20.txt",final.table.out = "test.binthresh.out.tsv"
             ,pidents.list = list(one=c(99,96,92,80),two=c(98,96,92,80)),SpeciesBL = "/home/tutorial/temp/disabled.sp.iran.txt",
            GenusBL = "/home/tutorial/temp/disabled.g.iran.txt",FamilyBL = "/home/tutorial/temp/disabled.f.iran.txt")

bin.thresh(
  blast.thresh.input = "/home/tutorial/temp/16S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.blast.filt.tempBLASTDBALL.tsv"
  ,tops=c(0),known_flags = "/home/tutorial/TOOLS/bastools/known_flags_downloaded_22-7-20.txt",final.table.out = "test.binthresh.out.tsv"
  ,pidents.list = list(one=c(0,0,0,0),two=c(100,100,100,100),three=c(99,99,92,80)))


plot.thresh(thresher.final.table = "test.binthresh.out.tsv",limit.plot.to.taxon = NULL,plot.at.level = "F")
# 

  