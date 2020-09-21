#post thresher blast

source("/home/tutorial/TOOLS/bastools/master_functions.R")
setwd("/home/tutorial/temp/")

ncbiTaxDir<-("/home/tutorial/TOOLS/metabinkit.install/db/")

bin.thresh(
blast.thresh.input = "/home/tutorial/temp/16S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.blast.filt.tempBLASTDBALL.tsv"
           ,tops=c(0,1,10,100),known_flags = "/home/tutorial/TOOLS/bastools/known_flags_downloaded_22-7-20.txt",
            final.table.out = "test.binthresh.out.tsv"
             ,pidents.list = list(one=c(99,97,94,85),two=c(98,96,93,80),three=c(97,95,92,80)),
            SpeciesBL = "/home/tutorial/temp/disabled.sp.iran.txt",
            GenusBL = "/home/tutorial/temp/disabled.g.iran.txt",
            FamilyBL = "/home/tutorial/temp/disabled.f.iran.txt",
            ncbiTaxDir = ncbiTaxDir)

plot.thresh(thresher.final.table = "test.binthresh.out.tsv",limit.plot.to.taxon = NULL,plot.at.level = "F")
