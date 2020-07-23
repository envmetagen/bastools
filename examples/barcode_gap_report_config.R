#Use straight self blast to find barcoding gaps
#e.g. find the lowest pident that will exclude other families

#step1
#starts with a blast file. This must has qseqid,staxid,saccver,pident,sseq. Theoretically must be a thorough blast (e.g. metabinkit default)
#creates a database using the sseq of the best hit of each staxid
#self-blasts the database, using thorough blast again
#step2
#from results it flags potentially spurious entries (see script for rules). Output can be inspected, and/or used by metabin.
#step3
#using corrected dtaabase, it plots the within taxon variation vs the among taxon variation (as ranges; at S, G and F level) 
#it plots for each taxon, plots can be limited by plot.limit.taxon
#it produces two summary plots 
#plots should help to choose thresholds at F, G, S level
#idea is that we want the lowest threshold (to keep lots of reads) that will still exclude most "among taxa" hits (to maximise resolution)
#choosing a threshold too low will allow all reads through and give poor resolution
#choosing a threshold too high will not allow many reads through
#however, using metabin, neither should give wrong results

bastoolsDir<-"/home/tutorial/TOOLS/bastools/"
input<-"/home/tutorial/temp/16S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.blast.filt.txt" #example
#input not required for step3, just selfblastout
ncbiTaxDir<-"/home/tutorial/TOOLS/metabinkit.install/db/"
KronaPath<-"/home/tutorial/TOOLS/Krona.install/bin/ktImportText"
outDir<-"/home/tutorial/temp/"
steps<-c("thresher") #options: "selfblast", "find_db_errors", "calc_bar_gap","metabin","thresher"
threshersteps<-c("threshplots") #"threshblast","threshbin","threshplots"

plot.limit.taxon<-c(";Mammalia;") # this is a grep, can be NULL #(step3) #can be multiple #suggest including ;xxx;
known_flags<-"/media/sf_Documents/WORK/CIBIO/temp/known_flags_downloaded_22-7-20.txt" #"/home/tutorial/SCRIPTS/PB_tests/known_flags.txt" #full path, can be NULL

TaxlevelTestall<-c("K","P","C","O","F") # not required for step3

#thresher
plot.at.level<-"F" #for thresher
limit.plot.to.taxon<-c("Mammalia","C") #for thresher, the taxon name and level can be NULL
counts.out<-"16S_PBthresh_counts.tsv"#for thresher
final.table.out<-"16S_PBthresh_final.table.tsv" #for thresher
divergence<-rep(5,length(TaxlevelTestall)) # not required for step3
use_flagged_accessions_bcg=F
use_flagged_accessions_mbk=T
#binning settings to loop through
tops<-c(0,1,100)
#order=S,G,F,AF
pidents.list<-list(one=c(99,97,95,90),two=c(98,94,92,88),three=c(93,85,75,60)) #can be more than three

#outfiles
krona_html_db="database.html" #(step1)
selfblastout="16S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.blast.filt.tempBLASTDB.tsv" #(step1)
flagged_accessions="16S_flagged.tsv" #(step2) #required step 3 only if use_flagged_accessions_bcg=T
flagged_error_detailed_table="16S_flagged_details.tsv" #(step2) #not required for step3
out_html="16S.test_barcode_gap_report.html" #full path
further_potential_errors="16S.test_non-flagged.potential.errors.tsv" #required step3

#knit
rmarkdown::render(input = paste0(bastoolsDir,"scripts/barcode_gap_report.Rmd"),output_file = paste0(outDir,out_html))
