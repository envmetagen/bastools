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
input<-"/home/tutorial/TOOLS/bastools/beta/PB_classifier/2019_August_002.VENE.lenFilt.trimmed.ids.SC4.pol.blast.txt" #example
ncbiTaxDir<-"/home/tutorial/TOOLS/metabinkit.install/db/"
#outDir<-"/home/tutorial/SCRIPTS/PB_tests/"
steps<-c("selfblast", "find_db_errors", "calc_bar_gap") #options: "selfblast", "find_db_errors", "calc_bar_gap"
plot.limit.taxon<-"Bivalvia" # this is a grep, can be NULL #step3
metabin_S<-98
metabin_G<-95
metabin_F<-85
use_flagged_accessions=T #use this file for metabinkit?
  
#outfiles
krona_html_db="/home/tutorial/SCRIPTS/PB_tests/database.html" #(step1)
selfblastout="/home/tutorial/SCRIPTS/PB_tests/tempBLASTDB.tsv" #(step1)
flagged_accessions="/home/tutorial/SCRIPTS/PB_tests/VENE_flagged.tsv" #(step2)
flagged_error_detailed_table="/home/tutorial/SCRIPTS/PB_tests/VENE_flagged_details.tsv" #(step2)
corrected_db="/home/tutorial/SCRIPTS/PB_tests/VENE_correct_database" #(step2, needed for plotting in step 3)
out_html="/home/tutorial/SCRIPTS/PB_tests/VENE_barcode_gap_report.html"
metabin_out="/home/tutorial/SCRIPTS/PB_tests/VENE_metabin" # .tsv will be added to this root

#knit
rmarkdown::render(input = paste0(bastoolsDir,"scripts/barcode_gap_report.Rmd"),output_file = out_html)