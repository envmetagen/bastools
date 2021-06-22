#taxatable filtering and report config

#usually will be for one primer/run, but if all settings can be applied globally then multiple tables can be provided
#It is intended that this is an iterative process. Run the inital script (examplez/example_initial_trial_POST_REBIN_SCRIPT_congif.R)
#with no filters, inspect, change settings, run again, repeat.

#############################################################################
#FILES

#change setting below as necessary
bastoolsDir<-"/home/tutorial/TOOLS/bastools/" #change to your bastools directory

#can be multiple...
taxatabs<-c("/media/sf_Documents/WORK/G-DRIVE/G-WORK/SHARED_FOLDERS/REPTILE/REPTILE_RUN/none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.taxatable.tf.txt")

# get master sheet first separately (better than including in pipe)
# source(paste0(bastoolsDir,"master_functions.R"))
# mastersheet<-google.overlord("https://docs.google.com/spreadsheets/d/1k1mAGogWq9rXcwBKDyxG9oZ0OrWreBRRcSEUw0RGwyk/edit?ts=5d776492#gid=1377121809")
# write.table(mastersheet,"MUSSELS_mastersheet.txt",append = F,quote = F,row.names = F,sep = "\t")

#path to pre-downloaded master sheet. Set to NULL if downloading from google. Tab - separated
master_sheet_file="/media/sf_Documents/WORK/G-DRIVE/G-WORK/SHARED_FOLDERS/REPTILE/master.sheet.postmbc.tsv"

#options for subsetting master sheet. This functions to sleect the samples you want to analyse.
#Each item in list is a column heading in master sheet and each character within the item should be what you want to include 
#(sample_type should always be lower case, even if it is not so on google)
subsetlist<-list(Replicate_Name="A")

#full path to this config file.
config="/home/tutorial/TOOLS/bastools/examples/taxatable.report.config.R"

#where do you want to save output? full path to target out html
out_html<-"/media/sf_Documents/WORK/G-DRIVE/G-WORK/SHARED_FOLDERS/REPTILE/RESULTS_AUG2020/test.html"

#full path to output the final taxatab, put to NULL if not desired
taxatab.out<-NULL

#################
#FILTERS

#problem taxa - specific taxa in high reads that are not targets, usually human or predator. This is a grep, so 
#form the taxa path accordingly and make sure you are only removing what you want (by checking output). 
#use problemTaxa<-c("NothingToAdd") to turn this off 
problemTaxa<-c("NothingToAdd") 
#Detection below fullpc % of overall taxatab read count will be removed (0.003=0.003%)
fullpc=0.003
#Detection below taxonpc % of taxon read count will be removed (0.1=0.1%)
taxonpc = 0
#Detection below samplepc % of sample read count will be removed (0.1=0.1%)
samplepc=0
#the absolute value for removing detections in pcr reps (maybe not necessary here, wait til after sumreps...?, as we consider anything in 2 reps to be true)
filter_dxn = 0
#the absolute value for removing detections after summing replicate level (not sure of level, could leave very low, as we consider anything in 2 reps to be true)
filter_dxn2 = 0
#sample_type used to describe your real samples (not negatives)
real = c("faecal_sample")
#negative types (as detailed in master sheet) and groups to which each one belongs (must be same order). Put neg.types=NULL and neg.groups=NULL 
#if no negtaives, or skipping this
neg.types=c("Extraction_Negative")
neg.groups=c("extraction_batch")
#Do first inspection of negatives prior to applying filters, useful, but long output if lots of negatives with reads
do.first.negative.inspection=T
#use.contamination.filter function?
use.contamination.filter=F
#remove contaminants, not just from batches/groups, but from the entire data set?
remove.entire.dataset=F
#only remove detections of a taxon if they are less than the number of reads for that taxon in the negative
rm.only.less.than=F
#grouping for removing detections in <1 rep
rep.rm<-"biomaterial"
#removing detections in <1 rep prior to rm.contaminants? 
rep.rm.first<-F
#removing detections in <1 rep post rm.contaminants? 
rep.rm.second<-F
#collapse all taxa at this level.  Put xLevel=NULL to skip this
xLevel<-NULL
#Keep only taxa at/below this level:
zLevel<-NULL
#collapse taxa pre-grouping?
aggregate.pre.grouping=F
#group taxa before filtering?
group.taxa.before.filtering=F
#sum reps by
sumrepsby<-"ss_sample_id"
#unwanted taxa - non-targets that are not required for final analysis. This is applied last. This is a grep, so 
#form the taxa path accordingly and make sure you are only removing what you want (by checking output).
#use problemTaxa<-c("NothingToAdd") to not use this 
unwantedTaxa<-c("NothingToAdd")
#grouping taxa. this should be a dataframe with the taxon to group in column 1 and taxon to group to in column 2. Done last.
#Can be made separately,e.g. in excel, this is just example. Ensure the order is correct. put taxa.to.group=NULL to skip this
taxa.to.group<-NULL

#################
#PLOTTING

#variables to use as groups for plots, set to NULL to skip making
plotting.vars<-c("biomaterial","predator_species") #these will be the x-axes, for each plotting.var a different plot is made
#the taxonomic levels to use for plotting. Set to NULL to use the taxatable as is (which will be printed regardless)
#for each plotting_level a different plot is made
plotting_levels<-c("species","genus","family") #applied only to barplots and pca plots, not heatmaps
#hide legend (for barplots, to avoid cowplot errors set to FALSE)
hidelegend=T
#hide legend for pca plots
hidelegend.pca=T
#split plots by facet, set to NULL if not required, usually for "below" plotting.vars (not implemeted for pca plots yet)
#facetcol is also used to create a comparison taxatable, which highlights taxa detected via one facet or another. This is only produced if the facet has two categories. 
facetcol=NULL
#add a colour bar to heatmap plots (only makes sense if "above" plotting.vars). 
#can be NULL, can be up to three - samples will be ordered according to the last one
colour.bar=c("predator_species")
#make krona plot of final taxatab, full path, put to F if not desired
krona.out<-F

plot.bars=F
plot.pca=F
show.pca.lines=F
plot.heatmaps=T

#knit
rmarkdown::render(input = paste0(bastoolsDir,"scripts/TAXATABLE_REPORT.Rmd"),output_file = out_html)
