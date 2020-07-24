#taxatable filtering and report config

#usually will be for one primer/run, but if all settings can be applied globally then multiple tables can be provided
#It is intended that this is an iterative process. Run the inital script (examplez/example_initial_trial_POST_REBIN_SCRIPT_congif.R)
#with no filters, inspect, change settings, run again, repeat.

#############################################################################
#FILES

#change setting below as necessary
bastoolsDir<-"/home/tutorial/TOOLS/bastools/" #change to your bastools directory

#can be multiple...
taxatabs<-c("/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/MUSSELS/Analysis_March_2020/2019_August_002.UNIO.lenFilt.trimmed.ids.SC1.pol.taxatab.tf.txt")

# get master sheet first separately (better than including in pipe)
# source(paste0(bastoolsDir,"master_functions.R"))
# mastersheet<-google.overlord("https://docs.google.com/spreadsheets/d/1k1mAGogWq9rXcwBKDyxG9oZ0OrWreBRRcSEUw0RGwyk/edit?ts=5d776492#gid=1377121809")
# write.table(mastersheet,"MUSSELS_mastersheet.txt",append = F,quote = F,row.names = F,sep = "\t")

#path to pre-downloaded master sheet. Set to NULL if downloading from google. Tab - separated
master_sheet_file="/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/MUSSELS/Analysis_March_2020/MUSSELS_mastersheet.txt"

#options for subsetting master sheet. This functions to sleect the samples you want to analyse.
#Each item in list is a column heading in master sheet and each character within the item should be what you want to include 
#(sample_type should always be lower case, even if it is not so on google)
subsetlist<-list(experiment_id=("2019_August_002"),primer_set=c("UNIO"))

#full path to this config file.
config="/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/MUSSELS/Analysis_March_2020/MUSSELS_ALL_PRS4.config.R"

#full path to target out html
out_html<-"/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/MUSSELS/Analysis_March_2020/MUSSELS_ALL_PRS4.html"

#full path to output the final taxatab, put to NULL if not desired
taxatab.out<-NULL

#################
#FILTERS

#problem taxa - specific taxa in high reads that are not targets, usually human or predator. This is a grep, so 
#form the taxa path accordingly and make sure you are only removing what you want (by checking output). 
#use problemTaxa<-c("NothingToAdd") to turn this off 
problemTaxa<-c("NothingToAdd") 
#Detection below taxonpc % of taxon read count will be removed (0.1=0.1%)
taxonpc = 0
#Detection below samplepc % of sample read count will be removed (0.1=0.1%)
samplepc=0
#the absolute value for removing detections in pcr reps (maybe not necessary here, wait til after sumreps...?, as we consider anything in 2 reps to be true)
filter_dxn = 0
#the absolute value for removing detections after summing replicate level (not sure of level, could leave very low, as we consider anything in 2 reps to be true)
filter_dxn2 = 0
#sample_type used to descirbe your real samples (not negatives)
real = c("water")
#negative types (as detailed in master sheet) and groups to which each one belongs (must be same order). Put neg.types=NULL and neg.groups=NULL 
#if no negtaives, or skipping this
neg.types=c("Field_Negative")
neg.groups=c("geographic location (region and locality)")
#Do first inspection of negatives prior to applying filters, useful, but long output if lots of negatives with reads
do.first.negative.inspection=F
#use.contamination.filter function?
use.contamination.filter=T
#remove contaminants, not just from batches/groups, but from the entire data set?
remove.entire.dataset=T
#only remove detections of a taxon if they are less than the number of reads for that taxon in the negative
rm.only.less.than=T
#grouping for removing detections in <1 rep
rep.rm<-"biomaterial"
#removing detections in <1 rep prior to rm.contaminants? 
rep.rm.first<-F
#removing detections in <1 rep post rm.contaminants? 
rep.rm.second<-F
#collapse all taxa at this level.  Put xLevel=NULL to skip this
xLevel<-NULL
#Keep only taxa at/below this level:
zLevel<-"family"
#collapse taxa pre-grouping?
aggregate.pre.grouping=F
#group taxa before filtering?
group.taxa.before.filtering=T
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
plotting.vars<-c("biomaterial","primer_set") #these will be the x-axes, for each plotting.var a different plot is made
#the taxonomic levels to use for plotting. Set to NULL to use the taxatable as is (which will be printed regardless)
#for each plotting_level a different plot is made
plotting_levels<-c("species","genus","family")
#hide legend (to avoid cowplot errors set to FALSE)
hidelegend=T
#split plots by facet, set to NULL if not required, usually for "below" plotting.vars
facetcol="pipe"
#add a colour bar to heatmap plots (only makes sense if "above" plotting.vars). can be NULL
colour.bar="geographic location (region and locality)"
#make krona plot of final taxatab, put to F if not desired
krona.out<-F

#knit
rmarkdown::render(input = paste0(bastoolsDir,"scripts/POST_REBIN_SCRIPT4.Rmd"),output_file = out_html)
