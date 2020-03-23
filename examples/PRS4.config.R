#taxatable filtering and report config

#usually will be for one primer/run, but if all settings can be applied globally then multiple tables can be provided

#It is intended that this is an iterative process. Run the inital script (examplez/example_initial_trial_POST_REBIN_SCRIPT_congif.R)
#with no filters, inspect, change settings, run again, repeat.

#STEPS
# Download and subset master sheet
# Import taxatabs (and merge if necessary) and subset according to subsetted master sheet
# Remove problem taxa - usually human, predator (along with NAs, nohits)
# Apply taxon and sample filters (doing this after problem taxa/NAs/no_hits to be more conservative in removals)
# Apply detection filter (this is an absolute value applied at the PCR replicate level)
# Produce a mini-report on taxa/reads in negatives
# Run the remove contaminants filter and output report **Note sure if this should be here or before taxon/sample/dxn filters**
# Remove detection occurring in one replicate only ***I tried this before contaminants filter but didnt like results**
# Sum replicates at chosen group (usually Sample_Name or biomaterial) 
# Apply 2nd detection filter (this is an absolute value applied at the summed rep level)
# Aggregate reads at chosen level, e.g. family
# Keep only taxa that reached chosen level
# Remove unwanted taxa (non-targets that we dont care about)
# Provide summary pre-grouping
# Do grouping (e.g. group a genus to a species because we know from local knowledge it couldnt be anything else in that genus)
# Provide summary post-grouping
# Produce some plots using desired variables (e.g. Site, biomaterial, date)
# Produce table/plot of step counts


#Potential improvements
#make re-ordering of steps an option

#############################################################################

#change setting below as necessary
bastoolsDir<-"/home/tutorial/TOOLS/bastools/" #change to your bastools directory
#before running for the very first time run this, then hash them out again:
# email="your.name@email.com"
# setwd(bastoolsDir) 
# googlesheets4::sheets_auth(email = email)

#Once settings are changed, open the Rmd script and change line 14 to match the path to this config file. Then Run All or knit

#can be multiple
taxatabs<-c("/media/sf_Documents/WORK/G-DRIVE/G-WORK/SHARED_FOLDERS/CRAYFISH/rebin/CRAY-HSJUN19BAS_COI.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.rebins.taxatable.tf.spliced.ALTEREDNAMES.txt")

# get sheet first separately (better than including in pipe)
# source(paste0(bastoolsDir,"master_functions.R"))
# mastersheet<-google.overlord("https://docs.google.com/spreadsheets/d/1k1mAGogWq9rXcwBKDyxG9oZ0OrWreBRRcSEUw0RGwyk/edit?ts=5d776492#gid=1377121809",email="basegeter@gmail.com")
# write.table(mastersheet,paste0(outDir,"mastersheet.txt"),append = F,quote = F,row.names = F,sep = "\t")

#use pre-downloaded master sheet instead. Set to NULL if downloading from google. Tab - seperated
master_sheet_file="/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/CRAYFISH/STATS/Stats_2020_Feb/master_sheet_all.tsv"

#options for subsetting master sheet. This functions to sleect the samples you want to analyse.
#Each item in list is a column heading in master sheet and each character within the item should be what you want to include 
#(sample_type should always be lower case, even if it is not so on google)
subsetlist<-list(experiment_id="HSJUN19BAS",primer_set="LERAY-XT",MPLX="N",sample_type=c("Extraction_Negative","GIT_contents","PCR_negative"),
                 Replicate_Name=c("ZYMO"),substudy="main")

#problem taxa - specific taxa in high reads that are not targets, usually human or predator. This is a grep, so 
#form the taxa path accordingly and make sure you are only removing what you want (by checking output). 
#use problemTaxa<-c("NothingToAdd") to turn this off 
problemTaxa<-c("Eukaryota;Arthropoda;Malacostraca;Decapoda;Cambaridae;","Eukaryota;Arthropoda;Malacostraca;Decapoda;Astacidae;","omycota;"
               ,"Eukaryota;unknown;unknown;unknown;NA;NA;NA") 
#Detection below taxonpc % of taxon read count will be removed (0.1=0.1%)
taxonpc = 0
#Detection below samplepc % of sample read count will be removed (0.1=0.1%)
samplepc=0
#the absolute value for removing detections in pcr reps (maybe not necessary here, wait til after sumreps...?, as we consider anything in 2 reps to be true)
filter_dxn = 300
#the absolute value for removing detections after summing replicate level (not sure of level, could leave very low, as we consider anything in 2 reps to be true)
filter_dxn2 = 0
#sample_type used to descirbe your real samples (not negatives)
real = c("GIT_contents")
#negative types (as detailed in master sheet) and groups to which each one belongs (must be same order). Put neg.types=NULL and neg.groups=NULL 
#if no negtaives, or skipping this
neg.types=c("PCR_negative", "Extraction_Negative")
neg.groups=c("Sample_Plate","extraction_batch")
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
#sum reps by
sumrepsby<-"biomaterial"
#removing detections in <1 rep post rm.contaminants? 
rep.rm.second<-F
#collapse all taxa at this level.  Put xLevel=NULL to skip this
xLevel<-NULL
#Keep only taxa at/below this level:
zLevel<-"species"
#collapse taxa pre-grouping?
aggregate.pre.grouping=F
#group taxa before filtering?
group.taxa.before.filtering=T

#unwanted taxa - non-targets that are not required for final analysis. This is applied last. This is a grep, so 
#form the taxa path accordingly and make sure you are only removing what you want (by checking output).
#use problemTaxa<-c("NothingToAdd") to not use this 
unwantedTaxa<-c("NothingToAdd")
#this should be a dataframe with the taxon to group in column 1 and taxon to group to in column 2. Done last.
#Can be made separately,e.g. in excel, this is just example. Ensure the order is correct. put taxa.to.group=NULL to skip this
taxa.to.group<-NULL
#variables to use as groups for plots, set to NULL to skip making barplots (useful if running into memory problems)
plotting.vars<-NULL #note that 'geographic location (country and/or sea)' and 'geographic location (region and locality)'
#should be entered as 'country' or 'locality' respectively


#the taxonomic levels to use for plotting. Set to NULL to use the taxatable as is (which will be printed regardless)
#Note that multiple values of this and plotting.vars will make a lot of barplots!
plotting_levels<-c("class")
#hide legend (to avoid cowplot errors set to FALSE)
hidelegend=T
#split barplots by facet, set to null if not required
facetcol=NULL

#full path to this config file.
config="/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/CRAYFISH/STATS/Stats_2020_Feb/PRS3.config.R"

#full path to target out html
out_html<-"/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/CRAYFISH/STATS/Stats_2020_Feb/PRS3_COI_MAIN.html"

#full path to taxatab out, put to NULL if not desired
taxatab.out<-"/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/CRAYFISH/STATS/Stats_2020_Feb/PRS3_COI_MAIN.txt"
#make krona plot, put to F if not desired
krona.out=T

#if your computer has trouble allocating memory set this to FALSE, it will exclude the inspection of negative prior to applying filters,
# the memory issue is caused by trying to store too many plots in memory (which happend when there is lots of contamination)
#of course if you are applying filters, you may longer care about the first inspection anyway, so F is better
do.first.negative.inspection=F
#another thing to help in this regard is to reduce the number of plotting factors (e.g. to one) or disable making barplots at all (see above)
#but Now I notice closing and restarting R / Rstudio is the biggest help

#knit
rmarkdown::render(input = paste0(bastoolsDir,"scripts/POST_REBIN_SCRIPT4.Rmd"),output_file = out_html)
