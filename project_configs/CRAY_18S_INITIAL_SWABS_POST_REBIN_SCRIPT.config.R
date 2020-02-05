#taxatable filtering and report config

#usually will be for one primer/run, but if all settings can be applied globally then multiple tables can be provided

#It is intended that this is an iterative process. Run the inital script using the settings in this script 
#(examples/example_initial_trial_POST_REBIN_SCRIPT_config.R)

#with no filters, then inspect, change settings, run again, repeat.

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

#############################################################################
#taxatables. can be multiple
taxatabs<-c("/media/sf_Documents/WORK/G-DRIVE/G-WORK/SHARED_FOLDERS/CRAYFISH/rebin/CRAY-HSJUN19BAS_18S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.rebins.taxatable.tf.spliced.ALTEREDNAMES.txt")

#datasheet url
ss_url<-"https://docs.google.com/spreadsheets/d/1KZLoXHTgtkD0btSWjyAmFiGJ_cPcYITyfFSlzehisRI/edit#gid=1531090624"

#options for subsetting master sheet. This functions to sleect the samples you want to analyse.
#Each item in list is a column heading in master sheet and each character within the item should be what you want to include 
#(sample_type should always be lower case, even if it is not so on google)
subsetlist<-list(experiment_id="HSJUN19BAS",primer_set="18S",MPLX="N",substudy="swabtest",Replicate_Name="ZYMO")

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
real = c("swab")
#negative types (as detailed in master sheet) and groups to which each type belongs (must be same order). Put neg.types=NULL and neg.groups=NULL 
#if no negtaives. 
neg.types=c("PCR_negative", "Extraction_Negative")
neg.groups=c("Sample_Plate","extraction_batch")
#use.contamination.filter function? 
use.contamination.filter=T
#remove contaminants, not just from batches/groups, but from the entire data set? Be careful setting this to FALSE, could bias results.
remove.entire.dataset=T
#grouping for removing detections in <1 rep
rep.rm<-"biomaterial"
#sum reps by
sumrepsby<-"biomaterial"
#collapse all taxa at this level. Also used for subsequently keeping only taxa that attain this level. Put xLevel=NULL to skip this
xLevel<-"family"
#unwanted taxa - non-targets that are not required for final analysis. This is applied last. This is a grep, so 
#form the taxa path accordingly and make sure you are only removing what you want (by checking output).
#use problemTaxa<-c("NothingToAdd") to not use this 
unwantedTaxa<-c("NothingToAdd")
#this should be a dataframe with the taxon to group in column 1 and taxon to group to in column 2. Done last.
#Ensure the order is correct. put taxa.to.group=NULL to skip this
taxa.to.group<-NULL
#variables to use as groups for plots
plotting.vars<-c("biomaterial")
