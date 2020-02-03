#taxatable filtering and report config

#usually will be for one primer/run, but if all settings can be applied globally then multiple tables can be provided

#It is intended that this is an iterative process. Run, inspect, change settings, run again, repeat.

#STEPS
# Download and subset master sheet
# Import taxatabs (and merge if necessary) and subset according to subsetted master sheet
# Remove problem taxa - usually human, predator (along with NAs, nohits)
# Apply taxon and sample filters (doing this after problem taxa/NAs/no_hits to be more conservative in removals)
# Apply detection filter (this is an absolute value applied at the PCR replicate level)
# Produce a mini-report on taxa/reads in negatives
# Run the remove contaminants filter and output report (is this the correct time to do this?)********Doesnt work post sum reps anyway
# Remove detection occurring in one replicate only (should probably go before contaminants filter)**********but if it does, then we shouldnt 
#   remove single PCR negs. would need to modify function to recognize negs. 
#   In fact this could go after taxon/sample filters. It is a major cull, but a strongly supported one.
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
#make steps optional
#make re-ordering of steps an option
#update single rep detection removal, see above
#make table of taxa at each step (already a table made up to aggregating taxa. Any way of plotting this?)

#############################################################################

#change setting below as necessary
bastoolsDir<-"/home/tutorial/TOOLS/bastools/" #change to your bastools directory
#before running for the very first time run this, then hash them out again:
# email="your.name@email.com"
# setwd(bastoolsDir) 
# googlesheets4::sheets_auth(email = email)

###For HISEQ data
#colnames do not match current for hiseq so do
# test<-data.table::fread(taxatabs,data.table = F)
# colnames(test)<-c("taxon",paste(colnames(test[,-1]),"-HSBAS",sep = ""))
# write.table(test,file="/media/sf_Documents/WORK/CIBIO/temp/MOCK-COI.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.taxatable.tf.spliced.ALTEREDNAMES.txt",append = F,row.names = F,quote = F,sep = "\t")


#Once settings are changed, open bastools/CRAY_POST_REBIN_SCRIPT.Rmd
#change line 14 to match the path to this config file. Run CRAY_POST_REBIN_SCRIPT.Rmd (via Run All or knit). Run all is better for seeing some objects being created 
#along the way. Knit is better for viewing overall r.eport

#can be multiple
taxatabs<-c("/media/sf_Documents/WORK/CIBIO/temp/MOCK-COI.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.taxatable.tf.spliced.ALTEREDNAMES.txt")

#datasheet url
ss_url<-"https://docs.google.com/spreadsheets/d/1zWLuRROPoS49gMAvxdaySXSzhQoeXZWkFcojZkSL6r0/edit#gid=898638112"

#options for subsetting master sheet. This functions to sleect the samples you want to analyse.
#Each item in list is a column heading in master sheet and 
#each character within the item should be what you want to include (sample_type should always be lower case,
#even if it is not so on google)
subsetlist<-list(experiment_id="HSJUN19BAS",primer_set="LERAY-XT")

#problem taxa - specific taxa in high reads that are not targets, usually human or predator. This is applied early on. This is a grep, so 
#form the taxa path accordingly and make sure you are only removing what you want (by checking output). 
problemTaxa<-c("NothingToAdd")
#Detection below taxonpc % of taxon read count will be removed (0.1=0.1%)
taxonpc = 0.1
#Detection below samplepc % of sample read count will be removed (0.1=0.1%)
samplepc=0.1
#the absolute value for removing detections in pcr reps (maybe not necessary here, wait til after sumreps...?, as we consider anything in 2 reps to be true)
filter_dxn = 0
#the absolute value for removing detections after summing replicate level (not sure of level, could leave very low, as we consider anything in 2 reps to be true)
filter_dxn2 = 4
#sample_type used to descirbe your real samples (not negatives)
real = c("Tissue")
#negative types (as detailed in master sheet) and groups to which each one belongs (must be same order)
neg.types=c("")
neg.groups=c("")
#grouping for removing detections in <1 rep
rep.rm<-"Sample_Name"
#sum reps by
sumrepsby<-"Sample_Name"
#collapse all taxa at this level. Also used for subsequently keeping only taxa that attain this level 
xLevel<-NULL
#unwanted taxa - non-targets that are not required for final analysis. This is applied last. This is a grep, so 
#form the taxa path accordingly and make sure you are only removing what you want (by checking output).
unwantedTaxa<-c("NothingToAdd")
#this should be a dataframe with the taxon to group in column 1 and taxon to group to in column 2. Done last.
#Can be made separately,e.g. in excel, this is just example. Ensure the order is correct
taxa.to.group<-NULL
#variables to use as groups for plots
plotting.vars<-c("biomaterial")
