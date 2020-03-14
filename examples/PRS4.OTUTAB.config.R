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
  
  #can be multiple
  taxatabs<-c("/media/sf_Documents/WORK/CIBIO/temp/12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.otutab.tsv")
  
  #blastfiltfile
  blastfilt<-"/media/sf_Documents/WORK/CIBIO/temp/12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.blast.filt.txt"

    #full path to this config file.
  config="/home/tutorial/TOOLS/bastools/examples/PRS4.OTUTAB.config.R"
  
  #full path to target out html
  out_html<-"/media/sf_Documents/WORK/CIBIO/temp/PRS4.OTUTAB.TEST.html"
  
  #full path to otutab out, put to NULL if not desired
  taxatab.out<-"/media/sf_Documents/WORK/CIBIO/temp/PRS4.OTUTAB.out.txt"
  
  #datasheet url, not required if master_sheet_file is not NULL
  #url<-"https://docs.google.com/spreadsheets/d/1FUSaeVaYzms2EOGUoCAB4jaRKzguD3AKTsC8lYwaKP4/edit?ts=5dae01be#gid=0"
  #email="basegeter@gmail.com"
  
  #to get sheet separately
  #mastersheet<-google.overlord("https://docs.google.com/spreadsheets/d/1Kw8ONKtlX1hGr0LM-FvkowXYElEJ_taFNu2ZSXqp_n8/edit#gid=0")
  #write.table(mastersheet,"/media/sf_Documents/WORK/G-DRIVE/G-WORK/SHARED_FOLDERS/IRAN-ILLUMINA_SCRIPT_OUTPUT/After_rebinning/mastersheet.txt",append = F,quote = F,row.names = F,sep = "\t")
  
  #use pre-downloaded master sheet instead. Set to NULL if downloading from google. Tab - seperated
  master_sheet_file="/media/sf_Documents/WORK/G-DRIVE/G-WORK/SHARED_FOLDERS/IRAN-ILLUMINA_SCRIPT_OUTPUT/After_rebinning/mastersheet.txt"
  
  #options for subsetting master sheet. This functions to sleect the samples you want to analyse.
  #Each item in list is a column heading in master sheet and each character within the item should be what you want to include 
  #(sample_type should always be lower case, even if it is not so on google)
  subsetlist<-list(sample_type=c("faecal_sample","PCR_negative"), primer_set="12SV5.1")
  
  #Detection below taxonpc % of taxon read count will be removed (0.1=0.1%)
  taxonpc = 0.1
  #Detection below samplepc % of sample read count will be removed (0.1=0.1%)
  samplepc=0.1
  #the absolute value for removing detections in pcr reps (maybe not necessary here, wait til after sumreps...?, as we consider anything in 2 reps to be true)
  filter_dxn = 60
  #the absolute value for removing detections after summing replicate level (not sure of level, could leave very low, as we consider anything in 2 reps to be true)
  filter_dxn2 = 0
  #sample_type used to descirbe your real samples (not negatives)
  real = c("faecal_sample")
  #negative types (as detailed in master sheet) and groups to which each one belongs (must be same order). Put neg.types=NULL and neg.groups=NULL 
  #if no negtaives, or skipping this
  neg.types=c("PCR_negative")
  neg.groups=c("Sample_Plate")
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
  zLevel<-NULL
  
  
  #variables to use as groups for plots, set to NULL to skip making barplots (useful if running into memory problems)
  plotting.vars<-"library_layout"
  
  #the taxonomic levels to use for plotting. Set to NULL to use the taxatable as is (which will be printed regardless)
  #Note that multiple values of this and plotting.vars will make a lot of barplots!
  plotting_levels<-NULL
  
  #hide legend (to avoid cowplot errors set to FALSE)
  hidelegend=T
  
  #split barplots by facet, set to null if not required
  facetcol="predator"
  
  #make krona plot, put to F if not desired
  krona.out<-T
  
  #if your computer has trouble allocating memory set this to FALSE, it will exclude the inspection of negative prior to applying filters,
  # the memory issue is caused by trying to store too many plots in memory (which happend when there is lots of contamination)
  #of course if you are applying filters, you may longer care about the first inspection anyway, so F is better
  do.first.negative.inspection=T
  #another thing to help in this regard is to reduce the number of plotting factors (e.g. to one) or disable making barplots at all (see above)
  #but Now I notice closing and restarting R / Rstudio is the biggest help
  
  #knit
  rmarkdown::render(input = paste0(bastoolsDir,"scripts/OTU_FILTER.Rmd"),output_file = out_html)
