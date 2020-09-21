# bastools
A bunch of scripts for processing read and primer data

## Main layout
1. all the main functions are stored in master_functions.R This contains almost 200 functions for manipulating, plotting and processing metabarcoding data, downloading data, building blast databases, blasting, testing primers etc.
2. bigger scripts are in bastools/scripts
3. example config files for bigger scripts are in bastools/examples
4. a (rough) install file is located at /bastools/install/install.rlibs.and.robitools.R

## Other dependencies (not part of install file)
1. [krona tools](https://github.com/marbl/Krona/wiki)
2. [metabinkit](https://github.com/envmetagen/metabinkit)
3. [taxonkit](https://bioinf.shenwei.me/taxonkit/usage/) #installed by metabinkit
4. [blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

## Bigger scripts
1. illumina_post_MBC_script3.R
   - A script that takes as input fastas and otu tables and performs the following: BLASTs, adds taxonomy, checks for database taxonomic errors, plots barcode gap reports, testing multiple binning thresholds, bins BLAST results, creates taxa tables, creates contributor tables (to inspect which species contributed to each taxonomic bin), makes heatmaps and krona plots.
   - See [the example config file](examples/illuminaScript3.config.R) for more info 
2. TAXATABLE_REPORT.Rmd
   - A detailed html report which takes a taxa table and master sheet as input. Options include heatmaps, removing specific taxa, various filters, a function for reporting and removing contaminations from relevant batches/sites, PCoA plots, grouping taxa by taxonomic levels & mergeing taxatables. 
   - See [the example config file](examples/taxatable.report.config.R) for more info 
3. barcode_gap_report.Rmd
   - A subset of illumina_post_MBC_script3.R, used for checking database taxonomic errors, plotting barcode gap reports, testing multiple binning thresholds by self-blasting and loop-blasting the inputs. Produces html.
   - See [the example config file](examples/barcode_gap_report_config.R) for more info
4. inspecting.thresher_blast.output.R
   - A small script for re-running a small portion of above script
5. bas.minion.pipeline2.R
   - A script for processing minion data, largely defunkt since development of [msi](https://github.com/nunofonseca/msi) 
6. ENA_prep_script.R    
   - A script for preparing files for upload to ENA
   - See [the example config file](examples/preparing_ENA_files2.R) for more info 
7. primer_script_FULLGB_Feb_2020.R
   - long script to assess primer bias
