source("/home/tutorial/TOOLS/bastools/master_functions.R")

library(taxa)
library(metacoder)

#make a proportion of reads plot for each group

a<-data.table::fread("/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/REPTILE/reptile.teresatestZ.madagascariensisvsZ.brygooi_familylvl.taxatab.post.report.tsv"
                     ,data.table = F)
a<-tidy.taxon(a,rm.preceeding.above = "phylum")

ms<-data.table::fread("/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/REPTILE/master.sheet.postmbc.tsv"
                      ,data.table = F)

ms<-subset_mastersheet(master_sheet = ms,list(predator_species=c("Zonosaurus brygooi","Zonosaurus madagascariensis")))

#ms for samples in taxatab
ms<-ms[ms$biomaterial %in% colnames(a),]

a<-sumreps(taxatab = a,ms_ss = ms,grouping = "predator_species",current.grouping = "biomaterial",discard = F)

obj <- taxa::parse_tax_data(a,
                      class_cols = "taxon", # the column that contains taxonomic information
                      class_sep = ";", # The character used to separate taxa in the classification
)

obj$data$tax_data <- calc_obs_props(obj, "tax_data")
obj$data$taxon_counts <- calc_taxon_abund(obj, data = "tax_data")


obj$data$taxon_counts$Zb<-obj$data$taxon_counts$`Zonosaurus brygooi`
obj$data$taxon_counts$Zm<-obj$data$taxon_counts$`Zonosaurus madagascariensis`

heat_tree(obj, node_label = taxon_names, 
          #node_size = n_obs, 
          node_color = Zb)

heat_tree(obj, node_label = taxon_names, 
          #node_size = n_obs, 
          node_color = Zm)

###############################################################
#other things i tried

obj$data$diff_table <- compare_groups(obj,
                                      data = "tax_data",
                                      cols = ms$biomaterial, # What columns of sample data to use
                                      groups = ms$predator_species) # What category each sample is assigned to

print(obj$data$diff_table)

heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
          node_color = , # A column from `obj$data$diff_table`
          node_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
          node_color_range = c("cyan", "gray", "tan") # The color palette used
          #node_size_axis_label = "OTU count",
          #node_color_axis_label = "Log 2 ratio of median proportions",
          #layout = "davidson-harel", # The primary layout algorithm
          #initial_layout = "reingold-tilford"
          ) # The layout algorithm that initializes node locations


obj$data$taxon_counts <- calc_taxon_abund(obj, data = "tax_data")
obj$data$taxon_counts$total <- rowSums(obj$data$taxon_counts[, -1]) # -1 = taxon_id column
heat_tree(obj, node_label = taxon_names, 
          node_size = n_obs, 
          node_color = total)



