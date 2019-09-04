read.ecopcrdb<-function(ecopcr.db,ecopcr.tab,obitaxoR){
  message("Converting to tab")
  f<-process$new(command = "obitab", args=c("-o",ecopcr.db), echo_cmd = T,stdout=ecopcr.tab)
  f$wait()
  message("Reading in ecopcrdb")
  ECOPCRDB<-as.data.frame(data.table::fread(file = ecopcr.tab,header = TRUE,sep = "\t"))
  message("Adding taxonomy path")
taxon.table=get.classic.taxonomy.Bas(ECOPCRDB,obitaxdb = obitaxoR)
taxon.table$class_name_ok<-as.character(taxon.table$class_name_ok)
taxon.table$order_name_ok<-as.character(taxon.table$order_name_ok)
ECOPCRDB$kingdom_name<-taxon.table$kingdom_name_ok
ECOPCRDB$phylum_name<-taxon.table$phylum_name_ok
ECOPCRDB$class_name<-taxon.table$class_name_ok
ECOPCRDB$order_name<-taxon.table$order_name_ok
ECOPCRDB$family_name<-taxon.table$family_name_ok
ECOPCRDB$genus_name<-taxon.table$genus_name_ok
ECOPCRDB$species_name<-taxon.table$species_name_ok
ECOPCRDB
}
