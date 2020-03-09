#' Bin reads into taxa based on blast results
#'
#' @param blastfile blast tabular output, e.g. -outfmt 6. Can have any columns, but must have:
#'     qcovs, evalue, staxid, qseqid, pident
#' @param headers A space-separated string of the header names of \code{blastfile}.
#'     Usually this is the same as the blast command used.
#' @param ncbiTaxDir The full path to the directory containing ncbi taxonomy. 
#' @param obitaxdb The full path to the obitools-formatted ncbi taxonomy. 
#' @param out specify file to output results
#' @param min_qcovs The minimum query coverage for all hits
#' @param max_evalue The maximum evalue for all hits
#' @param top The percentage pident to consider hits. i.e. For each query find the top hit pident (% identity),
#'     then find all pidents within "top" of this value, discard others - e.g. if top set at 1, and top hit for
#'     query A has pident=99.8, then any hits to query A with pident above 98.8 are kept
#' @param spident The minimum pident for binning hits at species level
#' @param gpident The minimum pident for binning hits at genus level
#' @param fpident The minimum pident for binning hits at family level
#' @param abspident An absolute pident for binning hits above family level. Hits below this threshold will be
#'     excluded completely
#'
#' @return A tab-separated file consisting of ESV name and associated taxonomic path
#' @note It works like this:
#' \itemize{
#'     \item 1. Apply min_qcovs (query coverage) to whole table
#'     \item 2. Apply max_evalue (blast evalue) to whole table
#'     \item 3. For each query find the top hit pident (% identity),
#'     then find all pidents within "top" of this value, discard others - e.g. if top set at 1,
#'     and top hit for query A has pident=99.8, then any hits with pident above 98.8 are kept.
#'     \item 4. Add taxon path to all remaining hits
#'     \item 5. For species-level binning: first remove hits that do not have species-level information
#'     (i.e. where database entries were only labelled at genus or higher), then remove hits with pident<spident.
#'     Then, for each query, find the lowest-common-ancestor (lca) of the remaining hits.
#'     \item 6. Repeat above for genus and family level binning
#'     \item 7. Do a final binning for hits that had taxa identified at higher-than-family level
#'     \item 6. Final report=if a query had species-level lca, use that. If not, use genus. If not genus, use family.
#'     If not family use absolute.
#'
#' @examples
#' blastfile<-"blast.results.txt"
#' headers<-"qseqid sacc sseqid sallseqid sallacc bitscore score staxid staxids sscinames scomnames salltitles length pident qcovs"
#' ncbiTaxDir<-"/Documents/TAXONOMIES/"
#' obitaxdb<-"/Documents/obitax_26-4-19"
#' bin.blast.bas(blastfile,headers,ncbiTaxDir,obitaxdb,out="blast_bins.txt",min_qcovs=70,max_evalue=0.001,top=1,spident=99,gpident=97,fpident=93)
#' @export
bin.blast2<-function(filtered_blastfile,ncbiTaxDir,
                     obitaxdb,out,spident=98,gpident=95,fpident=92,abspident=80,
                     disabledTaxaFiles=NULL,disabledTaxaOut=NULL,
                     force=F,full.force=F,consider_sp.=F){
  t1<-Sys.time()
  
  if(is.null(out)) stop("out not specified")
  if(is.null(ncbiTaxDir)) stop("ncbiTaxDir not specified")
  if(is.null(obitaxdb)) stop("obitaxdb not specified")
  
  require(treemap)
  
  ###################################################
  #read in obitools taxonomy
  obitaxdb2=ROBITaxonomy::read.taxonomy(obitaxdb)
  
  btab<-data.table::fread(filtered_blastfile,sep="\t",data.table = F)
  
  #preparing some things for final step
  total_hits<-length(btab$qseqid) #for info later
  total_queries<-length(unique(btab$qseqid))
  qseqids<-as.data.frame(unique(btab$qseqid))
  qseqids$qseqid<-qseqids$`unique(btab$qseqid)`
  qseqids$`unique(btab$qseqid)`=NULL
  
  #read and check disabled taxa file(s) 
  if(!is.null(disabledTaxaFiles)){
    
    disabledTaxaDf<-merge.and.check.disabled.taxa.files(disabledTaxaFiles,disabledTaxaOut,force = force,full.force = full.force)
    
    #get taxids at 7 levels
    disabledTaxaDf<-add.lineage.df(disabledTaxaDf,ncbiTaxDir,as.taxids = T)
    
    disabledSpecies<-disabledTaxaDf[disabledTaxaDf$disable_species==T,]
    disabledSpecies<-disabledSpecies[!is.na(disabledSpecies$taxids),]
    #get children of all disabled species
    disabledSpecies$taxids<-disabledSpecies$S
    if(nrow(disabledSpecies)>0)  {
      childrenS<-get.children.taxonkit(disabledSpecies) 
    } else {
      childrenS<-NULL
    } 
    
    disabledGenus<-disabledTaxaDf[disabledTaxaDf$disable_genus==T,]
    disabledGenus<-disabledGenus[!is.na(disabledGenus$taxids),]
    #get children of all disabled genera
    disabledGenus$taxids<-disabledGenus$G
    if(nrow(disabledGenus)>0)  {
      childrenG<-get.children.taxonkit(disabledGenus) 
    } else {
      childrenG<-NULL
    }
    
    
    disabledFamily<-disabledTaxaDf[disabledTaxaDf$disable_family==T,]
    disabledFamily<-disabledFamily[!is.na(disabledFamily$taxids),]
    #get children of all disabled families
    disabledFamily$taxids<-disabledFamily$F
    if(nrow(disabledFamily)>0)  {
      childrenF<-get.children.taxonkit(disabledFamily) 
    } else {
      childrenF<-NULL
    }
    
    message("The following taxa are disabled at species level")
     if(!is.null(childrenS)){
      childrenAlldf<-as.data.frame(unique(c(childrenS,childrenG,childrenF)))
      colnames(childrenAlldf)<-"taxids"
      childrenNames<-add.lineage.df(childrenAlldf,ncbiTaxDir)
      childrenNames$pathString<-paste("Family",childrenNames$F,childrenNames$G,childrenNames$S,sep = "/")
      childrenNames$pathString<-lapply(childrenNames$pathString, gsub, pattern = "unknown", replacement = "", fixed = TRUE)
      disabledtree <- data.tree::as.Node(childrenNames)
      print(disabledtree,limit = NULL)
    } else (message("No species disabled"))
    
    message("The following taxa are disabled at genus level")
      if(!is.null(childrenG)){
      childrenAlldf<-as.data.frame(unique(c(childrenG,childrenF)))
      colnames(childrenAlldf)<-"taxids"
      childrenNames<-add.lineage.df(childrenAlldf,ncbiTaxDir)
      childrenNames$pathString<-paste("Family",childrenNames$F,childrenNames$G,sep = "/")
      childrenNames$pathString<-lapply(childrenNames$pathString, gsub, pattern = "unknown", replacement = "", fixed = TRUE)
      disabledtree <- data.tree::as.Node(childrenNames)
      print(disabledtree,limit = NULL)
      } else (message("No genera disabled"))
  
    
    message("The following taxa are disabled at family level")
      if(!is.null(childrenF)){
      childrenAlldf<-as.data.frame(unique(c(childrenF)))
      colnames(childrenAlldf)<-"taxids"
      childrenNames<-add.lineage.df(childrenAlldf,ncbiTaxDir)
      childrenNames$pathString<-paste("Family",childrenNames$F,sep = "/")
      childrenNames$pathString<-lapply(childrenNames$pathString, gsub, pattern = "unknown", replacement = "", fixed = TRUE)
      disabledtree <- data.tree::as.Node(childrenNames)
      print(disabledtree,limit = NULL)
      } else (message("No families disabled"))
  }
  
  #species-level assignments
  message("binning at species level")
  
  btabsp<-btab[btab$S!="unknown",]
  
  if(!is.null(disabledTaxaFiles)){
    btabsp<-btabsp[!btabsp$taxids %in% unique(c(childrenS,childrenG,childrenF)),]
  }
  
  if(consider_sp.==F){
    message("Not considering species with 'sp.', numbers or more than one space")
    if(length(grep(" sp\\.",btabsp$S,ignore.case = T))>0) btabsp<-btabsp[-grep(" sp\\.",btabsp$S,ignore.case = T),]
    if(length(grep(" .* .*",btabsp$S,ignore.case = T))>0) btabsp<-btabsp[-grep(" .* .*",btabsp$S,ignore.case = T),]
    if(length(grep("[0-9]",btabsp$S))>0) btabsp<-btabsp[-grep("[0-9]",btabsp$S),]
  } else(message("Considering species with 'sp.', numbers or more than one space"))
  
  btabsp<-btabsp[btabsp$pident>spident,]
  if(length(btabsp$taxids)>0){
    lcasp = aggregate(btabsp$taxids, by=list(btabsp$qseqid),function(x) ROBITaxonomy::lowest.common.ancestor(obitaxdb2,x))
    
    #get lca names
    colnames(lcasp)<-gsub("x","taxids",colnames(lcasp))
    if(sum(is.na(lcasp$taxids))>0){
      message("************
              ERROR: Some taxids were not recognized by ROBITaxonomy::lowest.common.ancestor, probably need to update obitaxdb using NCBI2obitaxonomy
              *************")
      problem.contributors<-btabsp[!duplicated(btabsp[lcasp[is.na(lcasp$taxids),1] %in% btabsp$qseqid,c("taxids","K","P","C","O","F","G","S")]),
                                   c("taxids","K","P","C","O","F","G","S")]
      for(i in 1:length(problem.contributors$taxids)){
        problem.contributors$validate[i]<-ROBITaxonomy::validate(obitaxdb2,problem.contributors$taxids[i])
      }
      print(problem.contributors[is.na(problem.contributors$validate),])
    }
    lcasp<-add.lineage.df(df = lcasp,ncbiTaxDir)
    colnames(lcasp)<-gsub("Group.1","qseqid",colnames(lcasp))
  } else {
    lcasp<-data.frame(matrix(nrow=1,ncol = 10))
    colnames(lcasp)<-c("taxids","qseqid","old_taxids","K","P","C","O","F","G","S")
  }
  
  rm(btabsp)
  
  #genus-level assignments
  message("binning at genus level") 
  
  btabg<-btab[btab$G!="unknown",]  
  #reason - can have g=unknown and s=known (e.g. Ranidae isolate), these should be removed
  #can have g=unknown and s=unknown (e.g. Ranidae), these should be removed
  #can have g=known and s=unknown (e.g. Leiopelma), these should be kept
  
  if(!is.null(disabledTaxaFiles)){
    btabg<-btabg[!btabg$taxids %in% unique(c(childrenG,childrenF)),]
  }
  
  btabg<-btabg[btabg$pident>gpident,] ####line changed 
  if(length(btabg$taxids)>0){
  lcag = aggregate(btabg$taxids, by=list(btabg$qseqid), function(x) ROBITaxonomy::lowest.common.ancestor(obitaxdb2,x))
  
  #get lca names
  colnames(lcag)<-gsub("x","taxids",colnames(lcag))
  if(sum(is.na(lcag$taxids))>0){
    message("************
            ERROR: Some taxids were not recognized by ROBITaxonomy::lowest.common.ancestor, probably need to update obitaxdb using NCBI2obitaxonomy
            *************")
    problem.contributors<-btabg[!duplicated(btabg[lcag[is.na(lcag$taxids),1] %in% btabg$qseqid,c("taxids","K","P","C","O","F","G","S")]),
                                c("taxids","K","P","C","O","F","G","S")]
    for(i in 1:length(problem.contributors$taxids)){
      problem.contributors$validate[i]<-ROBITaxonomy::validate(obitaxdb2,problem.contributors$taxids[i])
    }
    print(problem.contributors[is.na(problem.contributors$validate),])
  }
  lcag<-add.lineage.df(df = lcag,ncbiTaxDir)
  colnames(lcag)<-gsub("Group.1","qseqid",colnames(lcag))
  } else {
    lcag<-data.frame(matrix(nrow=1,ncol = 10))
    colnames(lcag)<-c("taxids","qseqid","old_taxids","K","P","C","O","F","G","S")
  }
  
  rm(btabg)
  
  #family-level assignments
  message("binning at family level") 
  #can have f=known, g=unknown, s=unknown, these should be kept
  #can have f=unknown, g=known, s=known, these should be kept
  #can have f=unknown, g=known, s=unknown, these should be kept
  #can have f=known, g=known, s=unknown, these should be kept
  #can have f=known, g=known, s=known, these should be kept
  #can have f=unknown, g=known, s=unknown, these should be kept
  #can have f=known, g=unknown, s=known, these should be kept
  
  #can have f=unknown, g=unknown, s=known, these should be removed - 
  #assumes that this case would be a weird entry (e.g. Ranidae isolate)
  
  #can have f=unknown, g=unknown, s=unknown, these should be removed
  
  btabf<-btab[!(btab$F=="unknown" & btab$G=="unknown" & btab$S=="unknown"),]  ####line changed 
  btabf<-btabf[!(btabf$F=="unknown" & btabf$G=="unknown"),] ####line changed 
  
  #this is ok, but when ncbi taxonomy does not have family-level assignment, we only get to order, stupid, or rather we do get to 
  #family or subfamily as lca, but final report puts it to "unknown"
  
  if(!is.null(disabledTaxaFiles)){
    btabf<-btabf[!btabf$taxids %in% unique(c(childrenF)),]
  }
  
  btabf<-btabf[btabf$pident>fpident,]
  lcaf = aggregate(btabf$taxids, by=list(btabf$qseqid), function(x) ROBITaxonomy::lowest.common.ancestor(obitaxdb2,x))
  
  #get lca names
  colnames(lcaf)<-gsub("x","taxids",colnames(lcaf))
  if(sum(is.na(lcaf$taxids))>0){
    message("************
            ERROR: Some taxids were not recognized by ROBITaxonomy::lowest.common.ancestor, probably need to update obitaxdb using NCBI2obitaxonomy
            *************")
    problem.contributors<-btabf[!duplicated(btabf[lcaf[is.na(lcaf$taxids),1] %in% btabf$qseqid,c("taxids","K","P","C","O","F","G","S")]),
                                c("taxids","K","P","C","O","F","G","S")]
    for(i in 1:length(problem.contributors$taxids)){
      problem.contributors$validate[i]<-ROBITaxonomy::validate(obitaxdb2,problem.contributors$taxids[i])
    }
    print(problem.contributors[is.na(problem.contributors$validate),])
  }
  lcaf<-add.lineage.df(df = lcaf,ncbiTaxDir)
  colnames(lcaf)<-gsub("Group.1","qseqid",colnames(lcaf))
  
  rm(btabf)
  
  #higher-than-family-level assignments
  message("binning at higher-than-family level")
  btabhtf<-btab[btab$K!="unknown",]
  
  # if(!is.null(disabledTaxaFiles)){
  #   btabhtf<-btabhtf[!btabhtf$taxids %in% childrenAll,]
  # }
  
  btabhtf<-btabhtf[btabhtf$pident>abspident,]
  lcahtf = aggregate(btabhtf$taxids, by=list(btabhtf$qseqid),
                     function(x) ROBITaxonomy::lowest.common.ancestor(obitaxdb2,x))
  
  #get lca names
  colnames(lcahtf)<-gsub("x","taxids",colnames(lcahtf))
  if(sum(is.na(lcahtf$taxids))>0){
    message("************
            ERROR: Some taxids were not recognized by ROBITaxonomy::lowest.common.ancestor, probably need to update obitaxdb using NCBI2obitaxonomy
            *************")
    problem.contributors<-btabhtf[!duplicated(btabhtf[lcahtf[is.na(lcahtf$taxids),1] %in% btabhtf$qseqid,c("taxids","K","P","C","O","F","G","S")]),
                                  c("taxids","K","P","C","O","F","G","S")]
    for(i in 1:length(problem.contributors$taxids)){
      problem.contributors$validate[i]<-ROBITaxonomy::validate(obitaxdb2,problem.contributors$taxids[i])
    }
    print(problem.contributors[is.na(problem.contributors$validate),])
  }
  
  lcahtf<-add.lineage.df(df = lcahtf,ncbiTaxDir)
  colnames(lcahtf)<-gsub("Group.1","qseqid",colnames(lcahtf))
  
  rm(btabhtf)
  
  ###################################################
  #combine
  sp_level<-lcasp[lcasp$S!="unknown",]
  g_level<-lcag[lcag$G!="unknown",]
  if(length(g_level$taxids)>0) g_level$S<-NA
  g_level<-g_level[!g_level$qseqid %in% sp_level$qseqid,]
  f_level<-lcaf[lcaf$F!="unknown",]
  if(length(f_level$taxids)>0) f_level$G<-NA
  if(length(f_level$taxids)>0) f_level$S<-NA
  f_level<-f_level[!f_level$qseqid %in% sp_level$qseqid,]
  f_level<-f_level[!f_level$qseqid %in% g_level$qseqid,]
  
  abs_level<-lcahtf
  if(length(abs_level$taxids)>0) abs_level$G<-NA
  if(length(abs_level$taxids)>0) abs_level$S<-NA
  if(length(abs_level$taxids)>0) abs_level$F<-NA
  abs_level<-abs_level[!abs_level$qseqid %in% sp_level$qseqid,]
  abs_level<-abs_level[!abs_level$qseqid %in% g_level$qseqid,]
  abs_level<-abs_level[!abs_level$qseqid %in% f_level$qseqid,]
  
  com_level<-rbind(sp_level,g_level,f_level,abs_level)
  com_level<-merge(x=qseqids, y = com_level[,c(2,4:10)], by = "qseqid",all.x = T)
  
  #info
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  
  write.table(x = com_level,file = out,sep="\t",quote = F,row.names = F)
  
  message(c("Complete. ",total_hits, " hits from ", total_queries," queries processed in ",t3," mins."))
  
  message("Note: if all hits for a particular OTU are removed due to filters, 
        the results will be NA for all taxon levels.
        If the lca for a particular OTU is above kingdom, e.g. cellular organisms or root, 
        the results will be unknown for all taxon levels.")
}

check.low.res.results<-function(pathofinterest,bins,btab){
  #first query otus that contributed to pathofinterest
  bins<-bins[bins$path==pathofinterest,"qseqid"]
  
  #next, query otus in btab
  btab<-btab[btab$qseqid %in% bins,]
  
  #find average pident for each contributor taxon
  contributors<-do.call(data.frame,aggregate(x = btab$pident,by=list(btab$path),
                                             FUN=function(x) c(mn = mean(x), n = range(x) )))
  colnames(contributors)<-c("contributors","mean.pident","low.pident","high.pident")
  
  #add pathofinterest for reference
  contributors$pathofinterest<-pathofinterest
  
  #add taxid for later
  btab2<-btab[!duplicated(btab$path),]
  btab3<-btab2[,c("taxids","path")]
  contributors<-merge(contributors,btab3,by.x = "contributors",by.y = "path")
  
  return(contributors)
}

check.low.res.df<-function(filtered.taxatab,filtered_blastfile, binfile,disabledTaxaFile=NULL,
                           spident=NULL,gpident=NULL,fpident=NULL,abspident=NULL,rm.excess=T,out=T){
  
  bins<-data.table::fread(binfile,sep = "\t",data.table = F)
  bins$path<-paste0(bins$K,";",bins$P,";",bins$C,";",bins$O,";",bins$F,";",bins$G,";",bins$S)
  taxatab.tf<-data.table::fread(filtered.taxatab,sep = "\t",data.table = F)
  taxatab.tf<-taxatab.tf[taxatab.tf$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]
  taxatab.tf<-taxatab.tf[taxatab.tf$taxon!="No_hits",]
  taxatab.tf<-taxatab.tf[taxatab.tf$taxon!="NA;NA;NA;NA;NA;NA;NA",]
  btab<-data.table::fread(file = filtered_blastfile,sep = "\t",data.table = F)
  btab$path<-paste0(btab$K,";",btab$P,";",btab$C,";",btab$O,";",btab$F,";",btab$G,";",btab$S)
  
  contributorlist<-list()
  for(i in 1:length(taxatab.tf$taxon)){
    contributorlist[[i]]<-check.low.res.results(pathofinterest = taxatab.tf$taxon[i],bins = bins,btab = btab)
  }
  
  contributordf<-do.call(rbind,contributorlist)
  
  #add number of reads and no. of "pcrs" pathofinterest
  taxatab.tf$readcounts<-rowSums(taxatab.tf[,2:length(colnames(taxatab.tf))])
  taxatab.tf$n.samples<-rowSums(taxatab.tf[,2:(length(colnames(taxatab.tf))-1)]>0)
  contributordf<-merge(contributordf,taxatab.tf[,c("taxon","readcounts","n.samples")],
                       by.x = "pathofinterest",by.y = "taxon")
  
  #add disabled taxa columns
  if(!is.null(disabledTaxaFile)){
    disabledTaxaDf<-read.table(disabledTaxaFile, header=T,sep = "\t")
    if(!"taxids" %in% colnames(disabledTaxaDf)) stop("No column called 'taxids'")
    if(!"disable_species" %in% colnames(disabledTaxaDf)) stop("No column called 'disable_species'")
    if(!"disable_genus" %in% colnames(disabledTaxaDf)) stop("No column called 'disable_genus'")
    if(!"disable_family" %in% colnames(disabledTaxaDf)) stop("No column called 'disable_family'")
    
    disabledSpecies<-disabledTaxaDf[disabledTaxaDf$disable_species==T,"taxids"]
    disabledSpecies<-disabledSpecies[!is.na(disabledSpecies)]
    
    disabledGenus<-disabledTaxaDf[disabledTaxaDf$disable_genus==T,"taxids"]
    disabledGenus<-disabledGenus[!is.na(disabledGenus)]
    
    disabledFamily<-disabledTaxaDf[disabledTaxaDf$disable_family==T,"taxids"]
    disabledFamily<-disabledFamily[!is.na(disabledFamily)]
    
    contributordf$species_disabled<-contributordf$taxids %in% disabledSpecies
    contributordf$genus_disabled<-contributordf$taxids %in% disabledGenus
    contributordf$family_disabled<-contributordf$taxids %in% disabledFamily
    
  } else {contributordf$species_disabled<-"FALSE"
  contributordf$genus_disabled<-"FALSE"
  contributordf$family_disabled<-"FALSE"
  }
  
  #add rank
  temprank<-stringr::str_count(contributordf$pathofinterest,";NA")
  temprank<-gsub(0,"species",temprank)
  temprank<-gsub(1,"genus",temprank)
  temprank<-gsub(2,"family",temprank)
  temprank<-gsub(3,"above_family",temprank)
  
  contributordf$rank<-temprank
  
  if(!is.null(spident) & !is.null(gpident) & !is.null(fpident) & !is.null(abspident)){
    #add may_be_improved
    max_pidents<-aggregate(contributordf$high.pident,by=list(contributordf$pathofinterest),FUN=max)
    best_is_species<-max_pidents[max_pidents$x>spident,]
    if(length(best_is_species$Group.1)>0) best_is_species$best_possible<-"species"
    max_pidents<-max_pidents[!max_pidents$x>spident,]
    best_is_genus<-max_pidents[max_pidents$x>gpident,]
    if(length(best_is_genus$Group.1)>0) best_is_genus$best_possible<-"genus"
    max_pidents<-max_pidents[!max_pidents$x>gpident,]
    best_is_family<-max_pidents[max_pidents$x>fpident,]
    if(length(best_is_family$Group.1)>0) best_is_family$best_possible<-"family"
    max_pidents<-max_pidents[!max_pidents$x>fpident,]
    best_is_above_family<-max_pidents[max_pidents$x>abspident,]
    if(length(best_is_above_family$Group.1)>0) best_is_above_family$best_possible<-"above_family"
    
    best_all<-rbind(best_is_species,best_is_genus,best_is_family,best_is_above_family)
    colnames(best_all)[1]<-"pathofinterest"
    best_all$x=NULL
    
    contributordf<-merge(contributordf,best_all)
    contributordf$may_be_improved<-contributordf$rank!=contributordf$best_possible
    contributordf<-contributordf[contributordf$may_be_improved=="TRUE",]
    contributordf<-contributordf[!contributordf$high.pident<abspident,]
    
  } else {(message("Not adding 'May be improved' column as no pidents provided"))}
  
  if(rm.excess==T){
    #if pathofinterest at genus level, dont output contributors from diferent genera
    contributordf$contributors<-as.character(contributordf$contributors)
    
    #1. make column to see whether genera match
    contributordf$contrGenus<-do.call(rbind,stringr::str_split(contributordf$contributors,";"))[,6]
    contributordf$pathGenus<-do.call(rbind,stringr::str_split(contributordf$pathofinterest,";"))[,6]
    contributordf$contr.path.genus.match<-contributordf$contrGenus==contributordf$pathGenus
    
    #2. If rank=genus, and columns dont match, then remove
    contributordf<-contributordf[!(contributordf$rank=="genus" & contributordf$contr.path.genus.match==FALSE),]
    
    #if pathofinterest at family level, dont output contributors from different families
    
    #1. make column to see whether families match
    contributordf$contrFam<-do.call(rbind,stringr::str_split(contributordf$contributors,";"))[,5]
    contributordf$pathFam<-do.call(rbind,stringr::str_split(contributordf$pathofinterest,";"))[,5]
    contributordf$contr.path.fam.match<-contributordf$contrFam==contributordf$pathFam
    
    #2. If rank=family, and columns dont match, then remove
    contributordf<-contributordf[!(contributordf$rank=="family" & contributordf$contr.path.fam.match==FALSE),]
    
    contributordf<-contributordf[order(contributordf$pathofinterest,-contributordf$high.pident),1:14]
  } else {
    
    contributordf<-contributordf[order(contributordf$pathofinterest,-contributordf$high.pident),1:12]
    
  }
  
  if(out==T){
    write.table(x = contributordf,file=gsub("spliced.txt","spliced.contr.txt",filtered.taxatab),
              sep="\t",quote = F,row.names = F)
  } else (return(contributordf))
}

#############################decided to split function
filter.blast<-function(blastfile,headers="qseqid evalue staxid pident qcovs",ncbiTaxDir,out,min_qcovs=70,
                       max_evalue=0.001,top=1){
  
  if(length(grep("qcovs",headers))<1) stop("qcovs not in headers")
  if(length(grep("evalue",headers))<1) stop("evalue not in headers")
  if(length(grep("qseqid",headers))<1) stop("qseqid not in headers")
  if(length(grep("pident",headers))<1) stop("pident not in headers")
  if(length(grep("staxid",headers))<1) stop("staxid not in headers")
  
  if(is.null(ncbiTaxDir)) stop("ncbiTaxDir not specified")
  if(is.null(out)) stop("out not specified")
  
  message("reading blast results")
  btab<-as.data.frame(data.table::fread(file = blastfile,sep = "\t"))
  colnames(btab)<-strsplit(headers,split = " ")[[1]]
  
  message("applying global min_qcovs threshold")
  btab<-btab[btab$qcovs>min_qcovs,]
  message("applying global max_evalue threshold")
  btab<-btab[btab$evalue<max_evalue,]
  message("applying global top threshold")
  if(length(btab[,1])==0) stop("No hits passing min_qcovs and max_evalue thresholds")
  topdf<-aggregate(x = btab[,colnames(btab) %in% c("qseqid","pident")],by=list(btab$qseqid),FUN = max)
  topdf$min_pident<-topdf$pident-top
  btab<-merge(btab,topdf[,c(2,4)],by = "qseqid", all.y = T) #can definitely do this differently and faster
  btab<-btab[btab$pident>btab$min_pident,]
  
  #add lineage to results
  message("adding taxonomic lineages")
  btab$taxids<-btab$staxid #add.lineage.df requires this colname
  btab<-add.lineage.df(btab,ncbiTaxDir)
  
  #remove crappy hits 
  #1. btab$S contains uncultured
  message("Removing species containing the terms: uncultured, environmental, 
          unidentified,fungal, eukaryote or unclassified")
  if(length(grep("uncultured",btab$S,ignore.case = T))>0) {
    btab<-btab[-grep("uncultured",btab$S,ignore.case = T),]}
  #2. btab$S contains environmental
  if(length(grep("environmental",btab$S,ignore.case = T))>0) {
    btab<-btab[-grep("environmental",btab$S,ignore.case = T),]}
  #3. btab$S contains unclassified
  if(length(grep("unclassified",btab$S,ignore.case = T))>0) {
    btab<-btab[-grep("unclassified",btab$S,ignore.case = T),]}
  #4. btab$S contains "unidentified"
  if(length(grep("unidentified",btab$S,ignore.case = T))>0) {
    btab<-btab[-grep("unidentified",btab$S,ignore.case = T),]}
  #5. btab$S contains "fungal "
  if(length(grep("fungal ",btab$S,ignore.case = T))>0) {
    btab<-btab[-grep("fungal ",btab$S,ignore.case = T),]}
  #6. btab$S contains "eukaryote"
  if(length(grep("eukaryote",btab$S,ignore.case = T))>0) {
    btab<-btab[-grep("eukaryote",btab$S,ignore.case = T),]}
  
  write.table(x = btab,file = out,sep="\t",quote = F,row.names = F)
  
}


merge.and.check.disabled.taxa.files<-function(disabledTaxaFiles,disabledTaxaOut,force=F,full.force=F){
  
  message("Note: Use force=T to ignore any contributor entries where no levels were disabled when consistency checking.
                 Use force=F to do more thorough consistency checks")
  
  disabledTaxaDFList<-list()
  
  for(i in 1:length(disabledTaxaFiles)){
    
    disabledTaxaDFList[[i]]<-data.table::fread(disabledTaxaFiles[i], data.table = F,sep = "\t")
    
    if(!"taxids" %in% colnames(disabledTaxaDFList[[i]]))  stop(paste("No column called 'taxids' in", disabledTaxaFiles[i]))
    if(!"disable_species" %in% colnames(disabledTaxaDFList[[i]])) stop(paste("No column called 'disable_species' in", disabledTaxaFiles[i]))
    if(!"disable_genus" %in% colnames(disabledTaxaDFList[[i]])) stop(paste("No column called 'disable_genus' in", disabledTaxaFiles[i]))
    if(!"disable_family" %in% colnames(disabledTaxaDFList[[i]])) stop(paste("No column called 'disable_family' in", disabledTaxaFiles[i]))
    if(!"contributors" %in% colnames(disabledTaxaDFList[[i]])) stop(paste("No column called 'contributors' in", disabledTaxaFiles[i]))
    
    disabledTaxaDFList[[i]]<-disabledTaxaDFList[[i]][,c("contributors","taxids","disable_species","disable_genus","disable_family")]
    disabledTaxaDFList[[i]]$file<-disabledTaxaFiles[i]
  }
  
  disabledTaxaDF<-do.call(rbind,disabledTaxaDFList)
  
  #make empty cells FALSE
  disabledTaxaDF$disable_species<-as.logical(disabledTaxaDF$disable_species)
  disabledTaxaDF$disable_genus<-as.logical(disabledTaxaDF$disable_genus)
  disabledTaxaDF$disable_family<-as.logical(disabledTaxaDF$disable_family)
  disabledTaxaDF[,c(-1,-2)][is.na(disabledTaxaDF[,c(-1,-2)])]<-FALSE

  disabledTaxaDF<-add.lineage.df(disabledTaxaDF,ncbiTaxDir,as.taxids = T)
  
  if(full.force) force=T
  
  if(force) {
    message("Using force=T")
    disabledTaxaDF$trues<-rowSums(disabledTaxaDF[,c("disable_species","disable_genus","disable_family")])
    disabledTaxaDF<-disabledTaxaDF[disabledTaxaDF$trues>0,-match("trues",colnames(disabledTaxaDF))]
    
    if(full.force) {
      message("Using full.force=T")
      
      #change families marked as TRUE to TRUE 
      splitdf<-split(disabledTaxaDF, f = as.factor(disabledTaxaDF$F))
      
      for (i in 1:(length(splitdf))) {
        if(TRUE %in% splitdf[[i]]$disable_family) splitdf[[i]]$disable_family<-TRUE
      }
      
      disabledTaxaDF<-do.call(rbind, splitdf)
      
      #change genus marked as TRUE to TRUE 
      splitdf<-split(disabledTaxaDF, f = as.factor(disabledTaxaDF$G))
      
      for (i in 1:(length(splitdf))) {
        if(TRUE %in% splitdf[[i]]$disable_genus) splitdf[[i]]$disable_genus<-TRUE
      }
      
      disabledTaxaDF<-do.call(rbind, splitdf)
      
      #change species marked as TRUE to TRUE 
      splitdf<-split(disabledTaxaDF, f = as.factor(disabledTaxaDF$S))
      
      for (i in 1:(length(splitdf))) {
        if(TRUE %in% splitdf[[i]]$disable_species) splitdf[[i]]$disable_species<-TRUE
      }
      
      disabledTaxaDF<-do.call(rbind, splitdf)
    }
      
  } else { message("Using force=F")}
  
  
  #check that identical paths have not been treated differently
  shouldstop1<-list()
  for(i in 1:length(unique(disabledTaxaDF$contributors))){
    temp<-disabledTaxaDF[disabledTaxaDF$contributors==unique(disabledTaxaDF$contributors)[i],]
    if(length(temp$contributors)>1) if(sum(duplicated(temp[,c("disable_species","disable_genus","disable_family")]))!=
                                       length(temp$contributors)-1) {
      message("inconsistent taxid disabling detected")
      print(temp)
      shouldstop1[[i]]<-temp
    }
  }
  
  #check that identical families have not been treated differently
  shouldstop2<-list()
  for(i in 1:length(unique(disabledTaxaDF$F))){
    temp<-disabledTaxaDF[disabledTaxaDF$F==unique(disabledTaxaDF$F)[i],]
    if(length(temp$contributors)>1) if(sum(duplicated(temp[,"disable_family"]))!=length(temp$contributors)-1) {
      message("inconsistent family disabling detected:")
      print(temp)
      shouldstop2[[i]]<-temp
    }
  }
  
  #check that identical genera have not been treated differently
  shouldstop3<-list()
  for(i in 1:length(unique(disabledTaxaDF$G))){
    temp<-disabledTaxaDF[disabledTaxaDF$G==unique(disabledTaxaDF$G)[i],]
    if(length(temp$contributors)>1) if(sum(duplicated(temp[,"disable_genus"]))!=length(temp$contributors)-1) {
      message("inconsistent genus disabling detected")
      print(temp)
      shouldstop3[[i]]<-temp
    }
  }
  
  #check that identical species have not been treated differently
  shouldstop4<-list()
  for(i in 1:length(unique(disabledTaxaDF$S))){
    temp<-disabledTaxaDF[disabledTaxaDF$S==unique(disabledTaxaDF$S)[i],]
    if(length(temp$contributors)>1) if(sum(duplicated(temp[,"disable_species"]))!=length(temp$contributors)-1) {
      message("inconsistent species disabling detected")
      print(temp)
      shouldstop4[[i]]<-temp
    }
  }
  
  if(length(shouldstop1)>0) for(i in 1:length(shouldstop1)) if(!is.null(shouldstop1[[i]])) stop("fix inconsistencies")
  if(length(shouldstop2)>0) for(i in 1:length(shouldstop2)) if(!is.null(shouldstop2[[i]])) stop("fix inconsistencies")
  if(length(shouldstop3)>0) for(i in 1:length(shouldstop3)) if(!is.null(shouldstop3[[i]])) stop("fix inconsistencies")
  if(length(shouldstop4)>0) for(i in 1:length(shouldstop4)) if(!is.null(shouldstop4[[i]])) stop("fix inconsistencies")
  
  #check that if genus disabled, species within that genus also should be 
  temp<-disabledTaxaDF[(disabledTaxaDF$disable_genus==T & disabledTaxaDF$disable_species==F),"contributors"]
    if(length(temp)>0) {
      print(temp)
      stop("Genus was disabled for this contributor, but species was not")
    }
  
  #check that if family disabled, species within that family also should be 
  temp<-disabledTaxaDF[(disabledTaxaDF$disable_family==T & disabledTaxaDF$disable_species==F),"contributors"]
  if(length(temp)>0) {
    print(temp)
    stop("Family was disabled for this contributor, but species was not")
  }
  
  #check that if family disabled, genus within that family also should be 
  temp<-disabledTaxaDF[(disabledTaxaDF$disable_family==T & disabledTaxaDF$disable_genus==F),"contributors"]
  if(length(temp)>0) {
    print(temp)
    stop("Family was disabled for this contributor, but genus was not")
  }
    
  write.table(disabledTaxaDF,disabledTaxaOut,col.names = T,row.names = F,quote = F,sep = "\t")
  
  return(disabledTaxaDF)
}

bin.blast.lite<-function(filtered_blastfile,ncbiTaxDir,obitaxdb,out){
  
  t1<-Sys.time()
  
  if(is.null(out)) stop("out not specified")
  if(is.null(ncbiTaxDir)) stop("ncbiTaxDir not specified")
  if(is.null(obitaxdb)) stop("obitaxdb not specified")
  
  ###################################################
  #read in obitools taxonomy
  obitaxdb2=ROBITaxonomy::read.taxonomy(obitaxdb)
  
  btab<-data.table::fread(filtered_blastfile,sep="\t",data.table = F)
  
  #preparing some things for final step
  total_hits<-length(btab$qseqid) #for info later
  total_queries<-length(unique(btab$qseqid))
  qseqids<-as.data.frame(unique(btab$qseqid))
  qseqids$qseqid<-qseqids$`unique(btab$qseqid)`
  qseqids$`unique(btab$qseqid)`=NULL
  
  message("binning")
  btabsp<-btab[btab$S!="unknown",]
  
  message("Considering species with 'sp.', numbers or more than one space")
  
  if(length(btabsp$taxids)>0){
    lcasp = aggregate(btabsp$taxids, by=list(btabsp$qseqid),function(x) ROBITaxonomy::lowest.common.ancestor(obitaxdb2,x))
    
    #get lca names
    colnames(lcasp)<-gsub("x","taxids",colnames(lcasp))
    if(sum(is.na(lcasp$taxids))>0){
      message("************
              ERROR: Some taxids were not recognized by ROBITaxonomy::lowest.common.ancestor, probably need to update obitaxdb using NCBI2obitaxonomy
              *************")
      problem.contributors<-btabsp[!duplicated(btabsp[lcasp[is.na(lcasp$taxids),1] %in% btabsp$qseqid,c("taxids","K","P","C","O","F","G","S")]),
                                   c("taxids","K","P","C","O","F","G","S")]
      for(i in 1:length(problem.contributors$taxids)){
        problem.contributors$validate[i]<-ROBITaxonomy::validate(obitaxdb2,problem.contributors$taxids[i])
      }
      print(problem.contributors[is.na(problem.contributors$validate),])
    }
    lcasp<-add.lineage.df(df = lcasp,ncbiTaxDir)
    colnames(lcasp)<-gsub("Group.1","qseqid",colnames(lcasp))
  } else {
    lcasp<-data.frame(matrix(nrow=1,ncol = 10))
    colnames(lcasp)<-c("taxids","qseqid","old_taxids","K","P","C","O","F","G","S")
  }
  
  
  ###################################################
  
  com_level<-merge(x=qseqids, y = lcasp[,c(2,4:10)], by = "qseqid",all = T)
  
  com_level[com_level=="unknown"]<-NA
  
  #info
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  
  write.table(x = com_level,file = out,sep="\t",quote = F,row.names = F)
  
  message(c("Complete. ",total_hits, " hits from ", total_queries," queries processed in ",t3," mins."))
  
  message("Note: if all hits for a particular OTU are removed due to filters, 
          the results will be NA for all taxon levels.
          If the lca for a particular OTU is above kingdom, e.g. cellular organisms or root, 
          the results will be unknown for all taxon levels.")
}

