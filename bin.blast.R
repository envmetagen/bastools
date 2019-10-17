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
bin.blast<-function(blastfile,headers="qseqid evalue staxid pident qcovs",ncbiTaxDir,obitaxdb,out,min_qcovs=70,
                        max_evalue=0.001,top=1,
                        spident=98,gpident=95,fpident=92,abspident=80){
  
  t1<-Sys.time()

  if(length(grep("qcovs",headers))<1) stop("qcovs not in headers")
  if(length(grep("evalue",headers))<1) stop("evalue not in headers")
  if(length(grep("qseqid",headers))<1) stop("qseqid not in headers")
  if(length(grep("pident",headers))<1) stop("pident not in headers")
  if(length(grep("staxid",headers))<1) stop("staxid not in headers")
  
  if(is.null(ncbiTaxDir)) stop("ncbiTaxDir not specified")
  if(is.null(obitaxdb)) stop("obitaxdb not specified")
  if(is.null(out)) stop("out not specified")
  
  message("reading blast results")
  btab<-as.data.frame(data.table::fread(file = blastfile,sep = "\t"))
  colnames(btab)<-strsplit(headers,split = " ")[[1]]
  
  #preparing some things for final step
  total_hits<-length(btab$qseqid) #for info later
  total_queries<-length(unique(btab$qseqid))
  qseqids<-as.data.frame(unique(btab$qseqid))
  qseqids$qseqid<-qseqids$`unique(btab$qseqid)`
  qseqids$`unique(btab$qseqid)`=NULL

  message("applying global min_qcovs threshold")
  btab<-btab[btab$qcovs>min_qcovs,]
  message("applying global max_evalue threshold")
  btab<-btab[btab$evalue<max_evalue,]
  message("applying global top threshold")
  if(length(btab[,1])==0) stop("No hits passing min_qcovs and max_evalue thresholds")
  topdf<-aggregate(x = btab[,colnames(btab) %in% c("qseqid","pident")],by=list(btab$qseqid),FUN = max)
  topdf$min_pident<-topdf$pident-topdf$pident*top/100
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
  
  
  ###################################################
  #read in obitools taxonomy
  obitaxdb2=ROBITaxonomy::read.taxonomy(obitaxdb)
  
  #species-level assignments
  message("binning at species level")
  message("Removing species with sp., numbers or more than one space")
  btabsp<-btab[btab$S!="unknown",]
  if(length(grep(" sp\\.",btabsp$S,ignore.case = T))>0) btabsp<-btabsp[-grep(" sp\\.",btabsp$S,ignore.case = T),]
  if(length(grep(" .* .*",btabsp$S,ignore.case = T))>0) btabsp<-btabsp[-grep(" .* .*",btabsp$S,ignore.case = T),]
  if(length(grep("[0-9]",btabsp$S))>0) btabsp<-btabsp[-grep("[0-9]",btabsp$S),]
  btabsp<-btabsp[btabsp$pident>spident,]
  if(length(btabsp$taxids)>0){
  lcasp = aggregate(btabsp$taxids, by=list(btabsp$qseqid),
                    function(x) ROBITaxonomy::lowest.common.ancestor(obitaxdb2,x))
  rm(btabsp)
  #get lca names
  colnames(lcasp)<-gsub("x","taxids",colnames(lcasp))
  lcasp<-add.lineage.df(df = lcasp,ncbiTaxDir)
  colnames(lcasp)<-gsub("Group.1","qseqid",colnames(lcasp))
  } else {
    lcasp<-data.frame(matrix(nrow=1,ncol = 10))
    colnames(lcasp)<-c("taxids","qseqid","old_taxids","K","P","C","O","F","G","S")
    }
  #genus-level assignments
  message("binning at genus level") ######for all the next steps I should only process qseqids that were not yet
  ######successful, wasting time at the moment
  btabg<-btab[btab$G!="unknown" & btab$S!="unknown",]
  btabg<-btabg[btabg$pident>gpident,]
  lcag = aggregate(btabg$taxids, by=list(btabg$qseqid),
                    function(x) ROBITaxonomy::lowest.common.ancestor(obitaxdb2,x))
  rm(btabg)
  #get lca names
  colnames(lcag)<-gsub("x","taxids",colnames(lcag))
  lcag<-add.lineage.df(df = lcag,ncbiTaxDir)
  colnames(lcag)<-gsub("Group.1","qseqid",colnames(lcag))

  #family-level assignments
  message("binning at family level")
  btabf<-btab[btab$F!="unknown" & btab$G!="unknown" & btab$S!="unknown",]
  btabf<-btabf[btabf$pident>fpident,]
  lcaf = aggregate(btabf$taxids, by=list(btabf$qseqid),
                   function(x) ROBITaxonomy::lowest.common.ancestor(obitaxdb2,x))
  rm(btabf)
  #get lca names
  colnames(lcaf)<-gsub("x","taxids",colnames(lcaf))
  lcaf<-add.lineage.df(df = lcaf,ncbiTaxDir)
  colnames(lcaf)<-gsub("Group.1","qseqid",colnames(lcaf))

  #higher-than-family-level assignments
  message("binning at higher-than-family level")
  btabhtf<-btab[btab$K!="unknown",]
  btabhtf<-btabhtf[btabhtf$pident>abspident,]
  lcahtf = aggregate(btabhtf$taxids, by=list(btabhtf$qseqid),
                   function(x) ROBITaxonomy::lowest.common.ancestor(obitaxdb2,x))
  rm(btabhtf)
  #get lca names
  colnames(lcahtf)<-gsub("x","taxids",colnames(lcahtf))
  lcahtf<-add.lineage.df(df = lcahtf,ncbiTaxDir)
  colnames(lcahtf)<-gsub("Group.1","qseqid",colnames(lcahtf))

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
  com_level<-merge(x=qseqids, y = com_level[,c(2,4:10)], by = "qseqid",all = T)

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

check.low.res.df<-function(filtered.taxatab,filtered_blastfile, binfile){

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
write.table(x = contributordf,file=gsub("spliced.txt","spliced.contr.txt",filtered.taxatab),
            sep="\t",quote = F,row.names = F)
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
  if(is.null(obitaxdb)) stop("obitaxdb not specified")
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
  topdf$min_pident<-topdf$pident-topdf$pident*top/100
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

bin.blast2<-function(filtered_blastfile,headers="qseqid evalue staxid pident qcovs",ncbiTaxDir,
                     obitaxdb,out,spident=98,gpident=95,fpident=92,abspident=80){
t1<-Sys.time()
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

#species-level assignments
message("binning at species level")
message("Removing species with sp., numbers or more than one space")
btabsp<-btab[btab$S!="unknown",]
if(length(grep(" sp\\.",btabsp$S,ignore.case = T))>0) btabsp<-btabsp[-grep(" sp\\.",btabsp$S,ignore.case = T),]
if(length(grep(" .* .*",btabsp$S,ignore.case = T))>0) btabsp<-btabsp[-grep(" .* .*",btabsp$S,ignore.case = T),]
if(length(grep("[0-9]",btabsp$S))>0) btabsp<-btabsp[-grep("[0-9]",btabsp$S),]
btabsp<-btabsp[btabsp$pident>spident,]
if(length(btabsp$taxids)>0){
  lcasp = aggregate(btabsp$taxids, by=list(btabsp$qseqid),
                    function(x) ROBITaxonomy::lowest.common.ancestor(obitaxdb2,x))
  rm(btabsp)
  #get lca names
  colnames(lcasp)<-gsub("x","taxids",colnames(lcasp))
  lcasp<-add.lineage.df(df = lcasp,ncbiTaxDir)
  colnames(lcasp)<-gsub("Group.1","qseqid",colnames(lcasp))
} else {
  lcasp<-data.frame(matrix(nrow=1,ncol = 10))
  colnames(lcasp)<-c("taxids","qseqid","old_taxids","K","P","C","O","F","G","S")
}
#genus-level assignments
message("binning at genus level") ######for all the next steps I should only process qseqids that were not yet
######successful, wasting time at the moment
btabg<-btab[btab$G!="unknown" & btab$S!="unknown",]
btabg<-btabg[btabg$pident>gpident,]
lcag = aggregate(btabg$taxids, by=list(btabg$qseqid),
                 function(x) ROBITaxonomy::lowest.common.ancestor(obitaxdb2,x))
rm(btabg)
#get lca names
colnames(lcag)<-gsub("x","taxids",colnames(lcag))
lcag<-add.lineage.df(df = lcag,ncbiTaxDir)
colnames(lcag)<-gsub("Group.1","qseqid",colnames(lcag))

#family-level assignments
message("binning at family level")
btabf<-btab[btab$F!="unknown" & btab$G!="unknown" & btab$S!="unknown",]
btabf<-btabf[btabf$pident>fpident,]
lcaf = aggregate(btabf$taxids, by=list(btabf$qseqid),
                 function(x) ROBITaxonomy::lowest.common.ancestor(obitaxdb2,x))
rm(btabf)
#get lca names
colnames(lcaf)<-gsub("x","taxids",colnames(lcaf))
lcaf<-add.lineage.df(df = lcaf,ncbiTaxDir)
colnames(lcaf)<-gsub("Group.1","qseqid",colnames(lcaf))

#higher-than-family-level assignments
message("binning at higher-than-family level")
btabhtf<-btab[btab$K!="unknown",]
btabhtf<-btabhtf[btabhtf$pident>abspident,]
lcahtf = aggregate(btabhtf$taxids, by=list(btabhtf$qseqid),
                   function(x) ROBITaxonomy::lowest.common.ancestor(obitaxdb2,x))
rm(btabhtf)
#get lca names
colnames(lcahtf)<-gsub("x","taxids",colnames(lcahtf))
lcahtf<-add.lineage.df(df = lcahtf,ncbiTaxDir)
colnames(lcahtf)<-gsub("Group.1","qseqid",colnames(lcahtf))

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
com_level<-merge(x=qseqids, y = com_level[,c(2,4:10)], by = "qseqid",all = T)

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