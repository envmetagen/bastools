#' @import ROBITaxonomy 
#' @include 02_class_metabarcoding.data.R
NULL

#' Choose between databases for taxonomic classifications
#' 
#' Chooses a sequence taxonomic assignment in order of preference for the different 
#' reference databases that have been used when the assignment is above a certain threshold
#' 
#' 
#' @param x a \code{\link{metabarcoding.data}} object
#' @param taxonomy a \code{\linkS4class{taxonomy.obitools}} instance
#' @param dbrank string or vector indicating reference database names ranked by order of preference
#' @param thresh a best_identity threshold for applying priority. Default is \code{0.95}
#' 
#' @return returns a data.frame with the refined taxonomic assignement and classic taxonomy description.
#' 
#' @examples
#' 
#' data(termes)
#' 
#' taxo=default.taxonomy()
#' 
#' #create artificial taxonomic assignments
#' attr(termes, "motus")["best_identity:DB1"] = sample(seq(0.5,1,0.001),size=nrow(termes$motus), replace=T)
#' attr(termes, "motus")["best_identity:DB2"] = sample(seq(0.5,1,0.001),size=nrow(termes$motus), replace=T)
#' attr(termes, "motus")["best_identity:DB3"] = sample(seq(0.5,1,0.001),size=nrow(termes$motus), replace=T)
#' attr(termes, "motus")["taxid_by_db:DB1"] = termes$motus$taxid
#' attr(termes, "motus")["taxid_by_db:DB2"] = sample(termes$motus$taxid,size=nrow(termes$motus), replace=F)
#' attr(termes, "motus")["taxid_by_db:DB3"] = sample(termes$motus$taxid,size=nrow(termes$motus), replace=F)
#' 
#' #Run taxo.decider
#' termes.ok = taxo.decider(termes, taxo, "DB2", 0.95)
#' head(termes.ok$motus[union(grep("DB",  colnames(termes.ok$motus)), grep("_ok", colnames(termes.ok$motus)))])
#' 
#' termes.ok = taxo.decider(termes, taxo, c("DB3", "DB1"), 0.95)
#' head(termes.ok$motus[union(grep("DB",  colnames(termes.ok$motus)), grep("_ok", colnames(termes.ok$motus)))])
#' 
#' #Quick look at the enhancement in taxonomic assignements
#' par(mfrow=c(1,4))
#' for(i in grep("best_identity.", colnames(termes.ok$motus))){
#' hist(termes.ok$motus[,i], breaks=20, ylim=c(1,21), main=colnames(termes.ok$motus)[i], xlab="assignment score")
#' }
#'     
#' @seealso \code{\linkS4class{taxonomy.obitools}}, and methods \code{\link{species}},\code{\link{genus}}, \code{\link{family}},\code{\link{kingdom}},
#'          \code{\link{superkingdom}},\code{\link{taxonatrank}}, \code{\link{taxonmicank}}
#'
#' @author Lucie Zinger
#' @keywords taxonomy
#' 
#' @export
#' 

taxo.decider = function(x, taxonomy, dbrank, thresh=0.95) {
  
  noms = colnames(x$motus)
  best_ids_names = noms[grep("best_identity.", noms)]
  best_ids = x$motus[,best_ids_names]
  taxids = x$motus[, gsub("best_identity", "taxid_by_db", best_ids_names)]
  dbs = unlist(lapply(strsplit(best_ids_names, "\\:"), "[[", 2))
  
  
  #Set max indices
  ind = as.vector(t(apply(best_ids,1,function(y) order(rank(-y, ties.method="max"), match(dbrank, dbs))))[,1])
  
  #Set default vector: db, bestids, taxids with max score
  db_ok = dbs[ind]
  best_identity_ok = best_ids[cbind(1:length(ind), ind)]
  taxids_by_db_ok = taxids[cbind(1:length(ind), ind)]
  
  #Get vector of db index that should be used according to condition > thresh
  db_choice = taxo.decider.routine(dbrank, best_ids, dbs, thresh)
  
  #Replacing by right values according to db_ok
  for(i in 1:length(dbrank)){
    db_ok[which(db_choice==i)] = dbrank[i]
    best_identity_ok[which(db_choice==i)] = best_ids[which(db_choice==i),grep(dbrank[i], colnames(best_ids))]
    taxids_by_db_ok[which(db_choice==i)] = taxids[which(db_choice==i),grep(dbrank[i], colnames(taxids))]
  }
  
  decision = data.frame(db_ok, best_identity_ok, taxids_by_db_ok)
  
  coltaxid = colnames(decision)[grep("taxid", colnames(decision))]
  
  attr(x, "motus") = data.frame(x$motus, decision)
  new.tax = get.classic.taxonomy(x, taxonomy, coltaxid)
  
  attr(x, "motus") = data.frame(x$motus, new.tax)
  
  return(x)  
}


taxo.decider.routine = function(dbrank, best_ids, dbs, thresh) {
  #Setting mask 
  mask = matrix(NA,nrow(best_ids),length(dbrank))
  colnames(mask)=dbrank
  #For each DB, see if condition T/F
  for(i in dbrank){
    mask[,i] = best_ids[,which(dbs==i)]>thresh
  }
  #Get the first occurence of T in the table
  out = apply(mask, 1, function(x) which(x==T)[1])
  return(out)
}


