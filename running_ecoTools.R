#' Find primers to amplify a group of taxa
#' @title Find primers to amplify a group of taxa.
#' @param ecopcrdb An ecopcrdb containing the sequences to be tested
#' @param max_error Max number of mismatches allowed
#' @param min_length Minimum amplicon length (without primers - I think)
#' @param max_length Maximum amplicon length (without primers - I think)
#' @param strict_match_p Proportion of sequences with identical primer match
#' @param e_match_p Proportion of sequences with primer match within \code{e}
#' @return Default: A list containing the following dataframes
#' \itemize{
#'     \item \code{primer_table}: A table with the following columns
#'         \itemize{
#'             \item \code{PID} Primer ID
#'             \item \code{Pf} Forward primer sequence
#'             \item \code{Pr} Reverse primer seqeunce
#'             \item \code{TmPf} Melting temperature forward primer (without mismatches)
#'             \item \code{MinTmPf} Minimum melting temperature forward primer
#'             \item \code{TmPr} Melting temperature reverse primer (without mismatches)
#'             \item \code{MinTmPr} Minimum melting temperature reverse primer
#'             \item \code{CGPf} Number of Cs or Gs in forward primer
#'             \item \code{CGPr} Number of Cs or Gs in reverse primer
#'             \item \code{Specificity} GG (Good-Good) means that both primer are specific to the target dataset,
#'                 GB or BG (Good-Bad or Bad-Good) means that only one of the two primers is specific to the target dataset
#'             \item \code{Amp_Records} Number of records in the target dataset that are amplified
#'             \item \code{P_Records} Proportion of records in the target dataset that are amplified
#'             \item \code{Amp_taxa} Number of taxa in the target dataset that are amplified
#'             \item \code{Amp_NT_taxa} Number of taxa in the non-target dataset that are amplified
#'             \item \code{Taxa_unique} Number of taxa with unique amplicons
#'             \item \code{P_taxa_unique} Proportion of taxa with unique amplicons
#'             \item \code{Min_length} Minimum amplicon length (excluding primers)
#'             \item \code{Max_length} Maximum amplicon length (excluding primers)
#'             \item \code{Mean_length} Mean amplicon length (excluding primers)}
#'     \item \code{metadata}: Parameters used during ecoPrimers command}
#' @note Specificty and Amp_NT_taxa have meaningless results unless a non-target dataset was specifiec during ecoPrimers command
#´
#' @examples
#' a<-"/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/IRANVERTS/TESTING_PRIMERS_FOR_PREY/ovids_cytB/taxids.cytb_ovids_5-4-19.ecopcrdb"
#' b<-ecoPrimers.Bas(a,max_error = 1,min_length = 100,max_length = 300,strict_match_p = 0.7,e_match_p = 0.8)
#' @export
ecoPrimers.Bas<-function(ecopcrdb,max_error,min_length,max_length,strict_match_p,e_match_p,p3=T){
  cb <- function(line, proc) {cat(line, "\n")}
  a<-processx::run(command = "ecoPrimers",
                   args=c("-d",ecopcrdb,"-e",max_error,"-l",min_length,"-L", max_length,"-q",strict_match_p,"-s",e_match_p,"-3",p3),
                   echo=F,echo_cmd = T,stderr_line_callback = cb)

  f<-max(stringr::str_locate_all(string = a$stderr,pattern = ": ")[[1]][,2])
  g<-stringr::str_length(string = a$stderr)
  h<-as.numeric(substr(a$stderr,start = f+1,stop = g-1))
  message(paste("Primer search complete.",h,"primers found"))
  if(!h>0) stop("stop", call. = F)

  b<-suppressWarnings(data.table::fread(input = a$stdout,sep = "\t"))
  ecoprimeroutput<-as.data.frame(b[,1:21])
  colnames(ecoprimeroutput)<- c("PID","Pf","Pr","TmPf","MinTmPf","TmPr","MinTmPr","CGPf","CGPr","Specificity","Amp_Records","IGNORE","P_Records",
                                "Amp_taxa","IGNORE2","Amp_NT_taxa","Taxa_unique","P_taxa_unique","Min_length","Max_length","Mean_length")
  ecoprimeroutput<-ecoprimeroutput[,-grep("IGNORE",colnames(ecoprimeroutput))]
  max(stringr::str_locate_all(string = a$stdout,pattern = "#")[[1]][,1])
  d<-substr(x = a$stdout,start = 1,stop=max(stringr::str_locate_all(string = a$stdout,pattern = "#")[[1]][,1]))
  ecoprimeroutput_META<-data.table::fread(input = d,sep = "\n")
  output_list <- list(ecoprimeroutput_META, ecoprimeroutput)
  names(output_list)<-c("metadata","primer_table")
  return(output_list)
}

#' Perform in silico PCR
#' @title Perform in silico PCR.
#' @param ecopcrdb An ecopcrdb containing the sequences to be tested
#' @param max_error Max number of mismatches allowed
#' @param min_length Minimum amplicon length (without primers - I think)
#' @param max_length Maximum amplicon length (without primers - I think)
#' @return A dataframe with in silico PCR results
#'
#' @export
ecoPCR.Bas<-function(Pf,Pr,ecopcrdb,max_error,min_length,max_length,out,buffer=NULL){
  #cb <- function(line, proc) {cat(line, "\n")}
  
  if(length(grep("I",Pf))>0)(Pf<-gsub("I","N",Pf))
  if(length(grep("I",Pr))>0)(Pf<-gsub("I","N",Pr))
  
  if(is.null(buffer)){
  system2(command = "ecoPCR",args=c(Pf, Pr,"-d",ecopcrdb,"-e",max_error,"-l",min_length,"-L", max_length,"-k","-c"),
          stdout = out,wait = T)}
  
  if(!is.null(buffer)){
    system2(command = "ecoPCR",
            args=c(Pf, Pr,"-d",ecopcrdb,"-e",max_error,"-l",min_length,"-L", max_length,"-k","-c","-D",buffer),
            stdout = out,wait = T)}
  
  b<-suppressWarnings(data.table::fread(input = out,sep = "|"))
  b<-as.data.frame(b)
  if(length(colnames(b))==1) stop("No hits",call. = F)
  colnames(b)<- c("AC","seq_length","taxid","rank","species","species_name","genus",
                                "genus_name","family","family_name","superkingdom","superkingdom_name",
                                "strand","forward_match","forward_mismatch","forward_tm","reverse_match",
                             "reverse_mismatch","reverse_tm","amplicon_length","sequence","definition")
  write.table(b,file = out,quote = F,row.names = F,sep = "\t")
  }
