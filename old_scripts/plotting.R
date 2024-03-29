#' #' Some colours I like to use for plots (n=29)
#' #'@export
#' MyCols <- c("dodgerblue2","#E31A1C", # red
#'             "green4",
#'             "#6A3D9A", # purple
#'             "#FF7F00", # orange
#'             "black","gold1",
#'             "skyblue2","#FB9A99", # lt pink
#'             "palegreen2",
#'             "#CAB2D6", # lt purple
#'             "#FDBF6F", # lt orange
#'             "gray70", "khaki2",
#'             "maroon","orchid1","deeppink1","blue1","steelblue4",
#'             "darkturquoise","green1","yellow4","yellow3",
#'             "darkorange4","brown","chartreuse2","coral4","darkslategray1","mediumseagreen")

#
#' Plot output from function ecopcr.hit.table
#' @title Plot output from function ecopcr.hit.table
#' @param x A dataframe in the format corresponding to "HT" or "HT_PC", as output by \code{\link[bastools]{ecopcr.hit.table}}.
#' @param cols29 A logical on whether to use a preselected, easy to read, colour palette (29 colours maximum)
#' @return A ggplot of the number (if input dataframe is of type "HT") or percentage (if input dataframe is of type "HT_PC") of species amplified
#'     for for each taxonomic class at each number of primer mismatches.
#' @note Trying to plot over 100 taxa will result in a warning.
#' @examples
#' ecopcrTidyAll<-ecopcr.hit.table("C:/Users/basti/Documents/WORK/CIBIO/AA_PROJECTS/MINION/45F-63R_4Mis_ALLVERTS_REFSEQ.ecopcroutput",
#'    "C:/Users/basti/Documents/WORK/CIBIO/STATS_AND_CODE/OBITOOLS/ecoPCRoutput/USING 16S/taxdump/taxdump",
#'    "C:/Users/basti/Documents/WORK/CIBIO/AA_PROJECTS/MINION/ALL_VERTS_REFSEQ_MTDNA_1Kbp-COX1-1Kbp.taxids.ecopcrdb.tab",obitaxdbisinR=F)
#'
#'    plot.ecopcr.hit.table(ecopcrTidyAll$HT)
#'    plot.ecopcr.hit.table(ecopcrTidyAll$HT_pc)
#'@export
basplot.ecopcr.hit.table <- function(ecopcr.hit.table.output,cols29=F,simplify=T){
  if(length(colnames(ecopcr.hit.table.output))>100) {print("Trying to plot over 100 taxa could take a while, and may be difficult to read")}
  if(simplify==T){
    if(cols29==T){
      d<-ecopcr.hit.table.output[,c("DB","Total")]
      a<-reshape2::melt(d)
      b<-a[a$DB=="0 mismatches" | a$DB=="1 mismatch" | a$DB=="2 mismatches"| a$DB=="3 mismatches"| a$DB=="4 mismatches"| a$DB=="5 mismatches" | a$DB=="6 mismatches",]
      e<-ggplot2::ggplot(data=b, aes(y=value, x=DB, group=variable)) + geom_point(aes(colour=variable)) +  geom_line(aes(colour=variable))+
        labs(y="Species amplified", x="") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_colour_manual(values = MyCols)
    }
    if(cols29==F){
      d<-ecopcr.hit.table.output[,c("DB","Total")]
      a<-reshape2::melt(d)
      b<-a[a$DB=="0 mismatches" | a$DB=="1 mismatch" | a$DB=="2 mismatches"| a$DB=="3 mismatches"| a$DB=="4 mismatches"| a$DB=="5 mismatches" | a$DB=="6 mismatches",]
      e<-ggplot2::ggplot(data=b, aes(y=value, x=DB, group=variable)) + geom_point(aes(colour=variable)) +  geom_line(aes(colour=variable))+
        labs(y="Species amplified", x="") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  }
  if(simplify==F){
    if(cols29==T){
      a<-reshape2::melt(ecopcr.hit.table.output)
      b<-a[a$DB=="0 mismatches" | a$DB=="1 mismatch" | a$DB=="2 mismatches"| a$DB=="3 mismatches"| a$DB=="4 mismatches"| a$DB=="5 mismatches" | a$DB=="6 mismatches",]
      e<-ggplot2::ggplot(data=b, aes(y=value, x=DB, group=variable)) + geom_point(aes(colour=variable)) +  geom_line(aes(colour=variable))+
        labs(y="Species amplified", x="") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_colour_manual(values = MyCols)
    }

    if(cols29==F){
      a<-reshape2::melt(ecopcr.hit.table.output)
      b<-a[a$DB=="0 mismatches" | a$DB=="1 mismatch" | a$DB=="2 mismatches"| a$DB=="3 mismatches"| a$DB=="4 mismatches"| a$DB=="5 mismatches" | a$DB=="6 mismatches",]
      e<-ggplot2::ggplot(data=b, aes(y=value, x=DB, group=variable)) + geom_point(aes(colour=variable)) +  geom_line(aes(colour=variable))+
        labs(y="Species amplified", x="") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  }
  e
}

#' Plot the distribution of GC content for degenerate primers
#' @title Plot the distribution of GC content for degenerate primers
#' @param primer1 Primer sequence.
#' @param primer2 Optional. Second primer sequence (i.e. reverse primer)
#' @note This function calculates GC content for all possible permutations of the degnerate primer provided.
#'     A couple fo internal functions are taken straight from \pkg{primerTree}
#' @return A plot of the distribution of GC content for degenerate primers.
#´
#' @examples
#' plot.GC.degenerates("TARTTAACAGCTADRHRY","TTCRAHTCCYYCYTYYYT")
#'
#' @export
basplot.GC.degenerates<-function(primer1,primer2=NULL){

  #define IUPAC, and two internal functions

  iupac = list( "M" = list("A", "C"),
                "R" = list("A", "G"),
                "W" = list("A", "T"),
                "S" = list("C", "G"),
                "Y" = list("C", "T"),
                "K" = list("G", "T"),
                "V" = list("A", "C", "G"),
                "H" = list("A", "C", "T"),
                "D" = list("A", "G", "T"),
                "B" = list("C", "G", "T"),
                "N" = list("A", "C", "G", "T"),
                "I" = list("A", "T", "C"))

  enumerate_primers = function(forward, reverse){
    forward_primers = enumerate_ambiguity(forward)
    data.frame(forward=forward_primers,
               reverse=rep(enumerate_ambiguity(reverse),
                           each=length(forward_primers)),
               stringsAsFactors = FALSE)
  }

  enumerate_ambiguity = function(sequence){
    search_regex = paste(names(iupac), collapse='|')
    locs = stringr::str_locate_all(sequence, search_regex)
    sequences = list()
    count = 1
    for (i in seq_len(nrow(locs[[1]]))){
      loc = locs[[1]][i,]
      ambiguity = stringr::str_sub(sequence, loc[1], loc[2])
      for(type in iupac[[ambiguity]]){
        new_seq = sequence
        stringr::str_sub(new_seq, loc[1], loc[2]) <- type
        sequences[[count]] = enumerate_ambiguity(new_seq)
        count = count + 1
      }
      return(unlist(sequences))
    }
    return(sequence)
  }


  #enumerate all primer possibilities
  perms_degen_F<-as.data.frame(enumerate_ambiguity(primer1))
  extra<-matrix(NA,length(perms_degen_F$`enumerate_ambiguity(primer1)`),5)
  perms_degen_F<-cbind(perms_degen_F,extra)
  colnames(perms_degen_F)<-c("primer_seq","A","C","G","T")

  #calculate ACGT content for each primer possibility
  for (i in 1:length(perms_degen_F$primer_seq)){
    perms_degen_F[i,"A"]<- stringr::str_count(perms_degen_F[i,"primer_seq"],"A")/nchar(as.character(perms_degen_F[i,"primer_seq"]))
    perms_degen_F[i,"C"]<- stringr::str_count(perms_degen_F[i,"primer_seq"],"C")/nchar(as.character(perms_degen_F[i,"primer_seq"]))
    perms_degen_F[i,"G"]<- stringr::str_count(perms_degen_F[i,"primer_seq"],"G")/nchar(as.character(perms_degen_F[i,"primer_seq"]))
    perms_degen_F[i,"T"]<- stringr::str_count(perms_degen_F[i,"primer_seq"],"T")/nchar(as.character(perms_degen_F[i,"primer_seq"]))
  }

  #calculate GC content for each primer possibility
  perms_degen_F$gc.content<-rowSums(perms_degen_F[,c("C","G")])

  #if primer2 is supplied repeat for this primer
  if (!is.null(primer2)){
    perms_degen_R<-as.data.frame(enumerate_ambiguity(primer2))
    extra<-matrix(NA,length(perms_degen_R$`enumerate_ambiguity(primer2)`),5)
    perms_degen_R<-cbind(perms_degen_R,extra)
    colnames(perms_degen_R)<-c("primer_seq","A","C","G","T","Tm")

    for (i in 1:length(perms_degen_R$primer_seq)){
      perms_degen_R[i,"A"]<- stringr::str_count(perms_degen_R[i,"primer_seq"],"A")/nchar(as.character(perms_degen_R[i,"primer_seq"]))
      perms_degen_R[i,"C"]<- stringr::str_count(perms_degen_R[i,"primer_seq"],"C")/nchar(as.character(perms_degen_R[i,"primer_seq"]))
      perms_degen_R[i,"G"]<- stringr::str_count(perms_degen_R[i,"primer_seq"],"G")/nchar(as.character(perms_degen_R[i,"primer_seq"]))
      perms_degen_R[i,"T"]<- stringr::str_count(perms_degen_R[i,"primer_seq"],"T")/nchar(as.character(perms_degen_R[i,"primer_seq"]))
    }

    #calculate GC content for each primer possibility
    perms_degen_R$gc.content<-rowSums(perms_degen_R[,c("C","G")])

  }

  #plot GC content
  h<-hist(perms_degen_F$gc.content, plot = F)
  h$density = h$counts/sum(h$counts)*100
  plot(h,freq=FALSE,main = "GC content distribution",xlab = "Proportion GC", ylab = "% primer permutations",col=rgb(0,0,1,1/4),ylim = c(0,100))
  legend("topright", c("primer1"), fill = rgb(0,0,1,1/4))
  if (!is.null(primer2)){
    hr<-hist(perms_degen_R$gc.content,plot=F)
    hr$density = hr$counts/sum(hr$counts)*100
    plot(hr,freq=FALSE,add=T,col=rgb(1,0,0,1/4),ylim = c(0,100))
    legend("topright", c("primer1","primer2"), fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))
  }
}

#' Plot the distribution of melting temperature for degenerate primers
#' @title Plot the distribution of melting temperature for degenerate primers
#' @param primer1 Primer sequence.
#' @param primer2 Optional. Second primer sequence (i.e. reverse primer)
#' @note This function calculates the melting temperature for all possible permutations of the degnerate primer provided.
#'     A couple fo internal functions are taken straight from \pkg{primerTree}
#' @return A plot of the distribution of melting temperatures for degenerate primers.
#´
#' @examples
#' plot.Tm.degenerates("TARTTAACAGCTADRHRY","TTCRAHTCCYYCYTYYYT")
#'
#' @export
basplot.Tm.degenerates<-function(primer1,primer2=NULL){

  #define IUPAC, and two internal functions

  iupac = list( "M" = list("A", "C"),
                "R" = list("A", "G"),
                "W" = list("A", "T"),
                "S" = list("C", "G"),
                "Y" = list("C", "T"),
                "K" = list("G", "T"),
                "V" = list("A", "C", "G"),
                "H" = list("A", "C", "T"),
                "D" = list("A", "G", "T"),
                "B" = list("C", "G", "T"),
                "N" = list("A", "C", "G", "T"),
                "I" = list("A", "T", "C"))

  enumerate_primers = function(forward, reverse){
    forward_primers = enumerate_ambiguity(forward)
    data.frame(forward=forward_primers,
               reverse=rep(enumerate_ambiguity(reverse),
                           each=length(forward_primers)),
               stringsAsFactors = FALSE)
  }

  enumerate_ambiguity = function(sequence){
    search_regex = paste(names(iupac), collapse='|')
    locs = stringr::str_locate_all(sequence, search_regex)
    sequences = list()
    count = 1
    for (i in seq_len(nrow(locs[[1]]))){
      loc = locs[[1]][i,]
      ambiguity = stringr::str_sub(sequence, loc[1], loc[2])
      for(type in iupac[[ambiguity]]){
        new_seq = sequence
        stringr::str_sub(new_seq, loc[1], loc[2]) <- type
        sequences[[count]] = enumerate_ambiguity(new_seq)
        count = count + 1
      }
      return(unlist(sequences))
    }
    return(sequence)
  }


  #enumerate all primer possibilities
  perms_degen_F<-as.data.frame(enumerate_ambiguity(primer1))
  extra<-matrix(NA,length(perms_degen_F$`enumerate_ambiguity(primer1)`),1)
  perms_degen_F<-cbind(perms_degen_F,extra)
  colnames(perms_degen_F)<-c("primer_seq","Tm")

  #calculate Tm for each primer possibility
  for (i in 1:length(perms_degen_F$primer_seq)){
    perms_degen_F[i,"Tm"]<-2*(stringr::str_count(perms_degen_F[i,"primer_seq"],"A")+stringr::str_count(perms_degen_F[i,"primer_seq"],"T")) +
      4*(stringr::str_count(perms_degen_F[i,"primer_seq"],"G")+stringr::str_count(perms_degen_F[i,"primer_seq"],"C"))
  }


  #if primer2 is supplied repeat for this primer
  if (!is.null(primer2)){
    perms_degen_R<-as.data.frame(enumerate_ambiguity(primer2))
    extra<-matrix(NA,length(perms_degen_R$`enumerate_ambiguity(primer2)`),1)
    perms_degen_R<-cbind(perms_degen_R,extra)
    colnames(perms_degen_R)<-c("primer_seq","Tm")

    for (i in 1:length(perms_degen_R$primer_seq)){
      perms_degen_R[i,"Tm"]<-2*(stringr::str_count(perms_degen_R[i,"primer_seq"],"A")+stringr::str_count(perms_degen_R[i,"primer_seq"],"T")) +
        4*(stringr::str_count(perms_degen_R[i,"primer_seq"],"G")+stringr::str_count(perms_degen_R[i,"primer_seq"],"C"))
    }
  }
  #plot Tm
  perms_degen_F$ForR<-"Primer1"
  if (is.null(primer2)){
  final.plot<-ggplot(perms_degen_F) +
    geom_histogram(aes(x=Tm,y=..density..,color=ForR),binwidth=0.1,alpha=0.2,position = "identity")+
    ylab("density")+xlab("Celsius")+ggtitle("Tm distribution")+
    scale_x_continuous(breaks = pretty(perms_degen_F$Tm, n = 15))
  }
  if (!is.null(primer2)){
  perms_degen_R$ForR<-"Primer2"
  perms_degen_combo<-rbind(perms_degen_F,perms_degen_R)
  final.plot<-ggplot(perms_degen_combo) +
    geom_histogram(aes(x=Tm,y=..density..,color=ForR),binwidth=0.1,alpha=0.2,position = "identity")+
    ylab("density")+xlab("Celsius")+ggtitle("Tm distribution")+
    scale_x_continuous(breaks = pretty(perms_degen_combo$Tm, n = 15))
  }
  final.plot
}

#' Plot the taxonomic resolution of a metabarcode
#' @title Plot the taxonomic resolution of a metabarcode
#' @param ecopcroutput.res A dataframe resulting from \code{ecopcroutput.res.Bas}.
#' @param forceLegend Forces a legend even when there are more than 30 groups. Can take a long time if set to TRUE.
#' @return A plot of the taxonomic resolution of a metabarcode.
#' @examples
#' plot.ecopcroutput.res.Bas(ecopcroutput.res)
#' @export
basplot.ecopcroutput.res<-function(ecopcroutput.res,forceLegend=F){

if(forceLegend==F){

if(length(unique(ecopcroutput.res$taxon))>30) message("Over 30 taxa, supressing legend")

if(length(unique(ecopcroutput.res$taxon))<30){a<-ggplot2::ggplot(ecopcroutput.res, ggplot2::aes(x="", y=Freq, fill=taxon))+
    ggplot2::facet_grid(. ~ resolution) +ggplot2::geom_bar(width = 1, stat = "identity")+ ggplot2::coord_polar("y", start=0)+
    ggplot2::scale_fill_manual(values = MyCols)}

if(length(unique(ecopcroutput.res$taxon))>30){a<-ggplot2::ggplot(ecopcroutput.res, ggplot2::aes(x="", y=Freq, fill=taxon))+
    ggplot2::facet_grid(. ~ resolution) +ggplot2::geom_bar(width = 1, stat = "identity")+ ggplot2::coord_polar("y", start=0)+
  ggplot2::theme(legend.position = "none")
}
}
  if(forceLegend==T){

      if(length(unique(ecopcroutput.res$taxon))<30){a<-ggplot2::ggplot(ecopcroutput.res, ggplot2::aes(x="", y=Freq, fill=taxon))+
      ggplot2::facet_grid(. ~ resolution) +ggplot2::geom_bar(width = 1, stat = "identity")+ ggplot2::coord_polar("y", start=0)+
      ggplot2::scale_fill_manual(values = MyCols)}

    if(length(unique(ecopcroutput.res$taxon))>30){a<-ggplot2::ggplot(ecopcroutput.res, ggplot2::aes(x="", y=Freq, fill=taxon))+
      ggplot2::facet_grid(. ~ resolution) +ggplot2::geom_bar(width = 1, stat = "identity")+ ggplot2::coord_polar("y", start=0)
    }
  }

  a
}

#' Make a primer mismatch logo plot
#' @title Make a primer mismatch logo plot
#' @param ecopcroutput A data frame such as that produced by \code{ecoPCR.Bas} or \code{ecopcr.hit.table}
#' @param primer1 Forward primer sequence.
#' @param primer2 Reverse primer sequence
#' @return DNA logo plot.
#´
#' @examples
#' a2<-system.file("extdata", "45F-63R_4Mis_ALLVERTS_REFSEQ.ecopcroutput", package = "bastools")
#' d2<-ecopcr.hit.table(a2,obitaxdb = b2,taxLevel = "class")
#' f2<-d2$ecopcroutput_clean
#' plot.ecopcr.primer.logo(f2,"TARTTAACAGCTANRHRY","TTCGATTCCTTCCTTTCT")
#'
#' @export
basplot.ecopcr.primer.logo<-function(ecopcroutput,forward=NULL,reverse=NULL){
  forward.shannon.scores=ROBIBarcodes::ecopcr.forward.shanon(ecopcroutput)
  reverse.shannon.scores=ROBIBarcodes::ecopcr.reverse.shanon(ecopcroutput)

  par(mfrow =c(2,1))
  ROBIBarcodes::dnalogoplot(forward.shannon.scores, primer = forward, main='Forward')
  ROBIBarcodes::dnalogoplot(reverse.shannon.scores, primer = reverse, main='Reverse')
  par(mfrow =c(1,1))
}




