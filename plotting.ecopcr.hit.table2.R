#' Some colours I like to use for plots (n=29)
#'@export
MyCols <- c("dodgerblue2","#E31A1C", # red
            "green4",
            "#6A3D9A", # purple
            "#FF7F00", # orange
            "black","gold1",
            "skyblue2","#FB9A99", # lt pink
            "palegreen2",
            "#CAB2D6", # lt purple
            "#FDBF6F", # lt orange
            "gray70", "khaki2",
            "maroon","orchid1","deeppink1","blue1","steelblue4",
            "darkturquoise","green1","yellow4","yellow3",
            "darkorange4","brown","chartreuse2","coral4","darkslategray1","mediumseagreen")

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
basplot.ecopcr.hit.table2 <- function(ecopcr.hit.table.output,specify_taxon=NULL){

  #ideas
  #if nothing specified, then just plot totals - a line plot for each taxlevel and ggarranged
  #if a certain taxon is specified e.g. Anura - a line plot for each taxlevel (wihtin anura) and ggarranged

  if(is.null(specify_taxon)){

    plot.melted<-lapply(X = melted,FUN=function(x)ggplot2::ggplot(data=x,ggplot2::aes(y=value, x=Mismatches, group=variable) +
                            ggplot2::geom_point(ggplot2::aes(colour=variable)) +
                            ggplot2::geom_line(ggplot2::aes(colour=variable))+
                            ggplot2::labs(y="Species amplified", x="Number of primer mismatches") +
                            ggplot2::theme(axis.text.x = ggplot2::element_text(size=5,angle = 45, hjust = 1),legend.position = "none") +
                            #scale_colour_manual(values = MyCols)+
                            ggplot2::scale_x_discrete(limits=c("none","one","two","three","four","five","six"))))

  melted<-list()
  melted<-lapply(X = ecopcr.hit.table.output[1:5],FUN = function(x)reshape2::melt(x))
  plot.melted<-list()
  plot.melted<-lapply(X = melted,FUN=function(x)ggplot2::ggplot(data=x,
                                    ggplot2::aes(y=value, x=Mismatches, group=variable)) +
                        ggplot2::geom_point(ggplot2::aes(colour=variable)) +
                        ggplot2::geom_line(ggplot2::aes(colour=variable))+
                        ggplot2::labs(y="Species amplified", x="Number of primer mismatches") +
                        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),legend.position = "none") +
                        #scale_colour_manual(values = MyCols)+
                        ggplot2::scale_x_discrete(limits=c("none","one","two","three","four","five","six")))

}

a<-melted[[1]]
b<-a[a$variable!="Total",]
Number_of_primer_mismatches <- factor(b$Mismatches, level = c("none","one","two","three","four","five","six"))
plot1<-ggplot(data=b,ggplot2::aes(y=value, x=Number_of_primer_mismatches)) +
  geom_col(aes(fill = variable))+
  ggplot2::labs(y="Species amplified", x="Number of primer mismatches")+
  ggplot2::theme(legend.text=element_text(size=8),legend.key = element_rect(size = 8),
                 axis.text.x = ggplot2::element_text(size=8,angle = 45, hjust = 1),legend.position = "right")

if(sum(unique(b$variable))>20){
  message(c("Over 20 taxa in plot. To avoid a difficult legend, it has been excluded. The taxa are ", as.character(unique(b$variable))))

  plot2<-ggplot(data=b,ggplot2::aes(y=value, x=Number_of_primer_mismatches)) +
    geom_col(aes(fill = variable))+
    ggplot2::labs(y="Species amplified", x="Number of primer mismatches")+
    ggplot2::theme(legend.text=element_text(size=5),legend.key = element_rect(size = 3),
                   axis.text.x = ggplot2::element_text(size=8,angle = 45, hjust = 1),legend.position = "none")
}

}
