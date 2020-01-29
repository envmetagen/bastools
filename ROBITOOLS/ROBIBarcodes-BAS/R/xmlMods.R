#'@include ROBIBarcodes.R
#'@import XML
#'
NULL


#
# Checks that bibutils is installed on the system
#
hasBibUtils = !system("bib2xml -h",
                       intern=FALSE,
                       ignore.stdout = TRUE, 
                       ignore.stderr = TRUE)

if (!hasBibUtils) {
  message("\n============================================================\n",
          "Bibutils are not installed on your system\n",
          "or is not correctly setup\n",
          "Consider to visit: http://sourceforge.net/projects/bibutils/\n",
          "============================================================\n")
}


#'Qualify the elements of the `mods` elements with the \code{mods:} namespace.
#'
#' The \code{bibutils} programs generate XML file not qualified by a schema.
#' To respect the OBIBarcodes schema the mods elements must be qualified.
#' This function take a XML document produce by a \code{bibutils} program and
#' add the \code{mods:} namespace.
#' 
#' @note This function modifies the document past in argument and returns 
#'       nothing
#'
#' @note This is an internal function and consequently has not to be called by
#'       end users. 
#'       
#' @param   modsdoc a XMLInternalDocument instance corresponding to a 
#'                  modsCollectionDefinition element.
#'                  
#' @export
addmodsnamespace = function(modsdoc) {
  
  root = xmlRoot(modsdoc)
  xmlNamespaces(root,set=TRUE)=c(mods="http://www.loc.gov/mods/v3")
  hiden=xpathApply(modsdoc,
                   path='/.//*',
                   fun= function(n) xmlNamespace(n,set=TRUE)="mods")

}

patchmodsID = function(modsdoc) {
  hiden= xpathApply(modsdoc,'/.//*[attribute::ID]',
                    function(n) xmlAttrs(n)=list(ID=paste('BI.',
                                                  toupper(gsub('[^A-Za-z0-9_]',
                                                               '_',
                                                               xmlAttrs(n)['ID']
                                                               )
                                                          ),
                                                          sep=""
                                                  ))
                   )
}

#'@export
bib2mods = function(bibfile,bibutils='bib2xml') {
  
  tmp=tempfile()
  xmlerr=system(paste(bibutils,bibfile,'>',tmp,sep=' '),
                intern=FALSE,
                ignore.stderr=TRUE)
  
  if (xmlerr!=0)
    stop(paste("Cannot run ",bibutils))
  xml = paste(tmp,collapse='\n')
  xml = xmlParseDoc(tmp,asText=FALSE)
  file.remove(tmp)
  addmodsnamespace(xml)
  patchmodsID(xml)
  mods=getNodeSet(xml,'/.//mods:mods[attribute::ID]')
  
  return(mods)
}

if (!hasBibUtils) {
  bib2mods = function(bibfile,bibutils) {  
    stop("Bibutils not install visit: http://sourceforge.net/projects/bibutils/")
  }
}