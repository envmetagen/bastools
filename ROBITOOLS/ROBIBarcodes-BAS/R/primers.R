#'@include xmlMods.R
#'@import XML
#'
NULL


add.primer.barcodedb = function(barcodedb,
                               name,
                               sequence,
                               coding,
                               documentation) {
'  <primer ID="PR.XXX">
    <name>Xxx</name>
    <sequence>CGATCGATGCTAGCTAGCTGAT</sequence>
    <coding>false</coding>
    </primer>'
  
}