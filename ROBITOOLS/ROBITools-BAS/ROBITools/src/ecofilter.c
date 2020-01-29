#include "ecoPCR.h"

int eco_is_taxid_included(	ecotaxonomy_t *taxonomy, 
							int32_t *restricted_taxid, 
							int32_t tab_len, 
							int32_t taxid)
{
	int i;
	ecotx_t *taxon;
	
	taxon = eco_findtaxonbytaxid(taxonomy, taxid);
	
	if (taxon)
		for (i=0; i < tab_len; i++)
			if ( (taxon->taxid == restricted_taxid[i]) ||
				 (eco_isundertaxon(taxon, restricted_taxid[i])) )
				return 1;
	
	return 0;
}
