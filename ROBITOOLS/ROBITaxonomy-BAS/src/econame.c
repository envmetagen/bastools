#include "ecoPCR.h"
#include <string.h>
#include <stdlib.h>

static econame_t *readnext_econame(FILE *f,econame_t *name,ecotaxonomy_t *taxonomy);

econameidx_t *read_nameidx(const char *filename,ecotaxonomy_t *taxonomy)
{

	int32_t      		count;
	FILE         		*f;
	econameidx_t		*indexname;
	int32_t      		i;
	
	f = open_ecorecorddb(filename,&count,0);
	
	if (f==NULL)
		return NULL;

	indexname = (econameidx_t*) ECOMALLOC(sizeof(econameidx_t) + sizeof(econame_t) * (count-1),"Allocate names");
	
	indexname->count=count;
	                                    
	for (i=0; i < count; i++){
		readnext_econame(f,(indexname->names)+i,taxonomy);
	}

	return indexname;
}

int32_t delete_nameidx(econameidx_t *nameidx)
{
	size_t i;

	if (nameidx) {
		for (i=0; i < nameidx->count; i++) {
			if (nameidx->names[i].name)
				ECOFREE(nameidx->names[i].name,
						"Desallocate name");
			if (nameidx->names[i].classname)
				ECOFREE(nameidx->names[i].classname,
						"Desallocate classname");
		}

		ECOFREE(nameidx,"Desallocate name index");

		return 0;
	}

	return 1;
}

econame_t *readnext_econame(FILE *f,econame_t *name,ecotaxonomy_t *taxonomy)
{
	
	econameformat_t *raw;
	int32_t  rs;
	
	raw = read_ecorecord(f,&rs);
	
	if (!raw)
		return NULL;

	if (is_big_endian())
	{
		raw->is_scientificname 	= swap_int32_t(raw->is_scientificname);
		raw->namelength 	    = swap_int32_t(raw->namelength);
		raw->classlength        = swap_int32_t(raw->classlength);
		raw->taxid  	        = swap_int32_t(raw->taxid); 
	}
	
	name->is_scientificname=raw->is_scientificname;
	
	name->name   	= ECOMALLOC((raw->namelength+1) * sizeof(char),"Allocate name");
	strncpy(name->name,raw->names,raw->namelength);
	name->name[raw->namelength]=0;
	
	name->classname = ECOMALLOC((raw->classlength+1) * sizeof(char),"Allocate classname");
	strncpy(name->classname,(raw->names+raw->namelength),raw->classlength);
	name->classname[raw->classlength]=0;
	
	name->taxon = taxonomy->taxons->taxon + raw->taxid;

	return name;
}

