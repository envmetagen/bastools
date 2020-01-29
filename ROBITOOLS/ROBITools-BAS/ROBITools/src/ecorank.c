#include "ecoPCR.h"
#include <string.h>
#include <stdlib.h>

static int compareRankLabel(const void *label1, const void *label2);

ecorankidx_t     *read_rankidx(const char *filename)
{
	int32_t      count;
	FILE         *f;
	ecorankidx_t *index;
	int32_t      i;
	int32_t      rs;
	char         *buffer;
	
	f = open_ecorecorddb(filename,&count,0);
	
	if (f==NULL)
		return NULL;

	index = (ecorankidx_t*) ECOMALLOC(sizeof(ecorankidx_t) + sizeof(char*) * (count-1),
	                                  "Allocate rank index");
	 
	index->count=count;                                 
	          
	for (i=0; i < count; i++)
		{
			buffer = read_ecorecord(f,&rs);
			index->label[i]=(char*) ECOMALLOC(rs+1,
			                                  "Allocate rank label");
			strncpy(index->label[i],buffer,rs);
		}
		
	return index;
}

int32_t rank_index(const char* label,ecorankidx_t* ranks)
{
	char **rep;
	
	rep = bsearch(label,ranks->label,ranks->count,sizeof(char*),compareRankLabel);
	
	if (rep)
		return rep-ranks->label;
//	else
//		ECOERROR(ECO_NOTFOUND_ERROR,"Rank label not found");
		
	return -1;
}


int compareRankLabel(const void *label1, const void *label2)
{
	return strcmp((const char*)label1,*(const char**)label2);
}
