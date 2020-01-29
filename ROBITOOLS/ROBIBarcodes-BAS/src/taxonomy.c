/*
 * taxonomy.c
 *
 *  Created on: 17 janv. 2013
 *      Author: coissac
 */

#include "robitax.h"


static ecorankidx_t *new_bdbecorankidx(void);
static ecotxidx_t *new_bdbecotxidx(int nbtaxa, int rootrank);

static ecorankidx_t *new_bdbecorankidx(void)
{
  ecorankidx_t *index;
  int i;
  size_t size;
  char* rank[]={"class",      "family",         "forma",
                "genus",      "infraclass",     "infraorder",
                "kingdom",    "motu",           "no rank",
                "order",      "parvorder",      "phylum",
                "species",    "species group",  "species subgroup",
                "subclass",   "subfamily",      "subgenus",
                "subkingdom", "suborder",       "subphylum",
                "subspecies", "subtribe",       "superclass",
                "superfamily","superkingdom",   "superorder",
                "superphylum","tribe",          "varietas"};



  index = (ecorankidx_t*) ECOMALLOC(sizeof(ecorankidx_t) + sizeof(char*) * (29),
	                                  "Allocate rank index");

  index->count=30;
  for (i=0; i < index->count; i++)
		{
			size = strlen(rank[i]);
			index->label[i]=(char*) ECOMALLOC(size+1,
			                                  "Allocate rank label");
			strcpy(index->label[i],rank[i]);
		}

	return index;
}

static ecotxidx_t *new_bdbecotxidx(int nbtaxa, int rootrank) {

	  ecotxidx_t *index;
	  int rootnamelen = 4; // "rootname="root"

	  index = (ecotxidx_t*) ECOMALLOC(sizeof(ecotxidx_t) + sizeof(ecotx_t) * (nbtaxa),
	                                    "Allocate taxonomy");

	  //initialize the taxonomy with the root taxon
	  index->count=nbtaxa+1;
	  index->maxtaxid=1;
	  index->buffersize=nbtaxa+1;

	  // Create the root taxon
	  index->taxon[0].taxid=1;
	  index->taxon[0].rank=rootrank;
	  index->taxon[0].farest=0;
	  index->taxon[0].parent=(ecotx_t*)1;

	  index->taxon[0].name= (char*) ECOMALLOC(sizeof(char) * rootnamelen +1,
	                                          "Allocate taxonomy root name");

	  strcpy(index->taxon[0].name, "root");

	  return index;

}

static int cmptaxid(const void* t1, const void* t2) {
	int et1 = ((ecotx_t*)t1)->taxid;
	int et2 = ((ecotx_t*)t2)->taxid;

	if (et1 < et2)
		return -1;

	if (et2 < et1)
		return 1;

	return 0;
}

static int findparentcmp(const void* p, const void* t) {
	int parent = (int)p;
	int et = ((ecotx_t*)t)->taxid;

	if (parent < et)
		return -1;

	if (et < parent)
		return 1;

	return 0;
}

SEXP R_buildbarcodetaxo(SEXP rdata) {
	SEXP names;
	SEXP taxids;
	SEXP ranks;
	SEXP partofs;
	SEXP Rtax;

	int nbtaxa;
	int rootrank;
	int i;

	int taxid;
	char *name;
	int rank;
	int partof;

	ecotxidx_t* parent;

	taxids = getAttrib(rdata, R_RowNamesSymbol);
	names  = VECTOR_ELT(rdata, 0);
	ranks  = VECTOR_ELT(rdata, 1);
	partofs= VECTOR_ELT(rdata, 2);

	nbtaxa = GET_LENGTH(taxids);

	ecotaxonomy_t *tax;

	tax = ECOMALLOC(sizeof(ecotaxonomy_t),
					"Allocate taxonomy structure");

	tax->ranks =new_bdbecorankidx();

	rootrank = rank_index("no rank",tax->ranks);

	tax->taxons=new_bdbecotxidx(nbtaxa,rootrank);
	tax->names =NULL;


	for (i=0; i< nbtaxa; i++) {
		taxid = atol(CHAR(STRING_ELT(taxids, i)) + 3);
		name  = (char*) CHAR(STRING_ELT(names, i));
		rank  = rank_index(CHAR(STRING_ELT(ranks, i)),tax->ranks);
		partof = atol(CHAR(STRING_ELT(partofs, i)) + 3);

    if (taxid > tax->taxons->maxtaxid)
      tax->taxons->maxtaxid=taxid;

		tax->taxons->taxon[i+1].taxid=taxid;
		tax->taxons->taxon[i+1].rank=rootrank;
		tax->taxons->taxon[i+1].farest=0;
		tax->taxons->taxon[i+1].parent=(ecotx_t*)((size_t)partof);

		tax->taxons->taxon[i+1].name= (char*) ECOMALLOC(sizeof(char) * strlen(name) +1,
		                                          "Allocate taxonomy root name");

		strcpy(tax->taxons->taxon[i+1].name, name);
	}

	qsort((void*)tax->taxons->taxon,nbtaxa+1,sizeof(ecotx_t),cmptaxid);

	for (i=0; i< nbtaxa+1; i++) {
		parent = (ecotxidx_t*) bsearch((void*)(tax->taxons->taxon[i].parent),
				                       (void*)(tax->taxons->taxon),
				                       nbtaxa+1,
                               sizeof(ecotx_t),
				                       findparentcmp);
                               
           
		if (parent==NULL)
			error("Error during taxonomy indexing");

		tax->taxons->taxon[i].parent=(struct ecotxnode *)parent;
	} 

	Rtax = PROTECT(R_MakeExternalPtr(tax, mkString("ROBITools NCBI Taxonomy pointer"), R_NilValue));
	R_RegisterCFinalizerEx(Rtax, (R_CFinalizer_t)R_delete_taxonomy,TRUE);

	UNPROTECT(1);

	return Rtax;

}


SEXP R_delete_taxonomy(SEXP Rtaxonomy)
{
	ecotaxonomy_t *ptax;
	SEXP pointer;

    ptax = (ecotaxonomy_t *) R_ExternalPtrAddr(Rtaxonomy);

    (void) delete_ecotaxonomy(ptax);

    // Clear the external pointer
    R_ClearExternalPtr(Rtaxonomy);

    return R_NilValue;

}
