/*
 * robitax.c
 *
 *  Created on: 17 janv. 2013
 *      Author: coissac
 */

#include "robitax.h"
#include <unistd.h>
//#include <regex.h>
#include "slre.h"

/**
 * Return a pointeur to an obitools taxonomy C structure
 * from an R instance of taxonomy.obitools
 *
 * The function checks if the pointer stored in the R object is set
 * to NULL. In this case this means that we have to load the taxonomy
 * from the disk.
 *
 * @param taxonomy an R object
 * @type  taxonomy SEXP
 *
 * @return a pointer to the C structure
 * @rtype  ecotaxonomy_t *
 */

ecotaxonomy_t *getTaxPointer(SEXP Rtaxonomy)
{


	char *pwd;
	SEXP pointer;
	SEXP rclass;
	SEXP rdir;
	SEXP rfile;
	ecotaxonomy_t *ptax;
	const char *class;
	const char *file;
	const char *dir;

  int saved;

    if (!IS_S4_OBJECT(Rtaxonomy) )
        error("argument not obitools.taxonomy instance");

	// We get the class name and compare it to "taxonomy.obitools"
    rclass = getAttrib(Rtaxonomy, R_ClassSymbol);
    class = CHAR(asChar(rclass));

    if (strcmp(class,"obitools.taxonomy"))
        error("argument not obitools.taxonomy instance");

    pointer = R_do_slot(Rtaxonomy,mkString("pointer"));
    saved = LOGICAL(R_do_slot(Rtaxonomy,mkString("saved")))[0];
    ptax = (ecotaxonomy_t *) R_ExternalPtrAddr(pointer);

    // If the external pointer is set to NULL we have to load
    // the taxonomy from file
    if (ptax==NULL && saved)
    {
    	pwd = getcwd(NULL,0);

    	rfile = R_do_slot(Rtaxonomy,mkString("dbname"));
    	file  = CHAR(asChar(rfile));

    	rdir  = R_do_slot(Rtaxonomy,mkString("workingdir"));
    	dir   = CHAR(asChar(rdir));

    	chdir(dir);

    	ptax = read_taxonomy(file,1);

    	R_SetExternalPtrAddr(pointer,(void*)ptax);

    	chdir(pwd);
    	free(pwd);
    }
    
    if (ptax==NULL && ! saved)
      error("The taxonomy instance is no more valid and must be rebuilt");

    return ptax;
}

SEXP R_delete_taxonomy(SEXP Rtaxonomy)
{
	ecotaxonomy_t *ptax;
//	SEXP pointer;

    ptax = (ecotaxonomy_t *) R_ExternalPtrAddr(Rtaxonomy);

    (void) delete_ecotaxonomy(ptax);

    // Clear the external pointer
    R_ClearExternalPtr(Rtaxonomy);

    return R_NilValue;

}



SEXP R_read_taxonomy(SEXP filename, SEXP altenative)
{
	int   alt;
	const char* file;
	SEXP  Rtax;

    if (! isString(filename))
        error("filename not character");
    file = CHAR(STRING_ELT(filename, 0));

    if (! isLogical(altenative))
        error("altenative not logical");
    alt = LOGICAL(altenative)[0];

    ecotaxonomy_t *taxonomy = read_taxonomy(file,alt);

    if (! taxonomy)
    	error("Cannot open the taxonomy database");

	Rtax = PROTECT(R_MakeExternalPtr(taxonomy, mkString("ROBITools NCBI Taxonomy pointer"), R_NilValue));
	R_RegisterCFinalizerEx(Rtax, (R_CFinalizer_t)R_delete_taxonomy,TRUE);

    UNPROTECT(1);


	return Rtax;
}


SEXP R_get_scientific_name(SEXP Rtaxonomy,SEXP Rtaxid)
{
	ecotx_t *taxon;
	ecotaxonomy_t *ptax;
	int taxid;
//	SEXP scname;

    ptax = getTaxPointer(Rtaxonomy);

    if (! isInteger(Rtaxid))
        error("taxid not integer");

    taxid = *INTEGER(Rtaxid);

    if (! (taxid > 0))
        error("taxid not positive");


	taxon = eco_findtaxonbytaxid(ptax, taxid);

	if (!taxon)
		return ScalarString(R_NaString);
		// error("unkown taxid");

	return mkString(taxon->name);

}

SEXP R_get_rank(SEXP Rtaxonomy,SEXP Rtaxid)
{
	ecotx_t *taxon;
	ecotaxonomy_t *ptax;
	int *taxid;
  int ntaxid;
  int i;
  SEXP results;
//	SEXP scname;

    ptax = getTaxPointer(Rtaxonomy);

    if (! isInteger(Rtaxid))
        error("taxid not integer");
        
    ntaxid = length(Rtaxid);
    
    results = PROTECT(allocVector(STRSXP, ntaxid));

    taxid = INTEGER(Rtaxid);

    for (i=0; i < ntaxid; i++)
    {
      if (taxid[i]== NA_INTEGER || taxid[i] <= 0)
        SET_STRING_ELT(results, i, R_NaString);
      else {
        taxon = eco_findtaxonbytaxid(ptax, taxid[i]);
        if (!taxon)
          SET_STRING_ELT(results, i, R_NaString);
        else
          SET_STRING_ELT(results, i, mkChar(ptax->ranks->label[taxon->rank]));
      }
    }

  UNPROTECT(1);
  
	return results;

}

SEXP R_findtaxonatrank(SEXP Rtaxonomy,SEXP Rtaxid,SEXP Rrank, SEXP Rname)
{
	ecotx_t *taxon;
	ecotx_t *rep;
	ecotaxonomy_t *ptax;
	int taxid;
	int name;
	const char *rank;
	int   rankidx;
//	SEXP scname;

    ptax = getTaxPointer(Rtaxonomy);

    if (! isInteger(Rtaxid))
        error("taxid not integer");

    taxid = *INTEGER(Rtaxid);

    if (! (taxid > 0))
        error("taxid not positive");

    if (! isString(Rrank))
        error("rank not a string");

    rank=CHAR(STRING_ELT(Rrank,0));

    rankidx=rank_index(rank,ptax->ranks);

    if (rankidx < 0)
        error("unkown rank name");

    if (! isLogical(Rname))
        error("name not logical");
    name = LOGICAL(Rname)[0];

	taxon = eco_findtaxonbytaxid(ptax, taxid);

	if (!taxon)
  {
		if (name)
			return ScalarString(R_NaString);
		else
			return ScalarInteger(R_NaInt);
  }
	rep = eco_findtaxonatrank(taxon,rankidx);

	if (!rep)
  {
		if (name)
			return ScalarString(R_NaString);
		else
			return ScalarInteger(R_NaInt);
  }
  
	if (name)
		return mkString(rep->name);
	else
		return ScalarInteger(rep->taxid);
}


SEXP R_get_species(SEXP Rtaxonomy,SEXP Rtaxid,SEXP Rname)
{
	ecotx_t *taxon;
	ecotx_t *rep;
	ecotaxonomy_t *ptax;
	int taxid;
	int name;
//	SEXP scname;

    ptax = getTaxPointer(Rtaxonomy);

    if (! isInteger(Rtaxid))
        error("taxid not integer");

    taxid = *INTEGER(Rtaxid);

    if (! (taxid > 0))
        error("taxid not positive");

    if (! isLogical(Rname))
        error("name not logical");
    name = LOGICAL(Rname)[0];

	taxon = eco_findtaxonbytaxid(ptax, taxid);

	if (!taxon)
  {
		if (name)
			return ScalarString(R_NaString);
		else
			return ScalarInteger(R_NaInt);
  }
  
	rep = eco_getspecies(taxon,ptax);

	if (!rep)
  {
		if (name)
			return ScalarString(R_NaString);
		else
			return ScalarInteger(R_NaInt);
  }
  
	if (name)
		return mkString(rep->name);
	else
		return ScalarInteger(rep->taxid);
}

SEXP R_get_genus(SEXP Rtaxonomy,SEXP Rtaxid,SEXP Rname)
{
	ecotx_t *taxon;
	ecotx_t *rep;
	ecotaxonomy_t *ptax;
	int taxid;
	int name;
//	SEXP scname;

    ptax = getTaxPointer(Rtaxonomy);

    if (! isInteger(Rtaxid))
        error("taxid not integer");

    taxid = *INTEGER(Rtaxid);

    if (! (taxid > 0))
        error("taxid not positive");

    if (! isLogical(Rname))
        error("name not logical");
    name = LOGICAL(Rname)[0];

	taxon = eco_findtaxonbytaxid(ptax, taxid);

	if (!taxon)
  {
		if (name)
			return ScalarString(R_NaString);
		else
			return ScalarInteger(R_NaInt);
  }
  
	rep = eco_getgenus(taxon,ptax);

	if (!rep)
  {
		if (name)
			return ScalarString(R_NaString);
		else
			return ScalarInteger(R_NaInt);
  }
  
	if (name)
		return mkString(rep->name);
	else
		return ScalarInteger(rep->taxid);
}

SEXP R_get_family(SEXP Rtaxonomy,SEXP Rtaxid,SEXP Rname)
{
	ecotx_t *taxon;
	ecotx_t *rep;
	ecotaxonomy_t *ptax;
	int taxid;
	int name;
//	SEXP scname;

    ptax = getTaxPointer(Rtaxonomy);

    if (! isInteger(Rtaxid))
        error("taxid not integer");

    taxid = *INTEGER(Rtaxid);

    if (! (taxid > 0))
        error("taxid not positive");

    if (! isLogical(Rname))
        error("name not logical");
    name = LOGICAL(Rname)[0];

	taxon = eco_findtaxonbytaxid(ptax, taxid);

	if (!taxon)
  {
		if (name)
			return ScalarString(R_NaString);
		else
			return ScalarInteger(R_NaInt);
  }
  
	rep = eco_getfamily(taxon,ptax);

	if (!rep)
  {
		if (name)
			return ScalarString(R_NaString);
		else
			return ScalarInteger(R_NaInt);
  }
  
	if (name)
		return mkString(rep->name);
	else
		return ScalarInteger(rep->taxid);
}

SEXP R_get_kingdom(SEXP Rtaxonomy,SEXP Rtaxid,SEXP Rname)
{
	ecotx_t *taxon;
	ecotx_t *rep;
	ecotaxonomy_t *ptax;
	int taxid;
	int name;
//	SEXP scname;

    ptax = getTaxPointer(Rtaxonomy);

    if (! isInteger(Rtaxid))
        error("taxid not integer");

    taxid = *INTEGER(Rtaxid);

    if (! (taxid > 0))
        error("taxid not positive");

    if (! isLogical(Rname))
        error("name not logical");
    name = LOGICAL(Rname)[0];

	taxon = eco_findtaxonbytaxid(ptax, taxid);

	if (!taxon)
  {
		if (name)
			return ScalarString(R_NaString);
		else
			return ScalarInteger(R_NaInt);
  }
  
	rep = eco_getkingdom(taxon,ptax);

	if (!rep)
  {
		if (name)
			return ScalarString(R_NaString);
		else
			return ScalarInteger(R_NaInt);
  }
  
	if (name)
		return mkString(rep->name);
	else
		return ScalarInteger(rep->taxid);
}

SEXP R_get_superkingdom(SEXP Rtaxonomy,SEXP Rtaxid,SEXP Rname)
{
	ecotx_t *taxon;
	ecotx_t *rep;
	ecotaxonomy_t *ptax;
	int taxid;
	int name;
//	SEXP scname;

    ptax = getTaxPointer(Rtaxonomy);

    if (! isInteger(Rtaxid))
        error("taxid not integer");

    taxid = *INTEGER(Rtaxid);

    if (! (taxid > 0))
        error("taxid not positive");

    if (! isLogical(Rname))
        error("name not logical");
    name = LOGICAL(Rname)[0];

	taxon = eco_findtaxonbytaxid(ptax, taxid);

	if (!taxon)
  {
		if (name)
			return ScalarString(R_NaString);
		else
			return ScalarInteger(R_NaInt);
  }
  
	rep = eco_getsuperkingdom(taxon,ptax);

	if (!rep)
  {
		if (name)
			return ScalarString(R_NaString);
		else
			return ScalarInteger(R_NaInt);
  }
  
	if (name)
		return mkString(rep->name);
	else
		return ScalarInteger(rep->taxid);
}

SEXP R_get_parent(SEXP Rtaxonomy,SEXP Rtaxid,SEXP Rname)
{
	ecotx_t *taxon;
	ecotx_t *rep;
	ecotaxonomy_t *ptax;
	int taxid;
	int name;
//	SEXP scname;

    ptax = getTaxPointer(Rtaxonomy);

    if (! isInteger(Rtaxid))
        error("taxid not integer");

    taxid = *INTEGER(Rtaxid);

    if (! (taxid > 0))
        error("taxid not positive");

    if (! isLogical(Rname))
        error("name not logical");
    name = LOGICAL(Rname)[0];

	taxon = eco_findtaxonbytaxid(ptax, taxid);

	if (!taxon)
  {
		if (name)
			return ScalarString(R_NaString);
		else
			return ScalarInteger(R_NaInt);
  }
  
	rep = taxon->parent;

	if (rep->taxid==taxid)
  {
		if (name)
			return ScalarString(R_NaString);
		else
			return ScalarInteger(R_NaInt);
  }
  
	if (name)
		return mkString(rep->name);
	else
		return ScalarInteger(rep->taxid);
}


SEXP R_validate_taxid(SEXP Rtaxonomy,SEXP Rtaxid)
{
	ecotx_t *taxon;
	ecotaxonomy_t *ptax;
	int taxid;
//	int name;
//	SEXP scname;

    ptax = getTaxPointer(Rtaxonomy);

    if (! isInteger(Rtaxid))
        error("taxid not integer");

    taxid = *INTEGER(Rtaxid);

    if (! (taxid > 0))
    	return ScalarInteger(R_NaInt);

	taxon = eco_findtaxonbytaxid(ptax, taxid);

	if (!taxon)
		return ScalarInteger(R_NaInt);
	else
		return ScalarInteger(taxon->taxid);
}


SEXP R_is_under_taxon(SEXP Rtaxonomy, SEXP Rtaxid, SEXP Rparent)
{
	ecotx_t *taxon;
	ecotaxonomy_t *ptax;
	int taxid;
  int *ptaxid;
  int ntaxid;
	int parent;
	SEXP rep;
  int* prep;
  size_t i;
//	SEXP isunder;

    ptax = getTaxPointer(Rtaxonomy);

    if (! isInteger(Rparent))
        error("parent not integer");
        
    if (GET_LENGTH(Rparent) > 1)
        error("parent must have a length equal to one.");

    parent = *INTEGER(Rparent);
    
    if (! isInteger(Rtaxid))
        error("taxid not integer");

    ntaxid = GET_LENGTH(Rtaxid);
    rep = PROTECT(NEW_LOGICAL(ntaxid));
    prep=LOGICAL(rep);

    ptaxid = INTEGER_POINTER(Rtaxid); 


    if (parent > 0) {
        taxon = eco_findtaxonbytaxid(ptax, parent);
        if (taxon) {
          for (i=0; i < ntaxid; i++) {
            if (ptaxid[i] <= 0)
              prep[i]=NA_LOGICAL;
            else {
              taxon = eco_findtaxonbytaxid(ptax, ptaxid[i]);
              prep[i] = eco_isundertaxon(taxon, parent);            
            }  
          }          
        }
    }    
    else
    {
      for (i=0; i < ntaxid; i++)
        prep[i]=NA_LOGICAL;
    }

  unprotect(1);
	return rep;
}


SEXP R_longest_path(SEXP Rtaxonomy,SEXP Rtaxid)
{
	ecotx_t *taxon;
	ecotaxonomy_t *ptax;
	int taxid;
//	int name;
//	SEXP scname;

    ptax = getTaxPointer(Rtaxonomy);

    if (! isInteger(Rtaxid))
        error("taxid not integer");

    taxid = *INTEGER(Rtaxid);

    if (taxid <= 0)
    	return ScalarInteger(R_NaInt);

	taxon = eco_findtaxonbytaxid(ptax, taxid);

	if (!taxon)
		return ScalarInteger(R_NaInt);
	else
		return ScalarInteger(taxon->farest);
}

SEXP R_rank_list(SEXP Rtaxonomy)
{
	int nrank;
	int i;
	ecotaxonomy_t *ptax;
	SEXP rNames;

    ptax = getTaxPointer(Rtaxonomy);

    nrank = ptax->ranks->count;

    rNames = PROTECT(allocVector(STRSXP, nrank));

	for (i=0; i < nrank;i++)
		SET_STRING_ELT(rNames, i, mkChar(ptax->ranks->label[i]));

	UNPROTECT(1);

	return rNames;
}

SEXP R_taxid_list(SEXP Rtaxonomy)
{
	int ntaxid;
	int i;
	ecotaxonomy_t *ptax;
	SEXP rTaxids;

    ptax = getTaxPointer(Rtaxonomy);
    ntaxid  = ptax->taxons->count;
    rTaxids = PROTECT(allocVector(INTSXP, ntaxid));

	for (i=0; i < ntaxid;i++)
		INTEGER(rTaxids)[i]=ptax->taxons->taxon[i].taxid;

	UNPROTECT(1);

	return rTaxids;

}

SEXP R_max_taxid(SEXP Rtaxonomy)
{
//	int nrank;
//	int i;
	ecotaxonomy_t *ptax;
//	SEXP rNames;

    ptax = getTaxPointer(Rtaxonomy);

    return ScalarInteger(ptax->taxons->maxtaxid);
}

SEXP R_length_taxonomy(SEXP Rtaxonomy)
{
	ecotaxonomy_t *ptax;

    ptax = getTaxPointer(Rtaxonomy);

    return ScalarInteger(ptax->taxons->count);
}

SEXP R_ecofind(SEXP Rtaxonomy, SEXP Rpattern, SEXP Rrank, SEXP Ralternative)
{
  ecotaxonomy_t *ptax;
  econame_t		*name;
  char*     pattern=NULL;
	int				re_match;
  SEXP      taxids;
  int32_t*  buffer;
  int32_t		tax_count	= 0;
  size_t 		j 		= 0;
  int32_t		rankfilter 	= 1;
  int*      ptaxid;
  char			*rankname=NULL;
  int32_t			nummatch 	= 0;
  int32_t         alternative = 0;

  size_t    bsize;
	
  ptax = getTaxPointer(Rtaxonomy);
  tax_count = ptax->taxons->count;
	
  if (! isString(Rpattern))
      error("pattern not a string");

  pattern= (char*) CHAR(STRING_ELT(Rpattern,0));

  if (! isNull(Rrank))
  {
    if (! isString(Rrank))
      error("rank not a string");

    rankname= (char*) CHAR(STRING_ELT(Rrank,0));
  }
  
  if (! isLogical(Ralternative))
      error("rank not a logical");
      
  alternative = LOGICAL(Ralternative)[0];

		
	nummatch=0;
  buffer = (int32_t*) malloc(100 * sizeof(int32_t));
  bsize=100;

  if (alternative && ptax->names!=NULL)
	  for (j=0,name=ptax->names->names;
			  j < ptax->names->count;
			  name++,j++)
	  {
		  if(rankname)
			  rankfilter = !(strcmp(rankname,ptax->ranks->label[name->taxon->rank]));

  	  re_match = slre_match(pattern, name->name, 
                            strlen(name->name), 
                            NULL, 0, 
                            SLRE_IGNORE_CASE);

  	  if (re_match > 0 && rankfilter)
		  {
			  buffer[nummatch]=name->taxon->taxid;
			  nummatch++;
			  if (nummatch==bsize) {
				  bsize*=2;
				  buffer = (int32_t*) realloc(buffer, bsize * sizeof(int32_t));
				  if (buffer==0)
				  {
					  // regfree(&re_preg);
					  error("Cannot allocate memory for the taxid list");
				  }
			  }
		  }

	  }
  else
	  for (j=0; j < ptax->taxons->count;j++)
	  {
		  if(rankname)
			  rankfilter = !(strcmp(rankname,ptax->ranks->label[ptax->taxons->taxon[j].rank]));

//		  re_match = regexec (&re_preg, ptax->taxons->taxon[j].name, 0, NULL, 0);
      re_match = slre_match(pattern, ptax->taxons->taxon[j].name, 
                            strlen(ptax->taxons->taxon[j].name), 
                            NULL, 0, 
                            SLRE_IGNORE_CASE);


//  	  if (!re_match && rankfilter)
  	  if (re_match > 0 && rankfilter)
		  {
			  buffer[nummatch]=ptax->taxons->taxon[j].taxid;
			  nummatch++;
			  if (nummatch==bsize) {
				  bsize*=2;
				  buffer = (int32_t*) realloc(buffer, bsize * sizeof(int32_t));
				  if (buffer==0)
				  {
					  // regfree(&re_preg);
					  error("Cannot allocate memory for the taxid list");
				  }
			  }
		  }

	  }

   	//regfree(&re_preg);

    taxids = PROTECT(NEW_INTEGER(nummatch));
    ptaxid = INTEGER(taxids);
    
    for (j=0; j < nummatch; j++)
      ptaxid[j]=buffer[j];
      
    free(buffer);

    UNPROTECT(1);
    return taxids;
}
