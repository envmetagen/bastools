#ifndef ECOPCR_H_
#define ECOPCR_H_

#include <stdio.h>
#include <inttypes.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>


//#ifndef H_apat
//#include "../libapat/apat.h"
//#endif

/*****************************************************
 * 
 *  Data type declarations
 * 
 *****************************************************/

/*
 * 
 *  Sequence types
 * 
 */

typedef struct {
	
	int32_t  taxid;
	char     AC[20];
	int32_t  DE_length;
	int32_t  SQ_length;
	int32_t  CSQ_length;
	
	char     data[1];
	
} ecoseqformat_t;

typedef struct {
	int32_t taxid;
	int32_t SQ_length;
	char    *AC;
	char    *DE;
	char    *SQ;
} ecoseq_t;

/*
 * 
 * Taxonomy taxon types
 * 
 */


typedef struct {
	int32_t  taxid;
	int32_t  rank;
	int32_t	 parent;
	int32_t  namelength;
	char     name[1];
	
} ecotxformat_t;

typedef struct ecotxnode {
	int32_t           taxid;
	int32_t           rank;
	int32_t  		      farest;
	struct ecotxnode  *parent;
	char              *name;
} ecotx_t;

typedef struct {
	int32_t count;
	int32_t maxtaxid;
  int32_t buffersize;
	ecotx_t taxon[1];
} ecotxidx_t;
 
	
/*
 * 
 * Taxonomy rank types
 * 
 */
	
typedef struct {
	int32_t count;
	char*   label[1];
} ecorankidx_t;

/*
 * 
 * Taxonomy name types
 * 
 */

typedef struct {
 	int32_t is_scientificname;
	int32_t  namelength;
	int32_t  classlength;
	int32_t  taxid;
	char     names[1];	
} econameformat_t;
 
 
 typedef struct {
 	char 	*name;
 	char 	*classname;
 	int32_t is_scientificname;
 	struct ecotxnode  *taxon;
} econame_t;

 
typedef struct {
	int32_t count;
	econame_t   names[1];
} econameidx_t;


 typedef struct {
	ecorankidx_t *ranks;
	econameidx_t *names;
	ecotxidx_t   *taxons;
} ecotaxonomy_t;

 
/*****************************************************
 * 
 *  Function declarations
 * 
 *****************************************************/

/*
 * 
 * Low level system functions
 * 
 */

int32_t is_big_endian();
int32_t swap_int32_t(int32_t);

void   *eco_malloc(int32_t chunksize,
                   const char *error_message,
                   const char *filename,
                   int32_t    line);
                   
                   
void   *eco_realloc(void *chunk,
                    int32_t chunksize,
                    const char *error_message,
                    const char *filename,
                    int32_t    line);
                    
void    eco_free(void *chunk,
                 const char *error_message,
                 const char *filename,
                 int32_t    line);
                 
void    eco_trace_memory_allocation();
void    eco_untrace_memory_allocation();

#define ECOMALLOC(size,error_message) \
	    eco_malloc((size),(error_message),__FILE__,__LINE__)
	   
#define ECOREALLOC(chunk,size,error_message) \
        eco_realloc((chunk),(size),(error_message),__FILE__,__LINE__)
        
#define ECOFREE(chunk,error_message) \
        eco_free((chunk),(error_message),__FILE__,__LINE__)
        
        


/*
 * 
 * Error managment
 * 
 */
 
  
void ecoError(int32_t,const char*,const char *,int);

#define ECOERROR(code,message) ecoError((code),(message),__FILE__,__LINE__)

#define ECO_IO_ERROR       (1)
#define ECO_MEM_ERROR      (2)
#define ECO_ASSERT_ERROR   (3)
#define ECO_NOTFOUND_ERROR (4)


/*
 * 
 * Low level Disk access functions
 * 
 */

FILE *open_ecorecorddb(const char *filename,
                       int32_t    *sequencecount,
                       int32_t    abort_on_open_error);
                       
void *read_ecorecord(FILE *,int32_t *recordSize);



/* 
 *   Read function in internal binary format
 */

FILE             *open_ecoseqdb(const char *filename,
                                int32_t    *sequencecount);
                                                                
ecoseq_t         *readnext_ecoseq(FILE *);

ecorankidx_t     *read_rankidx(const char *filename);

econameidx_t     *read_nameidx(const char *filename,ecotaxonomy_t *taxonomy);



	/**
	 * Read taxonomy data as formated by the ecoPCRFormat.py script.
	 * 
	 * This function is normaly uses internaly by the read_taxonomy
	 * function and should not be called directly.
	 * 
	 * @arg filename  path to the *.tdx file of the reformated db
	 * 
	 * @return pointer to a taxonomy index structure
	 */
 
ecotxidx_t       *read_taxonomyidx(const char *filename,const char *filename2);

ecotaxonomy_t    *read_taxonomy(const char *prefix,int32_t readAlternativeName);

ecotx_t *eco_findtaxonatrank(ecotx_t *taxon, int32_t rankidx);

ecotx_t *eco_findtaxonbytaxid(ecotaxonomy_t *taxonomy, int32_t taxid);

int eco_isundertaxon(ecotx_t *taxon, int other_taxid);

ecoseq_t *ecoseq_iterator(const char *prefix);



ecoseq_t *new_ecoseq();
int32_t   delete_ecoseq(ecoseq_t *);
ecoseq_t *new_ecoseq_with_data( char *AC,
								char *DE,
								char *SQ,
								int32_t   taxid
								);


int32_t delete_taxon(ecotx_t *taxon);
int32_t delete_taxonomy(ecotxidx_t *index);
int32_t delete_ecotaxonomy(ecotaxonomy_t *taxonomy);


int32_t rank_index(const char* label,ecorankidx_t* ranks);

//int32_t  delete_apatseq(SeqPtr pseq);
//PatternPtr buildPattern(const char *pat, int32_t error_max);
//PatternPtr complementPattern(PatternPtr pat);
//
//SeqPtr ecoseq2apatseq(ecoseq_t *in,SeqPtr out,int32_t circular);

//char *ecoComplementPattern(char *nucAcSeq);
//char *ecoComplementSequence(char *nucAcSeq);
//char *getSubSequence(char* nucAcSeq,int32_t begin,int32_t end);

ecotx_t *eco_getspecies(ecotx_t *taxon,ecotaxonomy_t *taxonomy);
ecotx_t *eco_getgenus(ecotx_t *taxon,ecotaxonomy_t *taxonomy);
ecotx_t *eco_getfamily(ecotx_t *taxon,ecotaxonomy_t *taxonomy);
ecotx_t *eco_getkingdom(ecotx_t *taxon,ecotaxonomy_t *taxonomy);
ecotx_t *eco_getsuperkingdom(ecotx_t *taxon,ecotaxonomy_t *taxonomy);

//int eco_is_taxid_ignored(int32_t *ignored_taxid, int32_t tab_len, int32_t taxid);
//int eco_is_taxid_included(ecotaxonomy_t *taxonomy, int32_t *included_taxid, int32_t tab_len, int32_t taxid);


ecotaxonomy_t *getTaxPointer(SEXP Rtaxonomy);

#endif /*ECOPCR_H_*/
