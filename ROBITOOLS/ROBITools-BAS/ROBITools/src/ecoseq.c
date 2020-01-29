#include "ecoPCR.h"
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

static FILE *open_seqfile(const char *prefix,int32_t index);


ecoseq_t *new_ecoseq()
{
	void *tmp;

	tmp = ECOMALLOC(sizeof(ecoseq_t),"Allocate new ecoseq structure");

	return tmp;
}

int32_t delete_ecoseq(ecoseq_t * seq)
{

	if (seq)
	{
		if (seq->AC)
			ECOFREE(seq->AC,"Free sequence AC");

		if (seq->DE)
			ECOFREE(seq->DE,"Free sequence DE");

		if (seq->SQ)
			ECOFREE(seq->SQ,"Free sequence SQ");

		ECOFREE(seq,"Free sequence structure");

		return 0;

	}

	return 1;
}

ecoseq_t *new_ecoseq_with_data( char *AC,
								char *DE,
								char *SQ,
								int32_t   taxid_idx
								)
{
	ecoseq_t *tmp;
	int32_t lstr;
	tmp = new_ecoseq();

	tmp->taxid=taxid_idx;

	if (AC)
		{
			lstr =strlen(AC);
			tmp->AC=ECOMALLOC((lstr+1) * sizeof(char),
			                  "Allocate sequence accession");
			strcpy(tmp->AC,AC);
		}

	if (DE)
		{
			lstr =strlen(DE);
			tmp->DE=ECOMALLOC((lstr+1) * sizeof(char),
			                  "Allocate sequence definition");
			strcpy(tmp->DE,DE);
		}

	if (SQ)
		{
			lstr =strlen(SQ);
			tmp->SQ=ECOMALLOC((lstr+1) * sizeof(char),
			                  "Allocate sequence data");
			strcpy(tmp->SQ,SQ);
		}
	return tmp;

}

/**
 * ?? used ??
 **/
FILE *open_ecoseqdb(const char *filename,
                    int32_t    *sequencecount)
{
	return open_ecorecorddb(filename,sequencecount,1);
}

ecoseq_t *readnext_ecoseq(FILE *f)
{
	char     *compressed=NULL;

	ecoseqformat_t *raw;
	ecoseq_t *seq;
	int32_t  comp_status;
	unsigned long int seqlength;
	int32_t  rs;
	char *c;
	int32_t i;

	raw = read_ecorecord(f,&rs);

	if (!raw)
		return NULL;

	if (is_big_endian())
	{
		raw->CSQ_length = swap_int32_t(raw->CSQ_length);
		raw->DE_length  = swap_int32_t(raw->DE_length);
		raw->SQ_length  = swap_int32_t(raw->SQ_length);
		raw->taxid      = swap_int32_t(raw->taxid);
	}

	seq = new_ecoseq();

	seq->taxid = raw->taxid;

    seq->AC    = ECOMALLOC(strlen(raw->AC) +1,
                           "Allocate Sequence Accesion number");
    strncpy(seq->AC,raw->AC,strlen(raw->AC));


    seq->DE    = ECOMALLOC(raw->DE_length+1,
                           "Allocate Sequence definition");
    strncpy(seq->DE,raw->data,raw->DE_length);

	seqlength = seq->SQ_length = raw->SQ_length;

    compressed = raw->data + raw->DE_length;

    seq->SQ = ECOMALLOC(seqlength+1,
                        "Allocate sequence buffer");

//    comp_status = uncompress((unsigned char*)seq->SQ,
//                             &seqlength,
//                             (unsigned char*)compressed,
//                             raw->CSQ_length);
//
    if (comp_status != Z_OK)
    	ECOERROR(ECO_IO_ERROR,"I cannot uncompress sequence data");

    for (c=seq->SQ,i=0;i<seqlength;c++,i++)
    	*c=toupper(*c);


	return seq;
}

/**
 * Open the sequences database (.sdx file)
 * @param	prefix	name of the database (radical without extension)
 * @param	index 	integer
 *
 * @return	file object
 */
FILE *open_seqfile(const char *prefix,int32_t index)
{
	char           filename_buffer[1024];
	int32_t        filename_length;
	FILE           *input;
	int32_t        seqcount;

	filename_length = snprintf(filename_buffer,
								1023,
	                           "%s_%03d.sdx",
	                           prefix,
	                           index);

		//	fprintf(stderr,"# Coucou %s\n",filename_buffer);


	if (filename_length >= 1024)
		ECOERROR(ECO_ASSERT_ERROR,"file name is too long");

	filename_buffer[filename_length]=0;

	input=open_ecorecorddb(filename_buffer,&seqcount,0);

	if (input)
		fprintf(stderr,"# Reading file %s containing %d sequences...\n",
				filename_buffer,
				seqcount);

	return input;
}

ecoseq_t *ecoseq_iterator(const char *prefix)
{
	static FILE    *current_seq_file= NULL;
	static int32_t current_file_idx = 1;
	static char    current_prefix[1024];
	ecoseq_t       *seq;

	if (prefix)
	{
		current_file_idx = 1;

		if (current_seq_file)
			fclose(current_seq_file);

		strncpy(current_prefix,prefix,1023);
		current_prefix[1023]=0;

		current_seq_file = open_seqfile(current_prefix,
		 							    current_file_idx);

		if (!current_seq_file)
			return NULL;

	}

	seq = readnext_ecoseq(current_seq_file);

	if (!seq && feof(current_seq_file))
	{
		current_file_idx++;
		fclose(current_seq_file);
		current_seq_file = open_seqfile(current_prefix,
		 							    current_file_idx);


		if (current_seq_file)
			seq = readnext_ecoseq(current_seq_file);
	}

	return seq;
}
