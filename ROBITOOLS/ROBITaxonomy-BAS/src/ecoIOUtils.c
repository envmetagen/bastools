#include "ecoPCR.h"
#include <stdio.h>
#include <stdlib.h>

#define SWAPINT32(x)     ((((x) << 24) & 0xFF000000) | (((x) <<  8) & 0xFF0000) | \
                         (((x) >>  8) & 0xFF00)     | (((x) >> 24) & 0xFF))


int32_t is_big_endian()
{
	int32_t i=1;
	
	return (int32_t)((char*)&i)[0];
}




int32_t swap_int32_t(int32_t i)
{
	return SWAPINT32(i);
}


/**
 * Read part of the file
 * @param	*f	the database
 * @param	recordSize the size to be read
 * 
 * @return	buffer
 */
void *read_ecorecord(FILE *f,int32_t *recordSize)
{
	static void *buffer    =NULL;
	int32_t      buffersize=0;
	int32_t      read;
	
	if (!recordSize)
		ECOERROR(ECO_ASSERT_ERROR,
		         "recordSize cannot be NULL");
		
	read = fread(recordSize,
	      		 1,
	      		 sizeof(int32_t),
	             f);
	             
	if (feof(f))
		return NULL;
	             
	if (read != sizeof(int32_t))
		ECOERROR(ECO_IO_ERROR,"Reading record size error");
		
	if (is_big_endian())
		*recordSize=swap_int32_t(*recordSize);
		
	if (buffersize < *recordSize)
	{
		if (buffer)
			buffer = ECOREALLOC(buffer,*recordSize,
			                    "Increase size of record buffer");
		else
			buffer = ECOMALLOC(*recordSize,
			                    "Allocate record buffer");
	}
	
	read = fread(buffer,
	             1,
				 *recordSize,
				 f);
				 
	if (read != *recordSize)
		ECOERROR(ECO_IO_ERROR,"Reading record data error");
		
	return buffer;	 
}





/**
 * Open the database and check it's readable
 * @param 	filename 		name of the database (.sdx, .rdx, .tbx)
 * @param 	sequencecount	buffer - pointer to variable storing the number of occurence    
 * @param 	abort_on_open_error		 	boolean to define the behaviour in case of error 
 * 										while opening the database
 * @return 	FILE type
 **/
FILE *open_ecorecorddb(const char *filename,
                       int32_t    *sequencecount,
                       int32_t    abort_on_open_error)
{
    FILE        *f;
	int32_t      read;
	
	f = fopen(filename,"rb");
	
	if (!f)
		{
			if (abort_on_open_error)
		 		ECOERROR(ECO_IO_ERROR,"Cannot open file");
		 	else
		 	{
		 		*sequencecount=0;
		 		return NULL;
		 	}
		}
		
	read = fread(sequencecount,
	      		 1,
	      		 sizeof(int32_t),
	      		 f);
	             
	if (read != sizeof(int32_t))
		ECOERROR(ECO_IO_ERROR,"Reading record size error");

	if (is_big_endian())
		*sequencecount=swap_int32_t(*sequencecount);
		
	return f;                  
}

