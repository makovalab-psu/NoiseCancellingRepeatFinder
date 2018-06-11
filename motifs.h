#ifndef motifs_H				// (prevent multiple inclusion)
#define motifs_H

#include "utilities.h"

// establish ownership of global variables

#ifdef motifs_owner
#define global
#else
#define global extern
#endif

//----------
//
// datatypes--
//
//----------

// motifs

typedef struct motif
	{
	struct motif* next;			// next motif in linked list
	char*	name;				// identifying name for the motif
	u8*		forwardNucs;		// motif nucleotides
	u8*		reverseNucs;		// (may be NULL) reverse complement of motif
								// .. nucleotides
	} motif;

//----------
//
// prototypes for functions in this module--
//
//----------

motif* add_motif   (motif* motifList, char* arg);
void   free_motifs (motif* motifList);
void   dump_motifs (FILE* f, motif* motifList);


#undef global
#endif // motifs_H
