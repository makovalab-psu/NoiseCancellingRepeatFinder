// motifs.c-- operations on nucleotide motifs.

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "seq_ops.h"
#include "utilities.h"

#define  motifs_owner			// (make this the owner of its globals)
#include "motifs.h"				// interface to this module

//----------
//
// add_motif--
//	Add a motif to our motif list.
//
//----------
//
// Arguments:
//	motif*	motifList:	The list of motifs to add to.  This may be NULL.
//	char*	arg:		Command-line argument describing the motif.
//
// Returns:
//	A pointer to the motif list (which may have changed);  failures result in
//	fatality.
//
//----------

motif* add_motif
   (motif*	motifList,
	char*	arg)
	{
	char*	colon;
	char*	motifName;
	u8*		motifNucs;
	u8*		reverseNucs;
	motif*	m, *mScan;
	int		ix;

	// parse the argument

	colon = strchr(arg,':');
	if (colon == NULL)
		{
		motifName =       copy_string (arg);
		motifNucs = (u8*) copy_string (arg);
		}
	else if (colon == arg)
		{
		motifName =       copy_string (arg+1);
		motifNucs = (u8*) copy_string (arg+1);
		}
	else
		{
		motifName =       copy_prefix (arg, colon-arg);
		motifNucs = (u8*) copy_string (colon+1);
		}

	// verify that the motif is DNA (non-empty), and convert it to uppercase

	if (motifNucs[0] == 0)
		{
		fprintf (stderr, "Invalid motif (contains no DNA): \"%s\"\n", arg);
		exit (EXIT_FAILURE);
		}

	for (ix=0 ; motifNucs[ix]!=0 ; ix++)
		{
		if (isACGT (motifNucs[ix]))
			motifNucs[ix] = toUpperACGT (motifNucs[ix]);
		else
			{
			fprintf (stderr, "Invalid motif (not DNA): \"%s\"\n", motifNucs);
			exit (EXIT_FAILURE);
			}
		}

	// reverse complement to motif

	reverseNucs = (u8*) copy_string ((char*)motifNucs);
	reverse_complement (reverseNucs);
	if (strcmp ((char*)motifNucs, (char*)reverseNucs) == 0)
		{ free (reverseNucs);  reverseNucs = NULL; }

	// allocate a new record for this motif

	m = (motif*) malloc (sizeof(motif));
	if (m == NULL)
		{
		fprintf (stderr, "failed to allocate record for \"%s\"\n", motifNucs);
		exit (EXIT_FAILURE);
		}

	m->next        = NULL;
	m->name        = motifName;
	m->forwardNucs = motifNucs;
	m->reverseNucs = reverseNucs;

	// add the new record to the end of the list

	if (motifList == NULL)		// this record will be the head of the new list
		motifList = m;
	else
		{
		for (mScan=motifList ; mScan->next!=NULL ; mScan=mScan->next)
			; // do nothing
		mScan->next = m;
		}

	return motifList;
	}

//----------
//
// free_motifs--
//	Relinquish the memory allocated for our motif list.
//
//----------
//
// Arguments:
//	motif*	motifList:	The list of motifs to dispose of.
//
// Returns:
//	(nothing)
//
//----------

void free_motifs
   (motif*	motifList)
	{
	motif*	m, *mNext;

	for (m=motifList ; m!=NULL ; m=mNext)
		{
		mNext = m->next;
		if (m->name        != NULL) free (m->name);
		if (m->forwardNucs != NULL) free (m->forwardNucs);
		if (m->reverseNucs != NULL) free (m->reverseNucs);
		free (m);
		}
	}

//----------
//
// dump_motifs--
//	Write out our motif list.
//
//----------
//
// Arguments:
//	FILE*	f:			The file to write to.
//	motif*	motifList:	The list of motifs to write.
//
// Returns:
//	(nothing)
//
//----------

void dump_motifs
   (FILE*	f,
	motif*	motifList)
	{
	motif*	m;

	for (m=motifList ; m!=NULL ; m=m->next)
		{
		if (m->name        != NULL) fprintf (f, "%s: ", m->name);
		if (m->forwardNucs != NULL) fprintf (f, "%s",   m->forwardNucs);
		if (m->reverseNucs != NULL) fprintf (f, "/%s", m->reverseNucs);
		fprintf (f, "\n");
		}
	}

