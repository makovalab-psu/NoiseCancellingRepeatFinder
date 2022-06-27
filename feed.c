// feed.c-- operations on series of nucleotide sequences.

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "seq_ops.h"
#include "utilities.h"

#define  feed_owner				// (make this the owner of its globals)
#include "feed.h"				// interface to this module

// miscellany

#define defaultNameSize 100
#define defaultNtSize   (10*1000)

//----------
//
// prototypes for private functions
//
//----------

static int load_fasta (feed* src);

//----------
//
// open_fasta_feed--
//	Prepare to read sequences one-by-one from a fasta file.
//
//----------
//
// Arguments:
//	char*	filename:		The name of the file to read from.  If this is NULL
//							.. or an empty string, we'll read from stdin.
//	int		reportStrands:	+1 or '+' => report forward strand only
//							-1 or '-' => report reverse strand only
//							otherwise => report forward strand, then reverse
//	u32		ntBytes:		How much to allocate for the nt[] array.  Zero
//							.. means just use the internal default.
//	u32		nameBytes:		How much to allocate for the name[] array.  Zero
//							.. means just use the internal default.
//
// Returns:
//	A pointer to a feed record, allocated from the heap;  the caller is
//	responsible for deallocating this memory;  failure to allocate results in
//	program termination.
//
//----------

feed* open_fasta_feed
   (char*	filename,
	int		reportStrands,
	u32		ntBytes,
	u32		nameBytes)
	{
	FILE*	f;
	int		fIsStdin;
	feed*	src;
	u8*		src0;
	u32		filenameOffset;
	char*	name;
	u8*		nt;
	u32		numBytes;

	// tighten relaxed representation of reportStrands

	if      (reportStrands == '+')                        reportStrands =  1;
	else if (reportStrands == '-')                        reportStrands = -1;
	else if ((reportStrands < -1) || (reportStrands > 1)) reportStrands =  0;

	// open the file

	if ((filename == NULL) || (*filename == 0))
		{
		f        = stdin;
		fIsStdin = true;
		filename = "(stdin)";
		}
	else
		{
		f = fopen (filename, "rt");
		if (f == NULL) goto cant_open_file;
		fIsStdin = false;
		}

	// allocate the structure, including a copy of the filename

	numBytes =  sizeof(feed);

	numBytes =  round_up_16(numBytes);  filenameOffset = numBytes;
	numBytes += (strlen(filename)+1) * sizeof(char);

	src = (feed*) malloc (numBytes);
	if (src == NULL) goto cant_allocate_struture;
	src0 = (u8*) src;

	src->filename = (char*) (src0 + filenameOffset);

	strcpy (/*to*/ src->filename, /*from*/ filename);
	src->state    = feed_needInput;
	src->f        = f;
	src->fIsStdin = fIsStdin;

	// allocate the name

	if (nameBytes != 0) numBytes = nameBytes;
	               else numBytes = defaultNameSize+1;
	name = (char*) malloc (numBytes);
	if (name == NULL) goto cant_allocate_name;
	src->name     = name;
	src->nameSize = numBytes;
	name[0] = 0;

	// allocate the nucleotide array

	if (ntBytes != 0) numBytes = ntBytes;
	             else numBytes = defaultNtSize+1;
	nt = (u8*) malloc (numBytes);
	if (nt == NULL) goto cant_allocate_nukes;
	src->nt   = nt;
	src->size = numBytes;
	src->len  = 0;
	nt[0] = 0;

	// intialize remaining fields

	src->lineNumber     = 0;
	src->charNumber     = 0;
	src->sequenceNumber = 0;
	src->ungetBuffer    = 0;
	src->eofReached     = false;
	src->reportStrands  = reportStrands;
	src->isRevComp      = false;
	src->start          = 0;
	src->end            = 0;

	// success

	return src;

	//////////
	// failure exits
	//////////

cant_open_file:
	fprintf (stderr, "failed to open \"%s\"\n", filename);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

cant_allocate_struture:
	fprintf (stderr, "failed to allocate feed for \"%s\" (%s bytes)\n",
	                 filename, ucommatize(numBytes));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

cant_allocate_name:
	fprintf (stderr, "failed to allocate name for \"%s\" (%s bytes)\n",
	                 filename, ucommatize(numBytes));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

cant_allocate_nukes:
	fprintf (stderr, "failed to allocate nucleotide array for \"%s\" (%s bytes)\n",
	                 filename, ucommatize(numBytes));
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)
	}

//----------
//
// free_fasta_feed--
//	Deallocate a feed record.  All subordinate elements are deallocated and
//	closed.
//
//----------
//
// Arguments:
//	feed*	src:	The feed to deallocate.
//
// Returns:
//	(nothing)
//
//----------

void free_fasta_feed
   (feed*	src)
	{
	if (src == NULL) return;

	if ((!src->fIsStdin) && (src->f != NULL)) fclose (src->f);

	if (src->nt   != NULL) free (src->nt);
	if (src->name != NULL) free (src->name);
	free (src);
	}

//----------
//
// another_sequence_in_feed--
//	Ask for another sequence from a feed.
//
// This routine will return a single sequence.  Repeated calls will return 
// all sequences from the feed's file, as well as their reverse complements
// (depending on src->reportStrands).
//
// The returned text consists only of A,C,G,T, and N (upper or lower case).
// Note that any IUPAC ambiguous nucleotides are converted to N.
//
//----------
//
// Arguments:
//	feed*	src:	The feed providing sequences.
//
// Returns:
//	true if there was another sequence-- as a side effect, that sequence is
//	loaded into the feed record;  false if sequence input has been exhausted
//
//----------

int another_sequence_in_feed
   (feed*	src)
	{
	if (src == NULL) return false;
	if (src->state == feed_exhausted) return false;

	// if we provided the last sequence as reverse complement, then we need a
	// fresh sequence from the file;  or, if we're only doing a single strand
	// from each sequence, then we again need a fresh sequence from the file

	if (src->reportStrands == 0)
		{
		if (src->state == feed_providedAsReverseComp)
			src->state = feed_needInput;
		}
	else
		src->state = feed_needInput;

	// if we provided the last sequence as forward (and we're to report both
	// strands), then we can just provide it as reverse complement

	if (src->state == feed_providedAsForward)
		{
		reverse_complement (src->nt);
		src->isRevComp = true;
		src->state     = feed_providedAsReverseComp;
		return true;
		}

	// otherwise, we need a fresh sequence from the file

	if (!load_fasta (src))  // file contains no more sequences
		{
		src->name[0]   = 0;
		src->len       = 0;
		src->nt[0]     = 0;
		src->isRevComp = false;
		src->state     = feed_exhausted;
		return false;
		}

	if (src->reportStrands >= 0)
		src->state = feed_providedAsForward;
	else // if (src->reportStrands == -1)
		{
		reverse_complement (src->nt);
		src->isRevComp = true;
		src->state     = feed_providedAsReverseComp;
		}

	return true;
	}

//----------
//
// load_fasta--
//	Read the next sequence from a fasta file.
//
//----------
//
// Arguments:
//	feed*	src:	The feed providing sequences.
//
// Returns:
//	true if there was another sequence in the file;  false if the file has been
//	exhausted
//
//----------

// make_enough_name--
//	make sure there's enough room for name[nameIx]

#define make_enough_name                                                      \
	if (nameIx >= src->nameSize)                                              \
		{                                                                     \
		u32 numBytes = round_up_1K(src->nameSize + src->nameSize/10);         \
		name = realloc (name, numBytes);                                      \
		if (name == NULL) goto cant_reallocate_name;                          \
		src->nameSize = numBytes;                                             \
		src->name     = name;                                                 \
		}

// make_enough_nt--
//	make sure there's enough room for nt[pos]

#define make_enough_nt                                                       \
	if (pos >= src->size)                                                    \
		{                                                                    \
		u32 numBytes = round_up_1K(src->size + src->size/10);                \
		nt = realloc (nt, numBytes);                                         \
		if (nt == NULL) goto cant_reallocate_nukes;                          \
		src->size = numBytes;                                                \
		src->nt   = nt;                                                      \
		}

// parsing states

#define inName   1
#define inHeader 2

// text_from_fasta--

static int load_fasta
   (feed*	src)
	{
	char*	name = src->name;
	u8*		nt   = src->nt;
	FILE*	f    = src->f;
	int		ch, prevCh;
	int		state, isFirstLine;
	u8		nuke;
	u32		numBytes = 0;
	u32		pos, nameIx;

	//fprintf (stderr, "load_fasta\n");

	if (src->eofReached) return false;

	src->sequenceNumber++;

	// scan the file, collecting nukes

	pos    = 0;
	nameIx = 0;

	isFirstLine = true;
	state = 0;
	prevCh = '\n';
	while (true)
		{
		ch = src->ungetBuffer;
		if (ch != 0) src->ungetBuffer = 0;
		        else ch = getc (f);
		if (ch == EOF)
			{ src->eofReached = true;  break; }

		//fprintf (stderr, "  ch=%02X\n", ch);

		if (ch == '\n')
			{
			src->lineNumber++;
			src->charNumber = 0;

			if ((state == inName) || (state == inHeader))
				{
				//fprintf (stderr, "  nameIx=%u nameSize=%u for 00\n", nameIx, src->nameSize);
				make_enough_name;
				name[nameIx++] = 0;	// terminator on sequence's name
				}

			state  = 0;
			prevCh = ch;
			isFirstLine = false;
			continue;
			}

		src->charNumber++;

		if (ch == '>')
			{
			if (prevCh != '\n') goto bad_nucleotide;

			if (!isFirstLine)
				{ src->ungetBuffer = '>';  break; }

			isFirstLine = false;
			state       = inName;
			nameIx      = 0;

			prevCh = ch;
			continue;
			}

		if (state == inName)
			{
			if (!isspace(ch))
				{
				//fprintf (stderr, "  nameIx=%u nameSize=%u for %02X\n", nameIx, src->nameSize, ch);
				make_enough_name;
				name[nameIx++] = ch;
				}
			else
				state = inHeader;
			continue;
			}

		if (state == inHeader)
			continue;

		nuke = ch;
		if (!isACGTN(nuke))
			{
			if (!isIUPAC(nuke)) goto bad_nucleotide;
			if (isupper((char) nuke)) nuke = 'N';
			                     else nuke = 'n';
			}

		//fprintf (stderr, "  pos=%u size=%u\n", pos, src->size);
		make_enough_nt;
		nt[pos++] = nuke;
		}

	if (state != 0) goto internal_error_state;	// $$$ sequence empty
												// .. change how we deal with
												// .. empty sequences

	//fprintf (stderr, "  pos=%u size=%u (end)\n", pos, src->size);
	make_enough_nt;
	nt[pos] = 0;		// (string terminator)
	src->len = pos;

	src->isRevComp = false;
	src->start     = 0;
	src->end       = pos;

	// if the sequence had no name, give it one

	if (name[0] == 0)
		{
		char* nameFmt = "sequence%u";
		nameIx = 1 + snprintf (NULL, 0, nameFmt, src->sequenceNumber);
		make_enough_name;
		sprintf (name, nameFmt, src->sequenceNumber);
		}

	// success

	return true;

	//////////
	// failure exits
	//////////

cant_reallocate_name:
	fprintf (stderr, "failed to reallocate name for \"%s\" (%s bytes)\n",
	                 src->filename, ucommatize(numBytes));
	exit(EXIT_FAILURE);
	return false; // (never reaches here)

cant_reallocate_nukes:
	fprintf (stderr, "failed to reallocate nucleotide array for \"%s\" (%s bytes)\n",
	                 src->filename, ucommatize(numBytes));
	exit(EXIT_FAILURE);
	return false; // (never reaches here)

bad_nucleotide:
	if (isalnum(ch))
		fprintf (stderr, "bad nucleotide '%c' (line %u, column %u) in \"%s\"\n",
		                 ch, src->lineNumber, src->charNumber, src->filename);
	else
		fprintf (stderr, "bad nucleotide 0x%02X (line %u, column %u) in \"%s\"\n",
		                 ch, src->lineNumber, src->charNumber, src->filename);
	exit(EXIT_FAILURE);
	return false; // (never reaches here)

internal_error_state:
	fprintf (stderr, "internal error for %s, state is %d\n",
	                 src->filename, state);
	exit(EXIT_FAILURE);
	return false; // (never reaches here)
	}

//----------
//
// dump_feed_sequence--
//	(For debugging) Write a feed's current sequence.
//
//----------
//
// Arguments:
//	FILE*	f:			The file to write to.
//	feed*	src:		The feed providing sequences.
//	int		wrapLength: Number of nucleotides for each line.
//
// Returns:
//	(nothing)
//
//----------

void dump_feed_sequence
   (FILE*	f,
	feed*	src,
	int 	wrapLength)
	{
	u32		pos;

	fprintf (f, ">%s %c\n", src->name, src->isRevComp?'-':'+');
	if (src->len == 0) return;

	for (pos=0 ; pos<src->len ; pos++)
		{
		if (wrapLength != 0)
			{ if ((pos != 0) && (pos % wrapLength == 0)) fprintf (f, "\n"); }
		fprintf (f, "%c", src->nt[pos]);
		}
	fprintf (f, "\n");
	}
