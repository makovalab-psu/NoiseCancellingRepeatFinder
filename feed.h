#ifndef feed_H					// (prevent multiple inclusion)
#define feed_H

#include "utilities.h"

// establish ownership of global variables

#ifdef feed_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef feed_owner
int		dbgShowFeed   = false;
#else
global int dbgShowFeed;
#endif

//----------
//
// datatypes--
//
//----------

// a feed (a series of nucleotide sequences)

typedef struct feed
	{
	int			state;			// one of feed_needInput, etc.
	char*		filename;		// our copy of the filename;  unlike other
								// .. dynamic elements in this object, this is
								// .. allocated within this block
	u32			lineNumber;		// number of line being read from the file
	u32			charNumber;		// number of character on line being read
	u32			sequenceNumber;	// number of sequence (in file) being read
	char		ungetBuffer;	// single character buffer of characters we've
								// .. read from the file but did not want;  a
								// .. NUL indicates an empty buffer
	int			eofReached;		// true => reader has reached end-of-file
	FILE*		f;				// the file or stream the sequences are to be
								// .. read from
	int			fIsStdin;		// true => out input is from stdin
	int			reportStrands;	// which strands to report
								//   +1   => report forward strand only
								//   -1   => report reverse strand only
								//   o.w. => report forward strand, then reverse
								//           .. strand
	u32			size;			// number of bytes allocated for nt[]
	u32			len;			// number of charaters in nt[], not including
								// .. the zero terminator
	u8*			nt;				// zero-terminated array of ACGTN (upper and
								// .. lower case); this is allocated as a
								// .. separate block in the heap
	u32			nameSize;		// number of bytes allocated for name[]
	char*		name;			// name of the sequence (e.g. "chr1"); this is
								// .. allocated as a separate block in the heap
	int			isRevComp;		// true => piece is reverse complement;  note
								// that start and end are always relative to
								// to the forward strand
	u32			start;			// start of piece within that sequence;  origin
								// .. zero
	u32			end;			// end of piece within that sequence;  the
								// .. interval is half-open, so end is one
								// .. beyond the last position within the
								// .. piece
	} feed;

// feed states

enum
	{
	feed_needInput,				// next sequence must come from input
	feed_providedAsForward,		// last sequence was provided as forward strand
	feed_providedAsReverseComp,	// last sequence was provided as reverse comp
	feed_exhausted				// there are no more sequences
	};

//----------
//
// prototypes for functions in this module--
//
//----------

feed* open_fasta_feed          (char* filename, int reportStrands,
                                u32 ntBytes, u32 nameBytes);
void  free_fasta_feed          (feed* src);
int   another_sequence_in_feed (feed* src);
void  dump_feed_sequence       (FILE* f, feed* src, int wrapLength);


#undef global
#endif // feed_H
