#ifndef loop_aligner_H			// (prevent multiple inclusion)
#define loop_aligner_H

#include "utilities.h"

// establish ownership of global variables

#ifdef loop_aligner_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef loop_aligner_owner
int        dbgLoopAlign    = false;
int        dbgLAAllocation = false;
int        dbgLACell       = false;
int        dbgLAColumn     = false;
int        dbgLATraceback  = false;
#else
global int dbgLoopAlign;
global int dbgLAAllocation;
global int dbgLACell;
global int dbgLAColumn;
global int dbgLATraceback;
#endif

//----------
//
// datatypes--
//
//----------

// alignment scoring

typedef struct lascores
	{
	u32		match;				// reward when sequence matches motif
	u32		mismatch;			// penalty when sequence mismatches motif
	u32		iOpen;				// first penalty when sequence has nt, motif doesn't
	u32		iExtend;			// other penalty when sequence has nt, motif doesn't
	u32		dOpen;				// first penalty when motif has nt, sequence doesn't
	u32		dExtend;			// other penalty when motif has nt, sequence doesn't
	u32		minScore;			// minimium score for a "suitable" alignment
	} lascores;

// alignment control

typedef struct lacontrol
	{
	lascores scoring;
	size_t	allocCells;			// number of cells allocated for dp and tb
	void*	dp;					// block of memory for DP matrices;  this is
								// .. allocated from the heap
	void*	tb;					// block of memory for traceback matrices;
								// .. this is NOT allocated separately from the
								// .. heap (it's carved from the same block as
								// .. dp
	u32		textChars;			// number of characters allocated for seqText
								// .. and qryText (combined, including 
								// .. terminating zeros)
	char*	seqText;			// place to build alignment text, for sequence
	char*	qryText;			// place to build alignment text, for query
	} lacontrol;

// alignment results

typedef struct alignment
	{
	int		active;				// false => the rest of the record is invalid
	char	strand;
	s32		score;
	u32		seqBaseCount;
	u32		qryBaseCount;
	char*	seqText;
	char*	qryText;
	u32		seqStart;
	u32		seqEnd;
	u32		mCount;				// number of matches in the alignment
	u32		mmCount;			// number of mismatches
	u32		iCount;				// number of inserted bases
	u32		dCount;				// number of deleted bases
	float	matchRatio;			// ratio of m to m+mm+i+d
	} alignment;

#ifndef loop_aligner_owner
global const alignment noAlignment;
#else
const alignment noAlignment =
	{
	false,						// active
	'?',						// strand (? means strand not known)
	0,							// score
	0,							// seqBaseCount
	0,							// qryBaseCount
	NULL,						// seqText
	NULL,						// qryText
	0,							// seqStart
	0,							// seqEnd
	0,							// mCount
	0,							// mmCount
	0,							// iCount
	0,							// dCount
	0.0							// matchRatio
	};
#endif // if/not loop_aligner_owner

// miscellany

#define noScore ((s32) 0x80000000)

//----------
//
// prototypes for functions in this module--
//
//----------

void      init_loop_align    (lacontrol* control, lascores* scoring,
                              u32 allocBytes);
void      free_loop_align    (lacontrol* control);
void      free_alignment     (alignment a);
alignment copy_alignment     (alignment a);
alignment loop_align         (lacontrol* control, u8* seq, u8* motif);
alignment loop_align_segment (lacontrol* control, u8* seq, u32 seqLen, u8* motif);
void      rescore_alignment  (lascores* scoring, alignment* a);
u32       alignment_crc      (char* seqName, char* motif, alignment* a);

#undef global
#endif // loop_aligner_H
