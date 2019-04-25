// loop_aligner.c-- find high-scoring local alignments between a looped query
// motif and a sequence

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "seq_ops.h"
#include "utilities.h"

#define  loop_aligner_owner		// (make this the owner of its globals)
#include "loop_aligner.h"		// interface to this module

#define max(u,v) (((u)>=(v))?(u):(v))	// danger: possible side effects

//----------
//
// datatypes--
//
//----------

typedef s32 dpcell;				// DP matrix cells are just alignment scores

typedef u8  tbcell;				// traceback cells are packed bit fields:
								//   00000000 => empty cell
								//   xxxxxx01 => cell links to "M" cell
								//   xxxxxx10 => cell links to "I" cell
								//   xxxxxx11 => cell links to "D" cell
								//   xxxxx0xx => cell links to colIx-1
								//   xxxxx1xx => cell links to colIx
								//   xxx00xxx => cell links to rowIx-1
								//   xxx01xxx => cell links to rowIx
								//   xxx01xxx => cell links to row motifEndIx

#define tbEmpty   0

#define tbWhich   0x03
#define tbWhichM  0x01
#define tbWhichI  0x02
#define tbWhichD  0x03

#define tbCol     0x04
#define tbColPrev 0x00
#define tbColSame 0x04

#define tbRow     0x18
#define tbRowPrev 0x00
#define tbRowSame 0x08
#define tbRowLoop 0x10

//----------
//
// prototypes for private functions--
//
//----------

static void enough_cells (lacontrol* control, u32 queryLen, u32 sequenceLen);

//----------
//
// init_loop_align--
//	This should be called before the first call of loop_align().
//
//----------
//
// Arguments:
//	lacontrol*	control:	The aligner control block to initialize.
//	lascores*	scoring:	Alignment scoring definitions.
//	u32			allocBytes:	How much to allocate for the internal data
//							.. structures.  Zero means just use the internal
//							.. default.
//
// Returns:
//	(nothing)
//
//----------

void init_loop_align
   (lacontrol*	control,
	lascores*	scoring,
	u32			allocBytes)
	{
	int			queryLen    = 10;
	int			sequenceLen = 2000;

	if (allocBytes != 0)
		{
		queryLen = (allocBytes + sequenceLen) / (sequenceLen+1);
		if (queryLen < 5) queryLen = 5;
		sequenceLen = (allocBytes + queryLen) / (queryLen+1);
		}

	control->scoring = *scoring;

	control->allocCells = 0;
	control->dp         = NULL;
	control->tb         = NULL;

	control->textChars  = 0;
	control->seqText    = NULL;
	control->qryText    = NULL;

	enough_cells (control, queryLen, sequenceLen);
	}

//----------
//
// free_loop_align--
//	This should be called after the final call of loop_align().
//
//----------
//
// Arguments:
//	lacontrol*	control:	The aligner control block.
//
// Returns:
//	(nothing)
//
//----------

void free_loop_align
   (lacontrol*	control)
	{
	if (control->dp != NULL)
		{ free (control->dp);       control->dp = NULL;       control->allocCells = 0; }
	if (control->seqText != NULL)
		{ free (control->seqText);  control->seqText = NULL;  control->textChars = 0; }
	}

//----------
//
// free_alignment--
//	Dispose of any memory used by an alignment, presumably created by a call to
//	loop_align().
//
//----------
//
// Arguments:
//	alignment	a:	The alignment to dispose of.
//
// Returns:
//	(nothing)
//
//----------

void free_alignment
   (alignment	a)
	{
	if (a.seqText != NULL) free (a.seqText);
	if (a.qryText != NULL) free (a.qryText);
	}

//----------
//
// copy_alignment--
//	Make a copy of an alignment.
//
//----------
//
// Arguments:
//	alignment	a:	The alignment to copy.
//
// Returns:
//	A copy of the alignment;  noAlignment if no alignment scores high enough.
//	The caller is responsible for deallocating any memory allocated for this
//	result, by eventually calling free_alignment().
//
//----------

alignment copy_alignment
   (alignment	a)
	{
	alignment	aCopy;

	aCopy = a;
	if (a.seqText != NULL) aCopy.seqText = copy_string (a.seqText);
	if (a.qryText != NULL) aCopy.qryText = copy_string (a.qryText);

	return aCopy;
	}

//----------
//
// loop_align, loop_align_segment--
//	Find the *single* highest-scoring local alignment between the looped query
//	motif and a sequence.
//
//----------
//
// Arguments:
//	lacontrol*	control:	Alignment control (including scoring).
//	u8*			seq:		The sequence to search in.
//	u32			seqLen:		(loop_align_segment only) Length of the sequence;
//	u8*			motif:		The looped motif to search for.
//
// Returns:
//	The alignment;  noAlignment if no alignment scores high enough.  The caller
//	is responsible for deallocating any memory allocated for this result, by
//	eventually calling free_alignment().
//
//----------

static void print_debug_column (u8* motif,
                                dpcell* dpColM, dpcell* dpColI, dpcell* dpColD,
                                tbcell* tbColM, tbcell* tbColI, tbcell* tbColD,
                                char seqNuc, u32 colIx);
static char* debug_traceback_string (tbcell tbVal);
static char  debug_traceback_which  (tbcell tbVal);

// loop_align--

alignment loop_align
   (lacontrol*	control,
	u8*			seq,
	u8*			motif)
	{
	return loop_align_segment (control, seq, strlen((char*)seq), motif);
	}
	
// loop_align_segment--

alignment loop_align_segment
   (lacontrol*	control,
	u8*			seq,
	u32			seqLen,
	u8*			motif)
	{
	u32			m, mm, iop, iex, dop, dex;
	u32			qryLen, motifStartIx, motifEndIx;
	size_t		numCells;
	dpcell*		dpM, *dpI, *dpD;
	dpcell*		dpColM,  *dpColI,  *dpColD;
	dpcell*		dpPrevM, *dpPrevI, *dpPrevD;
	tbcell*		tbM, *tbI, *tbD;
	tbcell*		tbColM, *tbColI, *tbColD;
	u32			offsetColToCol;
	u32			colIx, rowIx, bestColIx, bestRowIx, predRowIx;
	u32			startColIx, endColIx, endRowIx, c, r, cellIx;
	u8			seqNuc, qryNuc;
	s32			bestScore;
	s32			mForCell;
	s32			mScore, mFromM, mFromI, mFromD, loopMFromM, loopMFromI, loopMFromD;
	s32			iScore, iFromM, iFromI;
	s32			dScore, dFromM, dFromD, prevDScore;
	tbcell		predRowTb, tbVal, whichTb;
	char*		seqTextHead, *seqTextScan, *qryTextHead, *qryTextScan;
	alignment	a;

	if (dbgLoopAlign != 0)
		{
		u32 motifLen = strlen((char*)motif);
		fprintf (stderr, "loop_align_segment(");
		if (seqLen   < 50) fprintf (stderr, "%s#%u", seq, seqLen);
		              else fprintf (stderr, "seqlen_%u", seqLen);
		if (motifLen < 50) fprintf (stderr, ",%s", motif);
		              else fprintf (stderr, ",motiflen_%u", motifLen);
		fprintf (stderr, ")\n");
		}

	m   = control->scoring.match;
	mm  = control->scoring.mismatch;
	iop = control->scoring.iOpen;
	iex = control->scoring.iExtend;
	dop = control->scoring.dOpen;
	dex = control->scoring.dExtend;

	qryLen = strlen ((char*) motif);
	motifStartIx = 1;
	motifEndIx   = qryLen;

	// carve up memory for the DP and traceback matrices

	enough_cells (control, qryLen, seqLen);

	numCells = (qryLen+1) * (size_t) (seqLen+1);

	dpM = (dpcell*) control->dp;
	dpI = dpM + numCells;
	dpD = dpI + numCells;
	tbM = (tbcell*) control->tb;
	tbI = tbM + numCells;
	tbD = tbI + numCells;

	offsetColToCol = qryLen + 1;
	dpColM = dpM;  dpColI = dpI;  dpColD = dpD;
	tbColM = tbM;  tbColI = tbI;  tbColD = tbD;

	if (dbgLAAllocation)
		{
		fprintf (stderr, "  numCells=%s 0x%016lX\n", ucommatize(numCells), numCells*sizeof(dpcell));
		fprintf (stderr, "  dpM=%p dpI=%p dpD=%p\n", dpM, dpI, dpD);
		fprintf (stderr, "  tbM=%p tbI=%p tbD=%p\n", tbM, tbI, tbD);
		}

	// perform the DP process, column-by-column

	bestScore = 0;
	bestColIx = bestRowIx = 0;

	if (dbgLoopAlign != 0)
		fprintf (stderr, "  column [0]\n");

	colIx = 0;
	for (rowIx=0 ; rowIx<=qryLen ; rowIx++)
		{
		dpColM[rowIx] = dpColI[rowIx] = dpColD[rowIx] = 0;
		tbColM[rowIx] = tbColI[rowIx] = tbColD[rowIx] = tbEmpty;
		}

	if (dbgLACell || dbgLAColumn)
		print_debug_column (motif,dpColM,dpColI,dpColD,tbColM,tbColI,tbColD,'$',0);

	for (colIx=1 ; colIx<=seqLen ; colIx++)
		{
		seqNuc = toUpperACGTN(seq[colIx-1]);

		if ((dbgLoopAlign != 0) && (colIx % dbgLoopAlign == 0))
			fprintf (stderr, "  column [%d] %c\n", colIx, seqNuc);
//		else if (dbgLAAllocation)
//			{
//			fprintf (stderr, "  column [%d] %c", colIx, seqNuc);
//			// $$$ limiting this to u32 could result in wraparound of these reported values
//			fprintf (stderr, "  dM=%d dI=%d dD=%d", (u32)(dpColM-dpM), (u32)(dpColI-dpI), (u32)(dpColD-dpD));
//			fprintf (stderr, "  tM=%d tI=%d tD=%d", (u32)(tbColM-tbM), (u32)(tbColI-tbI), (u32)(tbColD-tbD));
//			fprintf (stderr, "\n");
//			}

		dpPrevM = dpColM;  dpColM += offsetColToCol;  tbColM += offsetColToCol;
		dpPrevI = dpColI;  dpColI += offsetColToCol;  tbColI += offsetColToCol;
		dpPrevD = dpColD;  dpColD += offsetColToCol;  tbColD += offsetColToCol;

		dpColM[0] = dpColI[0] = dpColD[0] = 0;
		tbColM[0] = tbColI[0] = tbColD[0] = tbEmpty;

		// pass 1 over the motif -- (mis)match, insertion, deletion

		for (rowIx=1 ; rowIx<=qryLen ; rowIx++)
			{
			qryNuc = toUpperACGTN(motif[rowIx-1]);

			// (mis)match

			if (qryNuc == seqNuc) mForCell = m;
			                 else mForCell = -mm;
			mFromM = dpPrevM[rowIx-1] + mForCell;
			mFromI = dpPrevI[rowIx-1] + mForCell;
			mFromD = dpPrevD[rowIx-1] + mForCell;

			loopMFromM = loopMFromI = loopMFromD = noScore;
			if (rowIx == motifStartIx)  // (mis)match from loopback
				{
				loopMFromM = dpPrevM[motifEndIx] + mForCell;
				loopMFromI = dpPrevI[motifEndIx] + mForCell;
				loopMFromD = dpPrevD[motifEndIx] + mForCell;
				}

			mScore = max(0,mFromM);
			mScore = max(mScore,mFromI);
			mScore = max(mScore,mFromD);
			if (loopMFromM != noScore)
				{
				mScore = max(mScore,loopMFromM);
				mScore = max(mScore,loopMFromI);
				mScore = max(mScore,loopMFromD);
				}
			dpColM[rowIx] = mScore;
			if (mScore > bestScore)
				{ bestScore = mScore;  bestColIx = colIx;  bestRowIx = rowIx; }

			if      (mScore == mFromM)     tbColM[rowIx] = tbWhichM + tbColPrev + tbRowPrev;
			else if (mScore == mFromI)     tbColM[rowIx] = tbWhichI + tbColPrev + tbRowPrev;
			else if (mScore == mFromD)     tbColM[rowIx] = tbWhichD + tbColPrev + tbRowPrev;
			else if (mScore == loopMFromM) tbColM[rowIx] = tbWhichM + tbColPrev + tbRowLoop;
			else if (mScore == loopMFromI) tbColM[rowIx] = tbWhichI + tbColPrev + tbRowLoop;
			else if (mScore == loopMFromD) tbColM[rowIx] = tbWhichD + tbColPrev + tbRowLoop;
			else                           tbColM[rowIx] = tbEmpty;

			// insertion

			iFromM = dpPrevM[rowIx] - iop;
			iFromI = dpPrevI[rowIx] - iex;
			dpColI[rowIx] = iScore = max(iFromM,iFromI);
			if (iScore > bestScore)
				{ bestScore = iScore;  bestColIx = colIx;  bestRowIx = rowIx; }

			if (iScore == iFromM) tbColI[rowIx] = tbWhichM + tbColPrev + tbRowSame;
			                 else tbColI[rowIx] = tbWhichI + tbColPrev + tbRowSame;

			// deletion

			dFromM = dpColM[rowIx-1] - dop;
			dFromD = dpColD[rowIx-1] - dex;
			dpColD[rowIx] = dScore = max(dFromM,dFromD);
			if (dScore > bestScore)
				{ bestScore = dScore;  bestColIx = colIx;  bestRowIx = rowIx; }

			if (dScore == dFromM) tbColD[rowIx] = tbWhichM + tbColSame + tbRowPrev;
			                 else tbColD[rowIx] = tbWhichD + tbColSame + tbRowPrev;
			}

		// pass 2 over the motif (except the motif end)-- deletions only;  we'll
		// exit the loop early if the deletions aren't improving the score

		for (rowIx=motifStartIx ; rowIx<motifEndIx ; rowIx++)
			{
			qryNuc = toUpperACGTN(motif[rowIx-1]);
			prevDScore = dpColD[rowIx];

			if (rowIx == motifStartIx)
				{
				predRowIx = motifEndIx;  // deletion from loopback
				predRowTb = tbRowLoop;
				}
			else
				{
				predRowIx = rowIx-1;     // deletion from cell below
				predRowTb = tbRowPrev;
				}

			dFromM = dpColM[predRowIx] - dop;
			dFromD = dpColD[predRowIx] - dex;
			dScore = max(dFromM,dFromD);

			if (dScore <= prevDScore)    // no (more) improvement from loopback
				break;

			dpColD[rowIx] = dScore;
			if (dScore > bestScore)
				{ bestScore = dScore;  bestColIx = colIx;  bestRowIx = rowIx; }

			if (dScore == dFromM) tbColD[rowIx] = tbWhichM + tbColSame + predRowTb;
			                 else tbColD[rowIx] = tbWhichD + tbColSame + predRowTb;
			}

		if (dbgLACell || dbgLAColumn)
			print_debug_column (motif,dpColM,dpColI,dpColD,tbColM,tbColI,tbColD,seqNuc,colIx);
		}

	// if the best score isn't good enough, quit now;  otherwise we'll report
	// the alignment

	if ((u32) bestScore < control->scoring.minScore) return noAlignment;

	a = noAlignment;
	a.active = true;
	a.score  = bestScore;

	// otherwise, traceback *an* optimal alignment

	seqTextHead = (char*) &control->seqText[control->textChars-1];  *seqTextHead = 0;
	qryTextHead = (char*) &control->qryText[control->textChars-1];  *qryTextHead = 0;


	colIx = endColIx = bestColIx;  startColIx = 0;
	rowIx = endRowIx = bestRowIx;

	whichTb = tbWhichM;
	while (whichTb != 0)
		{
		startColIx = colIx;

		cellIx = (colIx*offsetColToCol) + rowIx;
		if      (whichTb == tbWhichI)  tbVal = tbI[cellIx];
		else if (whichTb == tbWhichD)  tbVal = tbD[cellIx];
		else   /*whichTb == tbWhichM*/ tbVal = tbM[cellIx];

		if (dbgLATraceback)
			fprintf (stderr, "%c[%d][%d] --> %s\n",
			                  debug_traceback_which(whichTb),colIx,rowIx,
			                  debug_traceback_string(tbVal));

		whichTb = tbVal & tbWhich;
		if (whichTb == tbEmpty) break;

		if      ((tbVal & tbCol) == tbColPrev) c = colIx-1;
		                                  else c = colIx;
		if      ((tbVal & tbRow) == tbRowPrev) r = rowIx-1;
		else if ((tbVal & tbRow) == tbRowLoop) r = motifEndIx;
		                                  else r = rowIx;

		if (c != colIx) *(--seqTextHead) = (char) toUpperACGTN(seq[colIx-1]);
		           else *(--seqTextHead) = '-';
		if (r != rowIx) *(--qryTextHead) = (char) toUpperACGTN(motif[rowIx-1]);
		           else *(--qryTextHead) = '-';

		colIx = c;  rowIx = r;
		}

	seqTextScan = seqTextHead;    // the following loop is similar to
	qryTextScan = qryTextHead;    // .. rescore_alignment()
	a.mCount = a.mmCount = a.iCount = a.dCount = 0;
	for ( ; *seqTextScan!=0 ; seqTextScan++,qryTextScan++)
		{
		if (*qryTextScan == '-')
			a.iCount++;
		else if (*seqTextScan == '-')
			a.dCount++;
		else if (*qryTextScan == *seqTextScan)
			a.mCount++;
		else
			{
			*qryTextScan = (char) toLowerACGTN((u8) *qryTextScan);
			a.mmCount++;
			}
		}

	if (dbgLATraceback)
		{
		fprintf (stderr, "  %s\n", seqTextHead);
		fprintf (stderr, "  %s\n", qryTextHead);
		}

	a.seqText      = copy_string (seqTextHead);
	a.qryText      = copy_string (qryTextHead);
	a.seqStart     = startColIx;
	a.seqEnd       = endColIx;
	a.seqBaseCount = endColIx - startColIx;
	a.qryBaseCount = a.mCount + a.mmCount + a.dCount;
	a.matchRatio   = ((float)a.mCount) / (a.mCount + a.mmCount + a.iCount + a.dCount);

	return a;
	}


// debug printers for loop_align()

static void print_debug_column
   (u8*		motif,
	dpcell*	dpColM,
	dpcell*	dpColI,
	dpcell*	dpColD,
	tbcell*	tbColM,
	tbcell*	tbColI,
	tbcell*	tbColD,
	char	seqNuc,
	u32		colIx)
	{
	u32		qryLen = strlen ((char*) motif);
	char	s[100];
	u8		qryNuc;
	u32		rowIx;

	fprintf (stderr, "\n");

	fprintf (stderr, "         ");
	for (rowIx=0 ; rowIx<=qryLen ; rowIx++)
		{
		sprintf (s, "[%d]", rowIx);
		fprintf (stderr, " %-4s", s);
		}
	fprintf (stderr, "\n");

	fprintf (stderr, "          ");
	for (rowIx=0 ; rowIx<=qryLen ; rowIx++)
		{
		if (rowIx == 0) qryNuc = '$';
		           else qryNuc = toUpperACGTN(motif[rowIx-1]);
		fprintf (stderr, " %-4c", qryNuc);
		}
	fprintf (stderr, "\n");

	sprintf (s, "[%d]", colIx);
	fprintf (stderr, "M%-5s %c: ", s, seqNuc);
	for (rowIx=0 ; rowIx<=qryLen ; rowIx++)
		{
		sprintf (s, "%d", dpColM[rowIx]);
		fprintf (stderr, " %-4s", s);
		}
	fprintf (stderr, "\n");

	sprintf (s, "[%d]", colIx);
	fprintf (stderr, "          ");
	for (rowIx=0 ; rowIx<=qryLen ; rowIx++)
		fprintf (stderr, " %-4s", debug_traceback_string(tbColM[rowIx]));
	fprintf (stderr, "\n");

	sprintf (s, "[%d]", colIx);
	fprintf (stderr, "I%-5s    ", s);
	for (rowIx=0 ; rowIx<=qryLen ; rowIx++)
		{
		sprintf (s, "%d", dpColI[rowIx]);
		fprintf (stderr, " %-4s", s);
		}
	fprintf (stderr, "\n");

	sprintf (s, "[%d]", colIx);
	fprintf (stderr, "          ");
	for (rowIx=0 ; rowIx<=qryLen ; rowIx++)
		fprintf (stderr, " %-4s", debug_traceback_string(tbColI[rowIx]));
	fprintf (stderr, "\n");

	sprintf (s, "[%d]", colIx);
	fprintf (stderr, "D%-5s    ", s);
	for (rowIx=0 ; rowIx<=qryLen ; rowIx++)
		{
		sprintf (s, "%d", dpColD[rowIx]);
		fprintf (stderr, " %-4s", s);
		}
	fprintf (stderr, "\n");

	sprintf (s, "[%d]", colIx);
	fprintf (stderr, "          ");
	for (rowIx=0 ; rowIx<=qryLen ; rowIx++)
		fprintf (stderr, " %-4s", debug_traceback_string(tbColD[rowIx]));
	fprintf (stderr, "\n");

	fprintf (stderr, "==============\n");
	}

static char* debug_traceback_string
   (tbcell	tbVal)
	{
	static char s[4];
	char	*ss;

	if ((tbVal & tbWhich) == 0)
		strcpy (s, "-");
	else
		{
		ss = s;

		if      ((tbVal & tbWhich) == tbWhichI)  *(ss++) = 'I';
		else if ((tbVal & tbWhich) == tbWhichD)  *(ss++) = 'D';
		else   /*(tbVal & tbWhich) == tbWhichM*/ *(ss++) = 'M';

		if      ((tbVal & tbRow)   == tbRowPrev) *(ss++) = 'S';
		else if ((tbVal & tbRow)   == tbRowLoop) *(ss++) = 'L';

		if      ((tbVal & tbCol)   == tbColPrev) *(ss++) = 'W';

		*ss = 0;
		}

	return s;
	}

static char debug_traceback_which
   (tbcell	tbVal)
	{
	if      ((tbVal & tbWhich) == tbWhichI)  return 'I';
	else if ((tbVal & tbWhich) == tbWhichD)  return 'D';
	else if ((tbVal & tbWhich) == tbWhichM)  return 'M';
	else                                     return '-';
	}

//----------
//
// rescore_alignment--
//	Recompute the score of an alignment.  
//	Count the number of matches, mismatches, insertions and deletions in an
//	alignment. 
//
//----------
//
// Arguments:
//	lascores*	scoring:	Alignment scoring definitions.
//	alignment*	a:			The alignment.
//
// Returns:
//	nothing;  the score, mCount, mmCount, iCount and dCount fields in the
//	alignment are modified.
//
//----------

void rescore_alignment
   (lascores*	scoring,
	alignment*	a)
	{
	u32			m, mm, iop, iex, dop, dex;
	char*		seqTextScan, *qryTextScan;
	char		prevEvent;

	m   = scoring->match;
	mm  = scoring->mismatch;
	iop = scoring->iOpen;
	iex = scoring->iExtend;
	dop = scoring->dOpen;
	dex = scoring->dExtend;

	seqTextScan = a->seqText;
	qryTextScan = a->qryText;

	a->mCount = a->mmCount = a->iCount = a->dCount = 0;
	a->score = 0;  prevEvent = 'M';
	for ( ; *seqTextScan!=0 ; seqTextScan++,qryTextScan++)
		{
		if (*qryTextScan == '-')
			{
			a->iCount++;
			if (prevEvent != 'I') a->score -= iop;
			                 else a->score -= iex;
			prevEvent = 'I';
			}
		else if (*seqTextScan == '-')
			{
			a->dCount++;
			if (prevEvent != 'D') a->score -= dop;
			                 else a->score -= dex;
			prevEvent = 'D';
			}
		else if (*qryTextScan == *seqTextScan)
			{
			a->mCount++;
			a->score += m;
			prevEvent = 'M';
			}
		else
			{
			a->mmCount++;
			a->score -= mm;
			prevEvent = 'M';
			}
		}

	}

//----------
//
// alignment_crc--
//	Compute the cyclic redundancy check value of an alignment.  
//
//----------
//
// Arguments:
//	char*		seqName:	Name of the aligned sequence.
//	char*		motif:		Aligned motif.
//	alignment*	a:			The alignment.
//
// Returns:
//	The alignment's crc value.
//
//----------

u32 alignment_crc
   (char*		seqName,
	char*		motif,
	alignment*	a)
	{
	u32			crc = 0;

	if (seqName != NULL) crc = update_crc_string(crc,seqName);
	crc = update_crc    (crc,(u8)a->strand);
	crc = update_crc_u32(crc,(int)a->seqStart);

	if (motif != NULL) crc = update_crc_string(crc,motif);

	crc = update_crc_string(crc,a->seqText);
	crc = update_crc_string(crc,a->qryText);

	return crc;
	}

//----------
//
// enough_cells--
//	Make sure we have enough cells allocated to support a given alignment
//	problem.
//
//----------
//
// Arguments:
//	lacontrol*	control:		The aligner control block.
//	u32			queryLen:		The number of rows in the matrices.
//	u32			sequenceLen:	The number of columns in the matrices.
//
// Returns:
//	(nothing)
//
//----------

static void enough_cells
   (lacontrol*	control,
	u32			queryLen,
	u32			sequenceLen)
	{
	size_t		numCells, growthCells;
	u32			numChars, growthChars;
	size_t		bytesNeeded;
	void*		dp;
	char*		seqText;

	// allocate a block of memory for DP matrices and traceback matrices

	numCells = (queryLen+1) * (size_t) (sequenceLen+1);
	if (numCells > control->allocCells)
		{
		growthCells = ((size_t) control->allocCells) + 100 + (control->allocCells/10);
		if (numCells < growthCells)
			numCells = growthCells;

		bytesNeeded = 3 * numCells * sizeof(dpcell)
		            + 3 * numCells * sizeof(tbcell);

		if (dbgLAAllocation)
			{
			fprintf (stderr, "for enough_cells(%s , %s)\n",ucommatize(queryLen),ucommatize(sequenceLen));
			fprintf (stderr, "  numCells=%s 0x%016lX\n",ucommatize(numCells),numCells);
			fprintf (stderr, "  bytesNeeded=%s %016lX\n",ucommatize(bytesNeeded),bytesNeeded);
			}

		if (control->dp == NULL)
			{
			dp = malloc (bytesNeeded);
			if (dp == NULL)
				{
				fprintf (stderr, "failed to allocate %s bytes for %s DP cells\n",
								 ucommatize(bytesNeeded), ucommatize(numCells));
				exit (EXIT_FAILURE);
				}
			}
		else
			{
			dp = realloc (control->dp, bytesNeeded);
			if (dp == NULL)
				{
				fprintf (stderr, "failed to re-allocate %s bytes for %s DP cells\n",
								 ucommatize(bytesNeeded), ucommatize(numCells));
				exit (EXIT_FAILURE);
				}
			}

		control->allocCells = numCells;
		control->dp = dp;
		control->tb = ((char*) dp) + (3 * numCells * sizeof(dpcell));

		if (dbgLAAllocation)
			{
			fprintf (stderr, "allocated %s bytes for %s DP cells\n",
							 ucommatize(bytesNeeded), ucommatize(numCells));
			fprintf (stderr, "  dp=%p\n", control->dp);
			fprintf (stderr, "  tb=%p\n", control->tb);
			}
		}

	// allocate a block of memory for alignment text for both sequence and
	// query;  the maximum length is twice times the sequence length, which
	// is enuough to allow an insertion and deletion at every position

	numChars = 2 * (sequenceLen+2);
	if (numChars > control->textChars)
		{
		growthChars = control->textChars + 100 + (control->textChars/10);
		if (numChars < growthChars)
			numChars = growthChars;

		bytesNeeded = 2 * numChars;
		if (control->seqText == NULL)
			{
			seqText = (char*) malloc (bytesNeeded);
			if (seqText == NULL)
				{
				fprintf (stderr, "failed to allocate %s bytes for alignment text\n",
								 ucommatize(bytesNeeded));
				exit (EXIT_FAILURE);
				}
			}
		else
			{
			seqText = (char*) realloc (control->seqText, bytesNeeded);
			if (seqText == NULL)
				{
				fprintf (stderr, "failed to allocate %s bytes for alignment text\n",
								 ucommatize(bytesNeeded));
				exit (EXIT_FAILURE);
				}
			}

		control->textChars = numChars;
		control->seqText   = seqText;
		control->qryText   = seqText + numChars;

		if (dbgLAAllocation)
			{
			fprintf (stderr, "allocated %s bytes for alignment text\n",
							 ucommatize(bytesNeeded));
			fprintf (stderr, "  seqText=%p\n", control->seqText);
			fprintf (stderr, "  qryText=%p\n", control->qryText);
			}
		}

	}

