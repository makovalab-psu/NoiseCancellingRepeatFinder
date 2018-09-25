// noise_cancelling_repeat_finder.c-- DNA aligner to align sequences to simple
// loops -- e.g. a microsatellite with a known core motif.

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "loop_aligner.h"
#include "feed.h"
#include "motifs.h"
#include "seq_ops.h"
#include "clock.h"
#include "utilities.h"

char* programName = "NCRF";
char* programVersionMajor    = VERSION_MAJOR;
char* programVersionMinor    = VERSION_MINOR;
char* programVersionSubMinor = VERSION_SUBMINOR;
char* programRevisionDate    = REVISION_DATE;

//----------
//
// global data and types--
//
//----------

// command line options

motif*	motifList = NULL;

lascores defaultScoring    = { 2,  5,	// match, mismatch
                               4,  11,	// insert-open, -extend
                               5,  5,	// delete-open, -extend
                               0 };		// minimum score
lascores pacbioScoring     = { 10, 35,	// match, mismatch
                               33, 21,	// insert-open, -extend
                               6,  28,	// delete-open, -extend
                               0 };		// minimum score
lascores pacbioScoringV2   = { 10, 51,	// match, mismatch
                               47, 59,	// insert-open, -extend
                               43, 23,	// delete-open, -extend
                               0 };		// minimum score
lascores pacbioScoringV1   = { 6,  14,	// match, mismatch
                               7,  22,	// insert-open, -extend
                               12, 12,	// delete-open, -extend
                               0 };		// minimum score
lascores nanoporeScoring   = { 10, 63,	// match, mismatch
                               51, 98,	// insert-open, -extend
                               27, 34,	// delete-open, -extend
                               0 };		// minimum score
lascores nanoporeScoringV2 = { 10, 36,	// match, mismatch
                               33, 37,	// insert-open, -extend
                               37, 44,	// delete-open, -extend
                               0 };		// minimum score
lascores nanoporeScoringV1 = { 10, 74,	// match, mismatch
                               70, 99,	// insert-open, -extend
                               45, 30,	// delete-open, -extend
                               0 };		// minimum score
lascores scoring;
lacontrol control;
int     reportScoring     = false;

#define defaultMinLength 500
u32		minLength         = defaultMinLength;	// minimum repeat length for a
                                // .. qualifying  alignment
float	minMRatio         = 0.0;// minimum ratio of matches to columns (a.k.a
                                // .. 1-maxNoise

int		debridge          = true;
u32		debridgeLength    = 0;	// (valid only if debridge is true) minimum
								// repeat length for qualifying alignment
								// .. outside of bridges
float	debridgeMRatio    = 0.0;// (valid only if debridge is true) minimum
								// ratio of matches to columns outside of
								// .. bridges
float	debridgeMRatioDiff= 0.02;//(valid only if debridge is true)
								// .. effective difference between
								// .. debridgeMRatio and 1-minMRatio

u32		clumpLength       = 0;	// (valid only if debridge is true) minimum
								// .. repeat length for DISqualifying clump in
								// .. an alignment
float	clumpMRatio       = 0.0;// (valid only if debridge is true) minimum
								// .. ratio of errors to columns inside a clump
float	clumpMRatioDiff   = 0.10;//(valid only if debridge is true)
								// .. effective difference between clumpMRatio
								// .. and 1-debridgeMRatio

u32		alignerMem        = 0;

int		sequenceStrands   = 0;
u32		minSequenceLength = 0;
u32		sequenceSkip      = 0;
u32		sequenceLimit     = 0;
u32		sequenceProgress  = 0;
u32		sequenceMem       = 0;
u32		sequenceNameMem   = 0;
u32		debridgerMem      = 0;
u32		clumpDetectMem    = 0;

int		showForward       = false;
int		showEvents        = false;
int		showEventsByPos   = false;
int		showCrc           = false;
u32		nameFieldW        = 1;
u32		lengthFieldW      = 1;
u32		countFieldW       = 1;
u32		rangeFieldW       = 1;
int     markEndOfFile     = true;

int		dbgAllocation     = false;
int		dbgNoSequence     = false;
int		dbgReportMotifs   = false;
int		dbgReportSequences = false;
int     dbgReportScoring  = false;
int		dbgReportAlignments = false;
int		dbgSegments       = false;
int		dbgSubalignment   = false;
int		dbgInhibitSubalignment = false;
int		dbgDebridger      = false;
int		dbgDebridgerDetail = false;
int		dbgDebridgerStack = false;
int		dbgClumpDetection = false;
int		dbgClumpSegments  = false;
int		dbgClumpCharacters = false;
int		dbgFlankLength    = 0;
char*	dbgTriggerName    = NULL;
u32		dbgTriggerCrc     = 0;

// sequence segment lists

typedef struct seglist
	{
	struct seglist* next;		// next segment in linked list
	u32		start;				// start and end of segment, 0-based half
	u32		end;				// .. open
	} seglist;

static seglist* segmentPool = NULL;

// by-position event counters

typedef struct evcount
	{
	u32		m;					// number of matches
	u32		mm;					// number of mismatches
	u32		i;					// number of inserted bases
	u32		d;					// number of deleted bases
	u32		mmA, mmC, mmG, mmT;	// number of mismatches, by substituted base
	} evcount;

static evcount* eventCounts = NULL;

// debridging

#define minDebridgeMaxStack     (10+1)
#define defaultDebridgeMaxStack ((10*1000)+1)

typedef struct
	{
	u32		lower;				// stack index of most recent interval whose
								// .. left score is no higher than this one's
	u32		left, right;
	double	leftScore, rightScore;
	} interval;

u32       debridgeMaxStack = 0;
interval* debridgeStack = NULL;	// (note that entry zero is not used)

// error clump detection/removal

#define minNumClumpSegments     (1000+1)
#define defaultNumClumpSegments (5000+1)

u32       numClumpSegments = 0;
double*   clumpMinSums  = NULL;
u32*      clumpMinWhere = NULL;

#define clumpEntrySize (sizeof(double)+sizeof(u32)+sizeof(interval))

// sequence flank strings (for deugging)

char* prefixFlank = NULL;
char* suffixFlank = NULL;

// reporting

u32 alignmentsReported;

//----------
//
// prototypes--
//
//----------

int main (int argc, char** argv);

// private functions

static void      parse_options              (int _argc, char** _argv);
static void      report_motif_matches       (FILE* f, feed* seq, motif* m);
static alignment best_alignment             (feed* seq, u32 start, u32 end, motif* m,
                                             u32 minLength);
static alignment longest_subalignment       (lascores* scoring, alignment a,
                                             u32 minLength, float minMRatio);
static alignment subalignment               (lascores* scoring, alignment a,
                                             u32 lftIx, u32 rgtIx);
static void      report_debridged_alignment (FILE* f, feed* seq, motif* m, alignment a);
static seglist*  find_error_clumps          (feed* seq, alignment a);
static void      report_alignment           (FILE* f, feed* seq, motif* m, alignment a);
static void      report_alignment_for_debug (FILE* f, feed* seq, motif* m, alignment a);
static void      count_positional_events    (motif* m, alignment a);
static void      free_segment_pool          (void);
static seglist*  allocate_segment           (void);
static seglist*  add_segment                (seglist** segmentList, u32 start, u32 end);
static seglist*  pop_segment                (seglist** segmentList);
static void      unuse_segment              (seglist* s);

//----------
//
// main program--
//
//----------

int main
   (int		argc,
	char**	argv)
	{
	feed*	sequenceFeed = NULL;
	u32		numSequences;
	motif*	m;

	parse_options (argc, argv);

	if (dbgReportMotifs)
		dump_motifs (stderr, motifList);

	if (reportScoring)
		{
		printf ("M=%d MM=%d IO=%d IX=%d DO=%d DX=%d\n",
		        scoring.match,scoring.mismatch,
		        scoring.iOpen,scoring.iExtend,
		        scoring.dOpen,scoring.dExtend);
		exit (EXIT_SUCCESS);
		}

	// derive debridging and error clump detection settings from minMRatio

	debridgeMRatio = 0.0;

	if (debridge)
		{
		debridgeMRatio = minMRatio;
		minMRatio      = minMRatio - debridgeMRatioDiff;
		}

	if (clumpMRatio < 0.0)
		clumpMRatio = 0.0;
	else if (debridgeMRatio > 0.0)
		{
		clumpMRatio = 1 - (debridgeMRatio - clumpMRatioDiff);
		if (clumpLength == 0) clumpLength = 100;
		}

	//fprintf (stderr, "minMRatio      = %.2f\n", minMRatio);
	//fprintf (stderr, "debridgeMRatio = %.2f\n", debridgeMRatio);
	//fprintf (stderr, "clumpMRatio    = %.2f\n", clumpMRatio);
	//fprintf (stderr, "minLength      = %u\n",   minLength);
	//fprintf (stderr, "debridgeLength = %u\n",   debridgeLength);
	//fprintf (stderr, "clumpLength    = %u\n",   clumpLength);

	// allocate an array to count positional events (if needed)

	if (showEventsByPos)
		{
		motif* m;
		int    longestMotif = 0;

		for (m=motifList ; m!=NULL ; m=m->next)
			{
			int motifLen = strlen ((char*) m->forwardNucs);
			if (motifLen > longestMotif) longestMotif = motifLen;
			}

		if (longestMotif == 0) longestMotif = 1;

		eventCounts = (evcount*) malloc (longestMotif * sizeof(evcount));
		if (eventCounts == NULL)
			{
			fprintf (stderr, "failed to allocate positional events array (for %d positions)\n",
			                 longestMotif);
			exit (EXIT_FAILURE);
			}
		if (dbgAllocation)
			{
			fprintf (stderr, "  allocated %ld bytes for positional events array (for %d positions)\n",
			                 longestMotif * sizeof(evcount), longestMotif);
			}
		}

	// allocate for debridging

	if (debridgeMRatio > 0.0)
		{
		if (debridgerMem == 0)
			debridgeMaxStack = defaultDebridgeMaxStack;
		else
			{
			debridgeMaxStack = debridgerMem / sizeof(interval);
			if (debridgeMaxStack < minDebridgeMaxStack)
				debridgeMaxStack = minDebridgeMaxStack;
			}

		debridgeStack = (interval*) malloc (debridgeMaxStack * sizeof(interval));
		if (debridgeStack == NULL)
			{
			fprintf (stderr, "failed to allocate debridger stack (for %d intervals)\n",
			                 debridgeMaxStack-1);
			exit (EXIT_FAILURE);
			}
		if (dbgAllocation)
			{
			fprintf (stderr, "  allocated %ld bytes for debridger stack (for %d intervals)\n",
			                 debridgeMaxStack * sizeof(interval), debridgeMaxStack-1);

			}
		}

	// allocate for error clump detection

	if (clumpMRatio > 0.0)
		{
		if (clumpDetectMem == 0)
			numClumpSegments = defaultNumClumpSegments;
		else
			{
			numClumpSegments = clumpDetectMem / clumpEntrySize;
			if (numClumpSegments < minNumClumpSegments)
				numClumpSegments = minNumClumpSegments;
			}

		clumpMinSums = (double*) malloc ((numClumpSegments+1) * sizeof(double));
		if (clumpMinSums == NULL)
			{
			fprintf (stderr, "failed to allocate clump detection (for %d entries)\n",
			                 numClumpSegments);
			exit (EXIT_FAILURE);
			}
		if (dbgAllocation)
			{
			fprintf (stderr, "  allocated %ld bytes for clump detection sums (for %d entries)\n",
			                 (numClumpSegments+1) * sizeof(double), numClumpSegments);

			}

		clumpMinWhere = (u32*) malloc ((numClumpSegments+1) * sizeof(u32));
		if (clumpMinWhere == NULL)
			{
			fprintf (stderr, "failed to allocate clump detection (for %d entries)\n",
			                 numClumpSegments);
			exit (EXIT_FAILURE);
			}
		if (dbgAllocation)
			{
			fprintf (stderr, "  allocated %ld bytes for clump detection wheres (for %d entries)\n",
			                 (numClumpSegments+1) * sizeof(u32), numClumpSegments);

			}
		}

	// allocate sequence flank strings

	if (dbgFlankLength > 0)
		{
		prefixFlank = (char*) malloc (dbgFlankLength+1);
		if (prefixFlank != NULL) suffixFlank = (char*) malloc (dbgFlankLength+1);
		if ((prefixFlank == NULL) || (suffixFlank == NULL))
			{
			fprintf (stderr, "failed to allocate %d bytes for flanking strings\n",
			                 2 * (dbgFlankLength+1));
			exit (EXIT_FAILURE);
			}
		if (dbgAllocation)
			{
			fprintf (stderr, "  allocated %d bytes for flanking strings)\n",
			                 2 * (dbgFlankLength+1));

			}
		}

	// process queries

	if (dbgNoSequence) goto no_sequences;

	init_loop_align (&control, &scoring, alignerMem);
	sequenceFeed = open_fasta_feed (NULL, /*forward only*/ 1, sequenceMem, sequenceNameMem);

	reset_progress_clock ();
	alignmentsReported = 0;
	numSequences = 0;
	while (another_sequence_in_feed (sequenceFeed))
		{
		numSequences++;
		if ((sequenceLimit != 0) && (numSequences > sequenceLimit))
			{
			fprintf (stderr, "limit of %s sequences reached\n", ucommatize(sequenceLimit));
			break;
			}
		if (sequenceSkip != 0)
			{
			if (numSequences < sequenceSkip)
				continue;
			if (numSequences == sequenceSkip)
				fprintf (stderr, "skipped %s sequences\n", ucommatize(sequenceSkip));
			}
		if (sequenceProgress != 0)
			{
			if ((numSequences == 1) || (numSequences % sequenceProgress == 0))
				{
				report_progress_clock (stderr, "({time})");
				fprintf (stderr, " processing sequence %s: %s%c (%u alignments reported so far)\n",
				                 ucommatize(numSequences),
				                 sequenceFeed->name, sequenceFeed->isRevComp?'-':'+',
				                 alignmentsReported);
				reset_progress_clock ();
				}
			}

		if (dbgTriggerName != NULL)
			{
			int triggerLen = strlen(dbgTriggerName);
			if (strncmp(sequenceFeed->name,dbgTriggerName,triggerLen) == 0)
				{
				free(dbgTriggerName);		// we've hit the trigger, so
				dbgTriggerName = NULL;		// .. disable it
				fprintf(stderr,"triggering on sequence \"%s\"\n",sequenceFeed->name);
				}
			}

		if (minSequenceLength > 0)
			{
			u32 seqLen = strlen ((char*) sequenceFeed->nt);
			if (seqLen < minSequenceLength)
				continue;
			}

		if (dbgReportSequences)
			dump_feed_sequence (stdout, sequenceFeed, 100);

		for (m=motifList ; m!=NULL ; m=m->next)
			report_motif_matches (stdout, sequenceFeed, m);
		}

	if (markEndOfFile)
		printf ("# ncrf end-of-file\n");

	fprintf (stderr, "(%u alignments reported)\n", alignmentsReported);

no_sequences:

	// relinquish allocated memory

	if (dbgTriggerName != NULL) { free (dbgTriggerName);           dbgTriggerName = NULL; }
	if (sequenceFeed   != NULL) { free_fasta_feed (sequenceFeed);  sequenceFeed   = NULL; }
	if (motifList      != NULL) { free_motifs (motifList);         motifList      = NULL; }
	if (eventCounts    != NULL) { free (eventCounts);              eventCounts    = NULL; }
	if (debridgeStack  != NULL) { free (debridgeStack);            debridgeStack  = NULL; }
	if (prefixFlank    != NULL) { free (prefixFlank);              prefixFlank    = NULL; }
	if (suffixFlank    != NULL) { free (suffixFlank);              suffixFlank    = NULL; }
	free_loop_align (&control);
	free_segment_pool ();

	// proclaim success

	return EXIT_SUCCESS;
	}

//----------
//
// option parsing--
//
//----------

static void  chastise         (const char* format, ...);
static void  usage            (char* message);
static void  usage_scoring    (void);
static void  usage_allocation (void);
static void  usage_other      (void);


static void chastise (const char* format, ...)
	{
	va_list	args;

	va_start (args, format);
	if (format != NULL)
		vfprintf (stderr, format, args);
	va_end (args);

	usage (NULL);
	}

static void usage (char* message)
	{
	if (message != NULL) fprintf (stderr, "%s\n", message);

	fprintf (stderr, "%s-- Noise Cancelling Repeat Finder, to find tandem repeats in noisy reads\n",
	                  programName);
	fprintf (stderr, "  (version %s.%s.%s %s)\n",
	                  programVersionMajor, programVersionMinor, programVersionSubMinor, programRevisionDate);
	fprintf (stderr, "usage: cat <fasta> | %s [options]\n", programName);
	fprintf (stderr, "\n");
	//                123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (stderr, "  <fasta>               fasta file containing sequences; read from stdin\n");
	fprintf (stderr, "  [<name>:]<motif>      dna repeat motif to search for\n");
	fprintf (stderr, "                        (there can be more than one motif)\n");
	fprintf (stderr, "  --minmratio=<ratio>   discard alignments with a low frequency of matches;\n");
	fprintf (stderr, "                        ratio can be between 0 and 1 (e.g. \"0.85\"), or can be\n");
	fprintf (stderr, "                        expressed as a percentage (e.g. \"85%%\")\n");
	fprintf (stderr, "  --maxnoise=<ratio>    (same as --minmratio but with 1-ratio)\n");
	fprintf (stderr, "  --minlength=<bp>      discard alignments that don't have long enough repeat\n");
	fprintf (stderr, "                        (default is %u)\n",defaultMinLength);
	fprintf (stderr, "  --minscore=<score>    discard alignments that don't score high enough\n");
	fprintf (stderr, "                        (default is zero)\n");
	fprintf (stderr, "  --stats=events        show match/mismatch/insert/delete counts\n");
	fprintf (stderr, "  --positionalevents    show match/mismatch/insert/delete counts by motif\n");
	fprintf (stderr, "                        position (independent of --stats=events)\n");
	fprintf (stderr, "                        this may be useful for detecting positional bias\n");
	fprintf (stderr, "  --help=scoring        show options relating to alignment scoring\n");
	fprintf (stderr, "  --help=allocation     show options relating to memory allocation\n");
	fprintf (stderr, "  --help=other          show other, less frequently used options\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "  The output is usually passed through a series of the ncrf* post-processing\n");
	fprintf (stderr, "  scripts.\n");

	exit (EXIT_FAILURE);
	}


static void usage_scoring (void)
	{
	//                123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (stderr, "  Default scoring is    M=%d MM=%d IO=%d IX=%d DO=%d DX=%d\n",
	                                          defaultScoring.match,defaultScoring.mismatch,
	                                          defaultScoring.iOpen,defaultScoring.iExtend,
	                                          defaultScoring.dOpen,defaultScoring.dExtend);
	fprintf (stderr, "  --scoring=pacbio      use scoring for pacbio reads\n");
	fprintf (stderr, "                        (M=%d MM=%d IO=%d IX=%d DO=%d DX=%d)\n",
	                                          pacbioScoring.match,pacbioScoring.mismatch,
	                                          pacbioScoring.iOpen,pacbioScoring.iExtend,
	                                          pacbioScoring.dOpen,pacbioScoring.dExtend);
	fprintf (stderr, "  --scoring=nanopore    use scoring for nanopore reads\n");
	fprintf (stderr, "                        (M=%d MM=%d IO=%d IX=%d DO=%d DX=%d)\n",
	                                          nanoporeScoring.match,nanoporeScoring.mismatch,
	                                          nanoporeScoring.iOpen,nanoporeScoring.iExtend,
	                                          nanoporeScoring.dOpen,nanoporeScoring.dExtend);
	fprintf (stderr, "  --match=<reward>      (M)  reward when sequence matches motif\n");
	fprintf (stderr, "  --mismatch=<penalty>  (MM) penalty when sequence mismatches motif\n");
	fprintf (stderr, "  --iopen=<penalty>     (IO) first penalty when sequence has nt, motif doesn't\n");
	fprintf (stderr, "  --iextend=<penalty>   (IX) other penalty when sequence has nt, motif doesn't\n");
	fprintf (stderr, "  --dopen=<penalty>     (DO) first penalty when motif has nt, sequence doesn't\n");
	fprintf (stderr, "  --dextend=<penalty>   (DX) other penalty when motif has nt, sequence doesn't\n");
	fprintf (stderr, "  --report:scoring      report scoring parameters and quit\n");
	fprintf (stderr, "\n");
	//                123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (stderr, "  Mathematically, scaling all scores, including --minscore, by a constant\n");
	fprintf (stderr, "  factor will produce equivalent alignment results. However, larger M shortens\n");
	fprintf (stderr, "  the length susceptible to overflow. Overflow may occur if an alignment is\n");
	fprintf (stderr, "  longer than 2^32/M. For M=100 this is about 40 million. An M larger than 100\n");
	fprintf (stderr, "  is unlikely to be necessary.\n");

	exit (EXIT_FAILURE);
	}


static void usage_allocation (void)
	{
	//                123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (stderr, "  --allocate:aligner=<bytes>       suggest space for aligner data structures\n");
	fprintf (stderr, "  --allocate:sequence=<bytes>      suggest space for sequence nucleotides\n");
	fprintf (stderr, "  --allocate:sequencename=<bytes>  suggest space for sequence name\n");
	fprintf (stderr, "  --allocate:debridger=<bytes>     suggest space for debridger data structures\n");
	fprintf (stderr, "                                   (each stack entry requires %lu bytes)\n", sizeof(interval));
	fprintf (stderr, "  --allocate:clump=<bytes>         suggest space for error clump data\n");
	fprintf (stderr, "                                   structures\n");

	exit (EXIT_FAILURE);
	}


static void usage_other (void)
	{
	//                123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (stderr, "  --strand=forward      limit processing to forward strand of sequence\n");
	fprintf (stderr, "  --strand=reverse      limit processing to reverse strand of sequence\n");
	fprintf (stderr, "  --showforward         always show the match in motif-forward orientation\n");
	fprintf (stderr, "                        (default: reverse-complement matches shown in reverse)\n");
	fprintf (stderr, "  --head=<number>       limit the number of input sequences\n");
	fprintf (stderr, "  --progress[=<n>]      report processing of every nth sequence\n");
	fprintf (stderr, "  --version             report the program version and quit\n");

	exit (EXIT_FAILURE);
	}


static void parse_options (int _argc, char** _argv)
	{
	int			argc;
	char**		argv;
	char*		arg, *argVal, *argVal2, *argVal3, *argVal4;

	// skip program name

	//programName = _argv[0];
	argv = _argv+1;  argc = _argc - 1;

	//////////
	// set defaults
	//////////

	scoring = defaultScoring;

	//////////
	// scan arguments
	//////////

	if (argc == 0)
		chastise (NULL);

	while (argc > 0)
		{
		arg = argv[0];
		if (arg[0] == 0) goto next_arg;
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --minmratio=<ratio>, --maxnoise=<ratio>

		if ((strcmp_prefix (arg, "--minmratio=") == 0)
		 || (strcmp_prefix (arg, "--minratio=")  == 0))
			{
			minMRatio = string_to_double (argVal);
			if ((minMRatio < 0.0) || (minMRatio > 1.0))
				chastise ("mratio has to be between 0 and 1, e.g. \"0.85\" or \"85%%\" (%s)\n", arg);
			goto next_arg;
			}

		if (strcmp_prefix (arg, "--maxnoise=") == 0)
			{
			minMRatio = 1.0 - string_to_double (argVal);
			if ((minMRatio < 0.0) || (minMRatio > 1.0))
				chastise ("noise has to be between 0 and 1, e.g. \"0.15\" or \"15%%\" (%s)\n", arg);
			goto next_arg;
			}

		// --minlength=<bp>

		if (strcmp_prefix (arg, "--minlength=") == 0)
			{
			minLength = string_to_u32 (argVal);
			if (minLength == 0) minLength = 1;
			goto next_arg;
			}

		// --minscore=<score>

		if (strcmp_prefix (arg, "--minscore=") == 0)
			{
			scoring.minScore = string_to_u32 (argVal);
			if ((s32) scoring.minScore < 0)
				chastise ("unable to represent minscore in 31 bits (%s)"
				          "\ntry something below 2 billion\n", arg);
			goto next_arg;
			}

		// --stats=events

		if (strcmp (arg, "--stats=events") == 0)
			{ showEvents = true;  goto next_arg; }

		// --positionalevents

		if ((strcmp (arg, "--positionalevents") == 0)
		 || (strcmp (arg, "--stats=positionalevents") == 0))
			{ showEventsByPos = true;  goto next_arg; }

		// --crc (unadvertised)

		if (strcmp (arg, "--crc") == 0)
			{ showCrc = true;  goto next_arg; }

		// --strand=<whatever>

		if ((strcmp (arg, "--strand=forward") == 0)
		 || (strcmp (arg, "--strand=plus"   ) == 0)
		 || (strcmp (arg, "--strand=+"      ) == 0)
		 || (strcmp (arg, "--forwardonly"   ) == 0))
			{ sequenceStrands = 1;  goto next_arg; }

		if ((strcmp (arg, "--strand=reverse"          ) == 0)
		 || (strcmp (arg, "--strand=backward"         ) == 0)
		 || (strcmp (arg, "--strand=complement"       ) == 0)
		 || (strcmp (arg, "--strand=minus"            ) == 0)
		 || (strcmp (arg, "--strand=revcomp"          ) == 0)
		 || (strcmp (arg, "--strand=reversecomplement") == 0)
		 || (strcmp (arg, "--strand=-"                ) == 0)
		 || (strcmp (arg, "--reverseonly"             ) == 0))
			{ sequenceStrands = -1;  goto next_arg; }

		// --showforward

		if (strcmp (arg, "--showforward") == 0)
			{ showForward = true;  goto next_arg; }

		// --scoring=pacbio
		//    or --scoring=pacbio.v3 (unadvertised)

		if ((strcmp (arg, "--scoring=pacbio")   == 0)
		 || (strcmp (arg, "--scoring=pacbio.v3") == 0))
			{
			u32	saveMinScore = scoring.minScore;
			scoring = pacbioScoring;
			scoring.minScore = saveMinScore;
			goto next_arg;
			}

		// --scoring=pacbio.v2 (unadvertised)

		if (strcmp (arg, "--scoring=pacbio.v2") == 0)
			{
			u32	saveMinScore = scoring.minScore;
			scoring = pacbioScoringV2;
			scoring.minScore = saveMinScore;
			goto next_arg;
			}

		// --scoring=pacbio.v1 (unadvertised)

		if (strcmp (arg, "--scoring=pacbio.v1") == 0)
			{
			u32	saveMinScore = scoring.minScore;
			scoring = pacbioScoringV1;
			scoring.minScore = saveMinScore;
			goto next_arg;
			}

		// --scoring=nanopore
		//    or --scoring=nanopore.v3 (unadvertised)

		if ((strcmp (arg, "--scoring=nanopore")   == 0)
		 || (strcmp (arg, "--scoring=nanopore.v3") == 0))
			{
			u32	saveMinScore = scoring.minScore;
			scoring = nanoporeScoring;
			scoring.minScore = saveMinScore;
			goto next_arg;
			}

		// --scoring=nanopore.v2 (unadvertised)

		if (strcmp (arg, "--scoring=nanopore.v2") == 0)
			{
			u32	saveMinScore = scoring.minScore;
			scoring = nanoporeScoringV2;
			scoring.minScore = saveMinScore;
			goto next_arg;
			}

		// --scoring=nanopore.v1 (unadvertised)

		if (strcmp (arg, "--scoring=nanopore.v1") == 0)
			{
			u32	saveMinScore = scoring.minScore;
			scoring = nanoporeScoringV1;
			scoring.minScore = saveMinScore;
			goto next_arg;
			}

		// --match=<reward>

		if ((strcmp_prefix (arg, "--match=") == 0)
		 || (strcmp_prefix (arg, "M=")       == 0))
			{
			scoring.match = string_to_u32 (argVal);
			if (scoring.match < 1)
				chastise ("invalid match reward: \"%s\"\n", argVal);
			goto next_arg;
			}

		// --mismatch=<penalty>

		if ((strcmp_prefix (arg, "--mismatch=") == 0)
		 || (strcmp_prefix (arg, "MM=")         == 0))
			{
			scoring.mismatch = string_to_u32 (argVal);
			if (scoring.mismatch < 1)
				chastise ("invalid mismatch penalty: \"%s\"\n", argVal);
			goto next_arg;
			}

		// --iopen=<penalty>

		if ((strcmp_prefix (arg, "--iopen=") == 0)
		 || (strcmp_prefix (arg, "--iOpen=") == 0)
		 || (strcmp_prefix (arg, "IO=")      == 0))
			{
			scoring.iOpen = string_to_u32 (argVal);
			if (scoring.iOpen < 1)
				chastise ("invalid insert-open penalty: \"%s\"\n", argVal);
			goto next_arg;
			}

		// --iextend=<penalty>

		if ((strcmp_prefix (arg, "--iextend=") == 0)
		 || (strcmp_prefix (arg, "--iExtend=") == 0)
		 || (strcmp_prefix (arg, "IX=")        == 0))
			{
			scoring.iExtend = string_to_u32 (argVal);
			if (scoring.iExtend < 1)
				chastise ("invalid insert-extend penalty: \"%s\"\n", argVal);
			goto next_arg;
			}

		// (undocumented) I=<penalty>[,<penalty>]

		if (strcmp_prefix (arg, "I=") == 0)
			{
			argVal2 = strchr(argVal,',');
			if (argVal2 != NULL) *(argVal2++) = 0;

			if (argVal2 == NULL)
				{
				scoring.iOpen = string_to_u32 (argVal);
				if (scoring.iOpen < 1)
					chastise ("invalid insert-open penalty: \"%s\"\n", argVal);
				scoring.iExtend = scoring.iOpen;
				}
			else
				{
				scoring.iOpen = string_to_u32 (argVal);
				if (scoring.iOpen < 1)
					chastise ("invalid insert-open penalty: \"%s\"\n", argVal);
				scoring.iExtend = string_to_u32 (argVal2);
				if (scoring.iExtend < 1)
					chastise ("invalid insert-extend penalty: \"%s\"\n", argVal2);
				}
			goto next_arg;
			}

		// --dopen=<penalty>

		if ((strcmp_prefix (arg, "--dopen=") == 0)
		 || (strcmp_prefix (arg, "--dOpen=") == 0)
		 || (strcmp_prefix (arg, "DO=")      == 0))
			{
			scoring.dOpen = string_to_u32 (argVal);
			if (scoring.dOpen < 1)
				chastise ("invalid delete-open penalty: \"%s\"\n", argVal);
			goto next_arg;
			}

		// --dextend=<penalty>

		if ((strcmp_prefix (arg, "--dextend=") == 0)
		 || (strcmp_prefix (arg, "--dExtend=") == 0)
		 || (strcmp_prefix (arg, "DX=")        == 0))
			{
			scoring.dExtend = string_to_u32 (argVal);
			if (scoring.dExtend < 1)
				chastise ("invalid delete-extend penalty: \"%s\"\n", argVal);
			goto next_arg;
			}

		if ((strcmp (arg, "--report:scoring") == 0)
		 || (strcmp (arg, "--report:scores")  == 0))
			{ reportScoring = true;  goto next_arg; }

		// (undocumented) D=<penalty>[,<penalty>]

		if (strcmp_prefix (arg, "D=") == 0)
			{
			argVal2 = strchr(argVal,',');
			if (argVal2 != NULL) *(argVal2++) = 0;

			if (argVal2 == NULL)
				{
				scoring.dOpen = string_to_u32 (argVal);
				if (scoring.dOpen < 1)
					chastise ("invalid delete-open penalty: \"%s\"\n", argVal);
				scoring.dExtend = scoring.dOpen;
				}
			else
				{
				scoring.dOpen = string_to_u32 (argVal);
				if (scoring.dOpen < 1)
					chastise ("invalid delete-open penalty: \"%s\"\n", argVal);
				scoring.dExtend = string_to_u32 (argVal2);
				if (scoring.dExtend < 1)
					chastise ("invalid delete-extend penalty: \"%s\"\n", argVal2);
				}
			goto next_arg;
			}

		// --debridgediff=<ratiodiff>, --debridgelength=<bp> and --nodebridge (unadvertized)

		if (strcmp_prefix (arg, "--debridgediff=") == 0)
			{
			debridgeMRatioDiff = string_to_double (argVal);
			if ((debridgeMRatio <= 0.0) || (debridgeMRatio >= 1.0))
				chastise ("bridgediff noise has to be strictly between 0 and 1 (e.g. 0.02 or 2%)\n", arg);
			goto next_arg;
			}

		if ((strcmp_prefix (arg, "--debridgelength=")    == 0)
		 || (strcmp_prefix (arg, "--debridgeminlength=") == 0))
			{ debridgeLength = string_to_u32 (argVal);  goto next_arg; }

		if ((strcmp (arg, "--nodebridge") == 0)
		 || (strcmp (arg, "--debridge=off") == 0)
		 || (strcmp (arg, "--debridge:off") == 0))
			{
			// inhibit debridging
			debridge = false;
			goto next_arg;
			}

		// --clumpdiff=<ratiodiff>, --minclump=<bp> and --nodeclump (unadvertized)

		if (strcmp_prefix (arg, "--clumpdiff=") == 0)
			{
			clumpMRatioDiff = string_to_double (argVal);
			if ((clumpMRatioDiff <= 0.0) || (clumpMRatioDiff >= 1.0))
				chastise ("clumpdiff noise has to be strictly between 0 and 1 (e.g. 0.10 or 10%)\n", arg);
			goto next_arg;
			}

		if (strcmp_prefix (arg, "--minclump=") == 0)
			{ clumpLength = string_to_u32 (argVal);  goto next_arg; }

		if ((strcmp (arg, "--nodeclump") == 0)
		 || (strcmp (arg, "--declump=off") == 0)
		 || (strcmp (arg, "--declump:off") == 0))
			{
			// inhibit error-clump removal
			clumpMRatio = -1.0;
			goto next_arg;
			}

		// field widths (unadvertised)

		if ((strcmp_prefix (arg, "--fields=") == 0)
		 || (strcmp_prefix (arg, "F=")        == 0))
			{
			argVal2 = strchr(argVal,',');
			if (argVal2 == NULL) chastise ("Can't understand \"%s\"\n", arg);
			argVal3 = strchr(argVal2+1,',');
			if (argVal3 == NULL) chastise ("Can't understand \"%s\"\n", arg);
			argVal4 = strchr(argVal3+1,',');
			if (argVal4 == NULL) chastise ("Can't understand \"%s\"\n", arg);
			*(argVal2++) = *(argVal3++) = *(argVal4++) = 0;
			nameFieldW = string_to_u32 (argVal);
			if (nameFieldW < 1) nameFieldW = 1;
			lengthFieldW = string_to_u32 (argVal2);
			if (lengthFieldW < 1) lengthFieldW = 1;
			countFieldW = string_to_u32 (argVal3);
			if (countFieldW < 1) countFieldW = 1;
			rangeFieldW = string_to_u32 (argVal4);
			if (rangeFieldW < 1) rangeFieldW = 1;
			goto next_arg;
			}

		if ((strcmp_prefix (arg, "--namefield=") == 0)
		 || (strcmp_prefix (arg, "F1=")          == 0))
			{
			nameFieldW = string_to_u32 (argVal);
			if (nameFieldW < 1) nameFieldW = 1;
			goto next_arg;
			}

		if ((strcmp_prefix (arg, "--lengthfield=") == 0)
		 || (strcmp_prefix (arg, "F2=")          == 0))
			{
			lengthFieldW = string_to_u32 (argVal);
			if (lengthFieldW < 1) lengthFieldW = 1;
			goto next_arg;
			}

		if ((strcmp_prefix (arg, "--countfield=") == 0)
		 || (strcmp_prefix (arg, "F3=")          == 0))
			{
			countFieldW = string_to_u32 (argVal);
			if (countFieldW < 1) countFieldW = 1;
			goto next_arg;
			}

		if ((strcmp_prefix (arg, "--intervalfield=") == 0)
		 || (strcmp_prefix (arg, "F4=")          == 0))
			{
			rangeFieldW = string_to_u32 (argVal);
			if (rangeFieldW < 1) rangeFieldW = 1;
			goto next_arg;
			}

		// --allocate:sequence=<bytes>

		if (strcmp_prefix (arg, "--allocate:aligner=") == 0)
			{
			alignerMem = string_to_unitized_u32 (argVal, false /*units of 1,024*/);
			goto next_arg;
			}

		// --allocate:sequence=<bytes>

		if (strcmp_prefix (arg, "--allocate:sequence=") == 0)
			{
			sequenceMem = string_to_unitized_u32 (argVal, false /*units of 1,024*/);
			goto next_arg;
			}

		// --allocate:sequencename=<bytes>

		if (strcmp_prefix (arg, "--allocate:sequencename=") == 0)
			{
			sequenceNameMem = string_to_unitized_u32 (argVal, false /*units of 1,024*/);
			goto next_arg;
			}

		// --allocate:debridger=<bytes>

		if (strcmp_prefix (arg, "--allocate:debridger=") == 0)
			{
			debridgerMem = string_to_unitized_u32 (argVal, false /*units of 1,024*/);
			goto next_arg;
			}

		// --allocate:clump=<bytes>

		if (strcmp_prefix (arg, "--allocate:clump=") == 0)
			{
			clumpDetectMem = string_to_unitized_u32 (argVal, false /*units of 1,024*/);
			goto next_arg;
			}

		// --minseqlength=<bp> (unadvertised)

		if ((strcmp_prefix (arg, "--minseqlength=") == 0)
		 || (strcmp_prefix (arg, "--seqminlength=") == 0))
			{ minSequenceLength = string_to_u32 (argVal);  goto next_arg; }

		// --noendmark, etc. (unadvertised)

		if ((strcmp (arg, "--noendmark") == 0)
		 || (strcmp (arg, "--noeof")     == 0)
		 || (strcmp (arg, "--nomark")    == 0))
			{ markEndOfFile = false;  goto next_arg; }

		// --skip=<number> (unadvertised)

		if ((strcmp_prefix (arg, "--skip=") == 0))
			{
			sequenceSkip = string_to_unitized_u32 (argVal, true /*units of 1,000*/);
			goto next_arg;
			}

		// --head=<number>

		if ((strcmp_prefix (arg, "--head=") == 0))
			{
			sequenceLimit = string_to_unitized_u32 (argVal, true /*units of 1,000*/);
			goto next_arg;
			}

		// --progress[=<n>]

		if (strcmp (arg, "--progress") == 0)
			{ sequenceProgress = 1;  goto next_arg; }

		if (strcmp_prefix (arg, "--progress=") == 0)
			{
			sequenceProgress = string_to_unitized_u32 (argVal, true /*units of 1,000*/);
			goto next_arg;
			}

		// --version

		if ((strcmp (arg, "--version") == 0)
		 || (strcmp (arg, "-v")        == 0)
		 || (strcmp (arg, "-version")  == 0))
			{
			fprintf (stderr, "%s (version %s.%s.%s %s)\n",
			                  programName,
			                  programVersionMajor, programVersionMinor, programVersionSubMinor, programRevisionDate);
#ifdef __GNUC__
			fprintf (stderr, "  built with gcc-%d.%d.%d \"%s\"\n",
			                  __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__, __VERSION__);
#endif
			exit (EXIT_FAILURE);
			}

		// --help

		if ((strcmp (arg, "--help") == 0)
		 || (strcmp (arg, "--h")    == 0)
		 || (strcmp (arg, "-help")  == 0)
		 || (strcmp (arg, "-h")     == 0))
			usage(NULL);

		if ((strcmp (arg, "--help=scoring") == 0)
		 || (strcmp (arg, "--help:scoring") == 0)
		 || (strcmp (arg, "--help=scores")  == 0)
		 || (strcmp (arg, "--help:scores")  == 0)
		 || (strcmp (arg, "--help=score")   == 0)
		 || (strcmp (arg, "--help:score")   == 0))
			usage_scoring();

		if ((strcmp (arg, "--help=allocation") == 0)
		 || (strcmp (arg, "--help:allocation") == 0)
		 || (strcmp (arg, "--help=alloc")      == 0)
		 || (strcmp (arg, "--help:alloc")      == 0)
		 || (strcmp (arg, "--help=memory")     == 0)
		 || (strcmp (arg, "--help:memory")     == 0))
			usage_allocation();

		if ((strcmp (arg, "--help=others")        == 0)
		 || (strcmp (arg, "--help:others")        == 0)
		 || (strcmp (arg, "--help=other")         == 0)
		 || (strcmp (arg, "--help:other")         == 0)
		 || (strcmp (arg, "--help=miscellaneous") == 0)
		 || (strcmp (arg, "--help:miscellaneous") == 0))
			usage_other();

		// --debug=<whatever>

		if ((strcmp (arg, "--debug=nosequence") == 0)
		 || (strcmp (arg, "--debug=notarget")   == 0))
			{ dbgNoSequence = true;  goto next_arg; }

		if (strcmp (arg, "--debug=report:motif") == 0)
			{ dbgReportMotifs = true;  goto next_arg; }

		if ((strcmp (arg, "--debug=report:sequence") == 0)
		 || (strcmp (arg, "--debug=report:target")   == 0))
			{ dbgReportSequences = true;  goto next_arg; }

		if (strcmp (arg, "--debug=report:scoring") == 0)
			{ dbgReportScoring = true;  goto next_arg; }

		if ((strcmp (arg, "--debug=report:alignments") == 0)
		 || (strcmp (arg, "--debug=report:alignment")  == 0)
		 || (strcmp (arg, "--debug=report:aligned")    == 0))
			{ dbgReportAlignments = true;  goto next_arg; }

		if (strcmp (arg, "--debug=segments") == 0)
			{ dbgSegments = true;  goto next_arg; }

		if (strcmp (arg, "--debug=subalignment") == 0)
			{ dbgSubalignment = true;  goto next_arg; }

		if (strcmp (arg, "--debug=inhibit:subalignment") == 0)
			{ dbgInhibitSubalignment = true;  goto next_arg; }

		if ((strcmp (arg, "--debug=debridge")  == 0)
		 || (strcmp (arg, "--debug=debridger") == 0))
			{ dbgDebridger = true;  goto next_arg; }

		if ((strcmp (arg, "--debug=debridge+")  == 0)
		 || (strcmp (arg, "--debug=debridger+") == 0))
			{ dbgDebridgerDetail = true;  goto next_arg; }

		if ((strcmp (arg, "--debug=debridge:stack")  == 0)
		 || (strcmp (arg, "--debug=debridger:stack") == 0))
			{ dbgDebridgerStack = true;  goto next_arg; }

		if ((strcmp (arg, "--debug=clump") == 0)
		 || (strcmp (arg, "--debug=clumps") == 0)
		 || (strcmp (arg, "--debug=clumpdetect") == 0))
			{ dbgClumpDetection = true;  goto next_arg; }

		if ((strcmp (arg, "--debug=clumpsegs") == 0)
		 || (strcmp (arg, "--debug=clumpsegments") == 0))
			{ dbgClumpSegments = true;  goto next_arg; }

		if ((strcmp (arg, "--debug=clumpchars") == 0)
		 || (strcmp (arg, "--debug=clumpcharacters") == 0))
			{ dbgClumpCharacters = true;  goto next_arg; }

		if (strcmp_prefix (arg, "--debug=flank:") == 0)
			{
			argVal2 = strchr(argVal,':') + 1;
			dbgFlankLength = string_to_u32 (argVal2);
			goto next_arg;
			}

		if (strcmp_prefix (arg, "--debug=triggername:") == 0)
			{
			argVal2 = strchr(argVal,':') + 1;
			if (dbgTriggerName != NULL) free(dbgTriggerName);
			dbgTriggerName = copy_string(argVal2);
			goto next_arg;
			}

		if (strcmp_prefix (arg, "--debug=triggercrc:") == 0)
			{
			argVal2 = strchr(argVal,':') + 1;
			dbgTriggerCrc = string_to_u32 (argVal2);
			goto next_arg;
			}

		if (strcmp (arg, "--debug=loop_align") == 0)
			{ dbgLoopAlign = true;  goto next_arg; }

		if (strcmp (arg, "--debug=allocation") == 0)
			{ dbgAllocation = dbgLAAllocation = true;  goto next_arg; }

		if (strcmp (arg, "--debug=cell") == 0)
			{ dbgLACell = true;  goto next_arg; }

		if (strcmp (arg, "--debug=column") == 0)
			{ dbgLAColumn = true;  goto next_arg; }

		if (strcmp (arg, "--debug=traceback") == 0)
			{ dbgLATraceback = true;  goto next_arg; }

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("Can't understand \"%s\"\n", arg);

		// [<name>:]<motif>

		motifList = add_motif (motifList, arg);

	next_arg:
		argv++;  argc--;
		continue;
		}

	//////////
	// sanity checks
	//////////

	// make sure we got at least one motif

	if ((motifList == NULL) && (!reportScoring))
		chastise ("you have to give me at least one motif to search for\n");

	// warn the user if her scoring prefers a back-to-back insert-delete to a
	// mismatch

	if (scoring.mismatch > scoring.iOpen + scoring.dOpen)
		{
		fprintf (stderr, "WARNING: your scoring prefers a back-to-back insert-delete to a mismatch,\n");
		fprintf (stderr, "         because MM (%u) is more than IO+DO (%u)\n",
		                 scoring.mismatch, scoring.iOpen+scoring.dOpen);
		}

	// copy debridging min length from regular min length, if necessary 

	if ((debridge) && (debridgeLength == 0))
		debridgeLength = minLength;

	return;
	}

//----------
//
// report_motif_matches--
//	Search for and report matches between a sequence and a motif loop.
//
// We'll report more than one alignment if they exist, but the reported
// alignments will not overlap in the sequence.
//
//----------
//
// Arguments:
//	FILE*	f:		The file to write to.
//	feed*	seq:	The sequence to search in.
//	motif*	m:		The motif to search for.
//
// Returns:
//	(nothing)
//
//----------
//
// Implementation notes
//
// (1) This implementation, as regards looking for more than one alignment, is
// computationaly wasteful.  The initial search finds the highest scoring
// alignment, after which the search is repeated in the leftover segments on
// left and right of that alignment.  The DP and traceback matrices for the
// left segment are still intact but we make no use of them.
//
// Our expectation is that a very small fraction of sequences will have any
// alignment, so the wasted computation will be negligible.
//
// (2) A better process would be to do a quick seeding search to identify
// segments likely to match the motif.
//
//----------

static void report_motif_matches
   (FILE*		f,
	feed*		seq,
	motif*		m)
	{
	seglist*	queue, *seg;
	alignment	a, aSub;
	u32			seqLen, minSegment;
	u32			segStart, segEnd, aStart, aEnd;

	seqLen = strlen ((char*) seq->nt);
	if (seqLen == 0) return;

	if ((dbgSegments) && (dbgTriggerName == NULL))
		{
		fprintf (stderr, "===\n");
		fprintf (stderr, "report_motif_matches(%s,%s)\n", seq->name, m->name);
		}

	if (minLength > 0) minSegment = minLength;
	              else minSegment = 1;

	queue = NULL;
	add_segment (&queue, 0, seqLen);

	while (queue != NULL)
		{
		// pop a segment from the queue and find the best alignment in it;  if
		// the segment contains a suitable alignment, report it, then dispose
		// of it;  note that we don't report alignments that are too noisy, but
		// we will still look for alignments in the subintervals
		// $$$ need to check whether any subinterval of a too-noisy alignment
		//     .. would be long enough and noise-free enough;  what we want is
		//     .. the longest subinterval that meets the noise threshold;  that
		//     .. can be found in one pass;  we need to reflect those new bounds
		//     .. in dividing the segment, because there may be more than one
		//     .. suitable sub-alignment

		seg = pop_segment (&queue);
		segStart = seg->start;
		segEnd   = seg->end;
		unuse_segment (seg);

		a = best_alignment (seq, segStart, segEnd, m, minLength);
		if ((dbgSegments) && (dbgTriggerName == NULL))
			{
			if (a.active)
				fprintf (stderr, "(aligned at %d,%d)\n", a.seqStart, a.seqEnd);
			else
				fprintf (stderr, "(no suitable alignment found)\n");
			}
		if (!a.active) continue;

		if ((minMRatio != 0.0) && (a.matchRatio < minMRatio) && (!dbgInhibitSubalignment))
			{
			if ((dbgSegments) && (dbgTriggerName == NULL))
				fprintf (stderr, "poorly aligned at %d,%d\n", a.seqStart, a.seqEnd);

			aSub = longest_subalignment (&scoring, a, minLength, minMRatio);
			if (aSub.active)
				{
				if ((dbgSegments) && (dbgTriggerName == NULL))
					fprintf (stderr, "sub-aligned at %d,%d\n", aSub.seqStart, aSub.seqEnd);
				free_alignment (a);
				a = aSub;
				}
			}

		aStart = a.seqStart;
		aEnd   = a.seqEnd;

		if ((minMRatio == 0.0) || (a.matchRatio >= minMRatio))
			{
			if (debridgeMRatio > 0.0)
				report_debridged_alignment (f, seq, m, a);
			else
				{
				if (a.score >= (s32) control.scoring.minScore)
					report_alignment(f, seq, m, a);
				}
			}

		free_alignment (a);

		if ((dbgSegments) && (dbgTriggerName == NULL))
			fprintf (stderr, "aligned at %d,%d\n", aStart, aEnd);

		// add left and right leftover segments back to the queue (unless they
		// are too short)

		if (segEnd - aEnd >= minSegment)            // right leftover
			add_segment (&queue, aEnd, segEnd);
		if (aStart - segStart >= minSegment)        // left leftover
			add_segment (&queue, segStart, aStart);
		}

	if ((dbgSegments) && (dbgTriggerName == NULL))
		fprintf (stderr, "(report_motif_matches done)\n");
	}

//----------
//
// best_alignment--
//	Find the best alignment between a segment of a sequence and a motif loop.
//
//----------
//
// Arguments:
//	feed*	seq:		The sequence to search in.
//	u32		start, end:	Start and end of segment, 0-based half open
//	motif*	m:			The motif to search for.
//	u32		minLength:	Minimum repeat length for qualifying alignment.
//
// Returns:
//	The best alignment;  noAlignment if no alignment is suitable.  The caller
//	is responsible for deallocating any memory allocated for this result, by
//	eventually calling free_alignment().
//
//----------

static alignment best_alignment
   (feed*		seq,
	u32			start,
	u32			end,
	motif*		m,
	u32			minLength)
	{
	alignment	aForward, aReverse, aBest;

	if ((dbgSegments) && (dbgTriggerName == NULL))
		fprintf (stderr, "best_alignment(%s[%u..%u],%s)\n",
		                 seq->name, start, end, m->name);

	if (dbgReportScoring)
		{
		fprintf (stderr, "M=%d MM=%d IO=%d IX=%d DO=%d DX=%d --minscore=%d\n",
		        control.scoring.match,control.scoring.mismatch,
		        control.scoring.iOpen,control.scoring.iExtend,
		        control.scoring.dOpen,control.scoring.dExtend,
		        control.scoring.minScore);
		dbgReportScoring = false;
		}

	// try to align to both strands (or whichever strand we're interested in)

	aForward = noAlignment;
	aReverse = noAlignment;
	if (sequenceStrands != -1)
		aForward = loop_align_segment (&control, seq->nt+start, end-start, m->forwardNucs);
	if ((sequenceStrands != 1) && (m->reverseNucs != NULL))
		aReverse = loop_align_segment (&control, seq->nt+start, end-start, m->reverseNucs);

	// discard too-short alignments

	if (minLength > 0)
		{
		if ((aForward.active) && (aForward.qryBaseCount < minLength))
			{ free_alignment (aForward);  aForward.active = false; }
		if ((aReverse.active) && (aReverse.qryBaseCount < minLength))
			{ free_alignment (aReverse);  aReverse.active = false; }
		}

	// if we have no suitable alignment, we're done

	if ((!aForward.active) && (!aReverse.active))
		return noAlignment;

	// decide which alignment is best

	if (aForward.active) aForward.strand = '+';
	if (aReverse.active) aReverse.strand = '-';

	if (!aForward.active)
		aBest = aReverse;
	else if (!aReverse.active)
		aBest = aForward;
	else if (aForward.score >= aReverse.score)
		{ aBest = aForward;  free_alignment (aReverse); }
	else
		{ aBest = aReverse;  free_alignment (aForward); }

	// adjust position to account for the subinterval

	aBest.seqStart += start;
	aBest.seqEnd   += start;

	return aBest;
	}

//----------
//
// longest_subalignment--
//	Find the longest subinterval of an alignment, such that it meets length
//	and match ratio criteria.
//
//----------
//
// Arguments:
//	lascores*	scoring:	Alignment scoring definitions.
//	alignment	a:			The alignment to search for a subinterval of.
//	u32			minLength:	Minimum repeat length for qualifying alignment.
//	float		minMRatio:	Minimum ratio of matches to columns.
//
// Returns:
//	The longest subalignment meeting the criteria;  noAlignment if no
//	subalignment is suitable.  The caller is responsible for deallocating any
//	memory allocated for this result, by eventually calling free_alignment().
//
//----------
//
// Implementation Notes
//
// (1) Derivation of test for an interval meeting the match ratio criterion.
//		matchRatio >= minMRatio
//		==> matches / (matches+errors) >= minMRatio
//		==> matches >= (matches+errors) * minMRatio
//		==> matches * (1-minMRatio) >= errors * minMRatio
//		==> matches * (1-minMRatio)/minMRatio >= errors
//
//----------

// prototypes--

static void print_current_window (char* seqText, char* qryText,
                                  u32 lftIx, u32 rgtIx,
                                  u32 qryBases, u32 matches, u32 errors,
                                  float minMRatio);

// longest_subalignment--

static alignment longest_subalignment
   (lascores*	scoring,
	alignment	a,
	u32			minLength,
	float		minMRatio)
	{
	float		testRatio;
	u32			lftIx, rgtIx, goodRgtIx, bestLftIx, bestRgtIx;
	u32			windowBases, matches, errors;

	// handle trivial cases

	if (minLength == 0)
		return noAlignment;

	if (a.qryBaseCount < minLength)
		return noAlignment;

	if (a.matchRatio >= minMRatio)
		return copy_alignment (a);

	// scan the first L=minLength query (motif) bases of the alignment,
	// counting events;  this gives us the right end of the first window of
	// length L (the left end is at zero);  for all of what follows, the
	// window is defined by lftIx,rgtIx (indexes into qryText) such that
	//   lftIx < rgtIx (strictly less than)
	//   qryText[lftIx]   is ACGTacgt
	//   qryText[rgtIx-1] is ACGTacgt
	//   qryText[rgtIx]   is ACGTacgt, '-', or the zero terminator

	if (dbgSubalignment)
		fprintf (stderr, "longest_subalignment([%u..%u])\n", a.seqStart, a.seqEnd);

	windowBases = matches = errors = 0;
	for (rgtIx=0 ; (windowBases<minLength)&&(a.seqText[rgtIx]!=0) ; rgtIx++)
		{
		if (a.qryText[rgtIx] == a.seqText[rgtIx]) matches++;
		                                     else errors++;
		if (a.qryText[rgtIx] != '-') windowBases++;
		}

	if (windowBases < minLength)     // this should not happen unless the
		return noAlignment;          // ..  alignment record is messed up

	lftIx = 0;
	for (lftIx=0 ; a.qryText[lftIx]=='-' ; lftIx++)
		errors--;

	// perform "slide-and-wide" algorithm;  consider, in turn, each window of
	// length L until we find one that meets the match criteria;  (L is counted
	// in query bases, and is initially equal to minLength);  note that the
	// subinterval we are searching for, if it exists, *must* contain such a
	// window;  whenever we find a suitable window, we lengthen it on the right
	// end as long as it meets the match criteria;  this new length becomes our
	// new L for the rest of the search

	testRatio = (1-minMRatio) / minMRatio;   // (see note 1)
	bestLftIx = bestRgtIx = (u32) -1;        // -1 indicates no window found

	while (true)
		{
		// "slide" the window until it meets the match criteria (or reaches the
		// end of the alignment

		while (true)
			{
			if (dbgSubalignment)
				print_current_window (a.seqText, a.qryText, lftIx, rgtIx,
				                      windowBases, matches, errors, minMRatio);

			if (a.seqText[rgtIx] == 0) goto search_done;
			if (matches * testRatio >= errors) break; // (see note 1)

			if (dbgSubalignment)
				fprintf (stderr, "-slide-\n");

			// move the right end one query base, maintaining this invariant
			// upon entrance to and exit from this block:
			//   qryText[rgtIx-1] is ACGTacgt
			//   qryText[rgtIx]   is ACGTacgt, '-', or the zero terminator

			while (true)
				{
				if (a.seqText[rgtIx] == 0) goto search_done;
				if (a.qryText[rgtIx] != '-') break;
				errors++;  rgtIx++;   // skip the '-' and count it as an error
				}

			if (a.qryText[rgtIx] == a.seqText[rgtIx]) matches++;
												 else errors++;
			rgtIx++;                  // move one position past the base

			// move the left end one query base, maintaining this invariant
			// upon entrance to and exit from this block:
			//   qryText[lftIx] is ACGTacgt

			if (a.qryText[lftIx] == a.seqText[lftIx]) matches--;
												 else errors--;
			lftIx++;

			while (a.qryText[lftIx] == '-')
				{
				errors--;  lftIx++;   // skip the '-' and uncount it as an error
				}
			}

		// "wide"; we have a window that meets the match criteria; lengthen it
		// until it doesn't, again maintaining this invariant upon entrance to
		// and exit from this block:
		//   qryText[rgtIx-1] is ACGTacgt
		//   qryText[rgtIx]   is ACGTacgt, '-', or the zero terminator
		// note that this loop will always lengthen the window unless we're
		// already at the terminator

		if (dbgSubalignment)
			fprintf (stderr, "-widen-\n");

		goodRgtIx = rgtIx;
		while (true)
			{
			if (a.seqText[rgtIx] == 0) break;
			if (a.qryText[rgtIx] == '-')
				{ errors++;  rgtIx++;  continue; }  // skip the '-'

			if (a.qryText[rgtIx] == a.seqText[rgtIx]) matches++;
												 else errors++;
			rgtIx++;  windowBases++;

			if (dbgSubalignment)
				print_current_window (a.seqText, a.qryText, lftIx, rgtIx,
				                      windowBases, matches, errors, minMRatio);

			if (matches * testRatio < errors) break;  // (see note 1)
			goodRgtIx = rgtIx;
			}

		if (dbgSubalignment)
			fprintf (stderr, "-widen done-\n");

		// record the new window as the longest

		bestLftIx = lftIx;
		bestRgtIx = goodRgtIx;
		}
search_done:

	// if we found no suitable window, say so

	if (bestLftIx == (u32) -1)
		return noAlignment;

	// we found a suitable window, create the subalignment record

	return subalignment (scoring, a, bestLftIx, bestRgtIx);
	}


// print_current_window--

static void print_current_window
   (char*	seqText,
	char*	qryText,
	u32		lftIx,
	u32		rgtIx,
	u32		qryBases,
	u32		matches,
	u32		errors,
	float	minMRatio)
	{
	float	testRatio = (1-minMRatio) / minMRatio;  // (see note 1)
	char*	isGood = (matches * testRatio >= errors)? " (good)" : "";
	u32		ix;
	
	fprintf (stderr, "[%u..%u] b=%u m=%u e=%u r=%.1f%%%s\n",
					 lftIx, rgtIx,
					 qryBases, matches, errors,
					 (100.0*matches)/(matches+errors), isGood);

	for (ix=lftIx ; ix<rgtIx ; ix++)
		fprintf (stderr, "%c", (seqText[ix]==qryText[ix])?'=':'x');
	fprintf (stderr, "\n");

	for (ix=lftIx ; ix<rgtIx ; ix++)
		fprintf (stderr, "%c", seqText[ix]);
	fprintf (stderr, "\n");

	for (ix=lftIx ; ix<rgtIx ; ix++)
		fprintf (stderr, "%c", qryText[ix]);
	fprintf (stderr, "\n");

	fprintf (stderr, "\n");
	}

//----------
//
// subalignment--
//	Extract a subinterval from an alignment.
//
//----------
//
// Arguments:
//	lascores*	scoring:	Alignment scoring definitions.
//	alignment	a:			The alignment to extract the subinterval of.
//	u32			lftIx:		The first column in the alignment subinterval.
//	u32			rgtIx:		The first column *past* the alignment subinterval.
//
// Returns:
//	The subalignment.  The caller is responsible for deallocating any memory
//	allocated for this result, by eventually calling free_alignment().
//
//----------

static alignment subalignment
   (lascores*	scoring,
	alignment	a,
	u32			lftIx,
	u32			rgtIx)
	{
	u32			ix;
	alignment	aSub;

	if ((dbgDebridgerDetail) && (dbgTriggerName == NULL) && (dbgTriggerCrc == 0))
		{
		fprintf (stderr, "subalignment(%d-%d)\n", lftIx, rgtIx);
		if ((rgtIx > (u32)strlen(a.seqText)) || (rgtIx > (u32)strlen(a.qryText)))
			{
			fprintf (stderr, "INTERNAL ERROR rgtIx=%d |seqText|=%u |qryText|=%u\n",
			                 rgtIx,(u32)strlen(a.seqText),(u32)strlen(a.qryText));
			exit (EXIT_FAILURE);
			}
		}

	aSub = noAlignment;
	aSub.active  = true;
	aSub.strand  = a.strand;
	aSub.seqText = copy_prefix (&a.seqText[lftIx], rgtIx-lftIx);
	aSub.qryText = copy_prefix (&a.qryText[lftIx], rgtIx-lftIx);
	rescore_alignment (scoring, &aSub);
	aSub.seqBaseCount = aSub.mCount + aSub.mmCount + aSub.iCount;
	aSub.qryBaseCount = aSub.mCount + aSub.mmCount + aSub.dCount;
	aSub.matchRatio   = ((float)aSub.mCount) / (aSub.mCount + aSub.mmCount + aSub.iCount + aSub.dCount);

	aSub.seqStart = a.seqStart;
	for (ix=0 ; ix<lftIx ; ix++)
		{ if (a.seqText[ix] != '-') aSub.seqStart++; }
	aSub.seqEnd = aSub.seqStart + aSub.seqBaseCount;

	return aSub;
	}

//----------
//
// report_debridged_alignment--
//	Given an alignment between a sequence and a motif loop, remove noisy
//	bridges and report the alignments between those bridges.
//
//----------
//
// Arguments:
//	FILE*		f:		The file to write to.
//	feed*		seq:	The sequence the alignment is in.
//	motif*		m:		The motif the alignment contains.
//	alignment	a:		The alignment to report (possibly in pieces).
//
// Control Globals:
//	debridgeMRatio:
//	debridgeLength:
//	clumpMRatio:
//	clumpLength:
//
// Returns:
//	(nothing)
//
//----------
//
// Bridge removal is a two stage process.  The implementation of stage 1 is
// in this function.  The implementation of stage 2 is in find_error_clumps(),
// which is *called* by this function.
//
// Stage 1 seeks to find intervals of the alignment for which the fraction of
// matches exceeds a threshold.  The intervals found are such that *all*
// prefixes and suffixes of an interval satisfy that threshold.
//
// The clustering algorithm used is adapted from Zhang, Berman, Miller, "Post-
// processing long pairwise alignments", Bioinformatics Vol. 15 no. 12 1999.
//
// An interval with M matches and X error sites satifies the threshold R if
// M/(M+X) >= R
//
//    M/(M+X) >= R
//    =>  M >= RM+RX
//    =>  (1-R)M - RX >= 0
//
// So we can identify such regions by keeping a running sum of (1-R)M - RX.
// The latter is accomplished by adding 1-R to the sum whenever we encounter a
// match, and subtracting R whenever we encounter an error.
//
//----------

static void report_debridged_alignment
   (FILE*		f,
	feed*		seq,
	motif*		m,
	alignment	a)
	{
	u32			top, sp;
	double		score, weight;
	u32			aIx, leftIx, rightIx;
	alignment	aSub, bSub;
	u32			subAlignmentLength;
	seglist*	clumps, *seg;
	u32			prevEnd;

	if (dbgDebridger)
		{
		u32 crc = alignment_crc(seq->name,(char*)m->forwardNucs,&a);
		if (dbgTriggerCrc != 0)
			{
			if (crc == dbgTriggerCrc)	// we've hit the trigger, so
				dbgTriggerCrc = 0;		// .. disable it
			fprintf(stderr,"triggering on crc %08X\n",crc);
			}
		if ((dbgTriggerName == NULL) && (dbgTriggerCrc == 0))
			{
			fprintf (stderr, "\n~~~~~ debridging this alignment ~~~~~\n");
			report_alignment_for_debug(stderr,seq,m,a);
			fprintf (stderr, "# crc %08X\n", crc);
			fprintf (stderr, "# |seqText|=%u |qryText|=%u\n",
			                 (u32)strlen(a.seqText), (u32)strlen(a.qryText));
			}
		}

	// build intervals

	top   = 0;
	score = 0.0;

	for (aIx=0 ; a.seqText[aIx]!=0 ; aIx++)
		{
		if (a.seqText[aIx] == a.qryText[aIx]) weight = 1-debridgeMRatio;
		                                 else weight =  -debridgeMRatio;
		score += weight;

		if ((dbgDebridgerStack) && (dbgTriggerName == NULL) && (dbgTriggerCrc == 0))
			{
			if (a.seqText[aIx] == a.qryText[aIx])
				fprintf (stderr, "%3d: = %5.2f", aIx, score);
			else
				fprintf (stderr, "%3d: x %5.2f", aIx, score);
			}

		if (weight < 0.0)
			{
			if ((dbgDebridgerStack) && (dbgTriggerName == NULL) && (dbgTriggerCrc == 0))
				fprintf (stderr, "\n");
			continue;
			}

		if ((top > 0) && (debridgeStack[top].right == aIx))
			{
			// add this site to interval on top of stack

			debridgeStack[top].right      = aIx+1;
			debridgeStack[top].rightScore = score;

			if ((dbgDebridgerStack) && (dbgTriggerName == NULL) && (dbgTriggerCrc == 0))
				fprintf (stderr, " extending [%d] %d-%d %4.1f %4.1f",
				        top,
				        debridgeStack[top].left,
				        debridgeStack[top].right,
				        debridgeStack[top].leftScore,
				        debridgeStack[top].rightScore);
			}
		else
			{
			// create a one site interval

			if (++top >= debridgeMaxStack) goto stack_overflow;

			debridgeStack[top].left       = aIx;
			debridgeStack[top].leftScore  = score - weight;
			debridgeStack[top].right      = aIx+1;
			debridgeStack[top].rightScore = score;
			debridgeStack[top].lower      = top - 1;

			while ((debridgeStack[top].lower > 0)
			    && (debridgeStack[debridgeStack[top].lower].leftScore > debridgeStack[top].leftScore))
				{ debridgeStack[top].lower = debridgeStack[debridgeStack[top].lower].lower; }

			if ((dbgDebridgerStack) && (dbgTriggerName == NULL) && (dbgTriggerCrc == 0))
				fprintf (stderr, " creating  [%d] %d-%d %4.1f %4.1f -> %d",
				        top,
				        debridgeStack[top].left,
				        debridgeStack[top].right,
				        debridgeStack[top].leftScore,
				        debridgeStack[top].rightScore,
						debridgeStack[top].lower);
			}

		// merge intervals;  if there is a previous interval with a no-higher
		// left score and no-higher right score, merge this interval (and all
		// intervening ones) into that one

		while ((top > 1)
		    && (debridgeStack[top].lower > 0)
			&& (debridgeStack[debridgeStack[top].lower].rightScore <= debridgeStack[top].rightScore))
			{
			debridgeStack[debridgeStack[top].lower].right      = debridgeStack[top].right;
			debridgeStack[debridgeStack[top].lower].rightScore = debridgeStack[top].rightScore;
			top = debridgeStack[top].lower;

			if ((dbgDebridgerStack) && (dbgTriggerName == NULL) && (dbgTriggerCrc == 0))
				fprintf (stderr, "\n%*s merging   [%d] %d-%d %4.1f %4.1f",
				        12, "", top,
				        debridgeStack[top].left,
				        debridgeStack[top].right,
				        debridgeStack[top].leftScore,
				        debridgeStack[top].rightScore);
			}

		if ((dbgDebridgerStack) && (dbgTriggerName == NULL) && (dbgTriggerCrc == 0))
			fprintf (stderr, "\n");
		}

	// report intervals from the stack

	for (sp=1 ; sp<=top ; sp++)
		{
		// convert interval to sub-alignment

		leftIx  = debridgeStack[sp].left;
		rightIx = debridgeStack[sp].right;

		aSub = subalignment (&control.scoring, a, leftIx, rightIx);

		if ((dbgDebridgerDetail) && (dbgTriggerName == NULL) && (dbgTriggerCrc == 0))
			fprintf (stderr, " stack     [%d] columns %d-%d"
			                 "  sequence interval %d..%d"
			                 "  length %d querybp %d\n",
			                 sp, leftIx, rightIx,
			                 aSub.seqStart, aSub.seqEnd,
			                 aSub.seqEnd - aSub.seqStart,
			                 aSub.qryBaseCount);

		// if the sub-alignment is too short, discard it and move on to the
		// next interval

		if (aSub.qryBaseCount < debridgeLength)
			{
			free_alignment (aSub);
			continue;
			}

		// otherwise, the sub-alignment is long enough; if we're not removing
		// error clumps, report it and move on to the next interval 

		if ((dbgDebridger) && (!dbgDebridgerDetail) && (dbgTriggerName == NULL) && (dbgTriggerCrc == 0))
			fprintf (stderr, " stack     [%d] columns %d-%d"
			                 "  sequence interval %d..%d"
			                 "  length %d querybp %d\n",
			                 sp, leftIx, rightIx,
			                 aSub.seqStart, aSub.seqEnd,
			                 aSub.seqEnd - aSub.seqStart,
			                 aSub.qryBaseCount);

		if (clumpMRatio == 0.0)
			{
			if (aSub.score >= (s32) control.scoring.minScore)
				report_alignment(f,seq,m,aSub);
			free_alignment (aSub);
			continue;
			}

		// remove error clumps and report any between-clump sub-sub-alignments
		// that are long enough

		clumps = find_error_clumps(seq,aSub);
		if (clumps == NULL)
			{
			// no clumps, so report the full alignment
			if (aSub.score >= (s32) control.scoring.minScore)
				report_alignment(f,seq,m,aSub);
			free_alignment (aSub);
			if (dbgClumpDetection) fprintf (stderr, "(no clumps removed)\n");
			continue;
			}

		subAlignmentLength = strlen(aSub.seqText);

		prevEnd = 0;
		for (seg=pop_segment(&clumps) ; seg!=NULL ; seg=pop_segment(&clumps))
			{
			if (seg->start > prevEnd)
				{
				bSub = subalignment (&control.scoring, aSub, prevEnd, seg->start);
				if (bSub.qryBaseCount >= debridgeLength)
					{
					if (bSub.score >= (s32) control.scoring.minScore)
						report_alignment(f,seq,m,bSub);
					}
				if (dbgClumpDetection)
					fprintf (stderr, "between clumps: %u-%u seq: %u-%u%s\n",
					                 prevEnd, seg->start,
					                 bSub.seqStart,bSub.seqEnd,
					                 (bSub.qryBaseCount>=debridgeLength)? "" : " (short, discarded)");
				free_alignment (bSub);
				}
			prevEnd = seg->end;
			unuse_segment (seg);
			}

		if (subAlignmentLength > prevEnd)
			{
			bSub = subalignment (&control.scoring, aSub, prevEnd, subAlignmentLength);
			if (bSub.qryBaseCount >= debridgeLength)
				{
				if (bSub.score >= (s32) control.scoring.minScore)
					report_alignment(f,seq,m,bSub);
				}
			if (dbgClumpDetection)
				fprintf (stderr, "between clumps: %u-%u seq: %u-%u%s\n",
				                 prevEnd, subAlignmentLength,
				                 bSub.seqStart,bSub.seqEnd,
				                 (bSub.qryBaseCount>=debridgeLength)? "" : " (short, discarded)");
			free_alignment (bSub);
			}

		free_alignment (aSub);
		}

	return;

	// failure exits

stack_overflow:
	fprintf (stderr, "debridger stack overflow (only supports %d intervals)\n"
	                 "consider using --allocate:debridger to increase the stack size\n",
	                 debridgeMaxStack-1);
	exit (EXIT_FAILURE);
	}

//----------
//
// find_error_clumps--
//	Find intervals in an alignment that have error density above some
//	threshold.
//
//----------
//
// Arguments:
//	feed*		seq:	The sequence the alignment is in.
//	alignment	a:		The alignment to detect clumps in (possibly in pieces).
//
// Control Globals:
//	clumpMRatio:
//	clumpLength:
//
// Returns:
//	A pointer to a list of segments to be removed from the alignment.  The
//	caller is responsible for returning these elements to the segment pool,
//	with unuse_segment().
//
//----------
//
// An interval with X error sites and M matches satifies the threshold R if
// X/(X+M) >= R.
//
//	   X/(X+M) >= R
//	   =>  X >= RX+RM
//	   =>  (1-R)X - RM >= 0
//
// So we can identify such regions by keeping a running sum of (1-R)X - RM.
// The latter is accomplished by adding 1-R to the sum whenever we encounter an
// error and subtracting R whenever we encounter a match.  Any interval for
// which the running sum is at least as high on the right as on the left
// satisfies the threshold.  Any such interval is a clump unless it is covered
// by a longer clump.
//
// Note that it *is* possible for two overlapping intervals to each be above
// the threshold while the combined interval is below the threshold.  This
// function returns such combined intervals, the rationale being that every
// position within the interval is part of *some* interval that's above the
// threshold.
//
//----------

// more_clump_segments-- (local support)

static void more_clump_segments (u32 atLeast);
static void more_clump_segments (u32 atLeast)
	{
	numClumpSegments += 1000 + (numClumpSegments/10);
	if (numClumpSegments < atLeast) numClumpSegments = atLeast;

	clumpMinSums = (double*) malloc ((numClumpSegments+1) * sizeof(double));
	if (clumpMinSums == NULL)
		{
		fprintf (stderr, "failed to re-allocate clump detection (for %d entries)\n",
						 numClumpSegments);
		exit (EXIT_FAILURE);
		}
	if (dbgAllocation)
		{
		fprintf (stderr, "  re-allocated %ld bytes for clump detection sums (for %d entries)\n",
						 (numClumpSegments+1) * sizeof(double), numClumpSegments);

		}

	clumpMinWhere = (u32*) malloc ((numClumpSegments+1) * sizeof(u32));
	if (clumpMinWhere == NULL)
		{
		fprintf (stderr, "failed to re-allocate clump detection (for %d entries)\n",
						 numClumpSegments);
		exit (EXIT_FAILURE);
		}
	if (dbgAllocation)
		{
		fprintf (stderr, "  re-allocated %ld bytes for clump detection wheres (for %d entries)\n",
						 (numClumpSegments+1) * sizeof(u32), numClumpSegments);

		}

	}


// find_error_clumps--

static seglist* find_error_clumps
   (feed*		seq,
	alignment	a)
	{
	u32			minLength = clumpLength;
	double		rate      = clumpMRatio;
	u32			alignmentLength, aIx;
	double		minSum, valSum, val;
	u32			numMinSums, minScan, minIx;
	u32			start, end, prevStart, prevEnd;
	seglist*	clumps = NULL, *seg = NULL;
	seglist*	collapsedClumps = NULL;
	int			hasErrors;

	if (dbgClumpDetection)
		fprintf (stderr, "\nlooking for clumps in %s %c %u-%u\n",
		                 seq->name,a.strand,a.seqStart,a.seqEnd);

	// prescan to see whether the alignment contains any errors at all; if it
	// doesn't, we can just return now; the main benefit of this is that such
	// an alignment is the worst case for this algorithm's memory allocation

	hasErrors = false;
	for (aIx=0 ; a.seqText[aIx]!=0 ; aIx++)
		{
		if (a.seqText[aIx] != a.qryText[aIx]) hasErrors = true;
		}

	alignmentLength = aIx;

	if (!hasErrors)
		{
		if (dbgClumpDetection)
			fprintf (stderr, "the alignment has no errors, hence it has no clumps\n");
		return NULL;
		}

	// search for clumps, intervals for which the average is above (or at) the
	// threshold;  this is equivalent to intervals in which the sum, minus the
	// length times the average, is positive or zero
	//
	// nota bene: this algorithm can hypothetically suffer from round off error
	//            in the sums, which can make the reported intervals less
	//            precise than we'd like

	minSum = valSum = 0.0;

	clumpMinSums [0] = valSum;
	clumpMinWhere[0] = (u32) -1;
	numMinSums  = 1;
	minScan     = 0;

	prevStart   = (u32) -1;
	prevEnd     = (u32) -1;

	for (aIx=0 ; aIx<alignmentLength ; aIx++)
		{
		// invariant: minSums[minScan] <= valSum < minSums[minScan-1]
		// with the implied assumption that minSums[-1] == infinity

		if (a.seqText[aIx] != a.qryText[aIx]) val = 1-rate;
		                                 else val =  -rate;

		valSum += val;
		if (valSum < minSum)
			{
			if (numMinSums >= numClumpSegments)
				more_clump_segments (alignmentLength/3);
			minSum = valSum;
			clumpMinSums [numMinSums] = valSum;
			clumpMinWhere[numMinSums] = aIx;
			numMinSums++;
			}

		// (re-establish the invariant)
		if (val < 0)								// the sum has decreased
			{
			while (clumpMinSums[minScan] > valSum)
				minScan++;
			}
		else if (val > 0)							// the sum has increased
			{
			while ((minScan > 0) && (clumpMinSums[minScan-1] <= valSum))
				minScan--;
			}

		if (dbgClumpCharacters)
			fprintf (stderr, " [%u] %s %.2f %.2f %u..\n",
			                 aIx,(val<0)?"=":"x",val,valSum,clumpMinWhere[minScan]);

		// minScan points at the earliest index with minSums[minScan] <= valSum,
		// so the interval (minWhere[minScan]+1 to aIx) has sum >= 0

		minIx = clumpMinWhere[minScan];
		if (aIx - minIx < minLength) continue;

		start = minIx+1;
		end   = aIx;

		if (dbgClumpSegments)
			fprintf (stderr, "      setting %u..%u %.2f %.2f %.2f\n",
			                  start,end+1,
			                  clumpMinSums[minScan],valSum,
			                  clumpMinSums[minScan]-valSum);

		if ((prevStart == (u32) -1)
		 || (start > prevEnd+1))
			{
			// no overlap with previous interval (or no previous interval)
			if (dbgClumpSegments)
				fprintf (stderr, "      adding (%d,%d)\n",start,end+1);
			seg       = add_segment (&clumps,start,end+1);
			prevStart = start;
			prevEnd   = end;
			}
		else if (start >= prevStart)
			{
			// new interval extends previous interval on right
			if (dbgClumpSegments)
				fprintf (stderr, "      replacing (%d,%d) with (%d,%d)\n",
				                 seg->start,seg->end,seg->start,end+1);
			seg->end = end+1;
			prevEnd  = end;
			}
		else
			{
			// new interval extends previous interval on left and right
			if (dbgClumpSegments)
				fprintf (stderr, "      replacing (%d,%d) with (%d,%d)\n",
				                 seg->start,seg->end,start,end+1);
			seg->start = start;
			seg->end   = end+1;
			prevStart  = start;
			prevEnd    = end;
			}
		}

	// if we found no clumps, we're done

	if (clumps == NULL)
		return NULL;

	// nota bene: at this point, clumps are in right-to-left order (the head
	// of the list contains the rightmost clump)

	//if (dbgClumpDetection)
	//	{
	//	for (seg=clumps ; seg!=NULL ; seg=seg->next)
	//		fprintf (stderr, "untrimmed clump: %u-%u\n", seg->start,seg->end);
	//	}

	// trim any non-error ends off of the clumps

	for (seg=clumps ; seg!=NULL ; seg=seg->next)
		{
		while ((seg->start < seg->end) && (a.seqText[seg->start] == a.qryText[seg->start]))
			seg->start++;
		while ((seg->start < seg->end) && (a.seqText[seg->end-1] == a.qryText[seg->end-1]))
			seg->end--;
		}

	//if (dbgClumpDetection)
	//	{
	//	for (seg=clumps ; seg!=NULL ; seg=seg->next)
	//		fprintf (stderr, "trimmed clump:   %u-%u\n", seg->start,seg->end);
	//	}

	// collapse overlapping clumps;  as a side effect, the resulting list is in
	// left-to-right order

	start = end = (u32) -1;
	for (seg=pop_segment(&clumps) ; seg!=NULL ; seg=pop_segment(&clumps))
		{
		if (end == (u32) -1)
			{ start = seg->start;  end = seg->end; }
		else if (seg->end >= start)
			{ if (seg->start < start) start = seg->start; }
		else
			{
			add_segment (&collapsedClumps,start,end);
			start = seg->start;  end = seg->end;
			}
		unuse_segment (seg);
		}

	if (end != (u32) -1)
		add_segment (&collapsedClumps,start,end);

	if (dbgClumpDetection)
		{
		for (seg=collapsedClumps ; seg!=NULL ; seg=seg->next)
			{
			u32 positiveCount = 0;
			u32 negativeCount = 0;
			for (aIx=seg->start ; aIx<seg->end ; aIx++)
				{
				if (a.seqText[aIx] != a.qryText[aIx]) negativeCount++;
				                                 else positiveCount++;
				}
			fprintf (stderr, "clump to remove: %u-%u m=%u x=%u mRatio=%.2f%%\n",
			                 seg->start,seg->end,positiveCount,negativeCount,
							 (100.0*positiveCount)/(positiveCount+negativeCount));
			}
		}

	return collapsedClumps;
	}

//----------
//
// report_alignment--
//	Report one alignment between a sequence and a motif loop.
//
//----------
//
// Arguments:
//	FILE*		f:		The file to write to.
//	feed*		seq:	The sequence the alignment is in.
//	motif*		m:		The motif the alignment contains.
//	alignment	a:		The alignment to report.
//
// Returns:
//	(nothing)
//
//----------

// glue functions

static void report_alignment_core (FILE* f, feed* seq, motif* m, alignment a, int countIt);

static void report_alignment(FILE* f, feed* seq, motif* m, alignment a)
	{ report_alignment_core(f,seq,m,a,false); }

static void report_alignment_for_debug(FILE* f, feed* seq, motif* m, alignment a)
	{ report_alignment_core(f,seq,m,a,true); }


// report_alignment_core

static void report_alignment_core
   (FILE*		f,
	feed*		seq,
	motif*		m,
	alignment	a,
	int			forDebug)
	{
	static		int firstReport = true;
	char		rangeStr[100];		// WARNING:  vulnerable to string overflow!
	char		scoreStr[100];
	char		motifWithStrand[1000];
	char		seqLengthStr[100];
	char		seqBaseCountStr[100], qryBaseCountStr[100];
	char		statsStr[300];
	u32			nameW, lengthW, countW, rangeW;
	u32			prefixLen = 0, suffixLen = 0;
	char*		prefix, *suffix, *tempPtr;
	char		prefixStr[1], suffixStr[1];

	if ((dbgReportAlignments) && (!forDebug))
		fprintf (stderr, "aligned %s %d-%d (%s bp)\n",
		                 seq->name, a.seqStart, a.seqEnd,
		                 ucommatize(a.seqEnd-a.seqStart));

	// fetch flanks (for debugging)

	if (dbgFlankLength == 0)
		{
		prefix = prefixStr;  prefixStr[0] = 0;
		suffix = suffixStr;  suffixStr[0] = 0;
		}
	else
		{
		prefixLen = suffixLen = dbgFlankLength;
		if (a.seqStart < prefixLen) prefixLen = a.seqStart;
		if (seq->len - a.seqEnd < suffixLen) suffixLen = seq->len - a.seqEnd;

		prefix = prefixFlank;
		suffix = suffixFlank;

		if (prefixLen == 0)
			prefix[0] = 0;
		else
			{
			strncpy (/*to*/ prefix, /*from*/ (char*) &seq->nt[a.seqStart-prefixLen], /*chars*/ prefixLen);
			string_to_lowercase (prefix);
			}

		if (suffixLen == 0)
			suffix[0] = 0;
		else
			{
			strncpy (/*to*/ suffix, /*from*/ (char*) &seq->nt[a.seqEnd],             /*chars*/ suffixLen);
			string_to_lowercase (suffix);
			}
		}

	// flip the alignment, if needed

	if ((showForward) && (a.strand == '-'))
		{
		u32 temp;
		reverse_complement ((u8*) a.seqText);
		reverse_complement ((u8*) a.qryText);
		temp       = a.seqStart;
		a.seqStart = a.seqEnd;
		a.seqEnd   = temp;

		if (dbgFlankLength > 0)
			{
			if (prefixLen > 0) reverse_complement ((u8*) prefix);
			if (suffixLen > 0) reverse_complement ((u8*) suffix);
			temp      = prefixLen;
			prefixLen = suffixLen;
			suffixLen = temp;
			tempPtr   = prefix;
			prefix    = suffix;
			suffix    = tempPtr;
			}
		}

	// determine field widths so that things will line up nicely

	if (!forDebug)
		{
		if (firstReport) firstReport = false;
		            else fprintf (f, "\n");
		}

	if (strlen(m->name)+1+1 > sizeof(motifWithStrand))
		{
		fprintf (stderr, "INTERNAL ERROR\n"
		                 "length of motif name exceeds 'motifWithStrand' buffer in report_alignment (%d>%d)\n"
		                 "give motif a short name (or increase buffer size and recompile)\n",
		                 (int) strlen(m->name)+1, (int) sizeof(motifWithStrand)-1);
		exit (EXIT_FAILURE);
		}

	sprintf (rangeStr, "%u-%u", a.seqStart, a.seqEnd);
	sprintf (scoreStr, "score=%d", a.score);
	sprintf (motifWithStrand, "%s%c", m->name, a.strand);
	sprintf (seqLengthStr,    "%d",   seq->len);
	sprintf (seqBaseCountStr, "%dbp", a.seqBaseCount);
	sprintf (qryBaseCountStr, "%dbp", a.qryBaseCount);

	nameW = nameFieldW;
	if (strlen(seq->name)       > nameW) nameW = strlen(seq->name);
	if (strlen(motifWithStrand) > nameW) nameW = strlen(motifWithStrand);

	lengthW = lengthFieldW;
	if (strlen(seqLengthStr) > lengthW) lengthW = strlen(seqLengthStr);

	countW = countFieldW;
	if (strlen(seqBaseCountStr) > countW) countW = strlen(seqBaseCountStr);
	if (strlen(qryBaseCountStr) > countW) countW = strlen(qryBaseCountStr);

	rangeW = rangeFieldW;
	if (strlen(rangeStr) > rangeW) rangeW = strlen(rangeStr);
	if (strlen(scoreStr) > rangeW) rangeW = strlen(scoreStr);

	// print optional alignment events line

	if (showEvents)
		{
		int   ix;

		sprintf (statsStr, "# score=%d querybp=%d mRatio=%.1f%% m=%d mm=%d i=%d d=%d",
		         a.score,a.qryBaseCount,100*a.matchRatio,
		         a.mCount,a.mmCount,a.iCount,a.dCount);
		if (strlen(statsStr) > nameW+lengthW+countW+rangeW+3)
			nameW = strlen(statsStr) - (lengthW+countW+rangeW+3);

		fprintf (f, "%-*s ", nameW+lengthW+countW+rangeW+3, statsStr);
		if (prefixLen > 0)
			{ for (ix=0 ; ix<(int)prefixLen ; ix++) fprintf (f, "."); }
		for (ix=0 ; a.seqText[ix]!=0 ; ix++)
			{
			if (a.seqText[ix] == a.qryText[ix]) fprintf (f, "=");
			                               else fprintf (f, "x");
			}
		if (suffixLen > 0)
			{ for (ix=0 ; ix<(int)suffixLen ; ix++) fprintf (f, "."); }
		fprintf (f, "\n");
		}

	// print the alignment

	fprintf (f, "%-*s %-*s %-*s %-*s %s%s%s\n",
	            nameW,   seq->name,
	            lengthW, seqLengthStr,
	            countW,  seqBaseCountStr,
	            rangeW,  rangeStr,
	            prefix, a.seqText, suffix);
	fprintf (f, "%-*s %-*s %-*s %-*s %*s%s\n",
	            nameW,   motifWithStrand,
	            lengthW, "",
	            countW,  qryBaseCountStr,
	            rangeW,  scoreStr,
	            prefixLen, "", a.qryText);

	if (!forDebug) alignmentsReported++;

	// print optional alignment positional events lines

	if (showEventsByPos)
		{
		int motifLen = strlen ((char*) m->forwardNucs);
		int motifPos;

		count_positional_events(m,a);

		for (motifPos=0 ; motifPos<motifLen ; motifPos++)
			{
			evcount* ev = &eventCounts[motifPos];
			double matchRatio = ((float)ev->m) / (ev->m + ev->mm + ev->i + ev->d);
			fprintf (f, "# position %d [%c] mRatio=%.1f%% m=%d mm=%d i=%d d=%d"
			            " mmA=%d mmC=%d mmG=%d mmT=%d"
			            " x=%d\n",
			            motifPos,  m->forwardNucs[motifPos],
			            100*matchRatio, ev->m, ev->mm, ev->i, ev->d,
			            ev->mmA, ev->mmC, ev->mmG, ev->mmT,
			            ev->mm + ev->i + ev->d);
			}
		}

	// print other optional alignment lines

	if (showCrc)
		{
		u32 crc = alignment_crc(seq->name,(char*)m->forwardNucs,&a);
		fprintf (f, "# crc %08X\n", crc);
		}

	}

//----------
//
// count_positional_events--
//	Count the error events for an alinment, separately for each position in the
//	motif.
//
//----------
//
// Arguments:
//	motif*		m:	The motif the alignment contains.
//	alignment*	a:	The alignment to inspect.
//
// Returns:
//	nothing;  eventCounts[ix] is modified, for each position ix in the motif
//
//----------

static void count_positional_events
   (motif*		m,
	alignment	a)
	{
	int			motifLen = strlen ((char*) m->forwardNucs);
	char*		seqTextScan, *qryTextScan;
	evcount*	ev;
	int			motifPos, motifDirection, motifStart, motifIx;
	int			isMatch;
	u8			qryNuc;

	// determine position in motif of the start of the alignment, and whether
	// we're using the forward or reverse motif

	motifPos = -1;
	motifDirection = 0;

	for (motifStart=0 ; motifStart<motifLen ; motifStart++)
		{
		motifIx = motifStart;
		isMatch = true;
		for (qryTextScan=a.qryText ; *qryTextScan!=0 ; qryTextScan++)
			{
			if (*qryTextScan == '-') continue;
			qryNuc = toUpperACGTN((u8) *qryTextScan);
			if (qryNuc != m->forwardNucs[motifIx]) { isMatch = false;  break; }
			motifIx++;
			if (motifIx == motifLen) motifIx = 0;
			}
		if (isMatch)
			{ motifPos = motifStart;  motifDirection = 1;  break; }
		}

	if ((motifPos == -1) && (m->reverseNucs != NULL))
		{
		for (motifStart=0 ; motifStart<motifLen ; motifStart++)
			{
			motifIx = motifStart;
			isMatch = true;
			for (qryTextScan=a.qryText ; *qryTextScan!=0 ; qryTextScan++)
				{
				if (*qryTextScan == '-') continue;
				qryNuc = toUpperACGTN((u8) *qryTextScan);
				if (qryNuc != m->reverseNucs[motifIx]) { isMatch = false;  break; }
				motifIx++;
				if (motifIx == motifLen) motifIx = 0;
				}
			if (isMatch)
				{ motifPos = motifLen-1 - motifStart;  motifDirection = -1;  break; }
			}
		}

	if (motifPos == -1)
		{ // no match found, just compute relative to position 0, forward
		motifPos = 0;
		motifDirection = 1;
		}

	// clear counters

	for (motifIx=0 ; motifIx<motifLen ; motifIx++)
		{
		ev = &eventCounts[motifIx];
		ev->m = ev->mm = ev->i = ev->d
		      = ev->mmA = ev->mmC = ev->mmG = ev->mmT = 0;
		}

	// count events

	seqTextScan = a.seqText;
	qryTextScan = a.qryText;

	for ( ; *seqTextScan!=0 ; seqTextScan++,qryTextScan++)
		{
		ev = &eventCounts[motifPos];
		if (*qryTextScan == '-')
			{ ev->i++; }
		else if (*seqTextScan == '-')
			{ ev->d++;  motifPos += motifDirection; }
		else if (*qryTextScan == *seqTextScan)
			{ ev->m++;  motifPos += motifDirection; }
		else
			{
			ev->mm++;
			if (a.strand == '-')
				{
				switch (*seqTextScan)
					{
					case 'A': ev->mmT++; break;
					case 'C': ev->mmG++; break;
					case 'G': ev->mmC++; break;
					case 'T': ev->mmA++; break;
					default:             break;
					}
				}
			else // if (a.strand == '+')
				{
				switch (*seqTextScan)
					{
					case 'A': ev->mmA++; break;
					case 'C': ev->mmC++; break;
					case 'G': ev->mmG++; break;
					case 'T': ev->mmT++; break;
					default:             break;
					}
				}
			motifPos += motifDirection;
			}

		if      (motifPos < 0)         motifPos = motifLen-1;
		else if (motifPos >= motifLen) motifPos = 0;
		}

	}

//----------
//
// segment list operators--
//
//----------

static void free_segment_pool
   (void)
	{
	seglist*	s, *sNext;

	for (s=segmentPool ; s!=NULL ; s=sNext)
		{ sNext = s->next;  free (s); }

	segmentPool = NULL;
	}


static seglist*  allocate_segment
   (void)
	{
	seglist*	s;

	// allocate a new record, or reclaim one from the pool

	if (segmentPool == NULL)
		{
		s = (seglist*) malloc (sizeof(seglist));
		if (s == NULL)
			{
			fprintf (stderr, "failed to allocate segment record\n");
			exit (EXIT_FAILURE);
			}
		if (dbgAllocation)
			{
			fprintf (stderr, "  allocated %ld bytes for a segment record\n",
			                 sizeof(seglist));

			}
		}
	else
		{
		s = segmentPool;
		segmentPool = segmentPool->next;
		}

	return s;
	}


static seglist* add_segment
   (seglist**	segmentList,
	u32			start,
	u32			end)
	{
	seglist*	s;

	if ((dbgSegments) && (dbgTriggerName == NULL))
		fprintf (stderr, "add_segment(%u,%u)\n", start, end);

	s = allocate_segment();

	// add the record to the head of the list

	s->next  = (*segmentList);
	s->start = start;
	s->end   = end;

	(*segmentList) = s;
	return s;
	}


static seglist* pop_segment
   (seglist**	segmentList)
	{
	seglist*	s;

	// if the list is empty, return NULL as the popped record

	if (*segmentList == NULL)
		return NULL;

	// detach the record from the list

	s = *segmentList;
	(*segmentList) = (*segmentList)->next;

	if ((dbgSegments) && (dbgTriggerName == NULL))
		fprintf (stderr, "pop_segment() --> %u,%u\n", s->start, s->end);

	return s;
	}


static void unuse_segment
   (seglist*	s)
	{
	// add the record to the head of the pool

	s->next     = segmentPool;
	segmentPool = s;
	}

