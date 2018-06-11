// clock.c-- support for crude time profiling.

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include "utilities.h"

#define  clock_owner			// (make this the owner of its globals)
#include "clock.h"				// interface to this module

//----------
//
// interface to underlying system clock operations--
//
//----------

//#define useStandardClock  (define at build time)

#ifdef useStandardClock
#define read_clock() clock()
#define clocksPerSec CLOCKS_PER_SEC
#endif // useStandardClock

#ifndef useStandardClock
#include <sys/time.h>
#define read_clock() microsec_clock()
#define clocksPerSec 1000000
static u64 microsec_clock (void);
static u64 microsec_clock (void)
	{
	static int		failed = false;
	struct timeval	time1;
	int				err;

	if (failed)	return 0;	// (previous call to gettimeofday has failed)
	err = gettimeofday (&time1, NULL);
	if (err != 0) { failed = true;  return 0; }

	return (((u64) time1.tv_sec) * 1000000) + time1.tv_usec;
	}
#endif // not useStandardClock

//----------
//
// reset_progress_clock--
//	Reset the clock.  Subsequent calls to report_progress_clock() will report
//	time measured from this moment.
//
//----------
//
// Arguments:
//	(none)
//
// Returns:
//	(nothing)
//
//----------

static s64 progressClock = 0;

void reset_progress_clock (void)
	{
	progressClock = -((s64) read_clock());
	}

//----------
//
// report_progress_clock--
//	Report (to the user) the elpased time since the last call to
//	reset_progress_clock().
//
//----------
//
// Arguments:
//	FILE*	f:		The file to write the report to.
//	char*	fmt:	A format string.  This is NOT a printf format string.  If
//					.. this contains the substring "{time}", the report will
//					.. replace that substring with a representation of the
//					.. elapsed time.  If that substring is absent, the report
//					.. will put the elpased time before the string.
//
// Returns:
//	(nothing)
//
//----------

void report_progress_clock
   (FILE*	f,
	char*	fmt)
	{
	char*	fmtSuffix;
	float	secs;
	int		hours, mins;

	fmtSuffix = strstr (fmt, "{time}");
	if (fmtSuffix != NULL)
		{
		while (fmt < fmtSuffix) fprintf (f, "%c", *(fmt++));
		fmt += strlen("{time}");
		}

	progressClock += (s64) read_clock();
	secs = ((float)(progressClock)) / clocksPerSec;

	if (secs < 60)
		fprintf (f, "%.3fs%s", secs, fmt);
	else if (secs < 3600)
		{
		mins =  secs / 60;
		secs -= 60 * mins;
		fprintf (f, "%dm%06.3fs%s", mins, secs, fmt);
		}
	else
		{
		mins  =  secs / 60;
		secs  -= 60 * mins;
		hours =  mins / 60;
		mins  -= 60 * hours;
		fprintf (f, "%dh%02dm%06.3fs%s", hours, mins, secs, fmt);
		}
	}

