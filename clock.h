#ifndef clock_H					// (prevent multiple inclusion)
#define clock_H

#include "utilities.h"

// establish ownership of global variables

#ifdef clock_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef clock_owner
int		dbgShowclock   = false;
#else
global int dbgShowclock;
#endif

//----------
//
// prototypes for functions in this module--
//
//----------

void reset_progress_clock  (void);
void report_progress_clock (FILE *f, char* fmt);

#undef global
#endif // clock_H
