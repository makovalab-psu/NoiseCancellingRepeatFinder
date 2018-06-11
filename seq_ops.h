#ifndef seq_ops_H				// (prevent multiple inclusion)
#define seq_ops_H

#include "utilities.h"

//----------
//
// prototypes for functions in this module--
//
//----------

void reverse_complement (u8* nt);
void unmask_sequence    (u8* nt);

#undef global
#endif // seq_ops_H
