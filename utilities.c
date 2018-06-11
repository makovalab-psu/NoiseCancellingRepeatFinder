// utilities.c-- miscellaneous utility functions.

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <inttypes.h>
#include <math.h>
#include <float.h>

#define  utilities_owner		// (make this the owner of its globals)
#include "utilities.h"			// interface to this module

//----------
//
// copy_string, copy_prefix--
//	Create (in the heap) a copy of a string or a prefix of a string.
//
//----------
//
// Arguments:
//	const char*	s:	The string to copy.
//	u32			n:	(copy_prefix only) the number of characters to copy.
//
// Returns:
//	A pointer to new string;  failures result in fatality.
//
//----------

char* copy_string
   (const char*	s)
	{
	char*		ss;

	if (s == NULL) return NULL;

	ss = malloc (strlen(s) + 1);
	if (ss == NULL) goto cant_allocate;
	return strcpy (/*to*/ ss, /*from*/ s);

cant_allocate:
	fprintf (stderr, "failed to allocate copy of string \"%s\"\n", s);
	exit (EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


char* copy_prefix
   (const char*	s,
	u32			n)
	{
	char*		ss;

	if (s == NULL) return NULL;

	ss = malloc (n + 1);
	if (ss == NULL) goto cant_allocate;
	memcpy (/*to*/ ss, /*from*/ s, /*how much*/ n);
	ss[n] = 0;
	return ss;

cant_allocate:
	fprintf (stderr, "failed to allocate %d bytes for copy of string \"%s\"\n", n, s);
	exit (EXIT_FAILURE);
	return NULL; // (never reaches here)
	}

//----------
//
// strcmp_prefix--
//	Determine if a string contains another as a prefix.
//
//----------
//
// Arguments:
//	const char*	str1:	The string.
//	const char*	str2:	The prefix string.
//
// Returns:
//	The same as strcmp(prefix1,str2) would, where prefix1 is str1 truncated
//	to be no longer than str2.
//
//----------

int strcmp_prefix
   (const char*	str1,
	const char*	str2)
	{
	return strncmp (str1, str2, strlen(str2));
	}

//----------
//
// string_to_lowercase--
//	Convert a string to lowercase.
//
//----------
//
// Arguments:
//	char*	s:	The string to convert.
//
// Returns:
//	(nothing)
//
//----------

void string_to_lowercase
   (char*	s)
	{
	char*	ss;

	for (ss=s ; *ss!=0 ; ss++)
		{ if ((*ss >= 'A') && (*ss <= 'Z')) *ss += 'a'-'A'; }
	}

//----------
//
// string_to_u32--
//	Parse a string for the unsigned integer value it contains.
//
//----------
//
// Arguments:
//	const char*	s:	The string to parse.
//
// Returns:
//	The integer value of the string.  Note that the string *must not* contain
//	anything other than a valid integer-- failures result in program
//	termination.
//
//----------

u32 string_to_u32
   (const char*	s)
	{
	char*		ss;
	u32			v;
	char		extra;

	// skip to first non-blank

	ss = (char*) s;
	while ((*ss == ' ') || (*ss == '\t') || (*ss == '\n'))
		ss++;
	if (*ss == 0) goto empty_string;

	// if it begins with 0x, parse it as hex

	if (strcmp_prefix(ss,"0x") == 0)
		{
		if (sscanf (ss, "%X%c", &v, &extra) != 1) goto not_an_integer;
		return v;
		}

	// convert to number

	if (*ss == '-') goto not_an_integer;
	if (sscanf (ss, "%u%c", &v, &extra) != 1) goto not_an_integer;

	return v;

	//////////
	// failure exits
	//////////

empty_string:
	fprintf (stderr, "an empty string is not an unsigned integer\n");
	exit (EXIT_FAILURE);

not_an_integer:
	fprintf (stderr, "\"%s\" is not an unsigned integer\n", s);
	exit (EXIT_FAILURE);

	return 0;
	}

//----------
//
// string_to_unitized_u32--
//	Parse a string for the integer value it contains, allowing K, M, and G
//	suffixes.
//
//----------
//
// Arguments:
//	const char*	s:		The string to parse.
//	int	byThousands:	true  => K means one thousand
//						false => K means 1,024.
//
// Returns:
//	The integer value of the string.  Note that the string *must not* contain
//	anything other than a valid integer-- failures result in program
//	termination.
//
//----------

u32 string_to_unitized_u32
   (const char*	s,
	int			byThousands)
	{
	char		ss[20];
	int			len = strlen (s);
	char*		parseMe;
	u32			v;
	float		vf;
	char		extra;
	u32			mult;
	int			isFloat;

	// convert to number

	mult = 1;

	if (len >= (int) sizeof (ss))
		parseMe = (char*) s;
	else
		{
		parseMe = ss;
		strcpy (ss, s);

		if (len > 0)
			{
			switch (ss[len-1])
				{
				case 'K': case 'k':
					mult = (byThousands)? 1000 : 1024;
					break;
				case 'M': case 'm':
					mult = (byThousands)? 1000000 : 1024L * 1024L;
					break;
				case 'G': case 'g':
					mult = (byThousands)? 1000000000 : 1024L * 1024L * 1024L;
					break;
				}

			if (mult != 1)
				ss[len-1] = 0;
			}
		}

	isFloat = false;
	if (sscanf (parseMe, "%u%c", &v, &extra) != 1)
		{
		if (sscanf (parseMe, "%f%c", &vf, &extra) != 1) goto bad;
		isFloat = true;
		}

	if (isFloat)
		{
		if  (vf < 0) goto bad;
		if ((vf > 0) && (vf*mult+.5 > u32Max)) goto overflow;
		v = vf*mult + .5;
		}
	else if (mult != 1)
		{
		if ((v > 0) && ( v > u32Max / mult)) goto overflow;
		v *= mult;
		}

	return v;

bad:
	fprintf (stderr, "\"%s\" is not an unsigned integer\n", s);
	exit (EXIT_FAILURE);

overflow:
	fprintf (stderr, "\"%s\" is out of range for an unsigned integer\n", s);
	exit (EXIT_FAILURE);
	}

//----------
//
// string_to_double--
//	Parse a string for the double floating point value it contains.
//
//----------
//
// Arguments:
//	const char*	s:	The string to parse.
//
// Returns:
//	The value of the string.  Note that the string *must not* contain anything
//	other than a valid number-- failures result in fatality.
//
//----------

double string_to_double
   (const char*	s)
	{
	char*		ss;
	double		v;
	char		extra, extra2;

	// skip to first non-blank

	ss = (char*) s;
	while ((*ss == ' ') || (*ss == '\t') || (*ss == '\n'))
		ss++;
	if (*ss == 0) goto empty_string;

	// check for named constants

	if (strcmp (s,    "inf") == 0) return  DBL_MAX;
	if (strcmp (s,   "+inf") == 0) return  DBL_MAX;
	if (strcmp (s,   "-inf") == 0) return -DBL_MAX;
	if (strcmp (s,  "1/inf") == 0) return  DBL_MIN;
	if (strcmp (s, "+1/inf") == 0) return  DBL_MIN;
	if (strcmp (s, "-1/inf") == 0) return -DBL_MIN;

	// convert to number

	if (s[strlen(s)-1] == '%')
		{
		if (sscanf (s, "%lf%c%c", &v, &extra, &extra2) != 2) goto not_a_number;
		v /= 100.0;
		}
	else
		{
		if (sscanf (s, "%lf%c", &v, &extra) != 1) goto not_a_number;
		}

	return v;

empty_string:
	fprintf (stderr, "an empty string is not a number\n");
	exit (EXIT_FAILURE);

not_a_number:
	fprintf (stderr, "\"%s\" is not a number\n", s);
	exit (EXIT_FAILURE);
	}

//----------
//
// ucommatize--
//	Convert an integer to a string, including commas.
//
//----------
//
// Arguments:
//	const u64 v:	The number to convert.
//
// Returns:
//	A string representing that number, including commas.  (see note 1)
//
//----------
//
// notes:
//
// (1)	The memory containing the returned string belongs to this routine, as
//		static memory.  There are only five such memory blocks, and they are
//		used on successive calls.  So when you make more than five calls, the
//		results of previous calls are clobbered.
//
//----------

char* ucommatize
   (const u64	v)
	{
	static char	 s1[52];// (big enough for 128-bit decimal value with commas,
	static char	 s2[52];//  .. the biggest being
	static char	 s3[52];//  .. 340,282,366,920,938,463,463,374,607,431,768,211,455)
	static char	 s4[52];
	static char	 s5[52];
	static char* s = s5;
	int		len, commas;
	char*	src, *dst;

	if      (s == s1) s = s2;	// (ping pong)
	else if (s == s2) s = s3;
	else if (s == s3) s = s4;
	else if (s == s4) s = s5;
	else              s = s1;

	sprintf (s, "%jd", (intmax_t) v);	// $$$ this could overflow the buffer
										// $$$ .. if int_max_t > 128 bits

	len = strlen (s);
	commas = (len-1) / 3;

	if (commas != 0)
		{
		src = s + len - 1;
		dst = s + len + commas;  *(dst--) = 0;

		while (dst > src)
			{
			*(dst--) = *(src--);
			*(dst--) = *(src--);
			*(dst--) = *(src--);
			*(dst--) = ',';
			}

		}

	return s;
	}

//----------
//
// update_crc--
//	Incorporate the next byte into a cyclic redundancy check.
// update_crc_string--
//	Incorporate the next string into a cyclic redundancy check.
// update_crc_u32--
//	Incorporate the next 32-bit value into a cyclic redundancy check.
//
//----------
//
// Arguments:
//	u32 crc:	The crc of the previous bytes
//	u8	ch:		The byte to "add" to the crc.
//
// Returns:
//	The new crc value.
//
//----------
//
// Notes:
//	(1)	This code is based on some borrowed ages from
//		  http://remus.rutgers.edu/~rhoads/Code/crc-32b.c
//	(2)	To compute a crc over a string of bytes, do something like this:
//		  crc = 0
//		  for ch in string:
//		    crc = update_crc(crc,ch)
//
//----------

static int crcTableInitialized = false;
static u32 crcTable[256];


static void generate_crc_table (void)
	{
	u32	crc, poly;
	int	i, j;

	if (crcTableInitialized) return;
	crcTableInitialized = true;
	
	poly = 0xEDB88320L;
	for (i=0 ; i<256 ; i++)
		{
		crc = i;
		for (j=8 ; j>0 ; j--)
			{
			if (crc & 1) crc = (crc >> 1) ^ poly;
			        else crc >>= 1;
			}
		crcTable[i] = crc;
		}
	}


u32 update_crc (u32 crc, u8 ch)
	{
	if (! crcTableInitialized) generate_crc_table();
	return (crc>>8) ^ crcTable[(crc^ch) & 0xFF];
	}


u32 update_crc_string (u32 crc, char* s)
	{
	for ( ; *s!=0 ; s++)
		crc = update_crc(crc,(u8)*s);
	return crc;
	}


u32 update_crc_u32 (u32 crc, u32 val)
	{
	crc = update_crc(crc,(u8)(val>>24));
	crc = update_crc(crc,(u8)(val>>16));
	crc = update_crc(crc,(u8)(val>>8 ));
	crc = update_crc(crc,(u8)(val    ));
	return crc;
	}

