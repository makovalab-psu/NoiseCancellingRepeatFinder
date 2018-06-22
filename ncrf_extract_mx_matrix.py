#!/usr/bin/env python
"""
Extract the positional match-and-error-counts matrix from Noise Cancelling
Repeat Finder alignments.
"""

from sys        import argv,stdin,stdout,stderr,exit
from ncrf_parse import alignments,int_with_unit


def usage(s=None):
	message = """
usage: cat <output_from_NCRF> | ncrf_extract_mx_matrix [options]
  --head=<number>         limit the number of input alignments

The output matrix has R rows and 2M+1 columns, where R is the number of input
alignments and M is the length of the aligned motif. (It is assumed that all
alignments are to the same motif, but this is NOT validated).

The first column is the line number of the alignment in the input file.  The
next M columns are the positional counts for matches ("m"), and the final M
columns are the positional counts for errors ("x").  The output is intended
to be suitable as input to R."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	headLimit  = None
	requireEof = True

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg in ["--noendmark]","--noeof","--nomark"]):   # (unadvertised)
			requireEof = False
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# process the alignments

	alignmentNum = 0
	for a in alignments(stdin,requireEof):
		alignmentNum +=1 

		if (headLimit != None) and (alignmentNum > headLimit):
			print >>stderr, "limit of %d alignments reached" % headLimit
			break

		positionalStats = a.positional_stats()

		numPositions = len(positionalStats)
		vec = [None] * (2*numPositions+1)
		vec[0] = a.lineNumber

		for (pos,stats) in enumerate(positionalStats):
			if ("m" not in stats):
				raise ValueError, \
				      "\"m\" missing from positional information for alignment at line %d" \
				    % a.lineNumber
			if ("x" not in stats):
				raise ValueError, \
				      "\"x\" missing from positional information for alignment at line %d" \
				    % a.lineNumber

			vec[1+pos]              = stats["m"]
			vec[1+numPositions+pos] = stats["x"]

		print "\t".join(map(str,vec))


if __name__ == "__main__": main()
