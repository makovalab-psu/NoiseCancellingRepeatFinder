#!/usr/bin/env python
"""
Convert the output of Noise Cancelling Repeat Finder to a summary, a tab-
delimited table with one line of error info per alignment, primarily to show
the relative positions of errors along the alignment.
"""

from sys        import argv,stdin,stdout,stderr,exit
from os         import path as os_path
from ncrf_parse import alignments,parse_noise_rate


def usage(s=None):
	message = """
usage: cat <output_from_NCRF> | ncrf_error_positions [options]
  --minmratio=<ratio>  discard alignments with a low frequency of matches;
                       ratio can be between 0 and 1 (e.g. "0.85"), or can be
                       expressed as a percentage (e.g. "85%")
  --maxnoise=<ratio>   (same as --minmratio but with 1-ratio)

Typical output:
  #line seq             strand start end   querybp mRatio nErrors errors
  1     FAB41174_065680 -      1568  3021  1461    82.6%  261     0.023 0.024 0.025 ... 0.981 0.983
  11    FAB41174_029197 -      3908  5077  1189    82.4%  215     0.021 0.029 0.032 ... 0.989 0.996
  21    FAB41174_005950 -      2312  3334  1060    81.1%  205     0.023 0.027 0.036 ... 0.979 0.995
   ...

Note that lines will usually have different numbers of fields, since they have
different error counts.  "errors" is a vector of the relative positions of
errors along the alignment, with positions ranging from zero to one."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	minMRatio  = None
	requireEof = True

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--minmratio=")):
			minMRatio = parse_noise_rate(argVal)
			if (not (0.0 <= minMRatio <= 1.0)):
				exit("%s: mratio has to be between 0 and 1 (e.g. 0.85 or 85%%)\n%s"
				   % (os_path.basename(argv[0]),arg))
		elif (arg.startswith("--maxnoise=")):
			minMRatio = 1 - parse_noise_rate(argVal)
			if (not (0.0 <= minMRatio <= 1.0)):
				exit("%s: noise has to be between 0 and 1 (e.g. 0.15 or 15%%)\n%s"
				   % (os_path.basename(argv[0]),arg))
		elif (arg in ["--noendmark","--noeof","--nomark"]):   # (unadvertised)
			requireEof = False
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# process the alignments

	print "\t".join(["#line","seq","strand","start","end","querybp","mRatio","nErrors","errors"])

	for a in alignments(stdin,requireEof):
		if (minMRatio != None) and (a.mRatio < minMRatio):
			continue

		errorPositions = []
		for (ix,ch) in enumerate(a.errorText):
			if (ch == "x"):
				errorPositions += [float(ix)/(len(a.errorText)-1)]

		print "\t".join([str(a.lineNumber),
		                 a.seqName,a.strand,str(a.start),str(a.end),
		                 str(a.motifBaseCount),
		                 "%.1f%%" % (100*a.mRatio),
		                 str(len(errorPositions))]
		              + map(lambda x:"%.3f"%x,errorPositions))


if __name__ == "__main__": main()
