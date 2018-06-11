#!/usr/bin/env python
"""
Convert the output of Noise Cancelling Repeat Finder to bed format.
"""

from sys        import argv,stdin,stdout,stderr,exit
from os         import path as os_path
from ncrf_parse import alignments,parse_noise_rate


def usage(s=None):
	message = """
usage: cat <output_from_NCRF> | ncrf_to_bed [options]
  --minmratio=<ratio>  discard alignments with a low frequency of matches;
                       ratio can be between 0 and 1 (e.g. "0.85"), or can be
                       expressed as a percentage (e.g. "85%")
  --maxnoise=<ratio>   (same as --minmratio but with 1-ratio)

Typical output is shown below.  The 6th column ("score" in the bed spec) is
the match ratio times 1000 (e.g. 826 is 82.6%).
  FAB41174_065680 1568 3021 . - 826
  FAB41174_029197 3908 5077 . - 824
  FAB41174_005950 2312 3334 . - 811
   ..."""

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
		elif (arg in ["--noendmark]","--noeof","--nomark"]):   # (unadvertised)
			requireEof = False
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# process the alignments

	for a in alignments(stdin,requireEof):
		if (minMRatio != None) and (a.mRatio < minMRatio):
			continue

		print "\t".join([a.seqName,str(a.start),str(a.end),
		                 ".",
		                 "%d" % (1000*a.mRatio),
		                 a.strand])


if __name__ == "__main__": main()
