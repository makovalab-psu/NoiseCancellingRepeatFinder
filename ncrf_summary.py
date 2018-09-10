#!/usr/bin/env python
"""
Convert the output of Noise Cancelling Repeat Finder to a summary, a
tab-delimited table with one line of stats per alignment.
"""

from sys        import argv,stdin,stdout,stderr,exit
from os         import path as os_path
from ncrf_parse import alignments,parse_noise_rate


def usage(s=None):
	message = """
usage: cat <output_from_NCRF> | ncrf_summary [options]
  --minmratio=<ratio>  discard alignments with a low frequency of matches;
                       ratio can be between 0 and 1 (e.g. "0.85"), or can be
                       expressed as a percentage (e.g. "85%")
  --maxnoise=<ratio>   (same as --minmratio but with 1-ratio)

Typical output:
  #line motif seq        start end  strand seqLen querybp mRatio m    mm  i  d
  1     GGAAT FAB41174_6 1568  3021 -      3352   1461    82.6%  1242 169 42 50
  11    GGAAT FAB41174_2 3908  5077 -      7347   1189    82.4%  1009 125 35 55
  21    GGAAT FAB41174_0 2312  3334 -      4223   1060    81.1%  881  115 26 64
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
		elif (arg in ["--noendmark","--noeof","--nomark"]):   # (unadvertised)
			requireEof = False
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# process the alignments

	print "\t".join(["#line","motif",
	                 "seq","start","end","strand",
	                 "seqLen","querybp",
	                 "mRatio","m","mm","i","d"])

	userHasntBeenWarned = True
	for a in alignments(stdin,requireEof):
		seqLenStr      = "NA"
		mRatioStr      = "NA"
		nMatchStr      = "NA"
		nMismatchStr   = "NA"
		nInsertionsStr = "NA"
		nDeletionsStr  = "NA"

		if (a.seqLen != None):
			seqLenStr = str(a.seqLen)

		if (hasattr(a,"mRatio")):
			if (minMRatio != None) and (a.mRatio < minMRatio):
				continue
			mRatioStr = "%.1f%%" % (100*a.mRatio)

		if (hasattr(a,"nMatch")):
			nMatchStr = str(a.nMatch)

		if (hasattr(a,"nMismatch")):
			nMismatchStr = str(a.nMismatch)

		if (hasattr(a,"nInsertions")):
			nInsertionsStr = str(a.nInsertions)

		if (hasattr(a,"nDeletions")):
			nDeletionsStr = str(a.nDeletions)

		if (mRatioStr == "NA"):
			if (userHasntBeenWarned):
				print >>stderr, \
                     ("WARNING: input alignments did not contain an event stats line"
                    + "\n         (NCRF --stats=events would create that line)")
				userHasntBeenWarned = False

		print "\t".join([str(a.lineNumber),
		                 a.motif,
		                 a.seqName,str(a.start),str(a.end),a.strand,
		                 seqLenStr,str(a.motifBaseCount),
		                 mRatioStr,
		                 nMatchStr,nMismatchStr,nInsertionsStr,nDeletionsStr])


if __name__ == "__main__": main()
