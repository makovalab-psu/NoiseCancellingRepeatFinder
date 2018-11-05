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
usage: ncrf_cat <output_from_NCRF> | ncrf_summary_with_consensus [options]
  --minmratio=<ratio>  discard alignments with a low frequency of matches;
                       ratio can be between 0 and 1 (e.g. "0.85"), or can be
                       expressed as a percentage (e.g. "85%")
  --maxnoise=<ratio>   (same as --minmratio but with 1-ratio)

The input alignments must include consensus information. This can be
accomplished by ncrf_consensus_filter's --consensusonly option.

The output is similar to that of ncrf_summary, but the error stats are replaced
by a "consensus" column.

Typical output:
  #line motif seq        start end  strand seqLen querybp consensus
  1     GGAAT FAB41174_6 1568  3021 -      3352   1461    GGAAT
  11    GGAAT FAB41174_2 3908  5077 -      7347   1189    GGAT
  21    GGAAT FAB41174_0 2312  3334 -      4223   1060    GGAAT
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
	                 "consensus"])

	userHasntBeenWarned = True
	for a in alignments(stdin,requireEof):
		seqLenStr = "NA"
		if (a.seqLen != None):
			seqLenStr = str(a.seqLen)

		if (hasattr(a,"mRatio")):
			if (minMRatio != None) and (a.mRatio < minMRatio):
				continue

		consensuses = []
		for line in a.lines:
			if (not line.startswith("# consensus")): continue
			line = line.split(None,2)
			consensuses += [line[2]]

		if (consensuses == []):
			if (userHasntBeenWarned):
				print >>stderr, \
                     ("WARNING: input alignments did not contain a consensus line"
                    + "\n         (ncrf_consensus_filter would create that line)")
				userHasntBeenWarned = False
			consensus = "(missing)"
		else:
			consensus = ",".join(consensuses)

		print "\t".join([str(a.lineNumber),
		                 a.motif,
		                 a.seqName,str(a.start),str(a.end),a.strand,
		                 seqLenStr,str(a.motifBaseCount),
		                 consensus])


if __name__ == "__main__": main()
