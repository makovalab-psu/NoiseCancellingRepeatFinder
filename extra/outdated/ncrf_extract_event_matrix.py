#!/usr/bin/env python
"""
Extract the (NON-positional) match-and-error event counts matrix from
Noise Cancelling Repeat Finder alignments.
"""

from sys        import argv,stdin,stdout,stderr,exit
from os         import path as os_path
from math       import log
from ncrf_parse import alignments,int_with_unit


def usage(s=None):
	message = """
usage: ncrf_cat <output_from_NCRF> | ncrf_extract_event_matrix [options]
  --withheader            include a header line in the output
  --sumonly               include only a summation line in the output
                          (by default, we output a separate line for each
                          alignment, and no sum)
  --head=<number>         limit the number of input alignments

The output matrix has R rows and 9 columns, where R is the number of input
alignments.

The first column is the line number of the alignment in the input file. The
second column is the motif, e.g. GGAAT.  The third column is the match ratio
("mRatio"). The remaining columns are, respectively, the counts for matches
("m"), mismatches ("mm"), insertion opens ("io"), insertion extensions
("ix"), deletion opens ("do"), and deletion extensions ("dx").

The output is intended to be suitable as input to R, and can be used as input
to infer_scoring."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	writeHeader = False
	writeWhat   = "per alignment"
	headLimit   = None
	requireEof  = True

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg in ["--withheader","--with=header","--with:header"]):
			writeHeader = True
		elif (arg in ["--sumonly","--sum=only","--sum:only"]):
			writeWhat = "sum only"
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg in ["--noendmark","--noeof","--nomark"]):   # (unadvertised)
			requireEof = False
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# process the alignments

	sum = {"m":0, "mm":0, "io":0, "ix":0, "do":0, "dx":0}

	alignmentNum = 0
	for a in alignments(stdin,requireEof):
		alignmentNum +=1 

		if (headLimit != None) and (alignmentNum > headLimit):
			print >>stderr, "limit of %d alignments reached" % headLimit
			break

		(nMatch,nMismatch,nInsO,nInsX,nDelO,nDelX) = extract_events(a)

		if (writeHeader):
			print "\t".join(["line","motif","mRatio","m","mm","io","ix","do","dx"])
			writeHeader = False

		if (writeWhat == "per alignment"):
			vec = [a.lineNumber,a.motif,a.mRatio,nMatch,nMismatch,nInsO,nInsX,nDelO,nDelX]
			print "\t".join(map(str,vec))

		sum["m"]  += nMatch
		sum["mm"] += nMismatch
		sum["io"] += nInsO
		sum["ix"] += nInsX
		sum["do"] += nDelO
		sum["dx"] += nDelX

	sum["events"] = (sum["m"] + sum["mm"] + sum["io"] + sum["ix"] + sum["do"] + sum["dx"])

	if (alignmentNum == 0):
		print >>stderr, "WARNING: input contained no alignments"
	elif (writeWhat == "sum only"):
		mRatio = float(sum["m"]) / sum["events"]
		mRatio = "%.3f" % mRatio
		vec = ["all",a.motif,mRatio,sum["m"],sum["mm"],sum["io"],sum["ix"],sum["do"],sum["dx"]]
		print "\t".join(map(str,vec))


def extract_events(a):
	if (len(a.seqText) != len(a.motifText)):
		exit(("%s: internal error, alignment text lengths aren't the same" \
		       + " for alignment at line %d")
		   % (os_path.basename(argv[0]),a.lineNumber))

	nMatch = nMismatch = nInsO = nInsX = nDelO = nDelX = 0

	prevEvent = None
	for (ix,(seqCh,motifCh)) in enumerate(zip(a.seqText,a.motifText)):
		if (seqCh == "-") and (motifCh == "-"):
			exit(("%s: internal error, gap aligned to gap in column %d"
			        + " of alignment at line %d")
			   % (os_path.basename(argv[0]),ix+1,a.lineNumber))

		if (seqCh == "-"):
			if (prevEvent != "d"):
				nDelO += 1
				prevEvent = "d"
			else:
				nDelX += 1
		elif (motifCh == "-"):
			if (prevEvent != "i"):
				nInsO += 1
				prevEvent = "i"
			else:
				nInsX += 1
		elif (seqCh == motifCh):
			nMatch += 1
			prevEvent = "m"
		else: # if (seqCh != motifCh):
			nMismatch += 1
			prevEvent = "mm"

	return (nMatch,nMismatch,nInsO,nInsX,nDelO,nDelX)


if __name__ == "__main__": main()
