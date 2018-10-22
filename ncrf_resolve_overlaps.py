#!/usr/bin/env python
"""
Resolve overlapping alignments of different motifs.
"""

from sys        import argv,stdin,stdout,stderr,exit
from ncrf_parse import int_with_unit


def usage(s=None):
	message = """
usage: ncrf_resolve_overlaps <alignment_summary..> [options]
  <alignment_summary>  (cululative) File(s) containing aligment summaries
                       for which overlaps are to be resolved
  --head=<number>      limit the number of input aligment summaries

The alignment summaries are usually the output from ncrf_summary. Any file may
contain alignments for more than one motif.

Typical input file:
  #line motif seq        start end  strand seqLen querybp mRatio m    mm  i  d
  1     GGAAT FAB41174_6 1568  3021 -      3352   1461    82.6%  1242 169 42 50
  11    GGAAT FAB41174_2 3908  5077 -      7347   1189    82.4%  1009 125 35 55
  21    GGAAT FAB41174_0 2312  3334 -      4223   1060    81.1%  881  115 26 64
   ..."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse the command line

	inputFilenames = []
	headLimit      = None
	debug          = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			if (arg not in inputFilenames): inputFilenames += [arg]

	if (inputFilenames == []):
		usage("you have to give me at least one summary file")

	# collect the alignments

	seqToSummaries = {}
	seqOrder = []

	summaryNum = 0
	for filename in inputFilenames:
		f = file(filename,"rt")
		for summary in read_summary(f,filename):
			summaryNum += 1 

			if (headLimit != None) and (summaryNum > headLimit):
				print >>stderr, "limit of %d summaries reached" % headLimit
				break

			if (summary.seq not in seqToSummaries):
				seqOrder += [summary.seq]
				seqToSummaries[summary.seq] = []
			seqToSummaries[summary.seq] += [summary]

		f.close()

	# partition the alignments into overlapping clumps

	seqToClumps = {}

	for seq in seqToSummaries:
		seqToClumps[seq] = overlapping_clumps(seqToSummaries[seq])

	if ("clumps" in debug):
		for seq in seqToSummaries:
			for clump in seqToClumps[seq]:
				print >>stderr, "==="
				for summary in clump:
					print >>stderr,summary.line


# read_summary--
#	Read alignments from an NCRF alignment summary file, one per line

class Summary: pass

def read_summary(f,filename=None):
	if (filename == None): filename = "input"

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		if (len(fields) != 13):
			exit("%s: wrong number of columns at line %d (%d, expected %d)\n%s"
			   % (os_path.basename(argv[0]),lineNumber,len(fields),13,line))

		s = Summary()
		s.filename   = filename
		s.lineNumber = lineNumber

		try:
			s.line    = line
			# (fields[0] intentionally ignored)
			s.motif   = fields[1]
			s.seq     = fields[2]
			s.start   = int(fields[3])
			s.end     = int(fields[4])
			s.strand  = fields[5]
			s.seqLen  = int(fields[6])
			s.querybp = int(fields[7])
			s.mRatio  = fields[8]
			s.m       = int(fields[9])
			s.mm      = int(fields[10])
			s.i       = int(fields[11])
			s.d       = int(fields[12])
		except ValueError:
			assert (False), "bad alignment summary line (%d in %s): %s" \
			              % (lineNumber,filename,line)

		try:
			if (s.end < s.start): raise ValueError
		except ValueError:
			assert (False), "bad alignment summary line (bad interval at %d in %s): %s" \
			              % (lineNumber,filename,line)

		try:
			if (s.strand not in ["+","-"]): raise ValueError
		except ValueError:
			assert (False), "bad alignment summary line (bad strand at %d in %s): %s" \
			              % (lineNumber,filename,line)

		try:
			if (not s.mRatio.endswith("%")): raise ValueError
			s.mRatio = float(s.mRatio[:-1]) / 100
			if (not 0 <= s.mRatio <= 1): raise ValueError
		except ValueError:
			assert (False), "bad alignment summary line (bad mRatio at %d in %s): %s" \
			              % (lineNumber,filename,line)

		yield s


# overlapping_clumps--
#   Combine alignment summaries into overlapping clumps

def overlapping_clumps(summaries):
	summaries = [(summary.start,summary.end,summary) for summary in summaries]
	summaries.sort()

	clumps = []

	currentClump = None
	start = end  = None
	for (s,e,summary) in summaries:
		if (start == None):
			currentClump = [summary]
			(start,end) = (s,e)
		elif (s >= end):
			clumps += [currentClump]
			currentClump = [summary]
			(start,end) = (s,e)
		elif (e <= end):
			currentClump += [summary]
		else: # if (e > end):
			currentClump += [summary]
			end = e

	if (currentClump != None):
		clumps += [currentClump]

	return clumps


if __name__ == "__main__": main()
