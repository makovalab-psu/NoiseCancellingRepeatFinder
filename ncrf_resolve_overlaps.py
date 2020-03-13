#!/usr/bin/env python
"""
Resolve overlapping alignments of different motifs.
"""

from sys        import argv,stdin,stdout,stderr,exit
from gzip       import open as gzip_open
from ncrf_parse import int_with_unit


def usage(s=None):
	message = """
usage: ncrf_resolve_overlaps <alignment_summary..> [options]
  <alignment_summary>    (cumulative) file(s) containing aligment summaries
                         for which overlaps are to be resolved
  --head=<number>        limit the number of input aligment summaries
  --out=<name_template>  file to write overlap groups to; see discussion of
                         name template below; if this option is absent, all
                         output is written to the console

The name template either names a single file or a collection of files. See
below for some examples.

The input alignment summaries are usually the output from ncrf_summary. Any
input file may contain alignments for more than one motif.

A typical input file is shown below. However, we do not interpret any columns
other than motif, seq, start, and end. This allows, for example, the output
from ncrf_summary_with_consensus.

  #line motif seq        start end  strand seqLen querybp mRatio m    mm  i  d
  1     GGAAT FAB41174_6 1568  3021 -      3352   1461    82.6%  1242 169 42 50
  11    GGAAT FAB41174_2 3908  5077 -      7347   1189    82.4%  1009 125 35 55
  21    GGAAT FAB41174_0 2312  3334 -      4223   1060    81.1%  881  115 26 64
   ...

If the output name template includes the substring "{motif}", this substring is
replaced by a motif name and any un-overlapped alignments to that motif are
written to that file. If the template name doesn't include "{motif}", all
un-overlapped alignments and overlapping groups are written to one file. Note
that "{motif}" is the word "motif" surrounded by two curly brackets.

Overlapping groups are either written to the console (if no name template is
given), to the same file with alignments (if the name template doesn't contain
"{motif}"), or to a a file separate from the alignments (with "{motif}"
replaced by "overlaps").

This is summarized in the table below. We assume for this that the input only
contains two motifs, GGAAT and CATATA.

  name_template    | output
  -----------------+----------------------------------------------------------
  (none)           | un-overlapped and overlap groups written to the console
  -----------------+----------------------------------------------------------
  filename         | un-overlapped and overlap groups written to filename
  -----------------+----------------------------------------------------------
  filename.{motif} | un-overlapped GGAAT written to filename.GGAAT
                   | un-overlapped CATATA written to filename.CATATA
                   | overlap groups written to filename.overlaps
  -----------------+----------------------------------------------------------

Overlap groups are separated by a single blank line, as shown below (note that
this is a contrived example). When un-overlapped alignments and overlapped ones
are in the same file, the un-overlapped ones are first, with a blank line
separating each alignment.

  #line motif  seq        start end  strand seqLen querybp mRatio m    mm  i  d
  1     GGAAT  FAB41174_6 1568  3021 -      3352   1461    82.6%  1242 169 42 50
  1     CATATA FAB41174_6 1621  2607 -      3352    ...
  (blank line)
  21    GGAAT  FAB41174_0 2312  3334 -      4223   1060    81.1%  881  115 26 64
  41    CATATA FAB41174_0 3276  4098 -      4223    ...
  31    GGAAT  FAB41174_0 3966  4271 +      4223    ...
  (blank line)
   ... (more groups)"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global summaryHeaderLine
	global debug

	summaryHeaderLine = None

	# parse the command line

	inputFilenames = []
	outTemplate    = None
	headLimit      = None
	debug          = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--out=")):
			outTemplate = argVal
		elif (arg.startswith("--head=")):
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

	writeSingletonsSeparately = (outTemplate != None) and ("{motif}" in outTemplate)

	# collect the alignments

	seqToSummaries = {}
	seqOrder = []

	summaryNum = 0
	for filename in inputFilenames:
		if (filename.endswith(".gz")) or (filename.endswith(".gzip")):
			f = gzip_open(filename,"rt")
		else:
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

	# partition the alignments into overlapping groups

	seqToGroups = {}

	for seq in seqOrder:
		seqToGroups[seq] = overlapping_groups(seqToSummaries[seq])

	if ("groups" in debug):
		for seq in seqOrder:
			for group in seqToGroups[seq]:
				print >>stderr, "==="
				for summary in group:
					print >>stderr,summary.line

	# collect groups by motif subset

	subsetToGroups = {}
	singletons = set()

	for seq in seqOrder:
		for group in seqToGroups[seq]:
			subset = set([summary.motif for summary in group])
			subset = list(subset)
			subset.sort()
			subset = tuple(subset)

			if (subset not in subsetToGroups):
				subsetToGroups[subset] = [group]
			else:
				subsetToGroups[subset] += [group]

			if (len(subset) == 1):
				singletons.add(subset[0])

	# if we're to report un-overlapped alignments separately, do so now (and
	# remove them from the groups)

	singletons = list(singletons)
	singletons.sort()

	if (writeSingletonsSeparately):
		for motif in singletons:
			subset = (motif,)

			motifFilename = outTemplate.replace("{motif}",motif)
			motifF = file(motifFilename,"wt")
			print >>stderr, "writing to \"%s\"" % motifFilename

			if (summaryHeaderLine != None):
				print >>motifF,summaryHeaderLine
			for group in subsetToGroups[subset]:
				for summary in group:
					print >>motifF,summary.line

			del subsetToGroups[subset]

	# report overlapping alignment groups (and un-overlapped groups if we
	# didn't report them already)

	if (outTemplate == None):
		outF = stdout
		if (list(subsetToGroups) == []):
			print >>stderr, "no alignments to write to console"
	elif ("{motif}" not in outTemplate):
		outF = file(outTemplate,"wt")
		if (list(subsetToGroups) == []):
			print >>stderr, "no alignments to write to \"%s\"" % outTemplate
		else:
			print >>stderr,"writing to \"%s\"" % outTemplate
	else:
		outFilename = outTemplate.replace("{motif}","overlaps")
		outF = file(outFilename,"wt")
		if (list(subsetToGroups) == []):
			print >>stderr, "no alignments to write to \"%s\"" % outFilename
		else:
			print >>stderr,"writing to \"%s\"" % outFilename

	motifCountToSubsets = {}
	for subset in subsetToGroups:
		motifCount = len(subset)
		if (motifCount not in motifCountToSubsets):
			motifCountToSubsets[motifCount] = [subset]
		else:
			motifCountToSubsets[motifCount] += [subset]

	motifCounts = list(motifCountToSubsets)
	motifCounts.sort()

	isFirstGroup = True
	for motifCount in motifCounts:
		subsets = motifCountToSubsets[motifCount]
		subsets.sort()
		for subset in subsets:
			for group in subsetToGroups[subset]:
				if (isFirstGroup):
					if (summaryHeaderLine != None):
						print >>outF,summaryHeaderLine
					isFirstGroup = False
				else:
					print >>outF  # (line to separate groups)

				for summary in group:
					print >>outF,summary.line

	if (outF != stdout):
		outF.close()


# read_summary--
#	Read alignments from an NCRF alignment summary file, one per line

class Summary: pass

def read_summary(f,filename=None):
	global summaryHeaderLine

	if (filename == None): filename = "input"

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")):
			if (summaryHeaderLine == None): summaryHeaderLine = line
			continue

		fields = line.split()
		if (len(fields) < 5):
			exit("%s: not enough columns at line %d (%d, expected at least %d)\n%s"
			   % (os_path.basename(argv[0]),lineNumber,len(fields),5,line))

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
			#s.strand  = fields[5]
			#s.seqLen  = int(fields[6])
			#s.querybp = int(fields[7])
			#s.mRatio  = fields[8]
			#s.m       = int(fields[9])
			#s.mm      = int(fields[10])
			#s.i       = int(fields[11])
			#s.d       = int(fields[12])
		except ValueError:
			assert (False), "bad alignment summary line (%d in %s): %s" \
			              % (lineNumber,filename,line)

		try:
			if (s.end < s.start): raise ValueError
		except ValueError:
			assert (False), "bad alignment summary line (bad interval at %d in %s): %s" \
			              % (lineNumber,filename,line)

		#try:
		#	if (s.strand not in ["+","-"]): raise ValueError
		#except ValueError:
		#	assert (False), "bad alignment summary line (bad strand at %d in %s): %s" \
		#	              % (lineNumber,filename,line)

		#try:
		#	if (not s.mRatio.endswith("%")): raise ValueError
		#	s.mRatio = float(s.mRatio[:-1]) / 100
		#	if (not 0 <= s.mRatio <= 1): raise ValueError
		#except ValueError:
		#	assert (False), "bad alignment summary line (bad mRatio at %d in %s): %s" \
		#	              % (lineNumber,filename,line)

		yield s


# overlapping_groups--
#   Combine alignment summaries into overlapping groups

def overlapping_groups(summaries):
	summaries = [(summary.start,summary.end,summary) for summary in summaries]
	summaries.sort()

	groups = []

	currentGroup = None
	start = end  = None
	for (s,e,summary) in summaries:
		if (start == None):
			currentGroup = [summary]
			(start,end) = (s,e)
		elif (s >= end):
			groups += [currentGroup]
			currentGroup = [summary]
			(start,end) = (s,e)
		elif (e <= end):
			currentGroup += [summary]
		else: # if (e > end):
			currentGroup += [summary]
			end = e

	if (currentGroup != None):
		groups += [currentGroup]

	return groups


if __name__ == "__main__": main()
