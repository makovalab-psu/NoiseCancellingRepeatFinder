#!/usr/bin/env python
"""
Compare observed events from alignments to known truth in a genome.

References:
  [1] https://en.wikipedia.org/wiki/Sensitivity_and_specificity
"""

from sys        import argv,stdin,stdout,stderr,exit
from sets       import Set
from ncrf_parse import parse_noise_rate


def usage(s=None):
	message = """
usage: cat <alignment_summary> | observed_vs_genome_truth <truth_catalog> [options]
  <truth_catalog>        File containing aligment "truth" (see below)
  --motif=<motif>        (cumulative) Motifs represented in the summary. Truth
                         intervals for other motifs are discarded. If this is
                         not provided, we use all truth intervals.
  --detection=<portion>  threshold for a truth intervals to be considered as
                         "suffiently" covered; 0<portion<=1, but can be
                         expressed with a % sign
                         (default is 95%)
  --detail=<file>        Report separate true positive rates for each observed
                         interval

We'll consider each base in the genome as an item to be classified. The
alignment summary tells us the classification of each base -- each base in an
alignment interval is classified as a positive, and any other base is
classified as a negative.

The truth catalog tells us the correct classification of each base -- each base
in a truth interval should be classified as a positive (by an ideal
classifier), and any other base should be classified as a negative.

We report
  true postive rate         (TPR) = TP/(TP+FN)
  false negative rate       (FNR) = FN/(TP+FN)
  positive predictive value (PPV) = TP/(TP+FP)   a.k.a. precision
  false discovery rate      (FDR) = FP/(TP+FP)
  detection rate            fraction of truth intervals "suffiently" covered by
                            alignments

Since we don't know the genome size, we can't report anything involving true
negatives.

Note that we AREN'T considering whether the aligner called the correct event
(mm, ins, or del) for a base; only whether it covered that base with an aligned
interval.

The alignment summary is usually the output from ncrf_summary. It has 13
columns but only the aligned intervals and motif are used here ("seq", "start",
"end", and "motif", columns 3, 4, 5, and 2).

The truth catalog is usually the output from the --catalog option of
mock_motif_read. It has 12 columns but only the intervals and motif are used
here ("chrom", "start", "end", and "motif", columns 1, 2, 3, and 4). Intervals
must be distinct; overlaps are not allowed."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	truthFilename   = None
	detectionThresh = 0.95
	detailFilename  = None
	motifs          = None

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--truth=")) or (arg.startswith("--catalog=")):
			if (truthFilename != None):
				usage("unrecognized option: %s" % arg)
			truthFilename = argVal
		elif (arg.startswith("--motif=")):
			if (motifs == None): motifs = Set()
			motifs.add(argVal)
		elif (arg.startswith("--detection=")):
			detectionThresh = parse_noise_rate(argVal)
			if (detectionThresh <= 0):
				usage("detection threshold must be positive (%s)" % arg)
			if (detectionThresh > 1):
				usage("detection threshold cannot be more than 100% (%s)" % arg)
		elif (arg.startswith("--detail=")):
			detailFilename = argVal
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (truthFilename == None):
			truthFilename = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (truthFilename == None):
		usage("you have to give me the truth file")

	# collect the truth

	chromOrder = []
	chromSeen  = Set()

	truth = {}

	f = file(truthFilename,"rt")
	for (chrom,start,end,motif) in read_intervals(f,(0,1,2,3),truthFilename):
		if (motifs != None):
			motif = motif.split(".")[0]
			if (motif not in motifs): continue
		if (chrom not in chromSeen):
			chromOrder += [chrom]
			chromSeen.add(chrom)
		if (chrom not in truth): truth[chrom] = []
		truth[chrom] += [(start,end)]
	f.close()

	for chrom in truth:
		truth[chrom].sort()

	for chrom in truth:
		(prevStart,prevEnd) = truth[chrom][0]
		for (s,e) in truth[chrom][1:]:
			assert (start >= prevEnd), \
			       "truth interval %d..%d overlaps %d..%d" \
			     % (prevStart,prevEnd,s,e)
			(prevStart,prevEnd) = (s,e)

	# collect the observations

	observed = {}

	for (chrom,start,end,motif) in read_intervals(stdin,(2,3,4,1)):
		if (motifs != None):
			if (motif not in motifs): continue
		if (chrom not in chromSeen):
			chromOrder += [chrom]
			chromSeen.add(chrom)
		if (chrom not in observed): observed[chrom] = []
		observed[chrom] += [(start,end)]

	for chrom in observed:
		observed[chrom] = merge_intervals(observed[chrom])

	# compute true positives and false positives for each observed interval
	#
	# See reference [1].

	pTotal = tpTotal = fpTotal = 0
	for chrom in chromOrder:
		if (chrom not in observed): continue

		truthOnChrom = None
		if (chrom in truth):
			truthOnChrom = truth[chrom]

		for (start,end) in observed[chrom]:
			length = end-start
			tp = 0
			if (truthOnChrom != None):
				tp = overlap_count(start,end,truthOnChrom)
			fp = length - tp
			pTotal  += length
			tpTotal += tp
			fpTotal += fp

	# compute false negatives, and per-observed-interval true positives

	perInterval = []

	fnTotal = truthTotal = 0
	detected = 0
	for chrom in chromOrder:
		if (chrom not in truth): continue

		observedOnChrom = None
		if (chrom in observed):
			observedOnChrom = observed[chrom]

		for (start,end) in truth[chrom]:
			length = end-start
			tp = 0
			if (observedOnChrom != None):
				tp = overlap_count(start,end,observedOnChrom)
			fn = length - tp
			fnTotal     += fn
			truthTotal  += length
			perInterval += [(chrom,start,end,tp)]
			if (tp >= length*detectionThresh):  # (tp/length >= detectionThresh)
				detected += 1

	# report

	if (tpTotal+fnTotal == 0):
		print "%s\t%s/%s\tNA" \
		    % ("TPR",tpTotal,tpTotal+fnTotal)
		print "%s\t%s/%s\tNA" \
		    % ("FNR",fnTotal,tpTotal+fnTotal)
	else:
		print "%s\t%s/%s\t%5.3f%%" \
		    % ("TPR",tpTotal,tpTotal+fnTotal,100.0*tpTotal/(tpTotal+fnTotal))
		print "%s\t%s/%s\t%5.3f%%" \
		    % ("FNR",fnTotal,tpTotal+fnTotal,100.0*fnTotal/(tpTotal+fnTotal))

	if (tpTotal+fpTotal == 0):
		print "%s\t%s/%s\tNA" \
		    % ("PPV",tpTotal,tpTotal+fpTotal)
		print "%s\t%s/%s\tNA" \
		    % ("FDR",fpTotal,tpTotal+fpTotal)
	else:
		print "%s\t%s/%s\t%5.3f%%" \
		    % ("PPV",tpTotal,tpTotal+fpTotal,100.0*tpTotal/(tpTotal+fpTotal))
		print "%s\t%s/%s\t%5.3f%%" \
		    % ("FDR",fpTotal,tpTotal+fpTotal,100.0*fpTotal/(tpTotal+fpTotal))

	if (len(perInterval) == 0):
		print "%s\t%s/%s\tNA" \
		    % ("DETECTED",detected,len(perInterval))
	else:
		print "%s\t%s/%s\t%5.3f%%" \
		    % ("DETECTED",detected,len(perInterval),100.0*detected/len(perInterval))

	if (detailFilename != None):
		f = file(detailFilename,"wt")
		print >>f, "#%s\t%s\t%s\t%s/%s\t%s" \
			     % ("chrom","start","end","TP","TP+FN","TPR")
		for (chrom,start,end,tp) in perInterval:
			length = end-start
			print >>f, "%s\t%s\t%s\t%s/%s\t%5.1f%%" \
				     % (chrom,start,end,tp,length,100.0*tp/length)
		f.close()


# read_intervals--
#	Read intervals from a file, one per line

def read_intervals(f,columns,filename=None):
	if (filename == None): filename = "input"

	if (len(columns) == 3):
		(chromCol,startCol,endCol) = columns
		motifCol = None
	else:
		(chromCol,startCol,endCol,motifCol) = columns

	minColumns = 1+max(columns)
	numColumns = None

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		if (numColumns == None):
			assert (len(fields) >= minColumns), \
			       "not enough fields at line %d in %s (%d, expected at least %d)" \
			     % (lineNumber,filename,len(fields),minColumns)
			numColumns = len(fields)
		else:
			assert (len(fields) == numColumns), \
			       "inconsistent number of fields at line %d in %s (%d, expected %d)" \
			     % (lineNumber,filename,len(fields),numColumns)

		try:
			chrom = fields[chromCol]
			start = int(fields[startCol])
			end   = int(fields[endCol])
			if (end < start): raise ValueError
			if (motifCol != None): motif = fields[motifCol]
		except ValueError:
			assert (False), "bad line (%d in %s): %s" % (lineNumber,filename,line)

		if (motifCol == None):
			yield (chrom,start,end)
		else:
			yield (chrom,start,end,motif)


# merge_intervals--
#   Sort intervals and merge overlaps

def merge_intervals(intervals):
	intervals = list(intervals)
	intervals.sort()

	newIntervals = []

	start = None
	for (s,e) in intervals:
		if (start == None):
			(start,end) = (s,e)
		elif (s > end):
			newIntervals += [(start,end)]
			(start,end) = (s,e)
		elif (e > end):
			end = e

	if (start != None):
		newIntervals += [(start,end)]

	return newIntervals


# overlap_count--
#   Report the number of and interval's bases that overlap a list of intervals
#
# We assume the intervals is sorted and contains no overlaps.

def overlap_count(start,end,intervals):
	baseCount = 0

	for (s,e) in intervals:
		if (e <= start): continue
		if (s >= end):   break

		s = max(s,start)
		e = min(e,end)
		baseCount += e-s

	return baseCount


if __name__ == "__main__": main()
