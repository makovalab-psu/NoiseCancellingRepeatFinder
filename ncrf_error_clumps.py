#!/usr/bin/env python
"""
Identify clumps of errors in Noise Cancelling Repeat Finder alignments.
"""

from sys        import argv,stdin,stdout,stderr,exit
from os         import path as os_path
from ncrf_parse import alignments,parse_noise_rate,int_with_unit


def usage(s=None):
	message = """
usage: cat <output_from_NCRF> | ncrf_error_clumps [options]
  --maxmratio=<ratio>     identify clumps with no more matches than given ratio; 
                          ratio can be between 0 and 1 (e.g. "0.85"), or can be
                          expressed as a percentage (e.g. "85%")
                          (default is 85%)
  --minnoise=<ratio>      (same as --maxmratio but with 1-ratio)
  --mincolumns=<columns>  identified clumps must have at least this many
                          alignment columns
                          (default is 10)
  --head=<number>         limit the number of input alignments
  --report:clumps         report clumps to stderr
                          (by default clumps are just marked on alignments)"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse the command line

	maxMRatio    = 0.85
	minColumns   = 10
	headLimit    = None
	reportClumps = False
	requireEof   = True
	debug        = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--maxMRatio=")):
			maxMRatio = parse_noise_rate(argVal)
			if (not (0.0 <= minMRatio <= 1.0)):
				exit("%s: mratio has to be between 0 and 1 (e.g. 0.85 or 85%%)\n%s"
				   % (os_path.basename(argv[0]),arg))
		elif (arg.startswith("--minnoise=")):
			maxMRatio = 1 - parse_noise_rate(argVal)
			if (not (0.0 <= minMRatio <= 1.0)):
				exit("%s: noise has to be between 0 and 1 (e.g. 0.15 or 15%%)\n%s"
				   % (os_path.basename(argv[0]),arg))
		elif (arg.startswith("--mincolumns=")) or (arg.startswith("--mindenom=")):
			minColumns =int(argVal)
			if (minColumns < 2):
				usage("minimum length has to be at least two columns\n%s" % arg)
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg == "--report:clumps") or (arg == "--report=clumps"):
			reportClumps = True
		elif (arg in ["--noendmark","--noeof","--nomark"]):   # (unadvertised)
			requireEof = False
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
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

		if (a.errorText == None):
			exit("%s: alignment at line %d doesn't include error text"
			   % (os_path.basename(argv[0]),a.lineNumber))

		if ("detail" in debug):
			print >>stderr, "\nlooking for clumps in %s %c %u-%u" \
			              % (a.seqName,a.strand,a.start,a.end)

		clumps = find_clumps(a.errorText,1-maxMRatio,minColumns,
		                     positiveCh='x',negativeCh='=')

		clumpText = ["-"] * len(a.errorText)
		for (start,end) in clumps:
			for ix in xrange(start,end): clumpText[ix] = "*"
		clumpText = "".join(clumpText)

		prefixLen = 1 + a.lines[0].find(" =")
		if (prefixLen < 0):
			prefixLen = 1 + a.lines[0].find(" x")

		if (alignmentNum > 1): print
		a.lines.insert(3,"# %-*s%s" % (prefixLen-2,"noise clumps",clumpText))
		print a

		if (reportClumps):
			for (start,end) in clumps:
				errorCount = matchCount = 0
				for ch in a.errorText[start:end]:
					if   (ch == 'x'): errorCount += 1
					elif (ch == '='): matchCount += 1
				print >>stderr, "line %d (%d,%d) m=%s x=%s mRatio: %.2f%%" \
							  % (a.lineNumber,
							     start,end,matchCount,errorCount,
								 (100.0*matchCount)/(matchCount+errorCount))

	if (requireEof):
		print "# ncrf end-of-file"


# find_clumps--
#	Find intervals in a string that have a density above some threshold.
#
# An interval with P positives and N negatives sites satifies the threshold R
# if P/(P+N) >= R.
#
#	   P/(P+N) >= R
#	   =>  P >= RP+RN
#	   =>  (1-R)P - RN >= 0
#
# So we can identify such regions by keeping a running sum of (1-R)P - RN.
# The latter is accomplished by adding 1-R to the sum whenever we encounter a
# positive and subtracting R whenever we encounter a negative.  Any interval
# for which the running sum is at least as high on the right as on the left
# satisfies the threshold.  Any such interval is a clump unless it is covered
# by a longer clump.
#
# Note that it *is* possible for two overlapping intervals to each be above
# the threshold while the combined interval is below the threshold.  This
# function returns such combined intervals, the rationale being that every
# position within the interval is part of *some* interval that's above the
# threshold.

def find_clumps(text,rate,minLength,positiveCh='+',negativeCh='-'):

	if (not 0 < rate < 1):
		raise ValueError, \
		      "%s is not a valid rate for find_clumps" % rate

	# perform a pre-scan to determine whether the sum in the algorithm would
	# be strictly decreasing over the entire vector;  in this case we know
	# there are no clump intervals and we can just quit

	allDecreasing = True
	for ch in text:
		if (ch == positiveCh):
			allDecreasing = False
			break

	if (allDecreasing): return []

	# search for clumps, intervals for which the average is above (or at) the
	# threshold;  this is equivalent to intervals in which the sum, minus the
	# length times the average, is positive (or zero)
	#
	# nota bene: this algorithm can suffer from round off error in the sums,
	#            which can make the reported intervals less precise than we'd
	#            like

	if ("detail" in debug):
		print >>stderr, "text: %s" % text
		print >>stderr, "rate: %.2f%%" % (100.0*rate)
		positiveCount = negativeCount = 0
		for ch in text:
			if   (ch == positiveCh): positiveCount += 1
			elif (ch == negativeCh): negativeCount += 1
		print >>stderr, "pos/neg: %s/%s ratio: %.2f%%" \
		              % (positiveCount,negativeCount,
		                 (100.0*positiveCount)/(positiveCount+negativeCount))

	minSum = valSum = 0.0

	minSums    = [valSum]   # minSums [0] = valSum
	minWhere   = [-1]       # minWhere[0] = -1
	numMinSums = 1
	minScan    = 0

	prevStart  = -1
	prevEnd    = -1
	clumps     = []

	for (ix,ch) in enumerate(text):
		# invariant: minSums[minScan] <= valSum < minSums[minScan-1]
		# with the implied assumption that minSums[-1] == infinity

		if   (ch == positiveCh): val = 1 - rate
		elif (ch == negativeCh): val = -rate
		else:
			raise ValueError, \
			      "\"%s\" is not a valid character (position %d in \"%s\")" \
			    % (ch,ix+1,text)

		valSum += val
		if (valSum < minSum):
			minSum     =  valSum
			minSums    += [valSum]   # minSums [numMinSums] = valSum
			minWhere   += [ix]       # minWhere[numMinSums] = ix
			numMinSums += 1

		# (re-establish the invariant)

		if (val < 0):								# the sum has decreased
			while (minSums[minScan] > valSum):
				minScan += 1
		elif (val > 0):								# the sum has increased
			while (minScan > 0) and (minSums[minScan-1] <= valSum):
				minScan -= 1

		if ("characters" in debug):
			print >>stderr, " [%u] %s %.2f %.2f %d.." \
			              % (ix,ch,val,valSum,minWhere[minScan])

		# minScan points at the earliest index with minSums[minScan] <= valSum,
		# so the interval (minWhere[minScan]+1 to ix) has sum >= 0

		minIx = minWhere[minScan]
		if (ix - minIx < minLength): continue

		start = minIx+1
		end   = ix

		if ("detail" in debug):
			print >>stderr, "      setting %u..%u %.2f %.2f %.2f" \
			              % (start,end+1,
			                 minSums[minScan],valSum,minSums[minScan]-valSum)

		if (prevStart == -1) or (start > prevEnd+1):
			# no overlap with previous interval (or no previous interval)
			if ("detail" in debug):
				print >>stderr, "      adding (%d,%d)" % (start,end+1)
			clumps += [(start,end+1)]
			prevStart = start
			prevEnd   = end
		elif (start >= prevStart):
			# new interval extends previous interval on right
			if ("detail" in debug):
				print >>stderr, "      replacing (%d,%d) with (%d,%d)" \
				              % (clumps[-1][0],clumps[-1][1],prevStart,end+1)
			clumps[-1] = (prevStart,end+1)
			prevEnd = end
		else:
			# new interval extends previous interval on left and right
			if ("detail" in debug):
				print >>stderr, "      replacing (%d,%d) with (%d,%d)" \
				              % (clumps[-1][0],clumps[-1][1],start,end+1)
			clumps[-1] = (start,end+1)
			prevStart = start
			prevEnd   = end

	# trim any negative ends off of the clumps

	if (not "notrim" in debug):
		for (ix,(start,end)) in enumerate(clumps):
			while (text[start] == negativeCh) and (start < end):
				start += 1
			while (text[end-1] == negativeCh) and (start < end):
				end -= 1
			clumps[ix] = (start,end)

	# collapse overlapping clumps

	if (not "nocollapse" in debug):
		collapsedClumps = []

		start = end = None
		for (s,e) in clumps:
			if (start == None):
				(start,end) = (s,e)
			elif (s < end):
				end = max(end,e)
			else:
				collapsedClumps += [(start,end)]
				(start,end) = (s,e)

		if (start != None):
			collapsedClumps += [(start,end)]

		clumps = collapsedClumps

	# report clumps to the caller

	if ("rate" in debug):
		for (start,end) in clumps:
			positiveCount = negativeCount = 0
			for ch in text[start:end]:
				if   (ch == positiveCh): positiveCount += 1
				elif (ch == negativeCh): negativeCount += 1
			print >>stderr, "(%d,%d) %s/%s ratio: %.2f%%" \
						  % (start,end,positiveCount,negativeCount,
							 (100.0*positiveCount)/(positiveCount+negativeCount))
	elif ("clumps" in debug):
		print >>stderr, "[%s]" % ",".join(["(%d,%d)"%(s,e) for (s,e) in clumps])

	return clumps


if __name__ == "__main__": main()
