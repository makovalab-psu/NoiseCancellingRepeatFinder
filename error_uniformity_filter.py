#!/usr/bin/env python
"""
Filter Noise Cancelling Repeat Finder alignments, discarding alignments that
don't (seem to) have uniformly distributed matches across the motif unit.
"""

# TODO: pass data to R via a temporary file, rather than on the Rscript
#       command line

import subprocess
from sys        import argv,stdin,stdout,stderr,exit
from os         import path as os_path
from random     import seed as random_seed,sample as random_sample,randint
from ncrf_parse import alignments,parse_probability,int_with_unit

defaultPrngSeed = "NCRF.error_uniformity_filter"


# Implementation Notes:
#
# [1] Justification for "matches-insertions" (as opposed to just "matches").
#     When NCRF counts matches, a deletion (a base deleted from the motif)
#     outcomes in a decrease of the count of matches at that position.  There is
#     no similar effect for an insertion (a base inserted into the motif).
#     Thus if we base our test solely on matches, insertions will not
#     contribute to the test.  Subtracting insertion counts from match counts
#     gives insertions the same effect as deletions.
#
# [2] Here we describe the format of the file produced by the unadvertised
#     option --report:matrix.  The reason this option isn't advertised is that
#     is should only be useful to developers testing implementation details of
#     the underlying chi-squared test.
#
#     The file has NO header line and one row for each alignment;  2M+2 columns
#     per row, where M is the motif length.  Column 1 is the line number of the
#     alignment in the incoming file.  Column 2 is the outcome of the test. 
#     Columns 3 thru M+2 are positional match counts for that alignment, and
#     columns M+3 thru 2M+2 are the positional error counts. We assume all rows
#     have the same number of columns, i.e. that the same motif length is
#     represented in all rows.

def usage(s=None):
	message = """
usage: cat <output_from_NCRF> | error_uniformity_filter [options]
  --method=chi-squared  judge alignments by a chi-squared test
                        (this is the default)
  --method=min-max      judge alignments by a "min-max" test  
  --effectsize=<value>  effect size for chi-squared test
                        (default is 0.3)
  --power=<probability> "power of test" for chi-squared test, 1 minus Type II
                        error probability
                        (default is 0.8)
  --trials=<number>     number of trials for min-max-test
                        (default is 100)
  --discard:good        discard the "good" alignments instead of the "bad" ones
                        (by default we discard the "bad" alignments) 
  --discard:none        don't discard any alignments, just report the test
                        outcomes (see below for how they are reported) 
  --test:matches        perform the test using match counts
                        (this is the default)
  --test:errors         perform the test using error counts; this may increase
                        the likelihood that the test cannot be performed for
                        some alignments
  --test:matches-insertions perform the test using match counts less insertions
  --warn:untested       report each untestable alignment as a warning
  --seed=<string>       set random seed; this is only applicable to the min-max
                        test; two special cases are "default" (use the built-in
                        default seed) and "none" (don't seed the random number
                        generator)
  --head=<number>       limit the number of input alignments
  --batch=<number>      number of input alignments processed by each call to R;
                        our ability to call R fails if the command line we pass
                        it is too long
                        (default is 30)

In a "true" alignment to a given motif unit, we expect the errors to be
distributed randomly and uniformly among the positions in the unit.  (That is
an underlying assumption, but might not be true.)  This program discards
alignments that fail a statistical test based on that assumption.

Since error counts may be too small for the statistical test, we use match
counts instead.

The input alignments must include position event information.  This can be
accomplished by using the --positionalevents option of Noise Cancelling Repeat
Finder.

When test outcomes are reported for --discard:none, one of the following lines
is added to each alignment:
# positional chi-squared: match uniformity rejected
# positional chi-squared: match uniformity not rejected
# positional chi-squared: untested
The "untested" case indicates that the chi-squared test could not be
performed, usually because one of the positional match counts is too small."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global batchSize
	global debug

	# parse the command line

	testMethod     = "chi-squared"
	effectSize     = 0.3
	power          = 0.8
	numTrials      = 100
	discardWhich   = "bad"
	testWhich      = "matches"
	warnOnUntested = False
	headLimit      = None
	batchSize      = 30
	reportAs       = "ncrf"
	requireEof     = True
	prngSeed       = defaultPrngSeed
	debug          = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg in ["--method=chi-squared","--method=chi-square"]):
			testMethod = "chi-squared"
		elif (arg == "--method=min-max"):
			testMethod = "min-max"
		elif (arg.startswith("--effectsize=")):
			effectSize = parse_probability(argVal)
		elif (arg.startswith("--power=")):
			power = parse_probability(argVal)
		elif (arg.startswith("--trials=")):
			numTrials = int(argVal)
		elif (arg in ["--discard:bad","--discard=bad"]):
			discardWhich = "bad"
		elif (arg in ["--discard:good","--discard=good"]):
			discardWhich = "good"
		elif (arg in ["--discard:none","--discard=none"]):
			discardWhich = "none"
		elif (arg in ["--test:matches","--test=matches"]):
			testWhich = "matches"
		elif (arg in ["--test:errors","--test=errors"]):
			testWhich = "errors"
		elif (arg in ["--test:matches-insertions","--test=matches-insertions",
		              "--test:m-i","--test=m-i"]):
			testWhich = "matches-insertions"
		elif (arg == "--warn:untested") or (arg == "--warn=matrix"):
			warnOnUntested = True
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--batch=")):
			batchSize = int(argVal)
		elif (arg == "--report:matrix") or (arg == "--report=matrix"):   # (unadvertised)
			reportAs = "matrix"
		elif (arg in ["--noendmark]","--noeof","--nomark"]):   # (unadvertised)
			requireEof = False
		elif (arg.startswith("--seed=")):
			seed = argVal
			if (seed in ["none","None","NONE"]):
				prngSeed = None
			elif (seed in ["default","Default","DEFAULT"]):
				prngSeed = defaultPrngSeed
			else:
				# nota bene: if the seed is a number, use it as a number, since
				#            string seeds can produce different sequences on
				#            different versions/builds of python
				try:
					seed = int(seed)
				except ValueError:
					try:               seed = float(seed)
					except ValueError: pass
				prngSeed = seed
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (reportAs == "matrix"):
		discardWhich = "none"

	# initialize the PRNG, if needed

	if (testMethod == "min-max"):
		if (prngSeed != None):
			random_seed(prngSeed)
	else:
		if (prngSeed not in [None,defaultPrngSeed]):
			print >>stderr, "WARNING: ignoring request to use PRNG with \"%s\"" \
			              % testMethod

	# make sure the shell commands we're gonna use have been installed

	if (testMethod == "chi-squared"):
		if (not shell_command_exists("Rscript")):
			exit(("%s: Unable to run the shell command \"Rscript\";"
			       + "\n  .. Either R hasn't been installed, or the command-line shell"
			       + " can't find it.")
			   % os_path.basename(argv[0]))

	# collect the alignments; we need to collect the positional info for all
	# alignments, to feed to R in batches (doing them one-by-one was incredibly
	# slow); hopefully this won't become a memory problem

	(unitLength,alignmentList,mxMatrix) \
	  = collect_alignments(stdin,testWhich,
	                       headLimit=headLimit,requireEof=requireEof)

	numAlignments = len(alignmentList)

	# assess the alignments, batch-by-batch

	accepted = []
	for batchStartIx in xrange(0,numAlignments,batchSize):
		batchEndIx = min(batchStartIx+batchSize,numAlignments)
		if ("batch" in debug):
			print >>stderr, "using R for alignments %d thru %d" \
			              % (batchStartIx+1,batchEndIx)

		mxBatch = mxMatrix[batchStartIx:batchEndIx]
		aBatch  = alignmentList[batchStartIx:batchEndIx]

		if (testMethod == "chi-squared"):
			resultBatch = mx_significance_tests(mxBatch,testWhich,effectSize,power)
			if (type(resultBatch) == str):
				exit(("%s: internal error: having trouble with R"
				        + " (with alignment batch %d..%d)"
				        + "\nHere's what R reported:\n%s")
				   % (os_path.basename(argv[0]),batchStartIx,batchEndIx,resultBatch))
		elif (testMethod == "min-max"):
			resultBatch = min_max_tests(aBatch,mxBatch,testWhich,numTrials)
			if (type(resultBatch) == str):
				exit(("%s: internal error: having trouble with min-max test"
				        + " (with alignment batch %d..%d)"
				        + "\nHere's what was reported:\n%s")
				   % (os_path.basename(argv[0]),batchStartIx,batchEndIx,resultBatch))
		else:
			exit("%s: internal error: unrecognized test method: \"%s\""
			   % (os_path.basename(argv[0]),testMethod))

		if (len(resultBatch) != batchEndIx-batchStartIx):
			exit(("%s: internal error: number of test outcomes reported by R (%d)"
			       + "\n  .. doesn't match the number of tests given to R (%d)")
			   % (os_path.basename(argv[0]),len(resultBatch),batchEndIx-batchStartIx))
		accepted += resultBatch

		if (warnOnUntested):
			for alignmentNum in xrange(batchStartIx,batchEndIx):
				testOutcome = accepted[alignmentNum]
				if (testOutcome == None):
					print >>stderr, "WARNING: alignment number %d (at line %d) could not be tested" \
					              % (alignmentNum,1+alignmentList[alignmentNum].lineNumber)

	outcomeCount = {True:0, False:0, None:0}
	for (alignmentNum,a) in enumerate(alignmentList):
		testOutcome = accepted[alignmentNum]
		outcomeCount[testOutcome] += 1

	# process the alignments and their assessments
	# $$$ untested alignments should be processed by some other test -- for
	#     example (if we're testing by error counts), a perfect alignment
	#     currently gets discarded because it can't be tested

	if (reportAs == "matrix"):
		outcomeMapping = {True:"not_rejected", False:"rejected", None:"untested"}
	else: # if (reportAs == "ncrf"):
		if (testWhich == "matches-insertions"):
			outcomeMapping = {True:  "match-insert uniformity not rejected",
			                  False: "match-insert uniformity rejected",
			                  None:  "untested"}
		elif (testWhich == "errors"):
			outcomeMapping = {True:  "error uniformity not rejected",
			                  False: "error uniformity rejected",
			                  None:  "untested"}
		else: # if (testWhich == "matches"):
			outcomeMapping = {True:  "match uniformity not rejected",
			                  False: "match uniformity rejected",
			                  None:  "untested"}

	outcomeNameW = max([len(outcomeMapping[testOutcome]) for testOutcome in outcomeMapping])
	for testOutcome in [True,False,None]:
		outcomeName = outcomeMapping[testOutcome]
		count       = outcomeCount[testOutcome]
		print >>stderr, "%-*s %d (%.2f%%)" \
					  % (outcomeNameW+1,"%s:" % outcomeName,count,100.0*count/numAlignments)

	if (reportAs == "matrix"):
		# see note [1] above for the format of the matrix file
		for (alignmentNum,a) in enumerate(alignmentList):
			testOutcome = accepted[alignmentNum]
			vec = [a.lineNumber,outcomeMapping[testOutcome]] + mxMatrix[alignmentNum]
			print "\t".join(map(str,vec))
	else: # if (reportAs == "ncrf"):
		numKept = 0
		isFirst = True
		for (alignmentNum,a) in enumerate(alignmentList):
			testOutcome = accepted[alignmentNum]
			if (discardWhich == "good"):
				if (testOutcome == True): continue
			elif (discardWhich == "bad"):
				if (testOutcome != True): continue

			if (discardWhich == "none"):
				chiSquaredInfo = "# positional chi-squared: %s" % outcomeMapping[testOutcome]
				(startIx,endIx) = a.positional_stats_indexes()
				a.lines.insert(endIx,chiSquaredInfo)

			if (isFirst): isFirst = False
			else:         print
			print a
			numKept += 1

		print >>stderr, "kept %d of %d alignments, %.2f%%" \
					  % (numKept,numAlignments,100.0*numKept/numAlignments)

		if (requireEof):
			print "# ncrf end-of-file"


# collect_alignments--

def collect_alignments(f,testWhich,headLimit=None,requireEof=True):
	alignmentList = []
	mxMatrix = []
	unitLength = None

	alignmentNum = 0
	for a in alignments(f,requireEof):
		alignmentNum +=1 

		if (headLimit != None) and (alignmentNum > headLimit):
			print >>stderr, "limit of %d alignments reached" % headLimit
			break

		if (testWhich == "matches-insertions"):
			mxRow = positional_error_vector(a,modified="m-i")
		else:
			mxRow = positional_error_vector(a)
		if (mxRow == None):
			raise ValueError, \
			      "alignment at line %d does not contain positional information" \
			    % a.lineNumber

		if (unitLength == None):
			unitLength = len(mxRow) / 2
		elif (len(mxRow) != 2*unitLength):
			raise ValueError, \
			      "alignments have different motif lengths, %d and %d (detected at line %d)" \
			    % (unitLength,len(mxRow)/2,a.lineNumber)

		alignmentList += [a]
		mxMatrix += [mxRow]

	return (unitLength,alignmentList,mxMatrix)


# positional_error_vector--

def positional_error_vector(a,modified=None):
	positionalStats = a.positional_stats()
	if (positionalStats == None):
		exit("%s: alignment at line %d has no positional stats"
		   % (os_path.basename(argv[0]),a.lineNumber))
	numPositions = len(positionalStats)
	vec = [None] * (2*numPositions)

	for (pos,stats) in enumerate(positionalStats):
		if ("m" not in stats):
			raise ValueError, \
				  "\"m\" missing from positional information for alignment at line %d" \
				% a.lineNumber
		if ("x" not in stats):
			raise ValueError, \
				  "\"x\" missing from positional information for alignment at line %d" \
				% a.lineNumber
		if (modified == "m-i") and ("i" not in stats):
			raise ValueError, \
				  "\"i\" missing from positional information for alignment at line %d" \
				% a.lineNumber

		if (modified == "m-i"): vec[pos] = stats["m"] - stats["i"]
		else:                   vec[pos] = stats["m"]

		vec[numPositions+pos] = stats["x"]

	return vec


# mx_significance_tests--
#	Run R in command-line mode to compute the significance test for a batch
#	of alignments.
#
# Input is an Nx2M matrix, where N is the number of alignments and M is the
# length of the aligned motif unit.  The first M columns are the positional
# counts for matches ("m"), and the final M columns are the positional counts
# for errors ("x").
#
# testWhich indicates whether to base the test on the matches or the errors.
# This is "matches", "errors", or "matches-insertions".
#
# If the process was unsuccessful, the return value is a string describing the
# the reason for the failure (hopefully).
#
# Otherwise, output is an N-element vector, containing acceptance/rejection for
# each row of the input.  Each entry is one of the following:
#   True  ==> don't reject null hypothesis; the counts are uniform
#   False ==> reject null hypothesis;       the counts are biased
#   None  ==> unable to run the test
#
# Each chi-squared test is performed on a 2xM matrix, corresponding to one row
# of our input matrix.

def mx_significance_tests(mxMatrix,testWhich,effectSize,power):
	# create a short one-line R program;  note that here each command is
	# a separate string in a list, and separating semi-colons will be added
	# when we go to run it

	testErrorCounts = "TRUE" if (testWhich == "errors") \
	             else "FALSE"  # this includes "matches" or "matches-insertions"

	n = len(mxMatrix)
	mxFlat = [v for rowVec in mxMatrix for v in rowVec]

	selfPath = os_path.dirname(os_path.realpath(__file__))

	rCommand =  []
	rCommand += ["source('%s/error_uniformity_filter.r')" % selfPath]
	rCommand += ["mxFlat = c(%s)" % ",".join(map(str,flatten(mxMatrix)))]
	rArgs = "%d,mxFlat,%s,%.10f,%.10f" % (n,testErrorCounts,effectSize,power)
	if ("rverbose" in debug): rArgs += ",verbose=T"
	rCommand += ["do_mx_significance_tests(%s)" % rArgs]

	# run that R program

	if ("r" in debug):
		print >>stderr, "calling R"

	(retCode,out,err) = run_rscript(rCommand)
	if (retCode != 0):
		return err

	lines = out.split("\n")
	if (lines[-1] == ""): lines = lines[:-1]
	mapping = {"TRUE":True, "FALSE":False, "NA":None}
	for line in lines:
		if (line not in mapping):
			exit("%s: internal error: having trouble with R; try reducing the --batch setting (%d)"
			   % (os_path.basename(argv[0]),batchSize))
	return map(lambda x:mapping[x],lines)


# run_rscript--

def run_rscript(rCommand):
	command = ["Rscript","--no-init-file","-e","; ".join(rCommand)]
	if ("rscript" in debug):
		print >>stderr, " ".join(command)
	process = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	(out,err) = process.communicate()
	retCode = process.wait()

	if ("r" in debug):
		print >>stderr, "=== return code ==="
		print >>stderr, retCode
		print >>stderr, "=== stdout ==="
		print >>stderr, out
		print >>stderr, "=== stderr ==="
		print >>stderr, err

	return (retCode,out,err)


# shell_command_exists--

def shell_command_exists(commandName):
	command = ["which",commandName]
	process = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	(out,err) = process.communicate()
	retCode = process.wait()

	# print >>stderr, "=== return code ==="
	# print >>stderr, retCode
	# print >>stderr, "=== stdout ==="
	# print >>stderr, out
	# print >>stderr, "=== stderr ==="
	# print >>stderr, err

	return (retCode == 0)


# min_max_tests--
#	"Judge" whether alignments in a batch are uniform or not.
#
# Input is an Nx2M matrix, where N is the number of alignments and M is the
# length of the aligned motif unit.  The first M columns are the positional
# counts for matches ("m"), and the final M columns are the positional counts
# for errors ("x").
#
# testWhich indicates whether to base the test on the matches or the errors.
# This is "matches", "errors", or "matches-insertions".
#
# If the process was unsuccessful, the return value is a string describing the
# the reason for the failure (hopefully).
#
# Otherwise, output is an N-element vector, containing acceptance/rejection for
# each row of the input.  Each entry is one of the following:
#   True  ==> don't reject null hypothesis; the counts are uniform
#   False ==> reject null hypothesis;       the counts are biased
#   None  ==> unable to run the test

def min_max_tests(aBatch,mxBatch,testWhich,numTrials):
	motifLen = len(mxBatch[0]) / 2

	results = [None] * len(mxBatch)
	for rowIx in xrange(len(mxBatch)):
		if (testWhich == "errors"):
			testCounts = mxBatch[rowIx][motifLen:2*motifLen]
		else:  # this includes "matches" or "matches-insertions"
			testCounts = mxBatch[rowIx][0:motifLen]

		minCount = min(testCounts)
		maxCount = max(testCounts)

		numEvents    = sum(testCounts)
		alignmentLen = len(aBatch[rowIx].seqText)

		minTrialCount = maxTrialCount = None
		for trialNum in xrange(numTrials):
			positionToCount = {pos:0 for pos in xrange(motifLen)}
			motifOffset = randint(0,motifLen-1)
			for aPos in random_sample(xrange(alignmentLen),numEvents):
				mPos = (aPos + motifOffset) % motifLen
				positionToCount[mPos] += 1

			trialMin = min([positionToCount[pos] for pos in xrange(motifLen)])
			trialMax = max([positionToCount[pos] for pos in xrange(motifLen)])

			if (minTrialCount == None) or (trialMin < minTrialCount):
				minTrialCount = trialMin
			if (maxTrialCount == None) or (trialMax < maxTrialCount):
				maxTrialCount = trialMax

		if (minCount >= minTrialCount) and (maxCount <= maxTrialCount):
			results[rowIx] = True
		else:
			results[rowIx] = False

		if ("min-max" in debug):
			print >>stderr, "[%d] %d/%d <%s> real:%d..%d random:%d..%d %s" \
			              % (rowIx,numEvents,alignmentLen,
			                 " ".join(map(str,testCounts)),
			                 minCount,maxCount,minTrialCount,maxTrialCount,
			                 "pass" if (results[rowIx]) else "fail")

	return results


# flatten--

def flatten(matrix):
	return [v for rowVec in matrix for v in rowVec]


if __name__ == "__main__": main()
