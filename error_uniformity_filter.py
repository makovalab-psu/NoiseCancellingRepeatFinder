#!/usr/bin/env python
"""
Filter Noise Cancelling Repeat Finder alignments, discarding alignments that
don't (seem to) have uniformly distributed matches across the motif unit.
"""

import subprocess
from sys        import argv,stdin,stdout,stderr,exit
from os         import path as os_path
from random     import seed as random_seed,sample as random_sample,randint
from ncrf_parse import alignments,parse_probability,int_with_unit,commatize

defaultPrngSeed = "NCRF.error_uniformity_filter"


# Implementation Notes:
#
# [1] Justification for "matches-insertions" (as opposed to just "matches").
#     When NCRF counts matches, a deletion (a base deleted from the motif)
#     outcomes in a decrease of the count of matches at that position. There is
#     no similar effect for an insertion (a base inserted into the motif).
#     Thus if we base our test solely on matches, insertions will not
#     contribute to the test. Subtracting insertion counts from match counts
#     gives insertions the same effect as deletions.
#
# [2] The min-max method introduces some non-determinism. It is expected that
#     this will only affect alignments near the pass/fail boundary. Whether
#     an alignment is passed or failed depends to some extent on the state of
#     the PRNG. This state is affected by the initial seed (which will be the
#     same for all runs unless the user sets it) and how many times the PRNG
#     has been called in processing previous alignments. Thus, even without
#     changing the seeding, an alignment processed in different runs with
#     different parameters (or appearing at a different relative position in
#     the input) is processed with a different PRNG state. It is likely that
#     alignments not close to the pass/fail boundary would not be affected.
#
# [3] Here we describe the format of the file produced by the unadvertised
#     option --report:matrix. The reason this option isn't advertised is that
#     is should only be useful to developers testing implementation details of
#     the underlying chi-squared test.
#
#     The file has NO header line and one row for each alignment;  2M+2 columns
#     per row, where M is the motif length. Column 1 is the line number of the
#     alignment in the incoming file. Column 2 is the outcome of the test. 
#     Columns 3 thru M+2 are positional match counts for that alignment, and
#     columns M+3 thru 2M+2 are the positional error counts. We assume all rows
#     have the same number of columns, i.e. that the same motif length is
#     represented in all rows.
#
# [4] The chi-squared test is no longer advertised, because it was found to be
#     flawed.

def usage(s=None):
	message = """
usage: cat <output_from_NCRF> | error_uniformity_filter [options]
  --method=min-max      judge alignments by a "min-max" test  
                        (this is the default)
  --trials=<number>     number of trials for the min-max test
                        (default is 10K)
  --trials=<number>/<number> numbers of successes needed and rials for the
                        min-max test; e.g. "5/10K" will do 10 thousand trials
                        and require at least five success to "pass"
  --discard:good        discard the "good" alignments instead of the "bad" ones
                        (by default we discard the "bad" alignments) 
  --discard:none        don't discard any alignments, just report the test
                        outcomes (see below for how they are reported) 
  --test:matches-insertions perform the test using match counts less insertions
                        (this is the default)
  --test:matches        perform the test using match counts
  --test:errors         perform the test using error counts; this may increase
                        the likelihood that the test cannot be performed for
                        some alignments, since the counts may be too low
  --warn:untested       report each untestable alignment as a warning
  --seed=<string>       set random seed; this is only applicable to the min-max
                        test; two special cases are "default" (use the built-in
                        default seed) and "none" (don't seed the random number
                        generator)
  --subsample=<k>/<n>   Conceptually split the alignments by into <n> groups,
                        and only process the <k>th group; <k> ranges from 1 to
                        <n>
  --head=<number>       limit the number of input alignments
  --progress=<number>   periodically report how many alignments we've tested
  --batch=<number>      number of input alignments processed by each call to R;
                        our ability to call R fails if the command line we pass
                        it is too long
                        (default is 30)

In a "true" alignment to a given motif unit, we expect the errors to be
distributed randomly and uniformly among the positions in the unit. (That is
an underlying assumption, possibly not true.)  This program discards alignments
that fail a test based on that assumption.

Error counts are often too small for the statistical test, so a test based on
match counts is used instead. However, an insertion error does not reduce the
match count at any position, so by default we decrease matches by the number of
insertions at that position.

The input alignments must include position event information. This can be
accomplished by using the --positionalevents option of Noise Cancelling Repeat
Finder.

When test outcomes are reported for --discard:none, one of the following lines
is added to each alignment
# positional min-max: match-insert uniformity rejected
# positional min-max: match-insert uniformity not rejected
# positional min-max: untested
The "untested" case indicates that the min-max test could not be performed,
usually because one of the positional match counts is too small.

The user should be aware that the results aren't necessarily determinstic.
When a PRNG is in play (as for min-max), the result for an alignment depends
on the state of the PRNG; and the state of the PRNG depends on all the
alignments that preceded the one being tested."""

# options no longer advertised:
#  --method=chi-squared  judge alignments by a chi-squared test
#  --effectsize=<value>  effect size for chi-squared test
#                        (default is 0.3)
#  --power=<probability> "power of test" for chi-squared test, 1 minus Type II
#                        error probability
#                        (default is 0.8)

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global reportProgress,batchSize
	global debug

	# parse the command line

	testMethod      = "min-max"
	numTrials       = 10*1000       # (only used for testMethod == "min-max")
	numNeededToPass = 1             # (only used for testMethod == "min-max")
	effectSize      = 0.3           # (only used for testMethod == "chi-square")
	power           = 0.8           # (only used for testMethod == "chi-square")
	discardWhich    = "bad"
	testWhich       = "matches-insertions"
	warnOnUntested  = False
	subsampleK      = None
	subsampleN      = None
	headLimit       = None
	batchSize       = None  # (will be replace by method-specific result)
	reportAs        = "ncrf"
	requireEof      = True
	prngSeed        = defaultPrngSeed
	reportProgress  = None
	debug           = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg == "--method=min-max"):
			testMethod = "min-max"
		elif (arg.startswith("--trials=")):
			if ("/" in argVal):
				(numNeededToPass,numTrials) = map(int_with_unit,argVal.split("/",1))
				if (numTrials < 1):
					usage("bad value in: %s (trials must be at least 1)" % arg)
				if (not 1 <= numNeededToPass <= numTrials):
					usage("bad value in: %s (num-in-bounds must be in range 1..trials)" % arg)
			else:
				(numNeededToPass,numTrials) = (1,int_with_unit(argVal))
				if (numTrials < 1):
					usage("bad value in: %s (trials must be at least 1)" % arg)
		elif (arg in ["--method=chi-squared","--method=chi-square"]):    # (unadvertised, see [4])
			testMethod = "chi-squared"
		elif (arg.startswith("--effectsize=")):                          # (unadvertised, see [4])
			effectSize = parse_probability(argVal)
		elif (arg.startswith("--power=")):                               # (unadvertised, see [4])
			power = parse_probability(argVal)
		elif (arg in ["--discard:bad","--discard=bad"]):
			discardWhich = "bad"
		elif (arg in ["--discard:good","--discard=good"]):
			discardWhich = "good"
		elif (arg in ["--discard:none","--discard=none"]):
			discardWhich = "none"
		elif (arg in ["--test:matches-insertions","--test=matches-insertions",
		              "--test:m-i","--test=m-i"]):
			testWhich = "matches-insertions"
		elif (arg in ["--test:matches","--test=matches"]):
			testWhich = "matches"
		elif (arg in ["--test:errors","--test=errors"]):
			testWhich = "errors"
		elif (arg == "--warn:untested") or (arg == "--warn=matrix"):
			warnOnUntested = True
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--subsample=")):
			(subsampleK,subsampleN) = map(int,argVal.split("/",2))
			if (not 0 < subsampleK <= subsampleN):
				usage("bad subsample description in %s" % arg)
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg.startswith("--batch=")):
			batchSize = int(argVal)
		elif (arg == "--report:matrix") or (arg == "--report=matrix"):   # (unadvertised)
			reportAs = "matrix"
		elif (arg == "--report:silent") or (arg == "--report=silent"):   # (unadvertised)
			reportAs = "silent"
		elif (arg in ["--noendmark","--noeof","--nomark"]):   # (unadvertised)
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

	if (reportAs in ["matrix","silent"]):
		discardWhich = "none"

	if (testMethod == "chi-squared"):
		testDescription = "positional chi-squared"
		if (batchSize == None): batchSize = 30
	elif (testMethod == "min-max"):
		testDescription = "positional min-max"
		if (batchSize == None): batchSize = 1
	else:
		exit("%s: internal error: unrecognized test method: \"%s\""
		   % (os_path.basename(argv[0]),testMethod))

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
	                       headLimit=headLimit,
	                       subsampleK=subsampleK,subsampleN=subsampleN,
	                       requireEof=requireEof)

	numAlignments = len(alignmentList)
	if (reportProgress != None):
		print >>stderr, "progress: read %s alignments" \
		              % (commatize(numAlignments))

	# assess the alignments, batch-by-batch

	if (reportProgress != None):
		progressReported = -1

	accepted = []
	outcomeCount = {True:0, False:0, None:0}
	for batchStartIx in xrange(0,numAlignments,batchSize):
		alignmentsTested = batchStartIx
		if (reportProgress != None):
			rBlock = (progressReported+1)/reportProgress
			aBlock = (alignmentsTested+1)/reportProgress
			if (alignmentsTested == 0) or (aBlock != rBlock):
				print >>stderr, "progress: testing alignment %s (%d uniform, %d non-uniform, %d untested)" \
				              % (commatize(1+alignmentsTested),
				                 outcomeCount[True],
				                 outcomeCount[False],
				                 outcomeCount[None])
				progressReported = alignmentsTested

		batchEndIx = min(batchStartIx+batchSize,numAlignments)
		if ("batch" in debug):
			print >>stderr, "using R for alignments %d thru %d" \
			              % (batchStartIx+1,batchEndIx)

		mxBatch = mxMatrix[batchStartIx:batchEndIx]
		aBatch  = alignmentList[batchStartIx:batchEndIx]

		if (testMethod == "chi-squared"):
			batchResult = mx_significance_tests(mxBatch,testWhich,effectSize,power)
			if (type(batchResult) == str):
				exit(("%s: internal error: having trouble with R"
				        + " (with alignment batch %d..%d)"
				        + "\nHere's what R reported:\n%s")
				   % (os_path.basename(argv[0]),batchStartIx,batchEndIx,batchResult))
		else:  # if (testMethod == "min-max"):
			batchResult = min_max_tests(aBatch,mxBatch,batchStartIx,testWhich,numTrials,numNeededToPass)
			if (type(batchResult) == str):
				exit(("%s: internal error: having trouble with min-max test"
				        + " (with alignment batch %d..%d)"
				        + "\nHere's what was reported:\n%s")
				   % (os_path.basename(argv[0]),batchStartIx,batchEndIx,batchResult))

		if (len(batchResult) != batchEndIx-batchStartIx):
			exit(("%s: internal error: number of test outcomes reported by R (%d)"
			       + "\n  .. doesn't match the number of tests given to R (%d)")
			   % (os_path.basename(argv[0]),len(batchResult),batchEndIx-batchStartIx))
		accepted += batchResult

		if (warnOnUntested):
			for alignmentNum in xrange(batchStartIx,batchEndIx):
				testOutcome = accepted[alignmentNum]
				if (testOutcome == None):
					print >>stderr, "WARNING: alignment number %d (at line %d) could not be tested" \
					              % (alignmentNum,1+alignmentList[alignmentNum].lineNumber)

		for alignmentNum in xrange(batchStartIx,batchEndIx):
			testOutcome = accepted[alignmentNum]
			outcomeCount[testOutcome] += 1

	# process the alignments and their assessments
	# $$$ untested alignments should be processed by some other test -- for
	#     example (if we're testing by error counts), a perfect alignment
	#     currently gets discarded because it can't be tested

	if (reportAs in ["matrix","silent"]):
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
		reportStr = "%-*s %d" % (outcomeNameW+1,"%s:" % outcomeName,count)
		if (numAlignments > 0):
			reportStr += " (%.2f%%)" % (100.0*count/numAlignments)
		print >>stderr, reportStr

	if (reportAs == "matrix"):
		# see note [3] above for the format of the matrix file
		for (alignmentNum,a) in enumerate(alignmentList):
			testOutcome = accepted[alignmentNum]
			vec = [a.lineNumber,outcomeMapping[testOutcome]] + mxMatrix[alignmentNum]
			print "\t".join(map(str,vec))
	elif (reportAs == "silent"):
		pass
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
				testInfo = "# %s: %s" % (testDescription,outcomeMapping[testOutcome])
				(startIx,endIx) = a.positional_stats_indexes()
				a.lines.insert(endIx,testInfo)

			if (isFirst): isFirst = False
			else:         print
			print a
			numKept += 1

		reportStr = "kept %d of %d alignments" % (numKept,numAlignments)
		if (numAlignments > 0):
			reportStr += ", %.2f%%" % (100.0*numKept/numAlignments)
		print >>stderr, reportStr

		if (requireEof):
			print "# ncrf end-of-file"


# collect_alignments--

def collect_alignments(f,testWhich,
                       headLimit=None,
                       subsampleK=None,subsampleN=None,
                       requireEof=True):
	alignmentList = []
	mxMatrix = []
	unitLength = None

	alignmentNum = 0
	for a in alignments(f,requireEof):
		alignmentNum += 1 
		if    (reportProgress != None) \
		  and ((alignmentNum == 1) or (alignmentNum % reportProgress == 0)):
			print >>stderr, "progress: reading alignment %s" \
			              % (commatize(alignmentNum))

		if (headLimit != None) and (alignmentNum > headLimit):
			print >>stderr, "limit of %d alignments reached" % headLimit
			break

		if (subsampleN != None):
			if ((alignmentNum-1) % subsampleN != (subsampleK-1)): continue

		if (testWhich == "matches-insertions"):
			# note [1]
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
# length of the aligned motif unit. The first M columns are the positional
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
# each row of the input. Each entry is one of the following:
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
#	"Judge" whether alignments in a batch are uniform or not. This test is
#	non-deterministic (see note [2]).
#
# Input is an Nx2M matrix, where N is the number of alignments and M is the
# length of the aligned motif unit. The first M columns are the positional
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
# each row of the input. Each entry is one of the following:
#   True  ==> don't reject null hypothesis; the counts are uniform
#   False ==> reject null hypothesis;       the counts are biased
#   None  ==> unable to run the test
#
# The following paragraph only applies when numNeededToPass is 1.
#
# Note that for a given alignment, we stop generating random count vectors as
# soon as we have a range of min and max that will contain the alignment's min
# and max. In many cases this happens in the first few trials if the alignment
# will pass (is uniform); alignments that fail will necessarily require the
# full complement of trials. In one test with motifLen=5 and alignments longer
# than 500 bp, the average number of trials before an alignment passed was just
# under 50.
#
# If the user has requested "min-max:complete" in debug, we always run the full
# number of trials.


def min_max_tests(aBatch,mxBatch,batchStartIx,testWhich,numTrials,numNeededToPass=1):
	motifLen = len(mxBatch[0]) / 2

	results = [None] * len(mxBatch)
	for rowIx in xrange(len(mxBatch)):
		if (testWhich == "errors"):
			testCounts = mxBatch[rowIx][motifLen:2*motifLen]
		else:  # this includes "matches" or "matches-insertions"
			testCounts = mxBatch[rowIx][0:motifLen]

		minRealCount = min(testCounts)
		maxRealCount = max(testCounts)

		numEvents    = sum(testCounts)
		alignmentLen = len(aBatch[rowIx].seqText)

		minProvenOk = maxProvenOk = 0
		realPasses = False  # (until proven otherwise)
		trialsMade = 0
		for trialNum in xrange(numTrials):
			positionToCount = [0] * motifLen  # (list performed faster than dict)

			motifOffset = randint(0,motifLen-1)
			for aPos in random_sample(xrange(alignmentLen),numEvents):
				mPos = (aPos + motifOffset) % motifLen
				positionToCount[mPos] += 1

			trialsMade += 1

			if (minProvenOk < numNeededToPass) or ("min-max:complete" in debug):
				trialMin = min([positionToCount[pos] for pos in xrange(motifLen)])
				if (minRealCount >= trialMin): minProvenOk += 1
			if (maxProvenOk < numNeededToPass) or ("min-max:complete" in debug):
				trialMax = max([positionToCount[pos] for pos in xrange(motifLen)])
				if (maxRealCount <= trialMax): maxProvenOk += 1

			if (minProvenOk >= numNeededToPass) and (maxProvenOk >= numNeededToPass):
				realPasses = True
				if ("min-max:complete" not in debug):
					break

		results[rowIx] = realPasses

		if   ("min-max" in debug) \
		  or (("min-max:fail" in debug) and (not results[rowIx])):
			print >>stderr, ("[%d] %s"
			               + " (%d trials, real min>=%d trials, real max<=%d trials)"
			               + " evRatio=%d/%d, real min..max=%d..%d, counts=<%s>") \
			              % (1+batchStartIx+rowIx,
			                 "pass" if (results[rowIx]) else "fail",
			                 trialsMade,minProvenOk,maxProvenOk,
			                 numEvents,alignmentLen,
			                 minRealCount,maxRealCount,
			                 ",".join(map(str,testCounts)))

	return results


# flatten--

def flatten(matrix):
	return [v for rowVec in matrix for v in rowVec]


if __name__ == "__main__": main()
