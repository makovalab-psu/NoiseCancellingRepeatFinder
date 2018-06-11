#!/usr/bin/env python
"""
Pick a random subset of fasta sequences, matching a given length distribution.
"""

from sys           import argv,stdin,stdout,stderr,exit
from math          import ceil
from random        import seed as random_seed,randint
from interval_dict import IntervalDict,Interval
from ncrf_parse    import int_with_unit,commatize


def usage(s=None):
	message = """

usage: cat <fasta_file> | fasta_match_length_distribution [options]
  <distribution_spec>   (required) file describing the length distribution
  --remainder=<file>    write unfulfilled length distribution to a file
  --wrap=<length>       number of nucleotides per line in output fasta
  --seed=<string>       random number generator seed
  --progress=<number>   periodically report how many sequences we've read

The distribution spec file consists of three columns-- a length (or interval),
the number of sequences of that length we should output, and the number of
sequences of that length we expect to see in input.  An interval is two numbers
connected by a hyphen, e.g. "100-199" indicates sequences at least 100 bp long
but shorter than 200 bp.

The distribution of input lengths allows us to easily select any read with the
correct probability, without having to pre-scan the input to collect that
information. It also allows us to operate over a collection of files, using the
remainder file to communicate from the processing of one file to the next."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global debug

	# parse the command line

	distributionFilename = None
	remainderFilename    = None
	wrapLength           = 100
	reportProgress       = None
	debug                = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--remainder=")):
			remainderFilename = argVal
		elif (arg.startswith("--wrap=")):
			wrapLength = int(argVal)
			if (wrapLength <= 0): wrapLength = None
		elif (arg.startswith("--seed=")):
			random_seed(argVal)
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (distributionFilename == None):
			distributionFilename = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (distributionFilename == None):
		usage("you must provide a length-distribution filename")

	# read the distribution

	intervals = IntervalDict()

	distribF = file(distributionFilename,"rt")
	for spec in read_distribution_spec(distribF,distributionFilename):
		(lineNumber,minLength,maxLength,outCount,inCount) = spec

		interval = intervals.add(minLength,maxLength)
		if (interval == None): # interval overlaps an existing interval
			interval = Interval(minLength,maxLength)
			previous = intervals.overlapper(minLength,maxLength)
			assert (False), \
			       "%s (line %d) overlaps %s (line %d)" \
			     % (interval,lineNumber,previous,previous.lineNumber)

		interval.lineNumber = lineNumber
		interval.outCount   = outCount
		interval.inCount    = inCount

	distribF.close ()

	if ("distribution" in debug):
		for interval in intervals:
			print >>stderr, "%s %d %d" \
			              % (interval,interval.outCount,interval.inCount)

	# process the reads
	#
	# this filters reads based on the length (on the interval containing the
	# length); if we expect to see E more sequences of this length (including
	# this one), and we are to output N of those, we output this sequence with
	# probability N/E; and we adjust N and E for this length accordingly

	inputCount = outputCount = inputBp = outputBp = 0
	for (name,seq) in read_fasta_sequences(stdin):
		seqLen = len(seq)
		inputCount += 1
		inputBp    += seqLen

		if (reportProgress != None):
			if (inputCount % reportProgress == 0):
				print >>stderr, "%s sequences read, %s written (%.1f%%); %s nts read, %s written" \
				              % (commatize(inputCount),commatize(outputCount),
				                 100.0*outputCount/inputCount,
				                 commatize(inputBp),commatize(outputBp))

		try: interval = intervals[seqLen]
		except KeyError: continue

		if (interval.inCount <= 0):
			print >>stderr, "ERROR: for length %d (%s), actual input exceeded expected input count" \
			              % (seqLen,interval)
			if (remainderFilename != None):
				print >>stderr, "      (writing remainders to %s)" % remainderFilename
				remainderF = file(remainderFilename,"wt")
				write_remainders(remainderF,intervals)
				remainderF.close ()
			assert (False)

		if (interval.outCount == 0):
			keepSeq = False
		else:
			keepSeq = (randint(1,interval.inCount) <= interval.outCount)

		interval.inCount  -= 1
		if (not keepSeq): continue

		interval.outCount -= 1
		outputCount += 1
		outputBp    += seqLen
		print ">%s" % name
		if (wrapLength == None):
			print seq
		else:
			for i in range(0,seqLen,wrapLength):
				print seq[i:i+wrapLength]

	# write the remainders

	if (remainderFilename != None):
		remainderF = file(remainderFilename,"wt")
		write_remainders(remainderF,intervals)
		remainderF.close ()


# read_distribution_spec--
#	yield (lineNumber,minLength,maxLength,outCount,inCount) tuples

def read_distribution_spec(f,filename=None):
	if (filename == None): filename = "(input)"

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		if (len(fields) != 3):
			assert (False), \
			       "bad length-distribution at %s line %d (expected exactly three fields)\n%s" \
			     % (filename,lineNumber,line)

		try:
			if ("-" in fields[0]):
				(minLength,maxLength) = map(int,fields[0].split("-",1))
				if (minLength > maxLength): raise ValueError
			else:
				minLength = maxLength = int(fields[0])
			outCount = int(fields[1])
			inCount  = int(fields[2])
			if (outCount > inCount): raise ValueError
		except ValueError:
			assert (False), \
			       "bad length-distribution at %s line %d\n%s" \
			     % (filename,lineNumber,line)

		yield (lineNumber,minLength,maxLength,outCount,inCount)


# read_fasta_sequences--
#	yield fasta (name,sequence) pairs

def read_fasta_sequences(f):
	name  = "(no_name)"
	lines = []

	lineNum = 0
	for line in f:
		lineNum += 1
		line = line.rstrip()

		if (line.startswith(">")):
			if (lines != []):
				seq = "".join(lines)
				yield (name,seq)
			name  = line[1:].strip()
			lines = []
		else:
			lines += [line.upper()]

	if (lines != []):
		seq = "".join(lines)
		yield (name,seq)


# write_remainders--

def write_remainders(f,intervals):
	for interval in intervals:
		print >>f, "%s %d %d" \
		         % (interval,interval.outCount,interval.inCount)


if __name__ == "__main__": main()
