#!/usr/bin/env python
"""
Given a length distribution and a target distibution, create a spec to be used
with fasta_match_length_distribution
"""

from sys import argv,stdin,stdout,stderr,exit


def usage(s=None):
	message = """

usage: make_length_distribution_spec <target> <distribution>
  <target>             (required) length distribution file for the desired
                       fasta output of fasta_match_length_distribution
  <distribution>       (required) length distribution file corresponding to
                       a fasta file (or a set of fasta files) that will be
                       input to fasta_match_length_distribution

The length distribution files consist of two columns -- length, count -- a
length or interval and the number of sequences of that length in that
distribution.

The output is a distribution spec suitable for fasta_match_length_distribution."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	# parse the command line

	targetFilename       = None
	distributionFilename = None

	for arg in argv[1:]:
		if (targetFilename == None):
			targetFilename = arg
		elif (distributionFilename == None):
			distributionFilename = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (targetFilename == None):
		usage("you must provide a target length-distribution filename")

	if (distributionFilename == None):
		usage("you must provide a length-distribution filename")

	# read the target distribution

	intervalsSeen  = set()
	intervalsOrder = []

	targetF = file(targetFilename,"rt")
	targetIntervals = {}

	for line in read_length_distribution(targetF,targetFilename):
		(lineNumber,intervalStr,count) = line

		assert (intervalStr not in targetIntervals), \
		       "(at line %d) interval %s appears more than once in %s" \
		     % (lineNumber,intervalStr,targetFilename)

		targetIntervals[intervalStr] = count

		if (intervalStr not in intervalsSeen):
			intervalsSeen.add(intervalStr)
			intervalsOrder += [intervalStr]

	targetF.close ()

	# read the distribution

	distribF = file(distributionFilename,"rt")
	distribIntervals = {}

	for line in read_length_distribution(distribF,distributionFilename):
		(lineNumber,intervalStr,count) = line

		assert (intervalStr not in distribIntervals), \
		       "(at line %d) interval %s appears more than once in %s" \
		     % (lineNumber,intervalStr,distributionFilename)

		if (intervalStr not in intervalsSeen): continue

		distribIntervals[intervalStr] = count
		assert (count >= targetIntervals[intervalStr]), \
		       "(at line %d in %s) interval %s's count (%d) is less than target's (%d)" \
		     % (lineNumber,distributionFilename,
		        intervalStr,count,targetIntervals[intervalStr])

	distribF.close ()

	# write the distribution spec

	for intervalStr in intervalsOrder:
		if (intervalStr not in distribIntervals): continue
		if (targetIntervals[intervalStr] == 0): continue

		print "%s\t%d\t%d" \
		    % (intervalStr,targetIntervals[intervalStr],distribIntervals[intervalStr])


# read_length_distribution--
#	yield (lineNumber,intervalStr,count) tuples

def read_length_distribution(f,filename=None):
	if (filename == None): filename = "(input)"

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		if (len(fields) != 2):
			assert (False), \
			       "bad length-distribution at %s line %d (expected exactly two fields)\n%s" \
			     % (filename,lineNumber,line)

		try:
			intervalStr = fields[0]
			count     = int(fields[1])
		except ValueError:
			assert (False), \
			       "bad length-distribution at %s line %d\n%s" \
			     % (filename,lineNumber,line)

		yield (lineNumber,intervalStr,count)


if __name__ == "__main__": main()
