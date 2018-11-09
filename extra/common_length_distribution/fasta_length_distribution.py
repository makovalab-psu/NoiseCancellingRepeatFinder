#!/usr/bin/env python
"""
Read a fasta file and report the distribution of sequence lengths.
"""

from sys        import argv,stdin,stdout,stderr,exit
from ncrf_parse import int_with_unit,commatize

def usage(s=None):
	message = """

usage: cat <fasta_file> | fasta_length_distribution [options]
  --progress=<number>   periodically report how many sequences we've read

The resulting file has two columns -- a sequence length and the number of times
that length was observed. The lines are sorted by increasing length."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():

	# parse the command line

	reportProgress = None

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# process the fasta sequences

	lengthToCount = {}

	inputCount = inputBp = 0
	for (seqLen) in read_fasta_lengths(stdin):
		inputCount += 1
		inputBp    += seqLen

		if (reportProgress != None):
			if (inputCount % reportProgress == 0):
				print >>stderr, "%s sequences read (%s nts, avg=%s)" \
				              % (commatize(inputCount),commatize(inputBp),
				                 commatize(int(round(float(inputBp)/inputCount))))

		if (seqLen not in lengthToCount): lengthToCount[seqLen] =  1
		else:                             lengthToCount[seqLen] += 1

	# report the distribution

	lengths = [length for length in lengthToCount]
	lengths.sort()

	print "\n".join(["%d\t%d" % (length,lengthToCount[length]) for length in lengths])


# read_fasta_lengths--

def read_fasta_lengths(f):
	seqLen = None

	for line in f:
		line = line.rstrip()

		if (line.startswith(">")):
			if (seqLen != None): yield seqLen
			seqLen = 0
		elif (seqLen == None):
			seqLen = len(line)
		else:
			seqLen += len(line)

	if (seqLen != None): yield seqLen


if __name__ == "__main__": main()
