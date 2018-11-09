#!/usr/bin/env python
"""
Compute the sum of several sequence length distribution files.
"""

from sys        import argv,stdin,stdout,stderr,exit
from ncrf_parse import commatize

def usage(s=None):
	message = """

usage: cat <length_distribution_files> | sum_length_distributions [options]
  --report:totals   report total bp, number of sequences, and averge bp per
                    sequence

Input files are the same as the output of fasta_length_distribution -- two
columns, a sequence length and the number of times that length was observed.
The lines are sorted by increasing length.

Output is the same format."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():

	# parse the command line

	reportTotals = False

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg in ["--report:totals","--report:total"]):
			reportTotals = True
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# accumulate the length distributions

	lengthToCount = {}

	for (length,count) in read_length_counts(stdin):
		if (length not in lengthToCount): lengthToCount[length] =  count
		else:                             lengthToCount[length] += count

	# report the total distribution

	lengths = [length for length in lengthToCount]
	lengths.sort()

	print "\n".join(["%d\t%d" % (length,lengthToCount[length]) for length in lengths])

	if (reportTotals):
		numSequences = sum([lengthToCount[length] for length in lengths])
		if (numSequences == 0):
			print >>stderr, "0 sequences / 0 bp total"
		else:
			totalBp     = sum([lengthToCount[length]*length for length in lengths])
			avgSequence = int(round(float(totalBp) / numSequences))
			print >>stderr, "%s sequences / %s bp total / %s bp average" \
			              % (commatize(numSequences),commatize(totalBp),commatize(avgSequence))


# read_length_counts--
#	yield (length,count) pairs

def read_length_counts(f,filename=None):
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
			length = int(fields[0])
			count  = int(fields[1])
		except ValueError:
			assert (False), \
			       "bad length,count at %s line %d\n%s" \
			     % (filename,lineNumber,line)

		yield (length,count)


if __name__ == "__main__": main()
