#!/usr/bin/env python
"""
Given several length distributions, report the maximum common distibution.
"""

from sys        import argv,stdin,stdout,stderr,exit
from math       import ceil
from ncrf_parse import parse_probability


def usage(s=None):
	message = """

usage: common_length_distribution [options]
  <distribution_file>  (cumulative, at least 2 required) length distribution
                       file for one of the input components
  --scale=<number>     factor to multiply the maximum distribution by
                       (0 < number <= 1; default is 1)

The length distribution files consist of two columns -- length, count -- a
length or interval and the number of sequences of that length in that component.

The maximum common distibution is the minimum, over each length, of the
component counts for that length (example below). We call this "maximum"
because it give the largest number that can be sampled from each component to
achieve a common length distribution.

  component 1    component 2    maximum common
  375-378 800    375-378 600    375-378 600
  379-382 900    379-382 750    379-382 750
  383-386 1000   383-386 1040   383-386 1000
  387-391 1050   387-391 1200   387-391 1050
  392-395 1100   392-395 1000   392-395 1000"""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():

	# parse the command line

	componentFilenames = []
	scale              = None

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--scale=")):
			try:
				scale = parse_probability(argVal)
				if (scale == 0.0): raise ValueError
			except ValueError:
				usage("%s is not a valid scale setting" % arg)
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			componentFilenames += [arg]

	if (len(componentFilenames) < 2):
		usage("you must provide at least two length distribution filenames")

	# read the component length distributions

	intervalsSeen       = set()
	intervalsOrder      = []
	filenameToIntervals = {}

	for filename in componentFilenames:
		assert (filename not in filenameToIntervals), \
		       "length distribution %s appears more than once" % filename

		componentF = file(filename,"rt")
		filenameToIntervals[filename] = intervals = {}

		for line in read_length_distribution(componentF,filename):
			(lineNumber,intervalStr,count) = line

			assert (intervalStr not in intervals), \
			       "(at line %d) interval %s appears more than once in %s" \
			     % (lineNumber,intervalStr,filename)

			intervals[intervalStr] = count

			if (intervalStr not in intervalsSeen):
				intervalsSeen.add(intervalStr)
				intervalsOrder += [intervalStr]

		componentF.close ()

	# determine the maximum common distibution

	commonIntervals = {}

	for intervalStr in intervalsOrder:
		inAllComponents = True
		minCount = None
		for filename in componentFilenames:
			if (intervalStr not in filenameToIntervals[filename]):
				inAllComponents = False
				break
			count = filenameToIntervals[filename][intervalStr]
			if (minCount == None) or (count < minCount):
				minCount = count

		if (not inAllComponents): continue

		if (scale != None): minCount = int(ceil(minCount*scale))
		if (minCount > 0): print "%s\t%d" % (intervalStr,minCount)


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
