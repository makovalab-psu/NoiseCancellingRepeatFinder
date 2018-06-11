#!/usr/bin/env python
"""
Spread a given length distribution over several component distributions.
"""

from sys    import argv,stdin,stdout,stderr,exit
from random import seed as random_seed,sample as random_sample

class Interval: pass


def usage(s=None):
	message = """

usage: cat <distribution_spec> | spread_length_distribution [options]
  <components_file>    (required) file listing the components
  --input=<filespec>   (required) length distribution filename spec, e.g.
                       "{component}.dat"
  --output=<filespec>  (required) distribution spec filename spec, e.g.
                       "{component}.spec"
  --nowarn             don't report a warning for intervals in components that
                       aren't in the distribution spec; such intervals are not
                       processed (whether we report them or not)
  --seed=<string>      random number generator seed

The distribution spec file consists of three columns -- length, outCount,
inCount -- a length or interval, the number of sequences of that length to be
output, and the number of sequences of that length we expect would be seen in
input. This is the same format used as input to fasta_match_length_distribution.

The components file lists component names, one per line. These are mapped to
length distribution files by the input filespec.  Each of the component length
distribution files consists of two columns -- length, inCount -- a length or
interval and the number of sequences of that length in that component.

It is assumed that the lengths in the distribution spec and all the components
are a perfect match, and that a length's inCount in the distribution spec
matches the sum of the length's inCount in all the components.  The goal of
this program is then to spread the outCount of the distribution spec to the
components, and thus create a distribution spec for each component.

For example, suppose we have five lengths/intervals and these three components
and dsitribution spec:

  component 1    component 2    component 3    distribution spec
  375-378 826    375-378 103    375-378 6804   375-378 4700 7733
  379-382 831    379-382 83     379-382 6893   379-382 2700 7807
  383-386 835    383-386 92     383-386 7043   383-386 4800 7970
  387-391 1026   387-391 131    387-391 8950   387-391 1900 10107
  392-395 820    392-395 90     392-395 7387   392-395 5100 8297

The goal, for the first row, is to spread outCount 4700 among the three
components, essentially randomly sampling from inCount 7733=826+103+6804 with
probabilities corresponding to each component. In the typical result shown
below, that split was 4700=509+66+4125.

  component 1        component 2      component 3
  375-378 509 826    375-378 66 103   375-378 4125 6804
  379-382 283 831    379-382 28 83    379-382 2389 6893
  383-386 507 835    383-386 53 92    383-386 4240 7043
  387-391 199 1026   387-391 25 131   387-391 1676 8950
  392-395 515 820    392-395 55 90    392-395 4530 7387

Note that we make no effort to maintain proportions, but the sampling process
will usually report in something close to proportionality."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global debug

	# parse the command line

	componentListFilename = None
	inputFilespec         = None
	outputFilespec        = None
	warnForExtraIntervals = True
	reportSpreadProgress  = False
	reportOutputWrites    = False
	debug                 = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--input=")):
			inputFilespec = argVal
			if ("{component}" not in inputFilespec):
				usage("\"%s\" doesn't contain \"{component}\"" % inputFilespec)
		elif (arg.startswith("--output=")):
			outputFilespec = argVal
			if ("{component}" not in inputFilespec):
				usage("\"%s\" doesn't contain \"{component}\"" % outputFilespec)
		elif (arg == "--nowarn"):
			warnForExtraIntervals = False
		elif (arg == "--progress=spread"):
			reportSpreadProgress = True
		elif (arg == "--progress=written"):
			reportOutputWrites = True
		elif (arg.startswith("--seed=")):
			random_seed(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (componentListFilename == None):
			componentListFilename = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (componentListFilename == None):
		usage("you must provide a component list filename")

	if (inputFilespec == None):
		usage("you must provide an input component filename spec")

	if (outputFilespec == None):
		usage("you must provide an output component filename spec")

	if (outputFilespec == inputFilespec):
		usage("input and output component filename specs can't be the same!")

	# read the distribution spec

	intervalsOrder    = []
	combinedIntervals = {}

	for spec in read_distribution_spec(stdin):
		(lineNumber,intervalStr,outCount,inCount) = spec

		assert (intervalStr not in combinedIntervals), \
		       "(at line %d) interval %s appears more than once in the distribution spec" \
		     % (lineNumber,intervalStr)

		intervalsOrder += [intervalStr]
		combinedIntervals[intervalStr] = interval = Interval()
		interval.outCount = outCount
		interval.inCount  = inCount

	# read the component length distributions

	componentsOrder     = []
	componentToIntervals = {}

	componentListF = file(componentListFilename,"rt")

	lineNumber = 0
	for line in componentListF:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue

		name = line
		assert (name not in componentToIntervals), \
		       "(at line %d) component %s appears more than once in %s" \
		     % (lineNumber,name,componentListFilename)

		componentFilename = inputFilespec.replace("{component}",name)
		componentF = file(componentFilename,"rt")
		componentsOrder += [name]
		componentToIntervals[name] = componentIntervals = {}

		for line in read_length_distribution(componentF,componentFilename):
			(lineNumber,intervalStr,inCount) = line

			assert (intervalStr not in componentIntervals), \
			       "(at line %d) interval %s appears more than once in %s" \
			     % (lineNumber,intervalStr,componentFilename)

			if (intervalStr not in combinedIntervals):
				if (warnForExtraIntervals):
					print >>stderr, "WARNING: interval %s appears at line %d in %s but not in the distribution spec" \
					              % (intervalStr,lineNumber,componentFilename)
				continue

			componentIntervals[intervalStr] = interval = Interval()
			interval.inCount = inCount

		componentF.close ()

	componentListF.close ()

	# zero-fill any missing length/intervals in the components

	for name in componentsOrder:
		componentIntervals = componentToIntervals[name]
		for intervalStr in intervalsOrder:
			if (intervalStr not in componentIntervals):
				componentIntervals[intervalStr] = interval = Interval()
				interval.inCount = 0

	# verify the sums for each length/interval

	for intervalStr in intervalsOrder:
		expectedSum = combinedIntervals[intervalStr].inCount
		componentsSum = sum([componentToIntervals[name][intervalStr].inCount
		                     for name in componentsOrder])
		assert (componentsSum == expectedSum), \
		       "sums for %s don't match; expected %d, but sum is %d" \
		     % (intervalStr,expectedSum,componentsSum)

	# spread output distributions over the components

	intervalNum = 0
	for intervalStr in intervalsOrder:
		intervalNum += 1
		if (reportSpreadProgress):
			print >>stderr, "(%d of %d) spreading %s" \
			              % (intervalNum,len(intervalsOrder),intervalStr)
		combinedInterval      = combinedIntervals[intervalStr]
		componentsForInterval = [componentToIntervals[name][intervalStr]
		                         for name in componentsOrder]
		spread_distribution(combinedInterval,componentsForInterval)

	# write distribution spec for each component

	componentNum = 0
	for name in componentsOrder:
		componentNum += 1
		outputFilename = outputFilespec.replace("{component}",name)
		outputF        = file(outputFilename,"wt")

		if (reportOutputWrites):
			print >>stderr, "(%d of %d) writing %s" \
			              % (componentNum,len(componentsOrder),outputFilename)

		for intervalStr in intervalsOrder:
			interval = componentToIntervals[name][intervalStr]
			if (interval.inCount == 0): continue
			print >>outputF, "%s\t%s\t%s" \
			         % (intervalStr,interval.outCount,interval.inCount)

		outputF.close ()


# spread_distribution--
#   combinedInterval      is an Interval with inCount and outCount
#   componentsForInterval is a list of Intervals with inCount

def spread_distribution(combinedInterval,componentsForInterval):
	population = []
	for component in componentsForInterval:
		population += [component] * component.inCount
		component.outCount = 0

	for component in random_sample(population,combinedInterval.outCount):
		component.outCount += 1


# read_distribution_spec--
#	yield (lineNumber,intervalStr,outCount,inCount) tuples

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
			intervalStr = fields[0]
			outCount    = int(fields[1])
			inCount     = int(fields[2])
			if (outCount > inCount): raise ValueError
		except ValueError:
			assert (False), \
			       "bad length-distribution at %s line %d\n%s" \
			     % (filename,lineNumber,line)

		yield (lineNumber,intervalStr,outCount,inCount)


# read_length_distribution--
#	yield (lineNumber,intervalStr,inCount) tuples

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
			inCount     = int(fields[1])
		except ValueError:
			assert (False), \
			       "bad length-distribution at %s line %d\n%s" \
			     % (filename,lineNumber,line)

		yield (lineNumber,intervalStr,inCount)


if __name__ == "__main__": main()
