#!/usr/bin/env python
"""
Accumulate per-position counts into binned intervals.

References:
  [1] Using eval() safely in python (lybniz2.sourceforge.net/safeeval.html)
"""

from sys           import argv,stdin,stdout,stderr,exit
from math          import *
from interval_dict import IntervalDict,Interval


def usage(s=None):
	message = """

usage: cat <counts_file> | bin_position_counts [options]
  bin=<function>     (required) function that maps a position to its bin;
                     positions that map to the same unit interval are in the
                     same bin; an example is "pos/10"
  pos=<function>     (required) function that maps a bin to its position; this
                     should be the mathematical inverse of the to-bin function;
                     an example is "10*bin"
  --minpos=<number>  positions lower than this are ignored
  --maxpos=<number>  positions higher than this are ignored

The counts file consists of two columns-- a position and a count for that
position.  For example, this set of position,count pairs
   375 15    381 26    387 28    393 21    399 23
   376 17    382 27    388 23    394 20    400 27
   377 21    383 26    389 18    395 21    401 24
   378 26    384 23    390 27    396 22    402 24
   379 21    385 14    391 30    397 20    403 28
   380 21    386 19    392 21    398 25    404 30
would be represented in 30 lines in the file:
   375 15
    ...
   404 30

With the binning function bin=pos/10 and its inverse pos=10*bin, these counts
are collected into bins like this. Note that the bins may cover some positions
not present in the input (but no empty bins will be output).
   370-379 100
   380-389 225
   390-399 230
   400-409 133

The reason that 370-379 forms a bin is that the function bin=pos/10 maps each
of them into the unit interval 37 <= bin < 38
   370  37.0
   371  37.1
    ...
   378  37.8
   379  37.9

The bin mapping function should be strictly increasing over the range of
positions presented. It can be more sophisticated than demonstrated above. For
example, bin=30.5*log(pos-250) and inverse pos=250+exp(bin/30.5) will group
small values (pos about 400) into bins of size 5 and large values (pos about
300,000) into much larger bins."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global toBinFunction,fromBinFunction
	global debug

	# parse the command line

	toBinFunction   = None
	fromBinFunction = None
	minPosition     = None
	maxPosition     = None
	debug           = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--tobin=")) or (arg.startswith("bin=")):
			toBinFunction = argVal
		elif (arg.startswith("--frombin=")) or (arg.startswith("pos=")):
			fromBinFunction = argVal
		elif (arg.startswith("--minpos=")) or (arg.startswith("--min=")):
			minPosition = int(argVal)
		elif (arg.startswith("--maxpos=")) or (arg.startswith("--max=")):
			maxPosition = int(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (toBinFunction == None):
		usage("you must provide a position-to-bin expression")
	if (fromBinFunction == None):
		usage("you must provide a bin-to-position expression")

	if ("pos" not in toBinFunction):
		usage("\"%s\" doesn't contain the variable \"pos\"" % toBinFunction)

	if ("bin" not in fromBinFunction):
		usage("\"%s\" doesn't contain the variable \"bin\"" % fromBinFunction)

	# accumulate the counts

	intervals = IntervalDict()
	functionsAreVerified = False

	for (pos,count) in read_position_counts(stdin):
		if (minPosition != None) and (pos < minPosition): continue
		if (maxPosition != None) and (pos > maxPosition): continue

		if ("tobin" in debug):
				context = dict(safeDict)
				context["pos"] = float(pos)
				bin = eval(toBinFunction,{"__builtins__":None},context)
				print >>stderr, "pos=%s --> bin=%s" % (pos,bin)

		if (not functionsAreVerified):
			assert (verify_functions(pos)), \
			       "\"%s\" doesn't appear to be the inverse of \"%s\"" \
			     % (fromBinFunction,toBinFunction)
			functionsAreVerified = True

		try:
			interval = intervals[pos]
			interval.count += count
		except KeyError:
			(bin,loPos,hiPos) = position_to_bin_interval(pos)
			if (minPosition != None) and (loPos < minPosition):
				loPos = minPosition
			if (maxPosition != None) and (hiPos > maxPosition):
				hiPos = maxPosition

			if (not loPos <= pos <= hiPos):
				interval = Interval(loPos,hiPos)
				assert (False), \
				       "(internal?) error: interval %s (bin %d) does not contain pos %d" \
				     % (interval,bin,pos)

			interval = intervals.add(loPos,hiPos)
			if (interval == None): # interval overlaps an existing interval
				interval = Interval(loPos,hiPos)
				previous = intervals.overlapper(loPos,hiPos)
				assert (False), \
				       "(internal?) error: interval %s (pos %d, bin %d) overlaps %s (pos %d, bin %d)" \
				     % (interval,pos,bin,previous,previous.pos,previous.bin)

			interval.pos   = pos
			interval.bin   = bin
			interval.count = count

	# report the counts

	for interval in intervals:
		print "%s %d" % (interval,interval.count)


# verify_functions--
#   verify that fromBinFunction is the inverse of toBinFunction

def verify_functions(pos,epsilon=1e-6):
	context = dict(safeDict)
	context["pos"] = float(pos)
	bin = eval(toBinFunction,{"__builtins__":None},context)

	context = dict(safeDict)
	context["bin"] = float(bin)
	computedPos = eval(fromBinFunction,{"__builtins__":None},context)

	return (abs(computedPos-pos) <= epsilon)


# position_to_bin_interval--

def position_to_bin_interval(pos):
	# map the position to its bin

	context = dict(safeDict)
	context["pos"] = float(pos)
	bin = eval("floor(%s)" % toBinFunction,{"__builtins__":None},context)

	# map the bin to its lower end and higher end

	context = dict(safeDict)
	context["bin"] = float(bin)
	loPos = eval("ceil(%s)" % fromBinFunction,{"__builtins__":None},context)

	context = dict(safeDict)
	context["bin"] = float(bin+1)
	hiPos = eval("floor(%s)" % fromBinFunction,{"__builtins__":None},context)

	# map the bin's higher end back to a bin; if it's not in the right bin
	# reduce the higher end
	#
	# this tests for the case where a bin boundary falls on an integer; the
	# math above would erroneously assign that integer to both intervals

	context = dict(safeDict)
	context["pos"] = float(hiPos)
	hiBin = eval("floor(%s)" % toBinFunction,{"__builtins__":None},context)
	if (hiBin > bin): hiPos -= 1

	return (bin,loPos,hiPos)


# read_position_counts--
#	yield (position,count) pairs

def read_position_counts(f,filename=None):
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
			pos   = int(fields[0])
			count = int(fields[1])
		except ValueError:
			assert (False), \
			       "bad position,count at %s line %d\n%s" \
			     % (filename,lineNumber,line)

		yield (pos,count)


# (see reference [1], lybniz2.sourceforge.net/safeeval.html)

def round_int(x): return int(round(x))

def floor_int(x): return int(floor(x))

def ceil_int(x):  return int(ceil(x))

def log2(x):      return log(x) / log(2.0)


safeList = ["e", "exp", "log", "log10", "pi", "pow", "sin", "sqrt"]
safeDict = dict([(k,locals().get(k,None)) for k in safeList])
safeDict["int"]     = int
safeDict["float"]   = float
safeDict["ceil"]    = ceil_int
safeDict["floor"]   = floor_int
safeDict["round"]   = round_int
safeDict["log2"]    = log2
safeDict["max"]     = max
safeDict["min"]     = min


if __name__ == "__main__": main()
