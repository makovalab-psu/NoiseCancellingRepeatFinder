#!/usr/bin/env python
"""
Parse alignments output by Noise Cancelling Repeat Finder.
"""

from sys    import argv,exit
from os     import path as os_path
from re     import compile as re_compile
from string import maketrans
from math   import ceil


class Alignment(object):

	def __init__(self):
		pass

	def set_field_widths(self,nameFieldW,lengthFieldW,countFieldW,rangeFieldW):
		self.nameFieldW   = nameFieldW
		self.lengthFieldW = lengthFieldW
		self.countFieldW  = countFieldW
		self.rangeFieldW  = rangeFieldW

	def __str__(self):
		# we assume that self.name, self.start, and self.end may have been
		# changed

		if (self.errorText != None):
			(statsLine,seqLine,qryLine) = self.lines[:3]
			extraLines = self.lines[3:]
		else:
			statsLine = None
			(seqLine,qryLine) = self.lines[:2]
			extraLines = self.lines[2:]

		(seqName,seqLength,seqBaseCount,seqRange,seqText) = seqLine.split()
		(motifName,qryBaseCount,score,qryText) = qryLine.split()

		seqName  = self.seqName
		seqRange = "%d-%d" % (self.start,self.end)
		seqLength = "?" if (self.seqLen == None) else str(self.seqLen)

		nameW   = max(len(seqName),len(motifName))
		lengthW = len(seqLength)
		countW  = max(len(seqBaseCount),len(qryBaseCount))
		rangeW  = max(len(seqRange),len(score))

		if (hasattr(self,"nameFieldW")):   nameW   = max(nameW,  self.nameFieldW)
		if (hasattr(self,"lengthFieldW")): lengthW = max(lengthW,self.lengthFieldW)
		if (hasattr(self,"countFieldW")):  countW  = max(countW, self.countFieldW)
		if (hasattr(self,"rangeFieldW")):  rangeW  = max(rangeW, self.rangeFieldW)

		if (statsLine != None):
			fields = statsLine.split()
			stats     = " ".join(fields[:8])
			statsText = fields[8]

			if (len(stats) > nameW+lengthW+countW+rangeW+3):
				nameW = len(stats) - (lengthW+countW+rangeW+3)

			statsW = nameW+lengthW+countW+rangeW+3

		lines = []
		if (statsLine != None):
			lines += ["%-*s %s" % (statsW,stats,statsText)]

		lines += ["%-*s %-*s %-*s %-*s %s" \
		        % (nameW,   seqName,
		           lengthW, seqLength,
		           countW,  seqBaseCount,
		           rangeW,  seqRange,
		           seqText)]

		lines += ["%-*s %-*s %-*s %-*s %s" \
		        % (nameW,   motifName,
		           lengthW, "",
		           countW,  qryBaseCount,
		           rangeW,  score,
		           qryText)]

		return "\n".join(lines+extraLines)

	def positional_stats_indexes(self):
		# identify a single range of lines, *all* starting with "# position "

		startIx = endIx = None
		for (ix,line) in enumerate(self.lines):
			if (line.startswith("# position ")):
				if (startIx == None):
					startIx = ix
				elif (endIx != None):
					raise ValueError, \
					      "non-consecutive positional information for alignment at line %d:\n%s" \
					    % (self.lineNumber,line)
			else:
				if (startIx != None):
					endIx = ix

		if (startIx == None):
			raise ValueError
		elif (endIx == None):
			endIx = len(self.lines)

		return (startIx,endIx)

	def positional_stats(self):
		try:
			(startIx,endIx) = self.positional_stats_indexes()
		except ValueError:
			return None

		# parse the lines

		positionalStats = [None] * (endIx-startIx)

		for ix in xrange(startIx,endIx):
			line = self.lines[ix]
			expectedPos = ix - startIx
			try:
				(pos,stats) = parse_positional_stats(line)
			except ValueError,ex:
				raise ValueError, \
				      "%s\ncan't parse positional stats at line %d:\n%s" \
				    % (str(ex),self.lineNumber,line)
			if (pos != expectedPos):
				raise ValueError, \
				      "positional information out of order for alignment at line %d:\n%s" \
				    % (self.lineNumber,line)
			positionalStats[pos] = stats

		return positionalStats


def alignments(f,requireEof=True):
	a = None

	eofMarkerSeen = False

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()

		if (eofMarkerSeen) and (line != ""):
			exit("%s: alignment input contains additional stuff after end marker (starting with \"%s\")"
			   % (os_path.basename(argv[0]),line[:10]))

		if (line == "# ncrf end-of-file"):
			eofMarkerSeen = True
			continue

		if (line == ""):
			if (a != None):
				if (lineInBlock <= 3):
					exit("%s: incomplete alignment for block starting at line %d\n%s"
					   % (os_path.basename(argv[0]),a.lineNumber,a.lines[-1]))
				yield a
				a = None
			continue

		if (line.startswith("#")):
			fields = line.split()
			if (len(fields) < 2) or (not fields[1].startswith("score=")):
				if (a != None):a.lines += [line]
				continue

		if (a == None):
			a = Alignment()
			a.lineNumber = lineNumber
			a.lines = [line]
			lineInBlock =  1
		else:
			a.lines += [line]

		fields = line.split()

		if (lineInBlock == 1):
			if (len(fields) == 9) and (fields[0] == "#"):
				a.mRatio      = parse_mRatio(fields[3])
				a.nMatch      = parse_nMatch(fields[4])
				a.nMismatch   = parse_nMismatch(fields[5])
				a.nInsertions = parse_nInsertions(fields[6])
				a.nDeletions  = parse_nDeletions(fields[7])
				a.errorText   = fields[8]
				lineInBlock   = 2
				continue
			else:
				a.errorText = None
				lineInBlock = 2

		if (lineInBlock == 2):
			if (len(fields) != 5):
				exit("%s: wrong number of columns for block line 2 at line %d (%d, expected %d)\n%s"
				   % (os_path.basename(argv[0]),lineNumber,len(fields),5,line))
			try:
				a.seqName       = fields[0]
				a.seqLen        = None if (fields[1] == "?") else int(fields[1])
				a.seqBaseCount  = parse_base_count(fields[2])
				(a.start,a.end) = parse_interval(fields[3])
				a.seqText       = fields[4]
			except ValueError:
				exit("%s: failed to parse alignment block line 2 at line %d\n%s"
				   % (os_path.basename(argv[0]),lineNumber,line))
			lineInBlock = 3
			continue

		if (lineInBlock == 3):
			if (len(fields) != 4):
				exit("%s: wrong number of columns for alignment block line 3 at line %d (%d, expected %d)\n%s"
				   % (os_path.basename(argv[0]),lineNumber,len(fields),4,line))
			try:
				(a.motif,a.strand) = parse_stranded_motif(fields[0])
				a.motifBaseCount   = parse_base_count(fields[1])
				a.score            = parse_score(fields[2])
				a.motifText        = fields[3]
			except ValueError:
				exit("%s: failed to parse alignment block line 3 at line %d\n%s"
				   % (os_path.basename(argv[0]),lineNumber,line))
			lineInBlock = 4
			continue

	if (a != None):
		if (lineInBlock <= 3):
			exit("%s: incomplete alignment for block starting at line %d"
			   % (os_path.basename(argv[0]),a.lineNumber))
		yield a

	if (requireEof) and (not eofMarkerSeen):
		exit("%s: alignment input may have been truncated (end marker is absent)"
		   % os_path.basename(argv[0]))


# field parsers--

def parse_stranded_motif(s):
	if (len(s) == 0): raise ValueError
	strand = s[-1]
	if (strand in ["+","-"]): s = s[:-1]
	else:                     strand = "+"
	return (s,strand)

def parse_base_count(s):
	baseCountRe = re_compile("^(?P<basecount>[0-9]+)bp$")
	m = baseCountRe.match(s)
	if (m == None): raise ValueError(s)
	return int(m.group("basecount"))

def parse_interval(s):
	intervalRe = re_compile("^(?P<start>[0-9]+)-(?P<end>[0-9]+)$")
	m = intervalRe.match(s)
	if (m == None): raise ValueError
	start = int(m.group("start"))
	end   = int(m.group("end"))
	return (start,end)

def parse_score(s):
	scoreRe = re_compile("^score=(?P<score>[-]?[0-9]+)$")
	m = scoreRe.match(s)
	if (m == None): raise ValueError(s)
	return int(m.group("score"))

def parse_mRatio(s):
	mRatioRe = re_compile("^mRatio=(?P<mRatio>[0-9.]+)%$")
	m = mRatioRe.match(s)
	if (m == None): raise ValueError(s)
	return float(m.group("mRatio")) / 100

def parse_nMatch(s):
	nMatchRe = re_compile("^m=(?P<nMatch>[0-9]+)$")
	m = nMatchRe.match(s)
	if (m == None): raise ValueError(s)
	return int(m.group("nMatch"))

def parse_nMismatch(s):
	nMismatchRe = re_compile("^mm=(?P<nMismatch>[0-9]+)$")
	m = nMismatchRe.match(s)
	if (m == None): raise ValueError(s)
	return int(m.group("nMismatch"))

def parse_nInsertions(s):
	nInsertionsRe = re_compile("^i=(?P<nInsertions>[0-9]+)$")
	m = nInsertionsRe.match(s)
	if (m == None): raise ValueError(s)
	return int(m.group("nInsertions"))

def parse_nDeletions(s):
	nDeletionsRe = re_compile("^d=(?P<nDeletions>[0-9]+)$")
	m = nDeletionsRe.match(s)
	if (m == None): raise ValueError(s)
	return int(m.group("nDeletions"))

def parse_positional_stats(s):
	# example input:
	#   # position 0 [G] mRatio=78.1% m=2305 mm=490 i=114 d=41 x=645
	positionalRe = re_compile("^# position +(?P<pos>[0-9]+) +\[.\](?P<stats>.*)$")
	statRe       = re_compile("^(?P<name>[^=]*)=(?P<val>.*)$")
	m = positionalRe.match(s)
	if (m == None): raise ValueError(s)
	pos = int(m.group("pos"))
	statsText = m.group("stats")
	stats = {}
	for field in statsText.split():
		m = statRe.match(field)
		if (m == None): raise ValueError("%s (field \"%s\")" % (s,field))
		name = m.group("name")
		val  = m.group("val")
		if (val.endswith("%")): val = float(val[:-1]) / 100.0
		else:                   val = int(val)
		stats[name] = val
	return (pos,stats)


# reverse_complement--

complementMap = maketrans("ACGTSWRYMKBDHVNacgtswrymkbdhvn",
                          "TGCASWYRKMVHDBNtgcaswyrkmvhdbn")

def reverse_complement(nukes):
	return nukes[::-1].translate(complementMap)


# parse_probability--
#	Parse a string as a probability

def parse_probability(s,strict=True):
	scale = 1.0
	if (s.endswith("%")):
		scale = 0.01
		s = s[:-1]

	try:
		p = float(s)
	except:
		try:
			(numer,denom) = s.split("/",1)
			p = float(numer)/float(denom)
		except:
			raise ValueError

	p *= scale

	if (strict) and (not 0.0 <= p <= 1.0):
		raise ValueError

	return p


# parse_noise_rate--
#	Parse a string as noise rate

def parse_noise_rate(s):
	scale = 1.0
	if (s.endswith("%")):
		scale = 0.01
		s = s[:-1]

	try:
		p = float(s)
	except:
		try:
			(numer,denom) = s.split("/",1)
			p = float(numer)/float(denom)
		except:
			raise ValueError

	return p*scale


# float_or_fraction--
#	Parse a string as a floating point value, allowing fractions

def float_or_fraction(s):
	if ("/" in s):
		(numer,denom) = s.split("/",1)
		return float(numer)/float(denom)
	else:
		return float(s)


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return          int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))


# commatize--
#	Convert a numeric string into one with commas.

def commatize(s):
	if (type(s) != str): s = str(s)
	(prefix,val,suffix) = ("",s,"")
	if (val.startswith("-")): (prefix,val) = ("-",val[1:])
	if ("." in val):
		(val,suffix) = val.split(".",1)
		suffix = "." + suffix

	try:    int(val)
	except: return s

	digits = len(val)
	if (digits > 3):
		leader = digits % 3
		chunks = []
		if (leader != 0):
			chunks += [val[:leader]]
		chunks += [val[ix:ix+3] for ix in xrange(leader,digits,3)]
		val = ",".join(chunks)

	return prefix + val + suffix
