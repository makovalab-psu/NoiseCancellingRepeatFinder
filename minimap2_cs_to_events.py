#!/usr/bin/env python
"""
Convert the cs tag in minimap2 output to ncrf-style event counts .
"""

from sys import argv,stdin,stdout,stderr,exit
from os  import path as os_path

class Alignment: pass


def usage(s=None):
	message = """
usage: cat <output_from_minimap2> | minimap2_cs_to_events [options]
  --minquality=<qual>  discard low quality alignments
  --withheader         include a header line in the output
  --remove:cs          remove the cs tag
  --remove:tags        remove the all tags

The minimap2 output should include the cs tag, i.e. minimap2 should have been
run with the "--cs=short" option."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	minQuality    = None
	writeHeader   = False
	removeCsTag   = False
	removeAllTags = False

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--minquality=")) or (arg.startswith("--minqual=")):
			minQuality = int(argVal)
		elif (arg in ["--withheader","--with=header","--with:header"]):
			writeHeader = True
		elif (arg == "--remove:cs"):
			removeCsTag = True
		elif (arg == "--remove:tags"):
			removeAllTags = True
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# process the alignments

	for a in minimap2_alignments(stdin):
		if (minQuality != None) and (a.mappingQuality < minQuality):
			continue

		if (writeHeader):
			fields = ["#qName","qLen","qStart","qEnd","s",
			          "tName","tLen","tStart","tEnd",
			          "matches","events","qual",
			          "mRatio","m","mm","io","ix","do","dx"]
			if (not removeAllTags): fields += ["tags"]
			print "\t".join(fields)
			writeHeader = False

		mRatio = m = mm = io = ix = do = dx = "NA"
		if ("cs" in a.tags):
			(s,cs) = a.tags["cs"]
			(m,mm,io,ix,do,dx) = cs_to_events(cs)
			mRatio = "%.1f%%" % (100.0*m / (m+mm+io+ix+do+dx))

		fields = [a.queryName,a.queryLen,a.queryStart,a.queryEnd,a.strand,
		          a.targetName,a.targetLen,a.targetStart,a.targetEnd,
		          a.nMatches,a.nEvents,a.mappingQuality,
		          mRatio,m,mm,io,ix,do,dx]

		if (not removeAllTags):
			for tag in a.tagOrder:
				if (removeCsTag) and (tag == "cs"): continue
				(s,val) = a.tags[tag]
				fields += [s]

		print "\t".join(map(str,fields))		


def minimap2_alignments(f):
	a = None

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == "") or (line.startswith("#")):
			continue

		fields = line.split()
		if (len(fields) < 12):
			exit("%s: wrong number of columns at line %d (%d, expected %d)\n%s"
			   % (os_path.basename(argv[0]),lineNumber,len(fields),12,line))

		a = Alignment()

		try:
			a.queryName      =     fields[0]
			a.queryLen       = int(fields[1])
			a.queryStart     = int(fields[2])
			a.queryEnd       = int(fields[3])
			a.strand         =     fields[4]
			a.targetName     =     fields[5]
			a.targetLen      = int(fields[6])
			a.targetStart    = int(fields[7])
			a.targetEnd      = int(fields[8])
			a.nMatches       = int(fields[9])
			a.nEvents        = int(fields[10])
			a.mappingQuality = int(fields[11])

			a.tagOrder = []
			a.tags = {}
			for s in fields[12:]:
				(tag,kind,val) = s.split(":",2)
				if (tag in a.tags): raise ValueError
				if   (kind == "i"): val = int(val)
				elif (kind == "f"): val = float(val)
				a.tagOrder += [tag]
				a.tags[tag] = (s,val)

		except ValueError:
			exit("%s: failed to parse alignment at line %d\n%s"
			   % (os_path.basename(argv[0]),lineNumber,line))

		yield a


def cs_to_events(cs):
	#print >>stderr, "\"%s\"" % cs
	m = mm = io = ix = do = dx = 0

	matchStart = state = None
	for (i,ch) in enumerate(cs):
		if (ch == ":"):
			if (state in ["io","do","*"]): raise ValueError
			if (matchStart != None):
				if (matchStart == i): raise ValueError
				m += int(cs[matchStart:i])
			matchStart = i+1
		elif (ch.isdigit()):
			pass
		elif (ch in "+-*"):
			if (matchStart != None):
				if (matchStart == i): raise ValueError
				m += int(cs[matchStart:i])
				matchStart = None
			if   (ch == "+"): state = "io"
			elif (ch == "-"): state = "do"
			else:             state = "*"
		elif (ch in "acgt"):
			if   (state == "io"): (io,state) = (io+1,"ix")
			elif (state == "ix"): ix += 1
			elif (state == "do"): (do,state) = (do+1,"dx")
			elif (state == "dx"): dx += 1
			elif (state == "*"):  state = "mm"
			elif (state == "mm"): (mm,state) = (mm+1,None)
		else:
			raise ValueError

	if (state in ["io","do","*"]): raise ValueError
	if (matchStart != None):
		if (matchStart == len(cs)): raise ValueError
		m += int(cs[matchStart:])

	return (m,mm,io,ix,do,dx)


if __name__ == "__main__": main()
