#!/usr/bin/env python
"""
Sort the alignments output by Noise Cancelling Repeat Finder.
"""

from sys        import argv,stdin,stdout,stderr,exit
from os         import path as os_path
from ncrf_parse import alignments


def usage(s=None):
	message = """
usage: cat <output_from_NCRF> | ncrf_sort [options]
   --sortby=mratio[-|+]    sort by decreasing or increasing match ratio
                           (by default we sort by decreasing match ratio)
   --sortby=score[-|+]     sort by decreasing or increasing alignment score
   --sortby=match[-|+]     sort by decreasing or increasing alignment match
                           count
   --sortby=length[-|+]    sort by decreasing or increasing length; length is
                           the number of sequence bases aligned
   --sortby=name           sort by sequence name (and position)
   --sortby=position[-|+]  sort by sequence name (and decreasing or increasing
                           position)"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	sortBy     = "mratio-"
	requireEof = True

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg == "--sortby=mratio-") or (arg == "--sortby=mratio"):
			sortBy = "mratio-"
		elif (arg == "--sortby=mratio+"):
			sortBy = "mratio+"
		elif (arg == "--sortby=score-") or (arg == "--sortby=score"):
			sortBy = "score-"
		elif (arg == "--sortby=score+"):
			sortBy = "score+"
		elif (arg == "--sortby=match-") or (arg == "--sortby=match"):
			sortBy = "match-"
		elif (arg == "--sortby=match+"):
			sortBy = "match+"
		elif (arg == "--sortby=length-") or (arg == "--sortby=length"):
			sortBy = "length-"
		elif (arg == "--sortby=length+"):
			sortBy = "length+"
		elif (arg == "--sortby=name") or (arg == "--sortby=name+"):
			sortBy = "name,pos+"
		elif (arg == "--sortby=position") or (arg == "--sortby=position+") \
		  or (arg == "--sortby=pos")      or (arg == "--sortby=pos+"):
			sortBy = "name,pos+"
		elif (arg == "--sortby=position-") or (arg == "--sortby=pos-"):
			sortBy = "name,pos-"
		elif (arg in ["--noendmark]","--noeof","--nomark"]):   # (unadvertised)
			requireEof = False
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# collect the alignments

	alignmentList = []

	for a in alignments(stdin,requireEof):
		if (sortBy == "mratio-"):
			key = -a.mRatio
		elif (sortBy == "mratio+"):
			key = a.mRatio
		elif (sortBy == "score-"):
			key = -a.score
		elif (sortBy == "score+"):
			key = a.score
		elif (sortBy == "match-"):
			key = -a.nMatch
		elif (sortBy == "match+"):
			key = a.nMatch
		elif (sortBy == "length-"):
			key = -a.seqBaseCount
		elif (sortBy == "length+"):
			key = a.seqBaseCount
		elif (sortBy == "name,pos+"):
			key = (name_particle(a.seqName),a.start,a.end)
		elif (sortBy == "name,pos-"):
			key = (name_particle(a.seqName),-a.start,-a.end)
		else:
			exit("%s: internal error: unknown key \"%s\""
			   % (os_path.basename(argv[0]),sortBy))
		alignmentList += [(key,a)]

	# sort and print them

	alignmentList.sort()

	isFirst = True
	for (_,a) in alignmentList:
		if (isFirst): isFirst = False
		else:         print
		print "\n".join(a.lines)

	if (requireEof):
		print "# ncrf end-of-file"


# name_particle--
#	Split a sequence name into parts, so that sorting will produce a saner
#	result when names have numeric parts.
#
#	For example, "SRR2036394.36267" is returned as ("SRR",2036394,".",36267).

digits  = "0123456789"
letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
letters += letters.lower()

def name_particle(s):
	particles = []
	particle = []
	pState = None
	for ch in s:
		if   (ch in digits):  chState = "number"
		elif (ch in letters): chState = "letters"
		else:                 chState = "puntuation"

		if (pState == chState):
			particle += [ch]
			continue

		if (particle != []):
			particle = "".join(particle)
			if (pState == "number"): particle = int(particle)
			particles += [particle]

		particle = [ch]
		pState = chState

	if (particle != []):
		particle = "".join(particle)
		if (pState == "number"): particle = int(particle)
		particles += [particle]

	return tuple(particles)


if __name__ == "__main__": main()
