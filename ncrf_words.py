#!/usr/bin/env python
"""
Look for prominent "wrong" motifs (words) in the output of Noise Cancelling
Repeat Finder.
"""

from sys         import argv,stdin,stdout,stderr,exit
from collections import Counter
from ncrf_parse  import alignments,reverse_complement,int_with_unit,float_or_fraction


def usage(s=None):
	message = """
usage: cat <output_from_NCRF> | ncrf_words [options]
  --minwordratio=<r>  only show words that have counts that are at least r
                      times the motif word's count (e.g. r=0.5 would show the
                      words that occur at least half as often as the motif)
                      (default is 1.0)
  --head=<number>     limit the number of input alignments"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse the command line

	countRatio = 1
	headLimit  = None
	requireEof = True
	debug      = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--minwordratio=")) or (arg.startswith("--ratio=")) or (arg.startswith("R=")):
			countRatio = float_or_fraction(argVal)
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg in ["--noendmark","--noeof","--nomark"]):   # (unadvertised)
			requireEof = False
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# process the alignments

	alignmentNum = 0
	for a in alignments(stdin,requireEof):
		alignmentNum += 1 

		if (headLimit != None) and (alignmentNum > headLimit):
			print >>stderr, "limit of %d alignments reached" % headLimit
			break

		if (alignmentNum > 1): print
		print "\n".join(a.lines)

		motifText = a.motifText
		seqText   = a.seqText
		if ("noflip" in debug):
			pass
		elif (a.strand == "-") and (a.start < a.end):
			# alignment was reported in reverse complement of motif, so flip it
			motifText = reverse_complement(motifText)
			seqText   = reverse_complement(seqText)

		(motifChunks,seqChunks) = chunkify(a.motif,motifText,seqText)

		wordCounts = Counter()
		for word in seqChunks:
			word = word.replace("-","")
			if (word != a.motif): word = word.lower()
			wordCounts[word] += 1

		if (a.motif in wordCounts): motifCount = wordCounts[a.motif]
		else:                       motifCount = 0

		wordCounts = [(wordCounts[word],abs(len(word)-len(a.motif)),word)
		                for word in wordCounts
		                if (wordCounts[word] >= motifCount * countRatio)]
		wordCounts.sort()
		wordCounts.reverse()
		print "# aligned words %s" % \
		      " ".join(["%s:%d"%(word,count) for (count,_,word) in wordCounts])

		if ("chunks" in debug):
			if ("noflip" in debug):
				seqChunks   = [reverse_complement(word) for word in seqChunks  [::-1]]
				motifChunks = [reverse_complement(word) for word in motifChunks[::-1]]
			print "# words: %s" % " ".join(seqChunks)
			print "# motif: %s" % " ".join(motifChunks)

	if (requireEof):
		print "# ncrf end-of-file"


def chunkify(motif,motifText,seqText):
	motifText = motifText.upper()
	(pos,direction) = motif_position(motif,motifText)

	chunks = []
	start = ix = pos
	motifIx = 0
	while (ix < len(motifText)):
		textNuc = motifText[ix]
		ix += 1
		if (textNuc == "-"): continue
		motifIx += 1
		if (motifIx == len(motif)):
			chunks += [(start,ix)]
			start = ix
			motifIx = 0

	motifChunks = []		
	for (start,end) in chunks:
		motifChunks += [motifText[start:end]]

	seqChunks = []		
	for (start,end) in chunks:
		seqChunks += [seqText[start:end]]

	return (motifChunks,seqChunks)


def motif_position(motif,text):
	for (start,startNuc) in enumerate(text):
		if (startNuc == '-'): continue
		pos = start
		motifIx = 0
		isMatch = True
		while (pos < len(text)):
			nuc = text[pos]
			pos += 1
			if (nuc == '-'): continue
			if (nuc != motif[motifIx]):
				isMatch = False
				break
			motifIx += 1
			if (motifIx == len(motif)):
				break
		if (isMatch):
			return (start,1)

	motif = reverse_complement(motif)
	for (start,startNuc) in enumerate(text):
		if (startNuc == '-'): continue
		pos = start
		motifIx = 0
		isMatch = True
		while (pos < len(text)):
			nuc = text[pos]
			pos += 1
			if (nuc == '-'): continue
			if (nuc != motif[motifIx]):
				isMatch = False
				break
			motifIx += 1
			if (motifIx == len(motif)):
				break
		if (isMatch):
			return (start,-1)

	# no match found, just report position 0, forward
	return (0,1)


if __name__ == "__main__": main()
