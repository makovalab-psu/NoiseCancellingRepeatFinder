#!/usr/bin/env python
"""
Filter Noise Cancelling Repeat Finder alignments, discarding alignments that
has a consensus different than the motif unit.
"""

from sys        import argv,stdin,stdout,stderr,exit
from math       import log,exp
from ncrf_parse import alignments,reverse_complement,canonical_motif, \
                       parse_probability,int_with_unit,commatize


def usage(s=None):
	message = """
usage: ncrf_cat <output_from_NCRF> | ncrf_consensus_filter [options]
  --consensusonly     just report the consensus motif(s) for each alignment,
                      instead of filtering; these are added to the alignment
                      file with a "# consensus" tag
  --head=<number>     limit the number of input alignments
  --progress=<number> periodically report how many alignments we've tested"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global headLimit,reportProgress,requireEof
	global winnerThreshold,filterToKeep,reportConsensus,reportMsa
	global canonicalizeConsensuses
	global debug

	canonicalizeConsensuses = True

	# parse the command line

	filterToKeep    = "consensus"
	reportConsensus = False
	reportMsa       = False
	winnerThreshold = 0.50  # (see derive_consensuses)
	sliceWidth      = None
	sliceStep       = None
	headLimit       = None
	reportProgress  = None
	requireEof      = True
	debug           = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg == "--nonconsensus"):             # (unadvertised)
			filterToKeep    = "non consensus"
			reportMsa       = False
			reportConsensus = True
		elif (arg == "--nonconsensus,msa"):       # (unadvertised)
			filterToKeep    = "non consensus"
			reportMsa       = True
			reportConsensus = True
		elif (arg == "--consensusonly"):
			filterToKeep    = "no filter"
			reportMsa       = False
			reportConsensus = True
		elif (arg == "--filter,consensus"):       # (unadvertised)
			filterToKeep    = "consensus"
			reportMsa       = False
			reportConsensus = True
		elif (arg == "--msa"):                    # (unadvertised)
			filterToKeep    = "no filter"
			reportMsa       = True
			reportConsensus = True
		elif (arg.startswith("--winner=")) or (arg.startswith("W=")):   # (unadvertised)
			winnerThreshold = parse_probability(argVal)
		elif (arg.startswith("--slice=")):        # (unadvertised)
			if ("by" in argVal):
				(sliceWidth,sliceStep) = map(int_with_unit,argVal.split("by",1))
			else:
				sliceWidth = sliceStep = int_with_unit(argVal)
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
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

	if (sliceWidth == None):
		simple_consensus_filter(stdin)
	else:
		sliced_consensus_filter(stdin,sliceWidth,sliceStep)


# simple_consensus_filter--

def simple_consensus_filter(f):
	alignmentNum = 0
	alignmentsWritten = 0
	for a in alignments(f,requireEof):
		alignmentNum += 1 

		if (reportProgress != None):
			if (alignmentNum == 1) or (alignmentNum % reportProgress == 0):
				print >>stderr, "progress: testing alignment %s" \
				              % commatize(alignmentNum)

		if (headLimit != None) and (alignmentNum > headLimit):
			print >>stderr, "limit of %d alignments reached" % headLimit
			break

		motifText = a.motifText
		seqText   = a.seqText
		if ("noflip" in debug):
			pass
		elif (a.strand == "-") and (a.start < a.end):
			# alignment was reported in reverse complement of motif, so flip it
			motifText = reverse_complement(motifText)
			seqText   = reverse_complement(seqText)

		# derive consensus(es)

		seqChunks = chunkify(a.motif,motifText,seqText)

		if ("consensus" in debug):
			print >>stderr
			print >>stderr, "%d score=%d" % (a.lineNumber,a.score)

		consensuses = derive_consensuses(seqChunks,winnerThreshold=winnerThreshold)
		consensuses = list(consensuses)

		# discard the alignment if it meets the filtering criterion (if there
		# is any such criterion)

		if (filterToKeep == "consensus"):
			if (a.motif not in consensuses): continue  # (discard it)
		elif (filterToKeep == "non consensus"):
			if (a.motif in consensuses): continue  # (discard it)
		else: # if (filterToKeep == "no filter"):
			pass

		# copy the (unfiltered) alignment to the output

		if (alignmentsWritten > 0): print
		alignmentsWritten += 1

		print "\n".join(a.lines)

		# report the consensus, if we're supposed to

		if (reportConsensus):
			if (consensuses == []):
				print "# consensus (none)"
			else:
				canonicalized = []
				for motif in consensuses:
					if (motif != a.motif) and (canonicalizeConsensuses):
						(motif,strand) = canonical_motif(motif)
					canonicalized += [motif]
				print "# consensus %s" % ",".join(canonicalized)

		# report the MSA from which the consensus was derived, if we're
		# supposed to

		if (reportMsa):
			motifLen = len(a.motif)
			positionLength = [1] * motifLen
			for chunk in seqChunks:
				for (motifIx,seqNucs) in enumerate(chunk):
					if (seqNucs == None): continue
					positionLength[motifIx] = max(positionLength[motifIx],len(seqNucs))

			line = []
			for (motifIx,motifNuc) in enumerate(a.motif):
				line += [motifNuc.ljust(positionLength[motifIx],".")]
			print "# msa.query %s" % "".join(line)

			for chunk in seqChunks:
				line = []
				for (motifIx,seqNucs) in enumerate(chunk):
					if (seqNucs == None):
						line += ["."*positionLength[motifIx]]
					elif (seqNucs == a.motif[motifIx]):
						line += ["="*positionLength[motifIx]]
					else:
						line += [seqNucs.ljust(positionLength[motifIx],".")]
				print "# msa.seq   %s" % "".join(line)

	if (requireEof):
		print "# ncrf end-of-file"


# sliced_consensus_filter--

userHasBeenWarned = False

def sliced_consensus_filter(f,sliceWidth,sliceStep):
	global userHasBeenWarned

	if (reportMsa) and (not userHasBeenWarned):
		print >>stderr, "WARNING: sliced consensus doesn't report MSA, ignoring that request"
		userHasBeenWarned = True

	alignmentNum = 0
	alignmentsWritten = 0
	for a in alignments(f,requireEof):
		alignmentNum += 1 

		if (reportProgress != None):
			if (alignmentNum == 1) or (alignmentNum % reportProgress == 0):
				print >>stderr, "progress: testing alignment %s" \
				              % commatize(alignmentNum)

		if (headLimit != None) and (alignmentNum > headLimit):
			print >>stderr, "limit of %d alignments reached" % headLimit
			break

		motifText = a.motifText
		seqText   = a.seqText
		if ("noflip" in debug):
			pass
		elif (a.strand == "-") and (a.start < a.end):
			# alignment was reported in reverse complement of motif, so flip it
			motifText = reverse_complement(motifText)
			seqText   = reverse_complement(seqText)

		# look for consensus over each slice, separately

		consensuses = set()

		numSlices = (len(motifText) + sliceStep-1) / sliceStep  # (an overestimate)
		minSlice  = 10*len(a.motif)

		for sliceNum in xrange(numSlices):
			sliceStart = sliceNum * sliceStep
			sliceEnd   = min(sliceStart+sliceWidth,len(motifText))
			if (sliceEnd - sliceStart < minSlice): break

			motifTextSlice = motifText[sliceStart:sliceEnd]
			seqTextSlice   = seqText  [sliceStart:sliceEnd]

			# derive consensus(es)

			seqChunks = chunkify(a.motif,motifTextSlice,seqTextSlice)

			if ("consensus" in debug):
				print >>stderr
				print >>stderr, "%d score=%d slice.start=%d slice.end=%d" \
				              % (a.lineNumber,a.score,sliceStart,sliceEnd)

			sliceConsensuses = derive_consensuses(seqChunks,winnerThreshold=winnerThreshold)
			sliceConsensuses = list(sliceConsensuses)
			if (sliceConsensuses == []):
				consensuses.add(None)
			else:
				for word in sliceConsensuses:
					consensuses.add(word)

				if ("consensus" in debug):
					for word in sliceConsensuses:
						print >>stderr, "consensus %s" % word

		consensuses = list(consensuses)

		# discard the alignment if it meets the filtering criterion (if there
		# is any such criterion)

		if (filterToKeep == "consensus"):
			if (a.motif not in consensuses): continue  # (discard it)
		elif (filterToKeep == "non consensus"):
			if (a.motif in consensuses): continue  # (discard it)
		else: # if (filterToKeep == "no filter"):
			pass

		# copy the (unfiltered) alignment to the output

		if (alignmentsWritten > 0): print
		alignmentsWritten += 1

		print "\n".join(a.lines)

		# report the consensus, if we're supposed to

		if (reportConsensus):
			if (consensuses == []):
				print "# consensus (none)"
			else:
				canonicalized = []
				for motif in consensuses:
					if (motif == None): continue
					if (motif != a.motif) and (canonicalizeConsensuses):
						(motif,strand) = canonical_motif(motif)
					canonicalized += [motif]
				if (None in consensuses):
					canonicalized += ["(none)"]
				print "# consensus %s" % ",".join(canonicalized)

	if (requireEof):
		print "# ncrf end-of-file"


# chunkify--
#	Returns a list of sublists. Each sublist is as long as the motif, consisting
#	of the string (the nt or nts) matched to each position in the motif.

def chunkify(motif,motifText,seqText):
	motifText = motifText.upper()
	motifLen  = len(motif)
	(motifPos,direction) = position_in_motif(motif,motifText)
	if ("motifpos" in debug):
		print >>stderr, "motifPos=%d%s" % (motifPos,direction)

	chunks = []
	motifIx = motifPos
	chunk = [None] * motifIx
	for (ix,textNuc) in enumerate(motifText):
		seqNuc = seqText[ix]
		if (textNuc == "-"):     # insertion, extra character in seqText
			if (chunk != []):    # add inserted character to latest token
				if (chunk[-1] != None): chunk[-1] += seqNuc
			elif (chunks != []): # add inserted character to final token in latest chunk
				latestChunk = chunks[-1]
				latestChunk[-1] += seqNuc
				chunks[-1] = latestChunk
			continue

		chunk += [seqNuc]        # match or deletion
		motifIx += 1
		if (motifIx == motifLen):
			chunks += [chunk]
			chunk  =  []
			motifIx = 0

	if (chunk != []):
		if (len(chunk) < motifLen): chunk += [None] * (motifLen - len(chunk))
		chunks += [chunk]

	if (direction == "-"):
		for chunk in chunks:
			chunk.reverse()

	if ("chunks" in debug):
		for (chunkIx,chunk) in enumerate(chunks):
			print >>stderr, "chunk[%d] = %s" % (chunkIx,str(chunk))

	return chunks


# position_in_motif--
#	Returns the index in 'motif' at which the first position of 'text" aligns

def position_in_motif(motif,text):
	motifLen = len(motif)

	for rotIx in xrange(motifLen):
		rotMotif = motif[rotIx:] + motif[:rotIx]
		pos = 0
		motifIx = 0
		isMatch = True
		while (pos < len(text)):
			nuc = text[pos]
			pos += 1
			if (nuc == '-'): continue
			if (nuc != rotMotif[motifIx]):
				isMatch = False
				break
			motifIx += 1
			if (motifIx == motifLen):
				break
		if (isMatch):
			return (rotIx,"+")

	motif = reverse_complement(motifLen)
	for rotIx in xrange(motif):
		rotMotif = motif[rotIx:] + motif[:rotIx]
		pos = 0
		motifIx = 0
		isMatch = True
		while (pos < len(text)):
			nuc = text[pos]
			pos += 1
			if (nuc == '-'): continue
			if (nuc != rotMotif[motifIx]):
				isMatch = False
				break
			motifIx += 1
			if (motifIx == motifLen):
				break
		if (isMatch):
			return (rotIx,"-")

	# no match found, just report position 0,forward

	return (0,"+")


# derive_consensuses--
#	Yields a series of potential consensus motifs for a given alignment. The
#	input is a list of sublists (as created by chunkify).

def derive_consensuses(seqChunks,winnerThreshold=0.50):
	if (seqChunks == []): return
	motifLen = len(seqChunks[0])

	# count the number of times each motif position is observed to match a
	# particular 'token'

	ixToTokens = {}
	for motifIx in xrange(motifLen):
		ixToTokens[motifIx] = {}

	for chunk in seqChunks:
		for (motifIx,seqNucs) in enumerate(chunk):
			if (seqNucs == None): seqNucs = ""
			if (seqNucs not in ixToTokens[motifIx]):
				ixToTokens[motifIx][seqNucs] = 1
			else:
				ixToTokens[motifIx][seqNucs] += 1

	# at each position, sort the tokens from most-observed to least-observed,
	# then reduce them to best (or almost best)

	ixToWinners = {}
	for motifIx in xrange(motifLen):
		tokensSeen = ixToTokens[motifIx]
		ixToTokens[motifIx] = [(ixToTokens[motifIx][seqNucs],seqNucs) for seqNucs in tokensSeen]
		ixToTokens[motifIx].sort()
		ixToTokens[motifIx].reverse()
		tokensCount = sum([count for (count,_) in ixToTokens[motifIx]])

		if ("consensus" in debug):
			(bestCount,_) = ixToTokens[motifIx][0]
			s = []
			for (count,seqNucs) in ixToTokens[motifIx]:
				s += ["%d:\"%s\"" % (count,seqNucs)]
			print >>stderr, "# [%d]%s tokensCount=%d %s" % \
			                (motifIx,
			                 " !!!" if (bestCount < winnerThreshold*tokensCount) else "",
			                 tokensCount," ".join(s))

		ixToWinners[motifIx] = [ixToTokens[motifIx][ix]
		                          for (ix,(count,_)) in enumerate(ixToTokens[motifIx])
		                          if (count >= winnerThreshold*tokensCount)]

	# generate the possible consensus motifs
	#
	# at each position i, we previously determined which token (or tokens) is
	# a clear winner; here we just catenate all these winners to form the
	# consensus; if any position fails to have a winner we report no consensus
	#
	# $$$ change this to include positions that have more than one winner,
	#     reporting all possible paths through the winners list

	motifsReported = set()

	if (winnerThreshold != None):
		motif = []
		for motifIx in xrange(motifLen):
			tokensSeen = ixToWinners[motifIx]
			if (tokensSeen == []):   # (no consensus can be formed, because nothing
				return               #  .. in this column is a clear winner)
			(_,bestSeqNucs) = tokensSeen[0]
			if (bestSeqNucs != "-"): motif += [bestSeqNucs]
		motif = "".join(motif)

		if (motif not in motifsReported):
			yield motif
			motifsReported.add(motif)


if __name__ == "__main__": main()
