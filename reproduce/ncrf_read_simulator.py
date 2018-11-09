#!/usr/bin/env python
"""
Convert a fasta file to sampled "reads" in fasta or fastq format.
"""

from sys        import argv,stdin,stdout,stderr,exit
from math       import log,floor,ceil
from random     import Random,seed as random_seed,random as unit_random,choice as random_choice
from re         import compile
from string     import maketrans
from gzip       import open as gzip_open
from time       import strftime
from prob_table import ProbabilityTable

BAM_FREVERSE = 16   # the read is mapped to the reverse strand

nucToSubstitutions = \
	{
	"A": [    "C","G","T"],
	"C": ["A"    ,"G","T"],
	"G": ["A","C"    ,"T"],
	"T": ["A","C","G"    ]
	}


# error profiles
#
# pacbio   (v3) mm=1.7%  i=8.9%  d=4.3%  is derived from GIAB HG002 blasr alignments
# pacbio   (v1) mm=1.30% i=6.37% d=3.62% is from Guiblet et al (submitted)
# pacbio        mm=1%    i=12%   d=2%    is from readsim-1.6 (ref [1] below)
# nanopore (v3) mm=4.60% i=3.78% d=7.73% is derived from GIAB HG002 minimap2 alignments
# nanopore (v1) mm=7.4%  i=2%    d=10.2% is derived from data from Jain et al (ref [2] below)
# nanopore      mm=3%    i=3%    d=3%    is from readsim-1.6 (ref [1] below)
#
# References
#   [1] readsim-1.6
#         https://sourceforge.net/p/readsim
#   [2] Miten Jain et al. "Nanopore sequencing and assembly of a human genome
#       with ultra-long reads."
#         https://www.biorxiv.org/content/early/2017/04/20/128835
#       Data at
#         https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md
#       Later version, peer-reviewed, at
#         https://www.nature.com/articles/nbt.4060

errorProfilePacbioV3        = {"mm":0.0170, "i":0.0890, "d":0.0430}
errorProfilePacbioV1        = {"mm":0.0130, "i":0.0637, "d":0.0362}
errorProfilePacbioReadsim   = {"mm":0.01,   "i":0.12,   "d":0.02  }
errorProfileNanoporeV3      = {"mm":0.0460, "i":0.0378, "d":0.0773}
errorProfileNanoporeV1      = {"mm":0.074,  "i":0.02,   "d":0.102 }
errorProfileNanoporeReadSim = {"mm":0.03,   "i":0.03,   "d":0.03  }


def usage(s=None):
	message = """
usage: [cat <fasta_file> |] ncrf_read_simulator [options]
  <num>x<length>[,<length>]  number and length of reads to generate; the second
                             length is used for paired reads when the mates
                             have different lengths; if <length> is "stdin",
                             lengths are read from stdin
  [<weight>:]<fasta_file>    (cumulative) sample reads from the specified
                             genome file; if no files are given, the genome is
                             read from stdin; <weights> can be used to
                             control the mixture of reads from different
                             genomes, and are internally scaled by the length
                             of the corresponding genome; if the weight is
                             absent, it is 1.0 by default
  --insert=<avg>,<stdev>     length of inserts for paired reads
                             (by default, reads are not paired)
  --orientation=<T2T|H2H|..> orientation for paired reads; the orientation
                             can be one of the following:
                               T2T, RF, PE (these are equivalent)
                               H2H, FR, MP (these are equivalent)
                               H2T, FF     (these are equivalent)
                               T2H, RR     (these are equivalent)
                             (default is tail-to-tail)
  --noise=<probability>      inject random sequencing errors (substitutions);
                             each base suffers a substitution error with the
                             given probability 
  --indel=<open>,<extend>    inject random sequencing errors (indels); <open>
                             is the probability of starting an indel at a
                             particular base; <extend> is the probability of
                             extended the gap after each base in the gap;
                             insertions and deletions are equally probable
  --errors=pacbio            simulate pacbio error profile
  --errors=nanopore          simulate nanopore error profile
  --errors=<spec>            simulate error profile with the given spec; <spec>
                             looks like this:
                                mm:1%,i:12%,d:2%
  --lower                    show substitutions in lowercase
  --prohibitn[=<length>]     don't sample reads containing a run of Ns; if
                             <length> is not given, no Ns are allowed
                             (by default, we do not check for Ns)
  --name=<template>          name for reads (see description below)
                             (default is FAKE_READ_[6], but with 6 replaced by
                             the smallest number sufficient)
  --width=<characters>       number of characters for each line of output dna;
                             this is only relevant for fasta output
                             (default is 100)
  --fastq                    output in fastq format
                             (default is fasta format)
  --intervals                output intervals, without sequence
  --output[+]=<filename>     write output to the specified file; if <filename>
                             contains {mate} or {zmate}, paired reads are
                             written to two files; if the plus sign is used,
                             we append to the file(s)
                             (default is stdout)
  --quality=<character>      set fastq qualities to a constant string of the
                             specified character (nothing better offered yet)
                             (default is J)
  --sam=<filename>           write a sam file equivalent of the generated reads
                             (currently not available for paired reads)
  --step=<number>            rather than sampling randomly, start at the
                             beginning and output a read starting at every Nth
                             position
  --strand=<+|->             only create reads from the specified strand
                             (by default we random choose strands)
  --cigars=<filename>        write cigar strings to the specified file
  --seed=<string>            set random seed
  --progress=<number>        periodically report how many reads we've generated

<template> is like this: BASE{4}_[6] where {4} is replaced by four random
letters/numbers and [6] is replaced by the read number (in six digits). [*]
is the read number in 'just enough' digits. Other recognized fields are
   {chrom}   the name of the sequence the read was drawn from
   {uchrom}  the name of the sequence the read was drawn from, in uppercase
   {start}   the starting position on that sequence (origin 1)
   {zstart}  the starting position on that sequence (origin 0)
   {end}     the ending position on that sequence
   {pstart}  the pair's starting position on that sequence (origin 1)
   {zpstart} the pair's starting position on that sequence (origin 0)
   {pend}    the pair's ending position on that sequence
   {strand}  the orientation on that sequence
   {mate}    for paired reads (1 or 2)
   {zmate}   for paired reads (0 or 1)"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global nucToSubstitutions,fixedStrand,qualityChar,cigarF
	global debug

	# parse the command line

	numReads       = None
	genomes        = []
	genomeToWeight = None
	readLength     = None
	insertLength   = None
	orientation    = None
	subProb        = None
	insOpenProb    = None
	insExtendProb  = None
	delOpenProb    = None
	delExtendProb  = None
	showAsLower    = False
	prohibitN      = None
	nameTemplate   = None
	lineWidth      = 100
	randomSeed     = None
	samplingStep   = None
	fixedStrand    = None
	cigarFilename  = None
	cigarSeparator = ""
	qualityChar    = None
	samFilename    = None
	outputAs       = None
	outFilename    = None
	appendToOutput = False
	reportProgress = None
	debug          = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--insert=")):
			if ("," in argVal):
				(ins,std) = argVal.split(",",1)
				ins = abs(float_with_unit(ins))
				std = abs(float_with_unit(std))
			else:
				ins = abs(float_with_unit(argVal))
				std = 0.0
			insertLength = (ins,std)
		elif (arg.startswith("--orientation=")):
			if   (argVal.upper() == "PE"):           orientation = "T2T"
			elif (argVal.upper() == "PAIREDEND"):    orientation = "T2T"
			elif (argVal.upper() == "T2T"):          orientation = "T2T"
			elif (argVal.upper() == "TAIL-TO-TAIL"): orientation = "T2T"
			elif (argVal.upper() == "RF"):           orientation = "T2T"
			elif (argVal.upper() == "MP"):           orientation = "H2H"
			elif (argVal.upper() == "MATEPAIR"):     orientation = "H2H"
			elif (argVal.upper() == "H2H"):          orientation = "H2H"
			elif (argVal.upper() == "HEAD-TO-HEAD"): orientation = "H2H"
			elif (argVal.upper() == "FR"):           orientation = "H2H"
			elif (argVal.upper() == "H2T"):          orientation = "H2T"
			elif (argVal.upper() == "FF"):           orientation = "H2T"
			elif (argVal.upper() == "T2H"):          orientation = "T2H"
			elif (argVal.upper() == "RR"):           orientation = "T2H"
			else: usage("unrecognized orientation: %s" % argVal)
		elif (arg.upper() in ["--PE","--PAIREDEND","--T2T","--RF"]):
			orientation = "T2T"
		elif (arg.upper() in ["--MP","--MATEPAIR","--H2H","--FR"]):
			orientation = "H2H"
		elif (arg.upper() in ["--H2T","--FF"]):
			orientation = "H2T"
		elif (arg.upper() in ["--T2H","--RR"]):
			orientation = "T2H"
		elif (arg.startswith("--noise=")):
			if (argVal == 0) or (argVal == "none"):
				subProb = None
			else:
				subProb = parse_probability(argVal)
		elif (arg.startswith("--indel=")):
			if (argVal == 0) or (argVal == "none"):
				insOpenProb = insExtendProb = None
				delOpenProb = delExtendProb = None
			else:
				(insOpenProb,insExtendProb) = argVal.split(",",1)
				insOpenProb   = parse_probability(insOpenProb) / 2
				insExtendProb = parse_probability(insExtendProb)
				assert (insExtendProb < 1)
				delOpenProb   = insOpenProb
				delExtendProb = insExtendProb
		elif (arg.startswith("--errors=")):
			errorProfile = None
			if (argVal in ["pacbio","pacbio.v3","pacbio.GIAB","pacbio.giab"]):
				errorProfile = errorProfilePacbioV3
			elif (argVal == "pacbio.v2"):  # for historical reasons, v2 is an alias for v3
				errorProfile = errorProfilePacbioV3
			elif (argVal in ["pacbio.v1","pacbio.Guiblet","pacbio.guiblet"]):
				errorProfile = errorProfilePacbioV1
			elif (argVal in ["pacbio.readsim"]):
				errorProfile = errorProfilePacbioReadsim
			elif (argVal in ["nanopore","nanopore.v3","nanopore.GIAB","nanopore.giab"]):
				errorProfile = errorProfileNanoporeV3
			elif (argVal == "nanopore.v2"):  # for historical reasons, v2 is an alias for v3
				errorProfile = errorProfileNanoporeV3
			elif (argVal in ["nanopore.v1","nanopore.Jain","nanopore.jain"]):
				errorProfile = errorProfileNanoporeV1
			elif (argVal in ["nanopore.readsim"]):
				errorProfile = errorProfileNanoporeReadSim
			elif (":" in argVal):
				try:
					errorProfile = parse_error_spec(argVal)
				except ValueError:
					pass
			if (errorProfile == None):
				usage("\"%s\" is not a valid error spec" % argVal)
			subProb       = errorProfile["mm"]
			insOpenProb   = errorProfile["i"]
			delOpenProb   = errorProfile["d"]
			insExtendProb = delExtendProb = 0.0
		elif (arg == "--lower"):
			showAsLower = True
		elif (arg in ["--prohibitn","--prohibitN"]):
			prohibitN = 1
		elif (arg.startswith("--prohibitn=")) or (arg.startswith("--prohibitN=")):
			if   (argVal.upper() == "NONE"): prohibitN = None
			elif (argVal == "0"):            prohibitN = None
			else:
				prohibitN = int(argVal)
				assert (prohibitN > 0)
		elif (arg.startswith("--name=")):
			nameTemplate = argVal
		elif (arg.startswith("--width=")):
			if (argVal.upper() == "NONE"): lineWidth = None
			else:                          lineWidth = abs(int(argVal))
		elif (arg.startswith("--out=")) or (arg.startswith("--output=")):
			outFilename    = argVal
			appendToOutput = False
		elif (arg.startswith("--out+=")) or (arg.startswith("--output+=")):
			outFilename    = argVal
			appendToOutput = True
		elif (arg.startswith("--quality=")):
			qualityChar = argVal
		elif (arg.startswith("--sam=")):
			samFilename = argVal
		elif (arg == "--fastq"):
			outputAs = "fastq"
		elif (arg == "--fasta"):
			outputAs = "fasta"
		elif (arg == "--intervals"):
			outputAs = "intervals"
		elif (arg.startswith("--step=")):
			samplingStep = int(argVal)
		elif (arg in ["--strand=+","--strand=plus","--strand=forward"]):
			fixedStrand = "+"
		elif (arg in ["--strand=-","--strand=minus","--strand=reverse","--strand=complement"]):
			fixedStrand = "-"
		elif (arg.startswith("--cigars=")) or (arg.startswith("--cigar=")):
			cigarFilename  = argVal
			cigarSeparator = ""
		elif (arg.startswith("--cigars+=")) or (arg.startswith("--cigar+=")):
			cigarFilename  = argVal
			cigarSeparator = " "
		elif (arg.startswith("--seed=")):
			if (argVal.lower() == "none"):
				randomSeed = None
			else:
				randomSeed = argVal
				randomSeed = randomSeed.replace("{date}",strftime("%d/%m/%Y"))
				randomSeed = randomSeed.replace("{time}",strftime("%I/%M/%S"))
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (parse_as_read_spec(arg) != None):
			(numReads,readLength) = parse_as_read_spec(arg)
		elif (parse_as_genome_spec(arg) != None):
			(weight,filename) = parse_as_genome_spec(arg)
			if (genomeToWeight == None): genomeToWeight = {}
			if (filename not in genomeToWeight):
				genomeToWeight[filename] =  weight
				genomes += [filename]
			else:
				genomeToWeight[filename] += weight
		else:
			usage("unrecognized option: %s" % arg)

	if (numReads == None) and (readLength != "stdin"):
		usage("you need to tell me how many reads to generate\n"
		    + "(e.g. 10Kx100 would be ten thousand reads, each 100 bp)")

	if (insertLength == None):
		if (type(readLength) == tuple):
			usage("you gave two read lengths but failed to use --insert")
		if (orientation != None):
			usage("you gave orientation=%s but failed to use --insert" % orientation)
	else: # if (insertLength != None):
		if (type(readLength) != tuple):
			readLength = (readLength,readLength)
		if (orientation == None):
			orientation = "T2T"

	if (randomSeed != None):
		#print >>stderr, "setting seed to \"%s\"" % randomSeed
		random_seed(randomSeed)

	readLengths = None
	if (readLength == "stdin"):
		readLengths = []
		for line in stdin:
			line = line.rstrip()
			for field in line.split():
				if (numReads == None) or (len(readLengths) < numReads):
					readLengths += [int(field)]
		numReads = len(readLengths)

	if (nameTemplate == None):
		numDigits = len("%d" % numReads)
		nameTemplate = "FAKE_READ_[%d]" % numDigits
		if (insertLength != None):
			nameTemplate += "#0/{mate}"
	elif ("[*]" in nameTemplate):
		numDigits = len("%d" % numReads)
		nameTemplate = nameTemplate.replace("[*]","[%d]" % numDigits)

	if (outputAs == None) and (outFilename != None):
		if   (outFilename.endswith(".fa")):    outputAs = "fasta"
		elif (outFilename.endswith(".fasta")): outputAs = "fasta"
		elif (outFilename.endswith(".fq")):    outputAs = "fastq"
		elif (outFilename.endswith(".fastq")): outputAs = "fastq"
	if (outputAs == None):
		outputAs = "fasta"

	if (appendToOutput): outputMode = "at"
	else:                outputMode = "wt"

	if (outputAs != "fastq") and (qualityChar != None):
		usage("quality character can only be used with fastq")
	elif (outputAs == "fastq") and (qualityChar == None):
		qualityChar = "J"

	nucsNeeded = (outputAs in ["fasta","fastq"])

	if (showAsLower):
		for nuc in nucToSubstitutions:
			subs = nucToSubstitutions[nuc]
			subs = [sub.lower() for sub in subs]
			nucToSubstitutions[nuc] = subs

	cigarF = None
	if (cigarFilename != None):
		cigarF = file(cigarFilename,"wt")

	# read the input sequence(s)

	genomeSelector = None

	if (samplingStep != None):
		# $$$ this should read from the file(s), too
		if (genomeToWeight != None):
			print >>stderr, "WARNING: "
		if ("input" in debug):
			print >>stderr, "reading stdin"
		assert (readLengths == None), \
		       "stdin can't provide both DNA and read lengths"
		sampler = DnaStepper(stdin,step=samplingStep,
		                     nucsNeeded=nucsNeeded,prohibitN=prohibitN,
			                 debug=debug)
	elif (genomeToWeight != None):
		genomeToSampler = {}
		genomeWeights   = []
		for filename in genomeToWeight:
			if ("input" in debug):
				print >>stderr, "reading %s" % filename
			if (filename == "/dev/stdin"):
				assert (readLengths == None), \
			       "stdin can't provide both DNA and read lengths"
			weight = genomeToWeight[filename]
			if (filename.endswith(".gz")) or (filename.endswith(".gzip")):
				f = gzip_open(filename,"rt")
			else:
				f = file(filename,"rt")
			if (randomSeed == None):
				sampler = DnaSampler(f,nucsNeeded=nucsNeeded,prohibitN=prohibitN,
				                     debug=debug)
			else:
				sampler = DnaSampler(f,nucsNeeded=nucsNeeded,prohibitN=prohibitN,
				                     seed=randomSeed+"_sampler"+"_"+filename,
				                     debug=debug)
			f.close()
			assert (len(sampler) != 0), \
		       "problem with \"%s\"; is it empty?" % filename
			genomeToSampler[filename] = sampler
			genomeWeights += [(weight*len(sampler),filename)]
		genomeSelector = ProbabilityTable(genomeWeights,scaleToOne=True)
		if ("weights" in debug):
			print >>stderr, genomeSelector.pMap
	elif (randomSeed == None):
		if ("input" in debug):
			print >>stderr, "reading stdin"
		assert (readLengths == None), \
		       "stdin can't provide both DNA and read lengths"
		sampler = DnaSampler(stdin,nucsNeeded=nucsNeeded,prohibitN=prohibitN,
		                     debug=debug)
	else:
		if ("input" in debug):
			print >>stderr, "reading stdin"
		assert (readLengths == None), \
		       "stdin can't provide both DNA and read lengths"
		sampler = DnaSampler(stdin,nucsNeeded=nucsNeeded,prohibitN=prohibitN,
		                     seed=randomSeed+"_sampler",debug=debug)

	if ("input" in debug):
		print >>stderr, "input complete"

	if (randomSeed == None):
		namer = RandomNamer(nameTemplate)
	else:
		namer = RandomNamer(nameTemplate,seed=randomSeed+"_namer")

	#=== generate single reads ===

	if (genomeSelector != None):
		genomeToCount = {}
		for genome in genomeToSampler:
			genomeToCount[genome] = 0

	if (insertLength == None):
		if (outFilename == None):
			outF = stdout
		else:
			assert ("{mate}" not in outFilename)
			assert ("{zmate}" not in outFilename)
			outF = file(outFilename,outputMode)

		samF = None
		if (samFilename != None):
			if (samFilename.endswith(".gz")) or (samFilename.endswith(".gzip")):
				samF = gzip_open(samFilename,"wt")
			else:
				samF = file(samFilename,"wt")
			write_sam_header(samF,sampler)

		for readNum in xrange(numReads):
			if (reportProgress != None):
				if (readNum+1 == 1) or ((readNum+1) % reportProgress == 0):
					print >>stderr, "(generating read number %s)" \
					              % commatize(readNum+1)

			if (readLengths != None):
				readLength = readLengths[readNum]

			name = namer.next()
			if (genomeSelector != None):
				genome = genomeSelector.choice()
				genomeToCount[genome] += 1
				sampler = genomeToSampler[genome]

			(basesNeeded,cigar) = random_cigar(readLength,insOpenProb,insExtendProb,delOpenProb,delExtendProb)
			read = sampler.sample(basesNeeded)
			if (read == None):
				print >>stderr, "exhausted input sequences after %d reads" \
				              % readNum
				break
			(seq,context) = read
			context["cigar"] = cigar

			if (subProb != None):
				seq = random_subs(seq,subProb)
			if (insOpenProb != None) or (delOpenProb != None):
				seq = random_indels(seq,cigar)

			modifiedName = modify_name(name,context)
			if (outputAs == "fasta"):
				quals = "*"
				write_fasta(outF,modifiedName,seq,width=lineWidth)
			elif (outputAs == "fastq"):
				quals = random_quals(len(seq))
				write_fastq(outF,modifiedName,seq,quals)
			elif (outputAs == "intervals"):
				write_interval(outF,modifiedName,context)
			else:
				assert (False)

			if (samF != None):
				write_sam(samF,modifiedName,seq,quals,context)

			if (cigarF != None):
				print >>cigarF, "%s %s%s %d %d %s" \
				              % (modifiedName,
				                 context["chrom"],context["strand"],
				                 context["start"],context["end"],
				                 cigar_to_string(cigar,cigarSeparator))
			elif ("cigars" in debug):
				print >>stderr, "%s %s" \
				               % (modifiedName,cigar_to_string(cigar," "))

		if (outF != stdout):
			outF.close()

		if (samF != None):
			samF.close()

	#=== generate paired reads ===

	else:
		if (outFilename == None):
			outF1 = outF2 = stdout
		elif ("{mate}" in outFilename):
			outF1 = file(outFilename.replace("{mate}","1"),outputMode)
			outF2 = file(outFilename.replace("{mate}","2"),outputMode)
		elif ("{zmate}" in outFilename):
			outF1 = file(outFilename.replace("{zmate}","0"),outputMode)
			outF2 = file(outFilename.replace("{zmate}","1"),outputMode)
		else:
			outF1 = outF2 = file(outFilename,outputMode)

		if (samFilename != None):
			assert (False), "sam support is not implemented for paired reads"

		for readNum in xrange(numReads):
			if (reportProgress != None) and ((readNum+1) % reportProgress == 0):
				print >>stderr, "(generating read pair number %s)" \
				              % commatize(readNum+1)

			name = namer.next()
			if (genomeSelector != None):
				genome = genomeSelector.choice()
				genomeToCount[genome] += 1
				sampler = genomeToSampler[genome]

			(basesNeeded1,cigar1) = random_cigar(readLength[0],insOpenProb,insExtendProb,delOpenProb,delExtendProb)
			(basesNeeded2,cigar2) = random_cigar(readLength[1],insOpenProb,insExtendProb,delOpenProb,delExtendProb)

			pair = sampler.sample_pair(insertLength,[basesNeeded1,basesNeeded2],orientation)
			if (pair == None):
				print >>stderr, "exhausted input sequences after %d pairs" \
				              % readNum
				break
			(seq1,context1,seq2,context2) = pair
			context1["cigar"] = cigar1
			context2["cigar"] = cigar2

			if (subProb != None):
				seq1 = random_subs(seq1,subProb)
				seq2 = random_subs(seq2,subProb)
			if (insOpenProb != None) or (delOpenProb != None):
				seq1 = random_indels(seq1,cigar1)
				seq2 = random_indels(seq2,cigar2)

			if (outputAs == "fasta"):
				write_fasta(outF1,modify_name(name,context1),seq1,width=lineWidth)
				write_fasta(outF2,modify_name(name,context2),seq2,width=lineWidth)
			elif (outputAs == "fastq"):
				write_fastq(outF1,modify_name(name,context1),seq1,random_quals(len(seq1)))
				write_fastq(outF2,modify_name(name,context2),seq2,random_quals(len(seq2)))
			elif (outputAs == "intervals"):
				write_interval(outF1,modify_name(name,context1),context1)
				write_interval(outF2,modify_name(name,context2),context2)
			else:
				assert (False)

			if (cigarF != None):
				print >>cigarF, "%s %s%s %d %d %s" \
				              % (modify_name(name,context1),
				                 context1["chrom"],context1["strand"],
				                 context1["start"],context1["end"],
				                 cigar_to_string(cigar1))
				print >>cigarF, "%s %s%s %d %d %s" \
				              % (modify_name(name,context2),
				                 context1["chrom"],context2["strand"],
				                 context1["start"],context2["end"],
				                 cigar_to_string(cigar2))
			elif ("cigars" in debug):
				print >>stderr, "%s %s" \
				              % (modify_name(name,context1),cigar_to_string(cigar1," "))
				print >>stderr, "%s %s" \
				              % (modify_name(name,context2),cigar_to_string(cigar2," "))

		if (outF1 != stdout):
			outF1.close()
		if (outF2 != stdout) and (outF2 != outF1):
			outF2.close()

	if ("choices" in debug) or ("weights" in debug):
		if (genomeSelector == None):
			print >>stderr, "can't tell you how many were chosen from which!"
		else:
			if (insertLength == None): readKind = "reads"
			else:                      readKind = "pairs"
			for genome in genomes:
				print >>stderr, "%d %s chosen from %s" \
				              % (genomeToCount[genome],readKind,genome)


# random_cigar--
#	Choose a random sequence of matches and indels;  we alternate between the
#	two types, and thus will not have back-to-back match segments or indels;
#	segment lengths are modeled by the geometric distribution.
#
#	We return the cigar operator list, as well as the number of reference bases
#	needed for a read with this cigar.
#
#	Note: a cigar starting or ending with a deletion would make no sense.
#	Downstream programs (like aligners) can't possibly distinguish between the
#	presence or absence of such an event.  We may conceptually generate a
#	deletion at the start, as our first event, but we don't store it.  As for
#	the end, the basesToGo>0 loop can only end on an event that decreases
#	basesToGo, which deletions don't do.
#
# $$$ the use of lambda functions for length distributions is in prep for
#     .. allowing other user-specified length distributions; ideally
#     .. insertionLenFunc, etc. would be passed as arguments instead of the
#     .. current implementation with globals

prevInsOpenProb   = None
prevInsExtendProb = None
prevDelOpenProb   = None
prevDelExtendProb = None
matchLenFunc      = None
insertionLenFunc  = None
deletionLenFunc   = None

def random_cigar(readLen,insOpenProb,insExtendProb,delOpenProb,delExtendProb):
	global prevInsOpenProb,prevInsExtendProb
	global prevDelOpenProb,prevDelExtendProb
	global matchLenFunc,insertionLenFunc,deletionLenFunc

	if (insOpenProb == None): insOpenProb = 0
	if (delOpenProb == None): delOpenProb = 0
	indelOpenProb = insOpenProb + delOpenProb

	if   (readLen       <= 0): return (readLen,[])
	elif (indelOpenProb == 0): return (readLen,[("M",readLen)])

	indelIsInsertProb = insOpenProb / indelOpenProb

	if (insOpenProb != prevInsOpenProb):
		prevInsOpenProb = insOpenProb
		matchLenFunc = geometric_distribution_func(1-indelOpenProb)
	if (insExtendProb != prevInsExtendProb):
		prevInsExtendProb = insExtendProb
		insertionLenFunc  = geometric_distribution_func(insExtendProb)
	if (delExtendProb != prevDelExtendProb):
		prevDelExtendProb = delExtendProb
		deletionLenFunc   = geometric_distribution_func(delExtendProb)

	basesToGo = seqNeeded = readLen
	cigar = []

	while (basesToGo > 0):
		# first pass through the loop might not have an indel; subsequent
		# passes will always have an indel; note that if we get a deletion
		# in this first indel, we don't bother to save it (see note above)

		if (basesToGo < readLen) or (unit_random() < indelOpenProb):
			if (unit_random() < indelIsInsertProb):
				runLen = min(basesToGo,insertionLenFunc())
				cigar += [("I",runLen)]
				basesToGo -= runLen
				seqNeeded -= runLen
			elif (len(cigar) > 0):
				runLen = deletionLenFunc()
				cigar += [("D",runLen)]
				seqNeeded += runLen

		if (basesToGo > 0):
			runLen = min(basesToGo,matchLenFunc())
			cigar += [("M",runLen)]
			basesToGo -= runLen

	return (seqNeeded,cigar)


# random_subs--

def random_subs(seq,prob):
	errors = [ix for ix in xrange(len(seq)) if (unit_random() < prob)]
	if (errors == []): return seq

	seq  = list(seq)
	subs = []
	for ix in errors:
		nuc = seq[ix]
		if (nuc in nucToSubstitutions):
			seq[ix] = random_choice(nucToSubstitutions[nuc])
			subs += [(ix,nuc,seq[ix])]

	return "".join(seq)


# random_indels--

def random_indels(seq,cigar):
	newSeq = []
	pos = 0
	for (op,opLen) in cigar:
		if (op == "I"):
			newSeq += [random_choice("ACGT") for _ in xrange(opLen)]
		elif (op == "D"):
			pos += opLen
		else: # if (op == "M"):
			newSeq += [seq[pos:pos+opLen]]
			pos += opLen

	assert (pos == len(seq)), "internal error"

	return "".join(newSeq)


# random_quals--

def random_quals(seqLen):
	# $$$ improve this!  A markov chain would be a good choice
	numCopies = (seqLen+len(qualityChar)-1) / len(qualityChar)
	return (qualityChar * numCopies)[:seqLen]


# modify_name--

def modify_name(name,context):
	if ("chrom" in context):
		name = name.replace("{chrom}",  context["chrom"])
		name = name.replace("{uchrom}", context["chrom"].upper())
		name = name.replace("{CHROM}",  context["chrom"].upper())

	if ("start" in context):
		name = name.replace("{start}",  str(context["start"]+1))
		name = name.replace("{zstart}", str(context["start"]))

	if ("end" in context):
		name = name.replace("{end}",    str(context["end"]))

	if ("pstart" in context):
		name = name.replace("{pstart}", str(context["pstart"]+1))
		name = name.replace("{zpstart}",str(context["pstart"]))
		name = name.replace("{pzstart}",str(context["pstart"]))

	if ("pend" in context):
		name = name.replace("{pend}",   str(context["pend"]))

	if ("chrom" in context):
		name = name.replace("{strand}","F" if (context["strand"] == "+") else "R")

	if ("mate" in context):
		name = name.replace("{mate}",  str(context["mate"]))
		name = name.replace("{zmate}", str(context["mate"]-1))

	return name


# write_fasta--

def write_fasta(f,name,seq,width=None):
	print >>f, ">%s" % name
	if (width == None):
		print >>f, seq
	else:
		for ix in range(0,len(seq),width):
			print >>f, "%s" % seq[ix:ix+width]


# write_fastq--

def write_fastq(f,name,seq,quals):
	print >>f, "@%s" % name
	print >>f, seq
	print >>f, "+"
	print >>f, quals


# write_interval--

def write_interval(f,name,context):
	chrom  = context["chrom"]
	start  = context["start"]
	end    = context["end"]
	strand = context["strand"]
	print >>f, "%s\t%s%s\t%d\t%d" % (name,chrom,strand,start,end)


# write_sam_header--

def write_sam_header(f,sampler):
	print >>f, "@HD\tVN:1.3\tSO:unsorted"
	for chrom in sampler.chromOrder:
		print >>f, "@SQ\tSN:%s\tLN:%d" % (chrom,sampler.chromToLen[chrom])
	print >>f, "@PG	ID:ncrf_read_simulator\tPN:ncrf_read_simulator"


# write_sam--

def write_sam(f,name,seq,quals,context):
	chrom  = context["chrom"]
	start  = context["start"] + 1
	strand = context["strand"]

	flags = BAM_FREVERSE if (strand == "-") else 0
	mapQ  = 60
	cigar = cigar_to_string(context["cigar"])

	print >>f, "%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s" \
	         % (name,flags,chrom,start,mapQ,cigar,seq,quals)


# cigar_to_string--

def cigar_to_string(cigar,separator=""):
	return separator.join(["%d%s"%(opLen,op) for (op,opLen) in cigar])


# DnaSampler--
#	Class to sample reads or pairs from a collection of DNA sequences.

class DnaSampler(object):

	def __init__(self,f,nucsNeeded=True,prohibitN=None,seed=None,debug=None):
		if (debug == None): self.debug = []
		else:               self.debug = list(debug)

		self.nucsNeeded = nucsNeeded
		self.prohibitN  = prohibitN

		self.prng = Random()
		self.prng.seed(seed)

		self.chromOrder = []
		self.chromToLen = {}

		if (self.nucsNeeded):
			self.segOrder = []
			self.segToSeq = {}
			self.segToLen = {}
			for (chrom,segStart,seq) in self.read_fasta_sequences(f):
				assert ((chrom,segStart) not in self.segToSeq), \
				       "\"%s\" appears more than once" % chrom
				self.segOrder += [(chrom,segStart,len(seq))]
				self.segToSeq[(chrom,segStart)] = seq
				self.segToLen[(chrom,segStart)] = len(seq)
		else:
			self.segOrder = []
			self.segToLen = {}
			for (chrom,length) in self.read_fasta_lengths(f):
				assert ((chrom,0) not in self.segToLen), \
				       "\"%s\" appears more than once" % chrom
				self.segOrder += [(chrom,0,length)]
				self.segToLen[(chrom,0)] = length

		self.numSequences = 0
		self.totalLen  = 0
		self.minLen    = None
		self.maxLen    = None
		for (chrom,segStart) in self.segToLen:
			seqLen = self.segToLen[(chrom,segStart)]
			self.numSequences += 1
			self.totalLen     += seqLen
			if (self.minLen == None) or (seqLen < self.minLen):
				self.minLen = seqLen
			if (self.maxLen == None) or (seqLen > self.maxLen):
				self.maxLen = seqLen

		if ("segments" in self.debug):
			keys = [(chrom,segStart) for (chrom,segStart) in self.segToLen]
			keys.sort()
			for (chrom,segStart) in keys:
				seqLen = self.segToLen[(chrom,segStart)]
				print >>stderr, "%s %d..%d" % (chrom,segStart,segStart+seqLen)


	def __len__(self):
		return self.totalLen

	def read_fasta_sequences(self,f):
		chrom = "(no_name)"
		lines = []

		lineNum = 0
		for line in f:
			lineNum += 1
			line = line.rstrip()

			if (line.startswith(">")):
				if (lines != []):
					seq = "".join(lines)
					self.chromOrder += [chrom]
					self.chromToLen[chrom] = len(seq)
					for (segStart,seq) in self.nucleotide_segments(seq):
						yield (chrom,segStart,seq)
				chrom  = line[1:].strip()
				lines = []
			else:
				lines += [line.upper()]

		if (lines != []):
			seq = "".join(lines)
			self.chromOrder += [chrom]
			self.chromToLen[chrom] = len(seq)
			for (segStart,seq) in self.nucleotide_segments(seq):
				yield (chrom,segStart,seq)


	def nucleotide_segments(self,seq):
		if (self.prohibitN == None):
			yield (0,seq)
			return

		nStart = None
		nRuns = []
		seqLen = len(seq)
		for (ix,nuc) in enumerate(seq):
			if (nuc == "N"):
				if (nStart == None): nStart = ix
			elif (nStart != None):
				if (ix - nStart >= self.prohibitN): nRuns += [(nStart,ix)]
				nStart = None
		if (nStart != None):
			if (seqLen - nStart >= self.prohibitN): nRuns += [(nStart,seqLen)]
		else:
			nRuns += [(seqLen,seqLen)]  # sentinel

		(nStart,nEnd) = nRuns[0]
		if (nStart != 0):
			segEnd = nStart
			yield (0,seq[0:segEnd])
		segStart = nEnd

		for (nStart,nEnd) in nRuns[1:]:
			segEnd = nStart
			yield (segStart,seq[segStart:segEnd])
			segStart = nEnd


	def read_fasta_lengths(self,f):
		chrom  = "(no_name)"
		length = -1

		lineNum = 0
		for line in f:
			lineNum += 1
			line = line.rstrip()

			if (line.startswith(">")):
				if (length >= 0): yield (chrom,length)
				chrom  = line[1:].strip()
				length = 0
			else:
				length += len(line.upper())

		if (length >= 0): yield (chrom,length)


	def sample(self,seqLen):
		piece = self.sample_position(seqLen)
		if (piece == None): return None
		(chrom,segStart,start) = piece

		if   (fixedStrand != None):     strand = fixedStrand
		elif (self.prng.random() < .5): strand = "+"
		else:                           strand = "-"

		context = {}
		context["chrom" ] = chrom
		context["start" ] = segStart + start
		context["end"   ] = segStart + start + seqLen
		context["strand"] = strand

		if (self.nucsNeeded):
			seq = self.segToSeq[(chrom,segStart)][start:start+seqLen]
			if (strand == "-"): seq = reverse_complement(seq)
		else:
			seq = None

		return (seq,context)


	def sample_pair(self,seqLenDistrib,readLengths,orientation="T2T"):
		(mu,sigma) = seqLenDistrib
		(lftLen,rgtLen) = readLengths
		minLen = min(lftLen,rgtLen)

		maxTries = 3
		while (maxTries > 0):
			insertLen = abs(int(floor(self.prng.gauss(mu,sigma))))
			if (minLen <= insertLen <= self.maxLen): break
			maxTries -= 1
		assert (insertLen <= self.maxLen)

		piece = self.sample_position(insertLen)
		if (piece == None): return None
		(chrom,segStart,start) = piece

		if (fixedStrand != None): forwardStrand = (fixedStrand == "+")
		else:                     forwardStrand = (self.prng.random() < .5)

		if (orientation == "T2T"):
			if (forwardStrand):
				(strand1,start1) = ("-",start)
				(strand2,start2) = ("+",start+insertLen-rgtLen)
			else:
				(strand1,start1) = ("+",start+insertLen-lftLen)
				(strand2,start2) = ("-",start)
		elif (orientation == "H2H"):
			if (forwardStrand):
				(strand1,start1) = ("+",start)
				(strand2,start2) = ("-",start+insertLen-rgtLen)
			else:
				(strand1,start1) = ("-",start+insertLen-lftLen)
				(strand2,start2) = ("+",start)
		elif (orientation == "H2T"):
			if (forwardStrand):
				(strand1,start1) = ("+",start)
				(strand2,start2) = ("+",start+insertLen-rgtLen)
			else:
				(strand1,start1) = ("-",start+insertLen-lftLen)
				(strand2,start2) = ("-",start)
		elif (orientation == "T2H"):
			if (forwardStrand):
				(strand1,start1) = ("-",start)
				(strand2,start2) = ("-",start+insertLen-rgtLen)
			else:
				(strand1,start1) = ("+",start+insertLen-lftLen)
				(strand2,start2) = ("+",start)
		else:
			raise ValueError

		context1 = {}
		context1["chrom" ] = chrom
		context1["start" ] = segStart + start1
		context1["end"   ] = segStart + start1 + lftLen
		context1["strand"] = strand1
		context1["mate"  ] = 1

		context2 = {}
		context2["chrom" ] = chrom
		context2["start" ] = segStart + start2
		context2["end"   ] = segStart + start2 + rgtLen
		context2["strand"] = strand2
		context2["mate"  ] = 2

		context1["pstart"] = context2["pstart"] = context1["start" ]
		context1["pend"  ] = context2["pend"  ] = context2["end"   ]

		if (self.nucsNeeded):
			seq1 = self.segToSeq[(chrom,segStart)][start1:start1+lftLen]
			seq2 = self.segToSeq[(chrom,segStart)][start2:start2+rgtLen]
			if (strand1 == "-"): seq1 = reverse_complement(seq1)
			if (strand2 == "-"): seq2 = reverse_complement(seq2)
		else:
			seq1 = seq2 = None

		return (seq1,context1,seq2,context2)


	def sample_position(self,seqLen):
		if (seqLen > self.maxLen): raise ValueError

		if (seqLen <= self.minLen):
			numPlacements = self.totalLen - self.numSequences*(seqLen-1)
		else:
			numPlacements = self.totalLen  # (may be an overestimate)

		r = self.prng.randint(0,numPlacements-1)
		realNumPlacements = 0
		for (chrom,segStart) in self.segToLen:
			seqPlacements = self.segToLen[(chrom,segStart)] - (seqLen-1)
			if (r < seqPlacements): return (chrom,segStart,r)
			r -= seqPlacements
			assert (r >= 0)
			realNumPlacements += seqPlacements

		r = self.prng.randint(0,realNumPlacements-1)
		for (chrom,segStart) in self.segToLen:
			seqPlacements = self.segToLen[(chrom,segStart)] - (seqLen-1)
			if (r < seqPlacements): return (chrom,segStart,r)
			r -= seqPlacements
			assert (r >= 0)

		assert (False)


# DnaStepper--
#	Class to sample at a uniform sampling step.

class DnaStepper(DnaSampler):

	def __init__(self,f,step=1,nucsNeeded=True,prohibitN=None,debug=None):
		super(DnaStepper,self).__init__(f,nucsNeeded=nucsNeeded,prohibitN=prohibitN,debug=debug)
		self.step  = step
		self.where = None

	def sample_position(self,seqLen):
		if (seqLen > self.maxLen): raise ValueError

		if (self.where == None): self.where = (0,0)
		(chromIx,positionOnSeg) = self.where

		while (True):
			if (chromIx >= len(self.segOrder)):
				return None      # exhausted input sequences
			(chrom,segStart,segLen) = self.segOrder[chromIx]
			if (positionOnSeg+seqLen <= segLen): break
			(chromIx,positionOnSeg) = (chromIx+1,0)

		self.where = (chromIx,positionOnSeg+self.step)
		return (chrom,segStart,positionOnSeg)


# RandomNamer--
#	Class to generate 'random' names according to a template.

class RandomNamer(object):

	def __init__(self,template,seed=None,duplicatesOK=False,attemptsLimit=10000):
		self.template      = template
		self.attemptsLimit = attemptsLimit

		self.prng = Random()
		self.prng.seed(seed)

		# determine which tokens are present in the template

		self.randomized = []
		for w in xrange(20,0,-1):
			brack = "{%d}" % w				# e.g. "{5}" means 5 random
			fields = template.split(brack)	# .. characters, "different" for
			if (len(fields) > 1):			# .. each read
				self.randomized += [(w,brack,"R",None)]

			brack = "{%dL}" % w				# e.g. "{5L}" means 5 random
			fields = template.split(brack)	# .. letters, "different" for each
			if (len(fields) > 1):			# .. read
				self.randomized += [(w,brack,"L",None)]

			brack = "{%dD}" % w				# e.g. "{5D}" means 5 random
			fields = template.split(brack)	# .. digits, "different" for each
			if (len(fields) > 1):			# .. read
				self.randomized += [(w,brack,"D",None)]

		self.counted     = None				# e.g. "[5]" means 5 digit count
		self.countOffset = 0				# or "[5@10001] means start at 10001
		for w in xrange(13,-1,-1):			# ... but only longest is accepted
			if (w > 0):
				brack   = "[%d]" % w
				brackRe = compile("\[%d@(?P<start>[1-9][0-9]*)\]" % w)
				fmt     = "%%0%dd" % w
			else:
				brack   = "[]"
				brackRe = compile("\[@(?P<start>[1-9][0-9]*)\]")
				fmt     = "%d"
			fields = template.split(brack,1)
			if (len(fields) > 1):
				self.counted = (w,brack,fmt)
				break
			m = brackRe.search(template)
			if (m != None):
				start = m.group("start")
				self.countOffset = int(start)-1
				match = m.group(0)
				fields = template.split(match,1)
				self.template = fields[0] + brack + fields[1]
				self.counted = (w,brack,fmt)
				break

		# if necessary, create a hash to prevent duplicate names

		if (duplicatesOK) or (self.randomized == []):
			self.nameUsed = None
		else:
			self.nameUsed = {}

		self.count = 0


	def next(self):
		self.count += 1
		name = None
		for attempt in xrange(self.attemptsLimit):
			cand = self.generate_name()
			if (self.nameUsed == None) or (cand not in self.nameUsed):
				name = cand
				break
		assert (name != None), "failed to generate name %d for %s" \
		                     % (self.count,self.template)
		if (self.nameUsed != None):
			self.nameUsed[name] = True

		return name


	def generate_name(self):
		name = self.template

		# substitute, e.g. "{5}" with 5 random characters

		for (w,brack,kind,val) in self.randomized:
			fields = name.split(brack)
			name = fields[0]
			if (kind == ""):
				for f in fields[1:]:
					name += val + f
			elif (kind == "L"):
				for f in fields[1:]:
					name += self.rand_letters(w) + f
			elif (kind == "D"):
				for f in fields[1:]:
					name += self.rand_digits(w) + f
			else: # if (kind == "R"):
				for f in fields[1:]:
					name += self.rand_chars(w) + f

		# substitute, e.g. "[5]" with "%05d" % count

		if (self.counted != None):
			(w,brack,fmt) = self.counted
			fields = name.split(brack,1)
			name = (fmt % (self.count+self.countOffset)).join(fields)

		return name


	def rand_chars(self,n):
		return "".join([self.prng.choice("ACDEFGHJKLMNPQRTUVWXYZ0123456789") 
		                                                  for i in range(n)])

	def rand_letters(self,n):
		return "".join([self.prng.choice("ACDEFGHJKLMNPQRTUVWXYZ") for i in range(n)])

	def rand_digits(self,n):
		return "".join([self.prng.choice("0123456789") for i in range(n)])


# geometric_distribution--

def geometric_distribution_func(pExtend):
	return lambda: geometric_distribution(pExtend)

def geometric_distribution(pExtend):
	if (pExtend == 0): return 1
	u = unit_random()
	return int(floor(1+log(1-u)/log(pExtend)))


# reverse_complement--

complementMap = maketrans("ACGTSWRYMKBDHVNacgtswrymkbdhvn",
                          "TGCASWYRKMVHDBNtgcaswyrkmvhdbn")

def reverse_complement(nukes):
	return nukes[::-1].translate(complementMap)


# parse_as_read_spec--
#	Parse a string for <num>x<length>[,<length>]

def parse_as_read_spec(s):
	if (s.startswith("--")): return None
	if (s == "stdin"): return (None,"stdin")
	if ("x" not in s): return None
	(numReads,readLength) = s.split("x",1)

	try:
		numReads = int_with_unit(numReads)
		assert (numReads >= 0)
		if (readLength == "stdin"):
			pass
		elif ("," not in readLength):
			readLength = int_with_unit(readLength)
			assert (readLength > 0)
		else:
			(lft,rgt) = readLength.split(",",1)
			lft = int_with_unit(lft)
			rgt = int_with_unit(rgt)
			assert (lft > 0) and (rgt > 0)
			readLength = (lft,rgt)
	except ValueError:
		return None

	return (numReads,readLength)


# parse_as_genome_spec--
#	Parse a string for [<weight>:]<fasta_file>

def parse_as_genome_spec(s):
	if (s.startswith("--")): return None
	if (":" not in s):
		weight = 1.0
	else:
		(weight,s) = s.split(":",1)
		weight = weight.strip()
		s      = s.strip()
		try:               weight = parse_probability(weight,strict=False)
		except ValueError: return None

	return (weight,s)


# parse_error_spec--
#	A typical spec looks like this: mm:1%,i:12%,d:2%

def parse_error_spec(spec):
	fields = spec.split(",")
	if (len(fields) != 3): raise ValueError

	errorProfile = {}
	for field in fields:
		field = field.strip()
		if (":" not in field): raise ValueError
		(event,p) = field.split(":",1)
		if (event not in ["mm","i","d"]): raise ValueError
		if (event in errorProfile): raise ValueError  # (event specified twice)
		errorProfile[event] = parse_probability(p)

	return errorProfile


# float_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def float_with_unit(s):
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

	return float(s) * multiplier


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


if __name__ == "__main__": main()
