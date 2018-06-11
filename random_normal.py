#!/usr/bin/env python
"""
Generate random numbers, sampled from a normal distribution.
"""

from sys        import argv,stdin,stderr,exit
from random     import seed as random_seed,gauss
from ncrf_parse import int_with_unit,float_or_fraction


def usage(s=None):
	message = """
usage: random_normal <num_values> [options]
  <num_values>             number of values to generate
  --mu=<value>             mean value         (default is 0.0)
  --sigma=<value>          standard deviation (default is 1.0)
  --round                  round to integers
  --floor                  "round down" to integers
  --ceiling                "round up" to integers
  --precision=<digits>     set number of digits after decimal
  --seed=<string>          set random seed"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	# parse the command line

	numValues = None
	mu        = 0.0
	sigma     = 1.0
	roundEm   = None
	precision = None

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--mu=")):
			mu = float_or_fraction(argVal)
		elif (arg.startswith("--sigma=")):
			sigma = float_or_fraction(argVal)
		elif (arg == "--round"):
			roundEm = "round"
		elif (arg == "--floor"):
			roundEm = "floor"
		elif (arg == "--ceiling"):
			roundEm = "ceiling"
		elif (arg.startswith("--precision=")):
			precision = int(argVal)
		elif (arg.startswith("--seed=")):
			# nota bene: if the seed is a number, use it as a number, since
			#            string seeds can produce different sequences on
			#            different versions/builds of python
			seed = argVal
			try:
				seed = int(seed)
			except ValueError:
				try:               seed = float(seed)
				except ValueError: pass
			random_seed(seed)
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (numValues == None):
			numValues = int_with_unit(arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (numValues == None):
		numValues = 1

	if (roundEm != None) and (precision != None):
		usage("can't use --precision with --%s" % roundEm)

	if   (precision == None): vFmt = "%s"
	elif (precision <= 0):    vFmt = "%d."
	else:                     vFmt = "%%.%df" % precision

	# generate the values

	for _ in xrange(numValues):
		v = gauss(mu,sigma)
		if   (roundEm == "round"):   v = int(round(v))
		elif (roundEm == "floor"):   v = int(floor(v))
		elif (roundEm == "ceiling"): v = int(ceil (v))

		print vFmt % v


if __name__ == "__main__": main()
