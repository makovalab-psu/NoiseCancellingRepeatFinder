#!/usr/bin/env python
"""
Read event counts from a truth catalog and report overall induced event
rates.
"""

from sys  import argv,stdin,stdout,stderr,exit

class Observations: pass


def usage(s=None):
	message = """
usage: cat <truth_catalog> | truth_to_rates [options]
  (currently, there are no options)

The truth catalog is usually the output from the --catalog option of
mock_motif_read. It has 12 columns but only the events are used here ("m",
"mm", "i", and "d", columns 9, 10, 11, and 12)."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	for arg in argv[1:]:
		if (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# accumulate input data

	obs = Observations()
	obs.m = obs.mm = obs.i = obs.d = 0

	lineNumber = 0
	for line in stdin:
		lineNumber +=1

		line = line.strip()
		if (line.startswith("#")):    continue
		if (line.startswith("line")): continue

		fields = line.split()
		assert (len(fields) != 7), \
		      ("wrong number of columns in line %d (%d, expected %d)" \
		       + "\nwas mock_motif_read used without the --errors= option?" \
		       + "\n%s") \
		    % (lineNumber,len(fields),12,line)
		assert (len(fields) == 12), \
		      "wrong number of columns in line %d (%d, expected %d)\n%s" \
		    % (lineNumber,len(fields),12,line)

		try:
			obs.m  += int(fields[8])
			obs.mm += int(fields[9])
			obs.i  += int(fields[10])
			obs.d  += int(fields[11])
		except ValueError:
			assert (False), \
				  "failed to parse line %d\n%s" \
				% (lineNumber,line)

	obs.events = obs.m + obs.mm + obs.i + obs.d

	if (obs.events == 0):
		exit("ERROR: input is empty")

	print "mm=%.3f%% i=%.3f%% d=%.3f%%" \
	    % (100.0*obs.mm/obs.events,
	       100.0*obs.i/obs.events,
	       100.0*obs.d/obs.events)


if __name__ == "__main__": main()
