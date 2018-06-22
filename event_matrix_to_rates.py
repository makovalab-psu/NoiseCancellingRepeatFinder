#!/usr/bin/env python
"""
Read match-and-error event counts and report overall event rates.
"""

from sys  import argv,stdin,stdout,stderr,exit

class Observations: pass


def usage(s=None):
	message = """
usage: cat <event_counts_table> | event_matrix_to_rates [options]
  (currently, there are no options)

The input is usually the output from ncrf_extract_event_matrix.  It has 9
columns, first three of which are ignored here. The remaining columns are,
respectively, counts for matches ("m"), mismatches ("mm"), insertion opens
("io"), insertion extensions ("ix"), deletion opens ("do"), and deletion
extensions ("dx")."""

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
	obs.m = obs.mm = obs.io = obs.ix = obs.do = obs.dx = 0

	lineNumber = 0
	for line in stdin:
		lineNumber +=1

		line = line.strip()
		if (line.startswith("#")):    continue
		if (line.startswith("line")): continue

		fields = line.split()
		assert (len(fields) == 9), \
		      "wrong number of columns in line %d (%d, expected %d)\n%s" \
		    % (lineNumber,len(fields),9,line)

		try:
			obs.m  += int(fields[3])
			obs.mm += int(fields[4])
			obs.io += int(fields[5])
			obs.ix += int(fields[6])
			obs.do += int(fields[7])
			obs.dx += int(fields[8])
		except ValueError:
			assert (False), \
				  "failed to parse line %d\n%s" \
				% (lineNumber,line)

	obs.events = obs.m + obs.mm + obs.io + obs.ix + obs.do + obs.dx

	if (obs.events == 0):
		exit("ERROR: input is empty")

	print "mm=%.3f%% io=%.3f%% ix=%.3f%% do=%.3f%% dx=%.3f%%" \
	    % (100.0*obs.mm/obs.events,
	       100.0*obs.io/obs.events,
	       100.0*obs.ix/obs.events,
	       100.0*obs.do/obs.events,
	       100.0*obs.dx/obs.events)


if __name__ == "__main__": main()
