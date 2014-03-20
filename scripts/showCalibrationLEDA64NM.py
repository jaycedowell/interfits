#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to displaying the contents of a (LEDA-specific) FITS IDI CALIBRATION table.
"""

import os
import re
import sys
import math
import pyfits


# Regular expression for trying to get the stand number out of an antenna
# name
_annameRE = re.compile('^.*?(?P<id>\d{1,3})$')


def main(args):
	filename = args[0]
	hdulist = pyfits.open(filename)
	
	try:
		calData = hdulist['CALIBRATION'].data
		antData = hdulist['ARRAY_GEOMETRY'].data
	except Exception, e:
		print str(e)
		raise RuntimeError("CALIBRATION table not found in %s" % os.path.basename(filename))
		
	print "Ant.  Stand  Delay X    Phase X  Delay Y    Phase Y"
	print "             [ns]       [deg]    [ns]       [deg]  "
	print "---------------------------------------------------"
	for i in xrange(len(calData['ANTENNA_NO'])):
		aid = calData['ANTENNA_NO'][i]
		stand = _annameRE.match(antData['ANNAME'][i])
		stand = int(stand.group('id'))
		delayX = calData['DELAY_A'][i]*1e9
		phaseX = calData['PHASE_A'][i]*180/math.pi
		delayY = calData['DELAY_B'][i]*1e9
		phaseY = calData['PHASE_B'][i]*180/math.pi
		
		print "%3i   %3i    %9.4f  %5.4f   %9.4f  %5.4f" % (aid, stand, delayX, phaseX, delayY, phaseY)
		


if __name__ == "__main__":
	main(sys.argv[1:])