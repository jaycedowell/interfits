#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to plot the autocorrelations in a LEDA64-NM FITS-IDI file.
"""

import os
import sys
import numpy
import getopt
import pyfits
from datetime import datetime

from lsl import astro
from lsl.imaging import utils
from lsl.common.progress import ProgressBar

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """autoLEDA64NM.py - Plot the XX and YY autocorrelations in a LEDA64-NM FITS-IDI file

Usage: autoLEDA64NM.py [OPTIONS] file

Options:
-h, --help             Display this help information

"""
	
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseConfig(args):
	config = {}
	# Command line flags - default values
	config['args'] = []
	
	# Read in and process the command line flags
	try:
		opts, arg = getopt.getopt(args, "h", ["help",])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		else:
			assert False
			
	# Add in arguments
	config['args'] = arg
	
	# Return configuration
	return config


def main(args):
	# Parse the command line
	config = parseConfig(args)
	filename = config['args'][0]
	
	idi = utils.CorrelatedData(filename)
	aa = idi.getAntennaArray()
	lo = idi.getObserver()
	lo.date = idi.dateObs.strftime("%Y/%m/%d %H:%M:%S")
	jd = lo.date + astro.DJD_OFFSET
	lst = str(lo.sidereal_time())
	
	nStand = len(idi.stands)
	nChan = len(idi.freq)
	freq = idi.freq
	
	print "Raw Stand Count: %i" % nStand
	print "Final Baseline Count: %i" % (nStand*(nStand-1)/2,)
	print "Spectra Coverage: %.3f to %.3f MHz in %i channels (%.2f kHz/channel)" % (freq[0]/1e6, freq[-1]/1e6, nChan, (freq[-1] - freq[0])/1e3/nChan)
	print "Polarization Products: %i starting with %i" % (len(idi.pols), idi.pols[0])
	print "JD: %.3f" % jd
	
	print "Reading in FITS IDI data"
	try:
		nSets = idi.integrationCount
	except AttributeError:
		hdulist = pyfits.open(filename)
		uvData = hdulist['UV_DATA']
		jd = uvData.data['DATE'] + uvData.data['TIME']
		nSets = len(numpy.unique(jd))
		hdulist.close()
		
	autocorrsXX = []
	autocorrsYY = []
	for set in range(1, nSets+1):
		print "Set #%i of %i" % (set, nSets)
		dataDict = idi.getDataSet(set, includeAuto=True)
		
		autocorrsXX.append([None for i in xrange(32)])
		autocorrsYY.append([None for i in xrange(32)])
		for i in xrange(len(dataDict['bls']['xx'])):
			stnd1, stnd2 = dataDict['bls']['xx'][i]
			if stnd1 != stnd2:
				continue
				
			vis = dataDict['vis']['xx'][i]
			autocorrsXX[-1][stnd1] = vis
			
		for i in xrange(len(dataDict['bls']['yy'])):
			stnd1, stnd2 = dataDict['bls']['yy'][i]
			if stnd1 != stnd2:
				continue
				
			vis = dataDict['vis']['yy'][i]
			autocorrsYY[-1][stnd1] = vis
			
	autocorrsXX = numpy.array(autocorrsXX)
	autocorrsYY = numpy.array(autocorrsYY)
	
	fig = plt.figure()
	for i in xrange(32):
		ax = fig.add_subplot(8, 8, 2*i+1)
		for j in xrange(autocorrsXX.shape[0]):
			try:
				ax.plot(freq/1e6, numpy.log10(numpy.abs(autocorrsXX[j,i])**2)*10)
			except TypeError:
				pass
			ax.set_title("%i X" % (idi.stands[i],))
			ax.set_ylim([140, 190])
			
		ax = fig.add_subplot(8, 8, 2*i+2)
		for j in xrange(autocorrsYY.shape[0]):
			try:
				ax.plot(freq/1e6, numpy.log10(numpy.abs(autocorrsYY[j,i])**2)*10)
			except TypeError:
				pass
			ax.set_title("%i Y" % (idi.stands[i],))
			ax.set_ylim([140, 190])
	plt.show()
	

if __name__ == "__main__":
	numpy.seterr(all='ignore')
	main(sys.argv[1:])
