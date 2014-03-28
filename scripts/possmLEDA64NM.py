#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script for creating POSSM-style plots from a LEDA64-NM FITS-IDI file.
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
	print """possmLEDA64NM.py - Generate POSSM-style plots from a LEDA64-NM FITS-IDI file

Usage: possmLEDA64NM.py [OPTIONS] file

Options:
-h, --help             Display this help information
-s, --dataset          Data set to show (Default = All)
-n, --baselines        Display the first N baselines or zero for all
                       (Default = 100)
-x, --show-xx          Display the XX data (Default)
-y, --show-xy          Display the YY data
-w, --show-xy          Display the XY data
-z, --show-yx          Display the YX data

"""
	
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseConfig(args):
	config = {}
	# Command line flags - default values
	config['dataset'] = 0
	config['nBaselines'] = 100
	config['pol'] = 'xx'
	config['args'] = []
	
	# Read in and process the command line flags
	try:
		opts, arg = getopt.getopt(args, "hs:n:xywz", ["help", "dataset=", "baselines=", "show-xx", "show-yy", "show-xy", "show-yx"])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-s', '--dataset'):
			config['dataset'] = int(value)
		elif opt in ('-n', '--baselines'):
			config['nBaselines'] = int(value)
		elif opt in ('-x', '--show-xx'):
			config['pol'] = 'xx'
		elif opt in ('-y', '--show-yy'):
			config['pol'] = 'yy'
		elif opt in ('-w', '--show-xy'):
			config['pol'] = 'xy'
		elif opt in ('-z', '--show-yx'):
			config['pol'] = 'yx'
		else:
			assert False
			
	# Add in arguments
	config['args'] = arg
	
	# Return configuration
	return config


def main(args):
	# Parse the command line
	config = parseConfig(args)
	
	# Grab the filename and open it
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
	print "Final Baseline Count: %i" % (nStand * (nStand - 1) / 2,)
	print "Spectra Coverage: %.3f to %.3f MHz in %i channels (%.2f kHz/channel)" % (freq[0] / 1e6, freq[-1] / 1e6, nChan, (freq[-1] - freq[0]) / 1e3 / nChan)
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
		
	for set in range(1, nSets + 1):
		if config['dataset'] != 0 and config['dataset'] != set:
			continue
			
		print "Set #%i of %i" % (set, nSets)
		dataDict = idi.getDataSet(set, includeAuto=True)
		pols = dataDict['jd'].keys()
		nBL = len(dataDict['bls'][config['pol']])
		
		if config['nBaselines'] != 0:
			if config['nBaselines'] < nBL:
				print "Plotting the first %i baselines of %s" % (config['nBaselines'], config['pol'].upper(),)
			else:
				config['nBaselines'] = nBL
		else:
			config['nBaselines'] = nBL
			
		nFigs = int(numpy.ceil(1.0*config['nBaselines']/25))
		
		pb = ProgressBar(max=config['nBaselines'])
		i = 0
		for k in xrange(nFigs):
			fig = plt.figure()
			
			for j in range(min([config['nBaselines'], 25])):
				try:
					stnd1, stnd2 = dataDict['bls'][config['pol']][i]
					vis = dataDict['vis'][config['pol']][i]
					i += 1
				except IndexError:
					plt.draw()
					break
					
				amp = numpy.log10(numpy.abs(vis)) * 10
				phs = numpy.angle(vis) * 180 / numpy.pi
				
				ax = fig.add_subplot(10, 5, 2 * (j / 5) * 5 + j % 5 + 1)
				if ((phs + 360) % 360).std() < phs.std():
					ax.plot((phs + 360) % 360, linestyle=' ', marker='x')
					ax.set_ylim([0, 360])
				else:
					ax.plot(phs, linestyle=' ', marker='x')
					ax.set_ylim([-180, 180])
				ax.set_title('%i-%i %s' % (idi.stands[stnd1], idi.stands[stnd2], config['pol'].upper()))
				ax.set_ylabel('Phs')
				
				ax = fig.add_subplot(10, 5, 2 * (j / 5) * 5 + j % 5 + 1 + 5)
				ax.plot(amp, linestyle=' ', marker='x', color='green')
				ax.set_title('%i-%i %s' % (idi.stands[stnd1], idi.stands[stnd2], config['pol'].upper()))
				ax.set_ylabel('Amp')
				
				pb.inc(amount=1)
				if pb.amount != 0 and pb.amount % 10 == 0:
					sys.stdout.write(pb.show() + '\r')
					sys.stdout.flush()
			plt.draw()
			
		sys.stdout.write(pb.show() + '\r')
		sys.stdout.write('\n')
		sys.stdout.flush()
		plt.show()


if __name__ == "__main__":
	numpy.seterr(all='ignore')
	main(sys.argv[1:])