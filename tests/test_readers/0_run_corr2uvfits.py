#! /usr/bin/env python
# encoding: utf-8
"""
pipeline.py:
Data reduction pipeline for LEDA64 Owen's Valley
created: 08 July 2013
"""

import re, sys, datetime, os, subprocess, shutil
from optparse import OptionParser

import ephem
import numpy as np


OFFSET_DELTA, INT_TIME, N_INT = 1044480000, 8.53333*26/25, 100
(latitude, longitude, elevation) = ('37.240391', '-118.28247', 1184)

def callSubprocess(args, test=False):
    """ Run a subprocess call with some printing and stuff """
    for arg in args: print arg,
    if not test: subprocess.call(args)

if __name__ == "__main__":
    
    fileroot = 'data/test_lalc'
    
    ######
    # Run corr2uvfits on output files
    ######
    print('\nRunning corr2uvfits')
    print('-------------------')
    
    #if options.lock_to_init: corr2uvfits_flags = '-l'
    #else: 
    corr2uvfits_flags = ''
    corr2uvfits_args = ['./corr2uvfits', 
        '-A', '%s,%s'%(longitude,latitude),
        '-a', os.path.join('./', fileroot+'.LA'),
        '-c', os.path.join('./', fileroot+'.LC'),
        '-o', os.path.join('./', fileroot+'.uvfits'), 
        '-H', os.path.join('uvfits_headers', 'header.tpl'),
        '-I', os.path.join('uvfits_headers', 'instr_config.tpl'),
        '-S', os.path.join('uvfits_headers', 'antenna_locations.tpl'), 
        '-f', str(1),
        corr2uvfits_flags
        ]
        
    callSubprocess(corr2uvfits_args)

    print "DONE!"