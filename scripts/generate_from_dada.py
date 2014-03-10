#! /usr/bin/env python
# encoding: utf-8
"""
generate_from_dada
==================

Generate FITS-IDI files from dada files. 
"""

from interfits.ledafits import *

def generate_fitsidi(filename_in, filename_out=None):
    """ Generate a fitsidi file from a dada file 
    
    filename_in (str): name of .dada file to open
    filename_out (str): output filename. Defaults to input with
                        .fitsidi extension appended
    """
    if filename_out is None:
        filename_out = os.path.split(filename_in)[0] + '.fitsidi'
    
    uvw = LedaFits(filename_in)
    #uvw.loadAntArr()
    #uvw.generateUVW(src='ZEN', use_stored=False, update_src=True)
    #uvw.apply_cable_delays()
    uvw.exportFitsidi(filename_out)
    
if __name__ == '__main__':
    import sys, os
    
    try:
        filename_in  = sys.argv[1]
        filename_out = sys.argv[2]
    except:
        print "ERROR: you must enter a filename"
        print "USAGE: python generate_fitsidi.py <filename_in> <filename_out>"
        exit()
    
    generate_fitsidi(filename_in, filename_out)

