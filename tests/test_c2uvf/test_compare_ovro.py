#! /usr/bin/env python
# encoding: utf-8
"""
test_compare_uvf_idi

Compare old method of creating FITS files to new method. That is, compare
lconvert + corr2uvfits to interfits.
"""

import re, sys, datetime, os, subprocess, shutil, time
from optparse import OptionParser
from colorama import Fore, Back, Style

from datetime import timedelta
import ephem
import numpy as np
import pylab as plt


from test_main import *
from dada2uvfits import callSubprocess, h1, printRed

def create_test_files(idi=True, uvfits=True):
    """ Run corr2uvfits and LedaFits on the same data 
    
    This produces a test.fitsidi, and a test.uvfits file
    """
    
    dada_filename = '/workdata/leda/data-ovro/ledaovro5_12/2014-02-26-21:44:32_0000000000000000.000000.dada'

    t1 = time.time()
    if idi:
        # Create fitsidi test file
        l = LedaFits(dada_filename)
        l.exportFitsidi('data/test_ovro.fitsidi')
    t2 = time.time()

    if uvfits:
        # Create UVFITS file
        callSubprocess(['../lconverter/dada2uvfits.py', dada_filename, '-l'])
        # Remove intermediate files
        callSubprocess(['mv', 'Zenith_b1_d20140226_utc214432/Zenith_b1_d20140226_utc214432.uvfits', 'data/test_ovro.uvfits'])
        callSubprocess(['rm', '-rf', 'Zenith_b1_d20140226_utc214432'])
    t3 = time.time()  
    
    print "INTERFITS: %2.3fs"%(t2-t1)
    print "CORR2UVFS: %2.3fs"%(t3-t2)
    
    
def clip_test_files():
    """ Extracts the first integration from both files.
    
    This creates test_uvf_clipped.fitsidi and test_idi_clipped.fitsidi,
    allowing tests to be done on a subset of the data (for speed).
    """
    
    l_idi = LedaFits('data/test_ovro.fitsidi')
    l_idi.extract_integrations(0, 1)
    l_idi.exportFitsidi('data/test_ovro_idi_clipped.fitsidi')
    
    l_idi = LedaFits('data/test_ovro.uvfits')
    l_idi.extract_integrations(0, 1)
    l_idi.exportFitsidi('data/test_ovro_uvf_clipped.fitsidi')

def compare_test_files(filename1, filename2):
    """ Compare two files to see how well they match. """
    
    l_idi = LedaFits(filename1)
    l_uvf = LedaFits(filename2)
    
    h1("Comparing files")
        
    try:
        xyz_idi = l_idi.d_array_geometry["STABXYZ"]
        xyz_uvf = l_uvf.d_array_geometry["STABXYZ"]
        assert np.allclose(xyz_uvf, xyz_idi, rtol=0.001)
        print "PASS: Station positions within tolerance"
    except AssertionError:
        for ii in range(xyz_idi.shape[0]):
            print xyz_idi[ii], xyz_uvf[ii]
        printRed("FAIL: Station positions do not match")
    
    try:
        assert np.allclose(l_idi.d_uv_data["UU"], l_uvf.d_uv_data["UU"], rtol=0.001)
        assert np.allclose(l_idi.d_uv_data["VV"], l_uvf.d_uv_data["VV"], rtol=0.001)
        assert np.allclose(l_idi.d_uv_data["WW"], l_uvf.d_uv_data["WW"], rtol=0.001)
        print "PASS: UVW coordinates within tolerance"
    except:
        printRed("FAIL: UVW coordinates do not match")  

    flux_idi = l_idi.d_uv_data["FLUX"]
    flux_uvf = l_uvf.d_uv_data["FLUX"]
            
    for ii in range(4):
        try:
            pow_idi   = np.sqrt(flux_idi[:, 2*ii::8]**2 + flux_idi[:, 2*ii+1::8]**2)
            pow_uvf   = np.sqrt(flux_uvf[:, 2*ii::8]**2 + flux_uvf[:, 2*ii+1::8]**2)
            assert np.allclose(pow_idi, pow_uvf, rtol=0.1)
            print "PASS %s of 4: Bandpass magnitude matches"%(ii + 1)

        except:
            #plt.plot(pow_idi[0])
            #plt.plot(pow_uvf[0])
            #print np.max(pow_idi[0]), np.max(pow_uvf[0])
            #plt.show()
            #raise
            printRed("FAIL %s of 4: Bandpass magnitude DOES NOT MATCH"%(ii + 1))
            if ii == 3:
                print "    NOTE: corr2uvfits does not write YX* for autocorrelations"
                print "    NOTE: Testing cross-correlations only"
                try:
                    for jj in range(flux_uvf.shape[0]):
                        bl_id     = l_idi.d_uv_data["BASELINE"]
                        num_ok = 0
                        if bl_id[jj] % 256 != bl_id[jj] / 256:
                            pow_idi   = np.sqrt(flux_idi[jj, 6::8]**2 + flux_idi[jj, 7::8]**2)
                            pow_uvf   = np.sqrt(flux_uvf[jj, 6::8]**2 + flux_uvf[jj, 7::8]**2)
                            assert np.allclose(pow_idi, pow_uvf, rtol=0.1)
                    print "PASS 4 of 4: Cross correlations bandpass magnitude match"
                except AssertionError:
                    printRed("FAIL 4 of 4: Bandpass magnitude DOES NOT MATCH for cross-correlations"%(ii + 1))


        
    test_compare_headers(l_idi, l_uvf)

def test_compare_headers(uvf, lalc):

    ok_count = 0
    # Check all header values
    h1("Testing header keywords")
    h2("Common")
    ok_count += compare_dicts(uvf.h_common, lalc.h_common)
    h2("Parameters")
    ok_count += compare_dicts(uvf.h_params, lalc.h_params)
    h2("Antenna")
    ok_count += compare_dicts(uvf.h_antenna, lalc.h_antenna)
    h2("Array Geometry")
    ok_count += compare_dicts(uvf.h_array_geometry, lalc.h_array_geometry)
    h2("Frequency")
    ok_count += compare_dicts(uvf.h_frequency, lalc.h_frequency)
    h2("Source")
    ok_count += compare_dicts(uvf.h_source, lalc.h_source)
    h2("UV DATA")
    ok_count += compare_dicts(uvf.h_uv_data, lalc.h_uv_data)

    #assert ok_count == 7
    
    h1("Testing data tables")
    h2("Antenna")
    ok_count += compare_dicts(uvf.d_antenna, lalc.d_antenna)
    h2("Array Geometry")
    ok_count += compare_dicts(uvf.d_array_geometry, lalc.d_array_geometry)
    h2("Frequency")
    ok_count += compare_dicts(uvf.d_frequency, lalc.d_frequency)
    h2("Source")
    ok_count += compare_dicts(uvf.d_source, lalc.d_source)
        
if __name__ == '__main__':
    
    do_create_test_files  = True
    do_clip_test_files    = True
    do_compare_test_files = True
    
    if do_create_test_files:
        create_test_files(idi=True, uvfits=True)
    if do_clip_test_files:
        clip_test_files()
    if do_compare_test_files:
        compare_test_files('data/test_ovro_idi_clipped.fitsidi', 'data/test_ovro_uvf_clipped.fitsidi')
        #compare_test_files('data/test_ovro.fitsidi', 'data/test_ovro.uvfits')
