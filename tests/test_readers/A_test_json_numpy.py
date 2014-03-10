from test_main import compare_dicts
from interfits.interfits import InterFits, PrintLog
from interfits.ledafits import LedaFits
from interfits.lib.json_numpy import *
import os
import numpy as np
import pylab as plt
import time

pp = PrintLog()
h1 = pp.h1
h2 = pp.h2

import pprint
ppp = pprint.PrettyPrinter(indent=4)
ppp = ppp.pprint

def test_json():    
    npDict_orig = {
        'ARRNAM'  : np.chararray([2]),
        'CHICKEN' : np.array([ii for ii in range(128)]),
        'STRING'  : 'qwertyuiop',
        'LIST'    : [1,7,2,9]
        }

    # Dump to file
    npDict_filename = 'data/jsontest.json'
    dump_json(npDict_orig, npDict_filename)

    # Load from file
    npDict_ff = load_json(npDict_filename)

    # Check the two dictionaries match
    ok_count = 0
    ok_count += compare_dicts(npDict_orig, npDict_ff)
    assert ok_count == 1

def create_json_uvfits():
    filename_uvf  = 'data/test_lalc.fitsidi'
    uvf = InterFits(filename_uvf)

    out_path = 'data/test_lalc_json'
    if not os.path.exists(out_path):
        os.mkdir('data/test_lalc_json')

    dump_json(uvf.h_antenna, os.path.join(out_path, 'h_antenna.json'))
    dump_json(uvf.h_array_geometry, os.path.join(out_path, 'h_array_geometry.json'))
    dump_json(uvf.h_common, os.path.join(out_path, 'h_common.json'))
    dump_json(uvf.h_frequency, os.path.join(out_path, 'h_frequency.json'))
    dump_json(uvf.h_params, os.path.join(out_path, 'h_params.json'))
    dump_json(uvf.h_uv_data, os.path.join(out_path, 'h_uv_data.json'))
    dump_json(uvf.h_source, os.path.join(out_path, 'h_source.json'))

    dump_json(uvf.d_antenna, os.path.join(out_path, 'd_antenna.json'))
    dump_json(uvf.d_array_geometry, os.path.join(out_path, 'd_array_geometry.json'))
    dump_json(uvf.d_frequency, os.path.join(out_path, 'd_frequency.json'))
    dump_json(uvf.d_source, os.path.join(out_path, 'd_source.json'))
    # dump_json(uvf.d_uv_data) # DONT DUMP UVDATA - VERY LARGE!

def compare_json_uvfits():
    filename_uvf  = 'data/test_lalc.fitsidi'
    filename_json = 'data/test_lalc_json/h_antenna.json'

    uvf = InterFits(filename_uvf)
    idi = InterFits(filename_json)

    ok_count = 0
    # Check all header values
    h1("Testing header keywords")
    h2("Common")
    ok_count += compare_dicts(uvf.h_common, idi.h_common)
    h2("Parameters")
    ok_count += compare_dicts(uvf.h_params, idi.h_params)
    h2("Antenna")
    ok_count += compare_dicts(uvf.h_antenna, idi.h_antenna)
    h2("Array Geometry")
    ok_count += compare_dicts(uvf.h_array_geometry, idi.h_array_geometry)
    h2("Frequency")
    ok_count += compare_dicts(uvf.h_frequency, idi.h_frequency)
    h2("Source")
    ok_count += compare_dicts(uvf.h_source, idi.h_source)
    h2("UV DATA")
    ok_count += compare_dicts(uvf.h_uv_data, idi.h_uv_data)

    h1("Testing data tables")
    h2("Antenna")
    ok_count += compare_dicts(uvf.d_antenna, idi.d_antenna)
    h2("Array Geometry")
    ok_count += compare_dicts(uvf.d_array_geometry, idi.d_array_geometry)
    h2("Frequency")
    ok_count += compare_dicts(uvf.d_frequency, idi.d_frequency)
    h2("Source")
    ok_count += compare_dicts(uvf.d_source, idi.d_source)
    #h2("UV DATA")
    #ok_count += compare_dicts(uvf.d_uv_data, idi.d_uv_data)
    assert ok_count == 11

    print "PASS: All header and table data match"
    try:
        assert str(uvf.__repr__()) == str(idi.__repr__())
        print "PASS: __repr__ match"
    except:
        print "ERROR: __repr__ do not match"
        print str(uvf.__repr__())
        print str(idi.__repr__())
        print uvf.h_uv_data
        print idi.h_uv_data
        print uvf.h_array_geometry
        print idi.h_array_geometry
        raise

    print "PASS: Comparison test passed"

if __name__ == '__main__':

    test_json()
    create_json_uvfits()
    compare_json_uvfits()