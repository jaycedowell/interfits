from test_main import compare_dicts
from interfits.interfits import InterFits, PrintLog
from interfits.ledafits import LedaFits
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

def create_json_uvfits():
    filename_uvf  = 'data2/CygA_109ch.fitsidi'
    uvf = InterFits(filename_uvf)

    out_path = 'data2/CygA_json'
    if not os.path.exists(out_path):
        os.mkdir('data2/CygA_json')

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
    filename_uvf  = 'data2/CygA_109ch.fitsidi'
    filename_json = 'data2/CygA_json/h_antenna.json'

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

#def remove_miriad_lalc():
#    uvf = LedaFits('data/test_lalc.fitsidi')
#    uvf.remove_miriad_baselines()
#    uvf.exportFitsidi('data/test_lalc_255.fitsidi', '../config/config.xml')

if __name__ == '__main__':

    create_json_uvfits()
    compare_json_uvfits()