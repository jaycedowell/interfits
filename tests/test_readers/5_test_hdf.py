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

def test_hdf():
    idi = LedaFits('data/test_lalc.fitsidi')
    idi.exportHdf5('data/test_lalc.hdf')

    hdf = LedaFits('data/test_lalc.hdf')

    ok_count = 0
    # Check all header values
    h1("Testing header keywords")
    h2("Common")
    ok_count += compare_dicts(idi.h_common, hdf.h_common)
    h2("Params")
    ok_count += compare_dicts(idi.h_params, hdf.h_params)
    h2("Antenna")
    ok_count += compare_dicts(idi.h_antenna, hdf.h_antenna)
    h2("Array Geometry")
    ok_count += compare_dicts(idi.h_array_geometry, hdf.h_array_geometry)
    h2("Frequency")
    ok_count += compare_dicts(idi.h_frequency, hdf.h_frequency)
    h2("Source")
    ok_count += compare_dicts(idi.h_source, hdf.h_source)
    h2("UV DATA")
    ok_count += compare_dicts(idi.h_uv_data, hdf.h_uv_data)

    h1("Testing data tables")
    h2("Antenna")
    ok_count += compare_dicts(idi.d_antenna, hdf.d_antenna)
    h2("Array Geometry")
    ok_count += compare_dicts(idi.d_array_geometry, hdf.d_array_geometry)
    h2("Frequency")
    ok_count += compare_dicts(idi.d_frequency, hdf.d_frequency)
    h2("Source")
    ok_count += compare_dicts(idi.d_source, hdf.d_source)
    h2("UV DATA")
    ok_count += compare_dicts(idi.d_uv_data, hdf.d_uv_data)
    assert ok_count == 12
    print "PASS: All header and table data match"

    try:
        assert str(hdf.__repr__()) == str(idi.__repr__())
        print "PASS: __repr__ match"
    except:
        print "ERROR: __repr__ do not match"
        print str(hdf.__repr__())
        print str(idi.__repr__())
        print hdf.h_uv_data
        print idi.h_uv_data
        print hdf.h_array_geometry
        print idi.h_array_geometry
        raise


if __name__ == '__main__':
    test_hdf()