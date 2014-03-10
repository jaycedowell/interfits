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


def test_lalc():
    """ Check that LA->interfits->LA->interfits creates identical data. """

    lalc = LedaFits()
    lalc.filename = 'data/test_lalc.LA'
    lalc.readLfile(n_ant=256, n_chans=109)

    lalc.exportFitsidi('data/test_lalc_direct.fitsidi')
    uvf  = LedaFits('data/test_lalc_direct.fitsidi')

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

    h1("Testing data tables")
    h2("Antenna")
    ok_count += compare_dicts(uvf.d_antenna, lalc.d_antenna)
    h2("Array Geometry")
    ok_count += compare_dicts(uvf.d_array_geometry, lalc.d_array_geometry)
    h2("Frequency")
    ok_count += compare_dicts(uvf.d_frequency, lalc.d_frequency)
    h2("Source")
    ok_count += compare_dicts(uvf.d_source, lalc.d_source)

    assert ok_count == 11

    try:
        assert repr(lalc) == repr(uvf)
    except AssertionError:
        print "ERROR: __repr__ outputs do not match"
        print repr(uvf)
        print repr(lalc)

    os.remove('data/test_lalc_direct.fitsidi')

def test_lalc_json():

    lalc = LedaFits()
    lalc.filename = 'data/test_lalc.LA'
    lalc.readLfile(n_ant=256, n_chans=109)

    lalc.readJson('data/test_lalc_json/d_antenna.json')
    lalc.generateUVW()
    uu = lalc.d_uv_data["UU"]
    vv = lalc.d_uv_data["VV"]
    ww = lalc.d_uv_data["WW"]

    #print uu
    #print len(uu)

    print "Checking data...",
    try:
        assert len(uu) == len(vv) == len(ww)
    except:
        print len(uu), len(vv), len(ww)
        raise
    try:
        assert len(uu) == len(lalc.d_uv_data["BASELINE"])
    except AssertionError:
        print len(uu), len(lalc.d_uv_data["BASELINE"])
        raise
    try:
        assert len(uu) == len(lalc.d_uv_data["FLUX"])
    except:
        print len(uu), len(lalc.d_uv_data["BASELINE"])
        raise
    print "OK"

    lalc.remove_miriad_baselines()

    h2("Exporting sewed data...")
    lalc.exportFitsidi('data/test_lalc_uvw.fitsidi', '../config/config.xml')

    #print "Plotting..."
    #import pylab as plt
    #plt.plot(uu[0:32896], vv[0:32896], 'ro')
    #plt.show()

def test_compare_idi_lalc_json():
    lalc = LedaFits()
    lalc.filename = 'data/test_lalc.LA'
    lalc.readLfile(n_ant=256, n_chans=109)
    lalc.readJson('data/test_lalc_json/d_antenna.json')
    lalc.generateUVW()

def test_compare_cyg_lalc():
    lalc = LedaFits('data/test_lalc_uvw.fitsidi')
    uvf  = LedaFits('data/CygA.fitsidi')

    ok_count = 0
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

    h1("Testing data tables")
    h2("Antenna")
    ok_count += compare_dicts(uvf.d_antenna, lalc.d_antenna)
    h2("Array Geometry")
    ok_count += compare_dicts(uvf.d_array_geometry, lalc.d_array_geometry)
    h2("Frequency")
    ok_count += compare_dicts(uvf.d_frequency, lalc.d_frequency)
    h2("Source")
    ok_count += compare_dicts(uvf.d_source, lalc.d_source)
    h2("UV DATA")
    ok_count += compare_dicts(uvf.d_uv_data, lalc.d_uv_data)

def remove_miriad_lalc():
    uvf = LedaFits('data/test_lalc.fitsidi')
    uvf.remove_miriad_baselines()
    uvf.exportFitsidi('data/test_lalc_255.fitsidi')

if __name__ == '__main__':

    remove_miriad_lalc()
    test_lalc()
    test_lalc_json()
    test_compare_idi_lalc_json()
    test_compare_cyg_lalc()