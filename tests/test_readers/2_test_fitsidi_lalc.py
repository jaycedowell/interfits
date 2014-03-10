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

def test_compare_flux():

    uvf  = LedaFits('data/test_lalc.fitsidi')

    lalc = LedaFits()
    lalc.filename = 'data/test_lalc.LA'
    lalc.readLfile(n_ant=256, n_chans=109)


    lalc_flux = lalc.d_uv_data['FLUX']
    uvf_flux  = uvf.d_uv_data['FLUX']

    lalc_bls = lalc.d_uv_data['BASELINE']
    uvf_bls  = uvf.d_uv_data['BASELINE']

    assert np.allclose(lalc_bls, uvf_bls)
    print "PASS: BASELINE IDS MATCH"

    assert lalc_flux.shape == uvf_flux.shape
    print "PASS: FLUX SHAPE MATCH"
    print lalc_flux.shape
    print uvf_flux.shape

    try:
        assert lalc_flux.dtype == uvf_flux.dtype
        print "PASS: FLUX DTYPE MATCH"
    except AssertionError:
        print "ERROR: DTYPES DO NOT MATCH"
    print lalc_flux.dtype, uvf_flux.dtype

    print "Testing flux data..."
    for row in range(0, lalc_flux.shape[0]):
        if not row % 1000:
            print "\t %i of %i"%(row, lalc_flux.shape[0])
        try:
            xxl = lalc_flux[row][::8]**2 + lalc_flux[row][1::8]**2
            yyl = lalc_flux[row][2::8]**2 + lalc_flux[row][3::8]**2
            xyl = lalc_flux[row][4::8]**2 + lalc_flux[row][5::8]**2
            yxl = lalc_flux[row][6::8]**2 + lalc_flux[row][7::8]**2

            xxu = uvf_flux[row][::8]**2  + uvf_flux[row][1::8]**2
            yyu = uvf_flux[row][2::8]**2 + uvf_flux[row][3::8]**2
            xyu = uvf_flux[row][4::8]**2 + uvf_flux[row][5::8]**2
            yxu = uvf_flux[row][6::8]**2 + uvf_flux[row][7::8]**2

            assert np.allclose(xxl, xxu)
            assert np.allclose(yyl, yyu)
            assert np.allclose(xyl, xyu)
            assert np.allclose(yxl, yxu)

        except AssertionError:
            print "ERROR: Flux values do not agree"
            print uvf_flux[row, 0:10]
            print lalc_flux[row, 0:10]
            raise

    print "PASS: FLUX DATA MATCH"

    print "Testing stokes generator..."
    xxl, yyl, xyl, yxl = lalc.formatStokes()
    xxu, yyu, xyu, yxu = uvf.formatStokes()

    try:
        assert np.allclose(np.abs(xxl), np.abs(xxu))
        assert np.allclose(np.abs(yyl), np.abs(yyu))
        assert np.allclose(np.abs(xyl), np.abs(xyu))
        assert np.allclose(np.abs(yxl), np.abs(yxu))

    except AssertionError:
        print "ERROR: Flux values do not agree"
        raise
    print "PASS: FLUX DATA STOKES FORMAT"

def test_compare_headers():

    uvf  = LedaFits('data/test_lalc.fitsidi')

    lalc = LedaFits()
    lalc.filename = 'data/test_lalc.LA'
    lalc.readLfile(n_ant=256, n_chans=109)

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

    assert ok_count == 7


if __name__ == '__main__':
    
    test_compare_flux()
    test_compare_headers()
