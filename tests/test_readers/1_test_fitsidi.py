from test_main import compare_dicts
from interfits.interfits import InterFits, PrintLog
import os
import numpy as np

pp = PrintLog()
h1 = pp.h1
h2 = pp.h2

def test_generate_fitsidi():
    """ Generate a FITS-IDI file
    """

    filename_uvf  = 'data/test_lalc.uvfits'
    filename_idi = filename_uvf.rstrip('.uvfits')+'.fitsidi'
    uvf = InterFits(filename_uvf)
    uvf.exportFitsidi(filename_idi)


def test_compare_uv2idi():
    """ Compare a uvfits file to a fits-idi file to confirm they are equivalent """

    filename_uvf  = 'data/test_lalc.uvfits'
    filename_idi = filename_uvf.rstrip('.uvfits')+'.fitsidi'
    uvf = InterFits(filename_uvf)
    idi = InterFits(filename_idi)

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
    h2("UV DATA")
    ok_count += compare_dicts(uvf.d_uv_data, idi.d_uv_data)

    assert ok_count == 12

def test_compare_idi_generated():
    """ Check that creating a fits-idi->interfits->fits-idi creates identical data. """

    idi = InterFits('data/test_lalc.fitsidi')
    idi.exportFitsidi('data/test_lalc2.fitsidi')
    idi2 = InterFits('data/test_lalc2.fitsidi')


    ok_count = 0
    # Check all header values
    h1("Testing header keywords")
    h2("Common")
    ok_count += compare_dicts(idi.h_common, idi2.h_common)
    h2("Parameters")
    ok_count += compare_dicts(idi.h_params, idi2.h_params)
    h2("Antenna")
    ok_count += compare_dicts(idi.h_antenna, idi2.h_antenna)
    h2("Array Geometry")
    ok_count += compare_dicts(idi.h_array_geometry, idi2.h_array_geometry)
    h2("Frequency")
    ok_count += compare_dicts(idi.h_frequency, idi2.h_frequency)
    h2("Source")
    ok_count += compare_dicts(idi.h_source, idi2.h_source)
    h2("UV DATA")
    ok_count += compare_dicts(idi.h_uv_data, idi2.h_uv_data)

    h1("Testing data tables")
    h2("Antenna")
    ok_count += compare_dicts(idi.d_antenna, idi2.d_antenna)
    h2("Array Geometry")
    ok_count += compare_dicts(idi.d_array_geometry, idi2.d_array_geometry)
    h2("Frequency")
    ok_count += compare_dicts(idi.d_frequency, idi2.d_frequency)
    h2("Source")
    ok_count += compare_dicts(idi.d_source, idi2.d_source)
    h2("UV DATA")
    ok_count += compare_dicts(idi.d_uv_data, idi2.d_uv_data)

    assert ok_count == 12
    os.remove('data/test_lalc2.xml')
    os.remove('data/test_lalc2.fitsidi')

if __name__ == '__main__':
    
    test_generate_fitsidi()
    test_compare_uv2idi()
    test_compare_idi_generated()


