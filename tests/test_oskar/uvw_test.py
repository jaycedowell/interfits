#! /usr/bin/env python
# encoding: utf-8
"""
generate_from_dada
==================

Generate FITS-IDI files from dada. Test. 
"""

from test_main import *
import pylab as plt

if __name__ == '__main__':
    #uvw = LedaFits('uvf-zen.fitsidi')
    uvw  = LedaFits('vis_00.uvfits')
    uvw.leda_set_value("ARRNAM", "LWA1")
    uvw.leda_set_value("TELESCOP", "LEDA64NM")
    
    ww = uvw.d_uv_data["WW"]
    
    uvw.phase_to_src("CAS")
    
    ww2 = uvw.d_uv_data["WW"]
    
    plt.plot(ww)
    plt.plot(ww2)
    plt.show()
    