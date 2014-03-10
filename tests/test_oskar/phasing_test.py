#! /usr/bin/env python
# encoding: utf-8
"""
generate_from_dada
==================

Generate FITS-IDI files from dada. Test. 
"""

from test_main import *
import pylab as plt

def extract_uvw(l):
    uu = l.d_uv_data["UU"]
    vv = l.d_uv_data["VV"]
    ww = l.d_uv_data["WW"]
    dd = l.formatStokes()

    return uu, vv, ww, dd

if __name__ == '__main__':
    #l = LedaFits('uvf-zen.fitsidi')
    l = LedaFits('vis.uvfits')
    l.leda_set_value("ARRNAM", "LWA1")
    l.leda_set_value("TELESCOP", "LEDA64NM")
    l.telescope = "LWA1"
    l._initialize_site()

    freqs = l.formatFreqs()
    w     = 2 * np.pi * freqs # Angular freq

    #l.apply_cable_delays()
    # Extract data for zenith phased
    uu, vv, ww, dd = extract_uvw(l)
    
    plt.plot(ww)
    plt.show()
    
    # Phase to Cygnus
    l.phase_to_src('CYG')
    uu2, vv2, ww2, dd2 = extract_uvw(l)

    # Phase back to ZEN
    l.phase_to_src('ZEN')
    uu3, vv3, ww3, dd3 = extract_uvw(l)

    # We should have ww == ww' != ww_cyg
    try:
        print "Testing phase routines for correctness...",
        #assert np.allclose(ww, ww3)
        assert not np.allclose(ww, ww2)

        assert not np.allclose(dd, dd2)
        #assert np.allclose(dd, dd3)
        #assert np.allclose(np.abs(dd), np.abs(dd2))
        #assert np.allclose(np.abs(dd), np.abs(dd3))
        print "OK"
    except AssertionError:
        print "ERROR: "
        raise

    print "BL    ANTS    WW0    WW2"
    # Extract data for baseline 13
    for aa in range(32):
        ant1, ant2 = 13, aa + 1
        if ant1 == ant2:
            pass
        else:
            if aa < ant1:
                bl_id = ant2 * 256 + ant1
            else:
                bl_id = ant1 * 256 + ant2
            bls = l.d_uv_data["BASELINE"]

            d = dd[0, bls == bl_id][0]
            d2 = dd2[0, bls == bl_id][0]

            tg = ww[bls == bl_id]
            tg2 = ww2[bls == bl_id]
            #tg3 = ww3[bls == bl_id]

            p = np.exp(1j * w * tg2)


            print "%04i  %i %i   %2.2f  %2.2f"%(bl_id, ant1, ant2, tg*1e9, tg2*1e9)
            plt.figure(1)
            plt.subplot(8, 4, aa+1)
            plt.plot(np.angle(d))
            plt.plot(np.angle(p))
            plt.title("%i %i"%(ant1, ant2))
            plt.suptitle('ZENITH PHASED')

            plt.figure(2)
            plt.subplot(8, 4, aa+1)
            plt.plot(np.angle(d2))
            plt.title("%i %i"%(ant1, ant2))
            plt.suptitle('%s PHASED'%l.source)

    plt.show()


    