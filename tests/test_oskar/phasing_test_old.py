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
    uvw = LedaFits('vis.uvfits')
    uvw.leda_set_value("ARRNAM", "LWA1")
    uvw.leda_set_value("TELESCOP", "LEDA64NM")
    p_vec = uvw._compute_pointing_vec("CAS", debug=True)
    
    #uvw.phase_to_src("CAS")
    www = uvw.d_uv_data["WW"]
    
    ant1, ant2 = 1, 3
    bl_id = ant1 * 256 + ant2
    bls = uvw.d_uv_data["BASELINE"]
    ant_locs = uvw.d_array_geometry["STABXYZ"]
    
    ws = 1
    bl_vec   = ant_locs[ant1-1] - ant_locs[ant2-1]
    tg = np.dot(bl_vec, p_vec) / ledafits_config.SPEED_OF_LIGHT
    freqs = uvw.formatFreqs()

    
    print "Pointing vector: ", p_vec
    #Extract data for baseline 1-13
    for aa in range(32):
        ant1, ant2 = 13, aa + 1
        if ant1 == ant2:
            pass
        else:
            if aa < ant1:
                ws = -1
                bl_id = ant2 * 256 + ant1
            else:
                ws = -1
                bl_id = ant1 * 256 + ant2
            print bl_id
            bls = uvw.d_uv_data["BASELINE"]
            
            #ant_locs = uvw.d_array_geometry["STABXYZ"]
            freqs = uvw.formatFreqs()
            w  = 2 * np.pi * freqs # Angular freq
            d = uvw.formatStokes()[0, bls == bl_id][0]
            
            # Compute compensating phase
            bl_vec   = ant_locs[ant1-1] - ant_locs[ant2-1]
            #tg = np.dot(bl_vec, p_vec) / ledafits_config.SPEED_OF_LIGHT
            tg = www[bls == bl_id][0]
            
            print "\nBL: %i Ants: %i, %i, tg: %2.3f"%(bl_id, ant1, ant2, tg * 1e9)
            print "tg: %2.2f, bl_vec: %s"%(tg * 1e9, str(bl_vec))
            phs = np.exp(1j * w * tg * ws)
            plt.figure(1)
            plt.subplot(8, 4, aa+1)
            plt.plot(np.angle(phs), c='#cc0000')
            plt.plot(np.angle(d), c='#333333')
            plt.title("%i %i"%(ant1, ant2))
            
            plt.figure(2)
            plt.subplot(8, 4, aa+1)
            phs = np.exp(-1j * w * tg * ws)
            plt.plot(np.angle(phs * d))
            plt.title("%i %i"%(ant1, ant2))
    plt.show()
    
    